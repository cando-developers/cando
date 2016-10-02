
(in-package :amber.off)


(defclass off-header ()
  ((name :initarg :name :accessor name)
   (type :initarg :type :accessor type)
   (fields :initarg :fields :accessor fields)))

(define-condition bad-header ()
  ((line :initarg :line :accessor line)))

(define-condition wrong-off-unit-name (error)
  ((header :initarg :header :accessor header)
   (unit-name :initarg :unit-name :accessor unit-name)
     (expected-unit-name :initarg :expected-unit-name :accessor expected-unit-name)))

(define-condition wrong-off-unit-part (error)
  ((header :initarg :header :accessor header)
   (unit-part :initarg :unit-part :accessor unit-part)
     (expected-part :initarg :expected-part :accessor expected-part)))

(defun parse-off-header (line)
  (unless (char= (elt line 0) #\!)
    (error 'bad-header :line line))
  (let* ((parts (core:split-at-white-space (subseq line 1)))
         (name (elt parts 0))
         (type (intern (string-upcase (nth 1 parts)) :keyword))
         (fields (cond
                   ((eq type :table)
                    (loop for index from 2 below (length parts) by 2
                       collect (cons (intern (string-upcase (elt parts index)) :keyword)
                                     (intern (string-upcase (elt parts (1+ index))) :keyword))))
                   ((eq type :array)
                    (list (cons (intern (string-upcase (elt parts 2)) :keyword)
                                :only)))
                   ((eq type :single)
                    (list (cons (intern (string-upcase (elt parts 2)) :keyword)
                                :only)))
                   (t (error "Unknown OFF file entry type: ~a" type)))))
    (make-instance 'off-header
                   :name name
                   :type type
                   :fields fields)))

(defun parse-off-data-line (line header)
  (let* ((data (make-array (length (fields header)) :element-type t)))
    (with-input-from-string (sin line)
      (loop for x from 0 below (length (fields header))
         do (setf (elt data x) (read sin))))
    data))


(defun read-off-data-line (fin header eof-error-p eof)
  (let ((first-char (read-char fin eof-error-p eof)))
    (cond
      ((eq first-char eof) eof)
      ((eq first-char #\!)
       (unread-char first-char fin)
       nil)
      (t (parse-off-data-line (read-line fin eof-error-p eof) header)))))

(defun read-off-data-block (fin &optional (eof-error-p t) eof)
  (let ((header-line (read-line fin eof-error-p eof)))
    (if (eq header-line eof)
        eof
        (let ((header (parse-off-header header-line)))
          (when (eq header eof)
            (return-from read-off-data-block eof))
          (let ((data (make-array 0 :element-type t :fill-pointer t)))
            (loop for x = (read-off-data-line fin header nil :eof)
               until (or (null x) (eq x :eof))
               do (vector-push-extend x data))
            (values header data))))))

(defun read-off-lib-index (fin)
  (multiple-value-bind (header body)
      (read-off-data-block fin nil nil)
    (if header
        (progn
          (unless (string= "!index" (name header))
            (error "Could not find !index in off file ~a" fin))
          (map 'list (lambda (x) (elt x 0)) body))
        (error "Could not read off file ~a" fin))))

(defun read-off-unit-part (fin &optional expected-unit-name expected-part (eof-error-p t) eof)
  "* Arguments
fin : A Stream
expected-unit-name : A string
expected-part : A string
eof-error-p : bool
eof : An object
* Description
Read the enntry of the unit from fin.  
If expected-unit-name and expected-part are provided then check that the read parts
match what is expected - otherwise signal an error.   Use eof-error-p and eof as in read."
  (flet ((my-split (string &key (delimiterp #'(lambda (x) (eq x #\.))))
           (loop :for beg = (position-if-not delimiterp string)
              :then (position-if-not delimiterp string :start (1+ end))
              :for end = (and beg (position-if delimiterp string :start beg))
              :when beg :collect (subseq string beg end)
              :while end)))
    (multiple-value-bind (header body)
        (read-off-data-block fin eof-error-p eof)
      (when header 
        (let* ((header-parts (my-split (name header)))
               (unit-name (second header-parts))
               (unit-part (fourth header-parts)))
          (when expected-unit-name
            (unless (equal unit-name expected-unit-name)
              (error 'wrong-off-unit-name
                     :header (name header)
                     :unit-name unit-name
                     :expected-unit-name expected-unit-name)))
          (when expected-part
            (unless (equal unit-part expected-part)
              (error 'wrong-off-unit-part
                     :header (name header)
                     :unit-part unit-part
                     :expected-part expected-part)))
          (values body unit-name unit-part))))))

(defun assemble-residues (read-residues)
  (let ((residues (make-array (length read-residues) :adjustable nil :fill-pointer 0)))
    (dotimes (ri (length read-residues))
      (let* ((read-residue (elt read-residues ri))
             (residue (chem:make-residue (intern (elt read-residue 0) :keyword))))
        (vector-push residue residues)))
    residues))

(defun assemble-atoms (read-atoms)
  (let ((atoms (make-array (length read-atoms) :adjustable nil :fill-pointer 0)))
    (dotimes (ai (length read-atoms))
      (let* ((read-atom (elt read-atoms ai))
             (atom-name (intern (elt read-atom 0) :keyword))
             (atom-type (intern (elt read-atom 1) :keyword))
             (atomic-number (elt read-atom 6))
             (atomic-charge (elt read-atom 7))
             (element (cond
                        ((string= atom-name "EP" :end1 2) :LP)
                        ((< atomic-number 0)
                         (chem:element-from-atom-name-string (string atom-name)))
                        (t (chem:element-for-atomic-number atomic-number))))
             (atom (chem:make-atom atom-name element)))
        (chem:set-type atom atom-type)
        (chem:set-charge atom atomic-charge)
        (vector-push atom atoms)))
    atoms))

(defun assemble-hierarchy (read-hierarchy residues atoms)
  (let ((unit (chem:make-aggregate)))
    (flet ((object-of-type (type index)
             (cond
               ((equal type "U") unit)
               ((equal type "M") (error "Implement way to get molecule"))
               ((equal type "R") (elt residues (1- index)))
               ((equal type "A") (elt atoms (1- index)))
               (t (error "Unknown object type ~a" type)))))
      (let ((unit (chem:make-aggregate)))
        (dotimes (hi (length read-hierarchy))
          (let* ((hierarchy (elt read-hierarchy hi))
                 (outer-type (elt hierarchy 0))
                 (outer-index (elt hierarchy 1))
                 (inner-type (elt hierarchy 2))
                 (inner-index (elt hierarchy 3)))
            ;; Destructure hierarchy like below types = "U", "M", "R", "A" ...
            ;; #( "U" 0 "R" 1 )
            ;; #( "R" 1 "A" 1 )
            ;; #( "R" 1 "A" 2 )
            (chem:add-matter (object-of-type outer-type outer-index)
                             (object-of-type inner-type inner-index))))
        unit))))

(defun form-bonds (read-connectivity atoms)
  (dotimes (bi (length read-connectivity))
    (let* ((connect (elt read-connectivity bi))
           (b1 (elt connect 0))
           (b2 (elt connect 1))
           (bo (elt connect 2))
           (a1 (elt atoms (1- b1))) ; For some reason I started counting at 1
           (a2 (elt atoms (1- b2))))
      (unless (= bo 1)
        (error "There is a bond order that is not 1"))
      (chem:bond-to a1 a2 :single-bond))))

(defun set-positions (read-positions atoms)
  (loop for atom across atoms
     for position across read-positions
     do (let ((pos-vec (geom:vec (elt position 0) (elt position 1) (elt position 2))))
          (chem:set-position atom pos-vec))))

(defun read-off-unit-with-name (fin expected-name &optional eof-error-p eof)
  "Read a leap UNIT from an OFF file.
Return (values chem:aggregate read-residueconnect read-residues).
The read-residueconnect and read-residues can be used to contruct a chem:topology
if the caller wants to do that."
  (let ((read-atoms (read-off-unit-part fin expected-name "atoms" eof-error-p eof))
        (read-atoms-pert-info (read-off-unit-part fin expected-name "atomspertinfo" eof-error-p eof))
        (read-bound-box (read-off-unit-part fin expected-name "boundbox" eof-error-p eof))
        (read-childsequence (read-off-unit-part fin expected-name "childsequence" eof-error-p eof))
        (read-connect (read-off-unit-part fin expected-name "connect" eof-error-p eof))
        (read-connectivity (read-off-unit-part fin expected-name "connectivity" eof-error-p eof))
        (read-hierarchy (read-off-unit-part fin expected-name "hierarchy" eof-error-p eof))
        (read-name (read-off-unit-part fin expected-name "name" eof-error-p eof))
        (read-positions (read-off-unit-part fin expected-name "positions" eof-error-p eof))
        (read-residueconnect (read-off-unit-part fin expected-name "residueconnect" eof-error-p eof))
        (read-residues (read-off-unit-part fin expected-name "residues" eof-error-p eof))
        (read-residuesPdbSequenceNumber (read-off-unit-part fin expected-name "residuesPdbSequenceNumber" eof-error-p eof))
        (read-solventcap (read-off-unit-part fin expected-name "solventcap" eof-error-p eof))
        (read-velocities (read-off-unit-part fin expected-name "velocities" eof-error-p eof)))
    (let ((atoms (assemble-atoms read-atoms)))
      (set-positions read-positions atoms)
      (form-bonds read-connectivity atoms)
      (let ((residues (assemble-residues read-residues)))
        ;; Assemble the hierarchy but ignore the returned aggregate
        ;; because it doesn't have molecules, only residues
        (assemble-hierarchy read-hierarchy residues atoms)
        (let ((atom-id-map (make-hash-table :test #'eql))
              (residue-id-map (make-hash-table :test #'eql)))
          (loop for atom across atoms
             for index from 0
             do (setf (gethash index atom-id-map) atom))
          (loop for residue across residues
             for index from 0
             do (setf (gethash index residue-id-map) residue))
          (let ((molecule-id-map (cando:build-molecules-from-atom-connectivity atom-id-map residue-id-map))
                molecule-ids)
            (maphash (lambda (id molecule)
                       (push id molecule-ids))
                     molecule-id-map)
            (let ((sorted-molecule-ids (sort molecule-ids #'<=))
                  (aggregate (chem:make-aggregate)))
              (mapc (lambda (id)
                      (chem:add-matter aggregate (gethash id molecule-id-map)))
                    sorted-molecule-ids)
              (values aggregate read-residueconnect read-residues))))))))

(defun translate-off-object (unit-name unit residueconnect residues)
  "If the OFF object is a unit containing a single residue then turn it into a topology.
Otherwise it's an Aggregate. Return the object."
  (declare (symbol unit-name)
           (chem:aggregate unit)
           (vector residue-connect residues))
  (if (eq (length residues) 1)
      (let* ((res (chem:content-at (chem:content-at unit 0) 0))
             (topology (cando::make-topology-from-residue res))
             (in-plug-idx (1- (elt (elt residueconnect 0) 0))) ;; connect atoms 1- index
             (out-plug-idx (1- (elt (elt residueconnect 0) 1)))) ;; connect atoms 1- index
        (when (>= in-plug-idx 0)
          (let* ((in-plug-atom (chem:content-at res in-plug-idx))
                 (in-plug (chem:make-in-plug :-default nil (chem:atom-name in-plug-atom) :single-bond)))
            (chem:add-plug topology (chem:get-name in-plug) in-plug)))
        (when (>= out-plug-idx 0)
          (let* ((out-plug-atom (chem:content-at res out-plug-idx))
                 (out-plug (chem:make-out-plug  :+default nil nil (chem:atom-name out-plug-atom) :single-bond)))
            (chem:add-plug topology (chem:get-name out-plug) out-plug)))
        (warn "translated off object ~a into a chem:topology")
        topology)
      unit))

(defun read-off-lib (fin)
  "Read an OFF library into a hash-table indexed by the name of each object in the OFF file.
Return the hash-table."
  (declare (stream fin))
  (let ((index (read-off-lib-index fin))
        (object-ht (make-hash-table)))
    (loop for object-name-str in index
       for object-name = (intern object-name-str :keyword)
       for (unit residueconnect residues) = (multiple-value-list (read-off-unit-with-name fin object-name-str nil :eof))
       until (eq unit :eof)
       do (let ((object (translate-off-object object-name unit residueconnect residues)))
            (setf (gethash object-name object-ht) object)))
    object-ht))

(defun load-off (filename)
  "* Arguments
- filename : Pathname
* Description
Load the OFF file containing forms into new-leap."
  (with-open-file (fin filename :direction :input)
    (let ((ht (amber.off:read-off-lib fin)))
      (maphash (lambda (name form)
                 (amber:register-variable name form))
               ht))))
