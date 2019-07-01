(in-package :leap.off)

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
  (let* ((parts (core:split (subseq line 1) " "))
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
      (t (let ((line (read-line fin eof-error-p eof)))
           (if (eq line eof)
               eof
               (parse-off-data-line line header)))))))

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

(defun read-off-unit-part (fin &optional (eof-error-p t) eof)
  "* Arguments
fin : A Stream
eof-error-p : bool
eof : An object
* Description
Read the enntry of the unit from fin.  
Use eof-error-p and eof as in read."
  (flet ((my-split (string &key (delimiterp #'(lambda (x) (eq x #\.))))
           (loop :for beg = (position-if-not delimiterp string)
              :then (position-if-not delimiterp string :start (1+ end))
              :for end = (and beg (position-if delimiterp string :start beg))
              :when beg :collect (subseq string beg end)
              :while end)))
    (multiple-value-bind (header body)
        (read-off-data-block fin eof-error-p eof)
      (if header 
        (let* ((header-parts (my-split (name header)))
               (unit-name (second header-parts))
               (unit-part (fourth header-parts)))
          (values body unit-name unit-part))
        (values) ; Return nothing if we hit the eof
        ))))

(defvar *residue-types* '(("w" . :solvent)
                          ("p" . :protein)
                          ("n" . :nucleic)
                          ("s" . :saccharide))
  "Residue types in OFF files - defined in leap residues.h:84 RESTYPExxx")

(defun assemble-residues (read-residues)
  (let ((residues (make-array (length read-residues) :adjustable nil :fill-pointer 0)))
    (dotimes (ri (length read-residues))
      (let* ((read-residue (elt read-residues ri))
             (residue (chem:make-residue (intern (elt read-residue 0) :keyword))))
        (let ((residue-assoc (assoc (elt read-residue 4) *residue-types* :test #'string-equal)))
          (when residue-assoc (chem:setf-residue-type residue (cdr residue-assoc))))
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

(defun read-entire-off-file (fin &optional eof-error-p eof)
  (let ((result (make-hash-table :test #'equal)))
    (loop for part = (multiple-value-list (read-off-unit-part fin eof-error-p eof))
          if part
            do (destructuring-bind (body unit-name unit-part)
                   part
                 (let ((part-hash-table (gethash unit-name result)))
                   (unless part-hash-table
                     (setf part-hash-table (make-hash-table :test #'equal)))
                   (setf (gethash unit-part part-hash-table) body)
                   (setf (gethash unit-name result) part-hash-table)))
          else
            do (return-from read-entire-off-file result))))

(defun assign-molecule-name-and-type (molecule)
  (if (= (chem:content-size molecule) 1)
      (let ((residue (chem:content-at molecule 0)))
        ;; Single residues molecules get their residue name
        ;; and their residue type
        (chem:set-name molecule (chem:get-name residue))
        (chem:setf-molecule-type molecule (chem:residue-type residue)))
      (progn
        ;; Now things get tricky - this is a multi- (or zero-) residue molecule
        ;; we don't know what to name it or what type to give it
        ;; Let's bail out for now
        )))

(defun read-off-unit-with-name (entry-parts expected-name)
  "Read a leap UNIT from an OFF file.
Return (values chem:aggregate read-residueconnect read-residues).
The read-residueconnect and read-residues can be used to contruct a chem:topology
if the caller wants to do that."
  (let ((read-atoms (gethash "atoms" entry-parts))
        (read-atoms-pert-info (gethash "atomspertinfo" entry-parts))
        (read-bound-box (gethash "boundbox" entry-parts))
        (read-childsequence (gethash "childsequence" entry-parts))
        (read-connect (gethash "connect" entry-parts))
        (read-connectivity (gethash "connectivity" entry-parts))
        (read-hierarchy (gethash "hierarchy" entry-parts))
        (read-name (gethash "name" entry-parts))
        (read-positions (gethash "positions" entry-parts))
        (read-residueconnect (gethash "residueconnect" entry-parts))
        (read-residues (gethash "residues" entry-parts))
        (read-residuesPdbSequenceNumber (gethash "residuesPdbSequenceNumber" entry-parts))
        (read-solventcap (gethash "solventcap" entry-parts))
        (read-velocities (gethash "velocities" entry-parts)))
    (let ((atoms (assemble-atoms read-atoms)))
      (set-positions read-positions atoms)
      (when read-connectivity
        (form-bonds read-connectivity atoms))
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
                      (let ((molecule (gethash id molecule-id-map)))
                        (assign-molecule-name-and-type molecule)
                        (chem:add-matter aggregate molecule)))
                    sorted-molecule-ids)
              ;; If there is a bounding box then use it
              (when (> (elt (elt read-bound-box 0) 0) 0)
                (let ((xbox (elt (elt read-bound-box 2) 0))
                      (ybox (elt (elt read-bound-box 3) 0))
                      (zbox (elt (elt read-bound-box 4) 0)))
                  (chem:set-property aggregate :bounding-box (list xbox ybox zbox))))
              (values aggregate read-connect read-residues))))))))

(defun translate-off-object (unit-name unit connect residues)
  "If the OFF object is a unit containing a single residue then turn it into a topology.
Otherwise it's an Aggregate. Return the object."
  (declare (symbol unit-name)
           (chem:aggregate unit)
           (vector connect residues))
  (if (eq (length residues) 1)
      (let ((res (let ((r (chem:content-at (chem:content-at unit 0) 0)))
                   (chem:set-name r unit-name)
                   r)))
        (multiple-value-bind (topology constitution)
            (cando:make-simple-topology-from-residue res)
          (let ((in-plug-idx (1- (elt (elt connect 0) 0))) ;; connect atoms 1- index
                (out-plug-idx (1- (elt (elt connect 1) 0)))) ;; connect atoms 1- index
            (when (>= in-plug-idx 0)
              (let* ((in-plug-atom (chem:content-at res in-plug-idx))
                     (in-plug (chem:make-in-plug :-default nil (chem:atom-name in-plug-atom) :single-bond)))
                (chem:add-plug topology (chem:get-name in-plug) in-plug)
                (chem:add-plug constitution (chem:get-name in-plug) in-plug)))
            (when (>= out-plug-idx 0)
              (let* ((out-plug-atom (chem:content-at res out-plug-idx))
                     (out-plug (chem:make-out-plug  :+default nil nil (chem:atom-name out-plug-atom) :single-bond)))
                (chem:add-plug topology (chem:get-name out-plug) out-plug)
                (chem:add-plug constitution (chem:get-name out-plug) out-plug))))
          topology))
      unit))

(defun read-off-lib (fin)
  "Read an OFF library into a hash-table indexed by the name of each object in the OFF file.
Return the hash-table."
  (declare (stream fin))
  (let ((index (read-off-lib-index fin))
        (entire-file (read-entire-off-file fin))
        (object-ht (make-hash-table)))
    (maphash (lambda (entry-name-string entry-parts)
               (let ((entry-name (intern entry-name-string *package*)))
                 (multiple-value-bind (unit residue-connect residues)
                     (read-off-unit-with-name entry-parts entry-name)
                   (let ((object (translate-off-object entry-name unit residue-connect residues)))
                     (setf (gethash entry-name object-ht) object)))))
             entire-file)
    object-ht))

(defun bug () (print "In bug"))

(defun load-off (filename)
  "* Arguments
- filename : Pathname
* Description
Load the OFF file containing forms into new-leap."
  (let ((absolute-filename (leap.core:search-path filename)))
    (with-open-file (fin absolute-filename :direction :input)
      (let ((ht (leap.off:read-off-lib fin))
            names)
        (maphash (lambda (name form)
                   (leap.core:register-variable name form)
                   (push name names))
                 ht)
        (nreverse names)))))
