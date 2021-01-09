(in-package :leap.pdb)

(defparameter *map-pdb-names-to-topology-names* (make-hash-table :test #'equal)
  "Map (residue name . terminus) from pdb files to form names")

(defun lookup-topology-using-pdb-name-and-context (name context)
  "context is one of :head, :tail or nil?"
  (multiple-value-bind (topology-name foundp)
      (gethash (cons name context) *map-pdb-names-to-topology-names*)
    (if foundp
        (cando:lookup-topology topology-name t)
        (let ((topology (cando:lookup-topology (intern (string name)) nil)))
          topology))))

(defun add-pdb-res-map (mappings)
  "* Arguments
- mappings :: A list of two or three element lists.
* Description
The Name Map is used to map residue names read from PDB files to variable
names within LEaP.  The LIST is a LIST of LISTs:
      { {sublist} {sublist} ... }
Each sublist contains two or three entries to add to the Name Map:
      { [terminalflag] PDBName LEaPVar }
where the PDBName will be mapped to the LEaPVar. The terminalflag indicates
the special cases of terminal residues: allowable values are 0 or :head for beginning
residues (N-terminal for proteins, 5' for nucleic acids) and 1 or :tail for ending
residues (C-terminal for proteins, 3' for nucleic acids).  If the
terminalflag is given, the PDBName->LEaPVar name map will only be applied
for the appropriate terminal residue.  The `leaprc' file included with
the distribution contains default mappings."
  (unless (core:proper-list-p mappings)
    (error "You must provide a list of lists of mappings"))
  (labels ((ensure-symbol (name &optional (package :keyword))
             (if (symbolp name) name (intern name package)))
           (do-add-map (res-name term-sym var-name)
             (let* ((res-sym (ensure-symbol res-name))
                    (var-sym (ensure-symbol var-name *package*))
                    (key (cons res-sym term-sym))
                    (var-old (gethash key *map-pdb-names-to-topology-names*)))
               (when (and var-old (not (eq var-old var-sym)))
                 (warn "The PDB residue name/terminalflag ~a/~a was changed to ~a from ~a"
                       res-sym term-sym var-sym var-old))
               (setf (gethash key *map-pdb-names-to-topology-names*) var-sym))))
    (dolist (mapping mappings)
      (cond
        ((= (length mapping) 3)
         (destructuring-bind (terminus res-name var-name)
             mapping
           (let ((term-sym (case terminus
                             (0 :head)
                             (1 :tail)
                             (:head :head)
                             (:tail :tail)
                             (otherwise (error "Illegal terminalflag ~a" terminus)))))
             (do-add-map res-name term-sym var-name))))
        ((= (length mapping) 2)
         (destructuring-bind (res-name var-name)
             mapping
           (do-add-map res-name :main var-name)))
        ((null mapping) #|Do nothing|#)
        (t (error "Illegal argument ~a" mapping))))))
                                             
(defparameter *pdb-atom-map* (make-hash-table :test #'eq)
  "Old PDB files may contain atom names that have been renamed in the common convention.
For instance old DNA ribose rings were named C1*, C2*, C3*, C4* etc.
The PDB renamed them to C1', C2', C3', C4' respectively.
This hash table maps old names to new names.
There can only be one mapping from an old name and if a name is
ever replaced then a warning is generated.")

(defun add-pdb-atom-map (name-list)
  "* Arguments
- name-list :: A list of two element lists.
* Description
The atom Name Map is used to try to map atom names read from PDB files
to atoms within residues when the atom name in the PDB file does
not match a atom in the residue.  This enables PDB files to be read
in without extensive editing of atom names.  The LIST is a LIST of LISTs:
      { {sublist} {sublist} ... }
where each sublist is of the form
      { \"OddAtomName\" \"LibAtomName\" }
Many `odd' atom names can map to one `standard' atom name, but any single
odd atom name maps only to the last standard atom name it was mapped to."
  (unless (core:proper-list-p name-list)
    (error "You must provide a list of lists of names"))
  (dolist (old-new name-list)
    (destructuring-bind (old-name new-name)
        old-new
      (unless (and (= (length old-new) 2) old-name new-name)
        (error "Illegal argument - you must provide (\"old-name\" \"new-name\") pairs - you gave: ~a" old-new))
      (flet ((ensure-symbol (name)
               (if (symbolp name) name (intern name :keyword))))
        (let* ((old-sym (ensure-symbol old-name))
               (new-sym (ensure-symbol new-name))
               (old-map (gethash old-sym *pdb-atom-map*)))
          (when (and old-map (not (eq old-map new-sym)))
            (warn "The PDB atom name ~a was mapped to ~a and is now being mapped to ~a"
                  old-sym old-map new-sym))
          (setf (gethash old-sym *pdb-atom-map*) new-sym))))))

(defstruct pdb-id name context)
           
(defclass pdb-residue ()
  ((chain-id :initarg :chain-id :accessor chain-id)
   (res-seq :initarg :res-seq :accessor res-seq)
   (i-code :initarg :i-code :accessor i-code)
   (topology :initform nil :accessor topology)
   (name :initarg :name :accessor name)
   (residue-name :initarg :residue-name :accessor residue-name)
   (atom-names :initarg :atom-names :accessor atom-names)
   (context :initform nil :initarg :context :accessor context)
   (atom-serial-first :initform nil :initarg :atom-serial-first :accessor atom-serial-first)
   (atom-serial-last :initform nil :accessor atom-serial-last))
  (:documentation
     "Keep track of residue and context (:head :main :tail)"))

  
(defmethod print-object ((obj pdb-residue) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a -> ~a :atoms ~a" (name obj) (residue-name obj) (atom-names obj))))

(defclass disulphide ()
  ((chain-id1 :initarg :chain-id1 :reader chain-id1)
   (res-seq1 :initarg :res-seq1 :reader res-seq1)
   (chain-id2 :initarg :chain-id2 :reader chain-id2)
   (res-seq2 :initarg :res-seq2 :reader res-seq2)))

(defclass connect ()
  ((from :initarg :from :reader from)
   (to :initarg :to :reader to)))

(defun find-pdb-residue (scanner res-seq chain-id)
  (declare (optimize (debug 3)))
  (block found
    (loop for seq in (sequences scanner)
       do (loop for res in seq
             when (and (= res-seq (res-seq res))
                       (eq chain-id (chain-id res)))
             do (return-from found res)))))

(defclass pdb-reader ()
  ((stream :initarg :stream :reader stream)
   (current-residue-number :initform nil :initarg :current-residue-number :accessor current-residue-number)
   (current-residue :initform nil :initarg :current-residue :accessor current-residue)
   (current-i-code :initform nil :initarg :current-i-code :accessor current-i-code)))

(defclass scanner ()
  ((big-z :accessor big-z) ; If the pdb file has an extra big z coordinate
   (seen-residues :initform (make-hash-table) :accessor seen-residues)
   (sequences :initform nil :accessor sequences)
   (matrices :initform nil  :accessor matrices)
   (disulphides :initform nil :initarg :disulphides :accessor disulphides)
   (unknown-residues :initform (make-hash-table :test #'equal) :accessor unknown-residues)
   (connects :initform nil :accessor connects)))

(defun scanner-unknown-residue-key (res)
  (cons (name res) (context res)))

(defclass pdb-scanner (pdb-reader)
  ((previous-atom-serial :initform nil :accessor previous-atom-serial)
   (previous-residue :initform nil :accessor previous-residue)
   (current-residue-name :initform nil :accessor current-residue-name)
   (current-reverse-sequence :initform nil :accessor current-reverse-sequence
                             :documentation "Accumulate the current sequence in reverse order")
   (reversed-sequences :initform nil :accessor reversed-sequences
                       :documentation "A list of sequences (proper order) in pushed in reverse order")
   (scanner :initarg :scanner :reader scanner))
  (:documentation " Keep track of sequence and matrix info while scanning a PDB file"))

(defclass pdb-atom-reader (pdb-reader)
  ((sequences-index :initform 0 :accessor sequences-index)
   (sequence-index :initform 0 :accessor sequence-index)
   (sequences :initarg :sequences :accessor sequences)
   (previous-residue :initarg :previous-residue :accessor previous-residue)
   (previous-topology :initarg :previous-topology :accessor previous-topology)
   (aggregate :initform (chem:make-aggregate) :reader aggregate)
   (molecules :initform nil :initarg :molecules :accessor molecules)
   (finish-molecule :initform nil :accessor finish-molecule)
   (molecule :initform nil :initarg :molecule :accessor molecule)
   (current-residue :initform nil :initarg :current-residue :accessor current-residue)
   (current-topology :initform nil :initarg :current-topology :accessor current-topology)
   (connect-atoms :initarg :connect-atoms :reader connect-atoms)
   (current-residue-name :initform nil :initarg :current-residue-name :accessor current-residue-name)
   (serial-to-atom :initarg :serial-to-atom :reader serial-to-atom)
   ))

(defun ensure-molecule (reader chain-id)
  (unless (molecule reader)
;;;    (format t "Getting molecule at ~a~%" (sequences-index reader))
    (setf (molecule reader) (elt (molecules reader) (sequences-index reader)))
    #+(or)(format *debug-io* "Added molecule ~a to aggregate~%" (molecule reader)))
  (molecule reader))

(defun new-residue-p (residue-name residue-number i-code pdb-atom-reader)
  (if (or (not (eq residue-number (current-residue-number pdb-atom-reader)))
          (not (eq i-code (current-i-code pdb-atom-reader))))
      t
      (if (not (eq residue-name (current-residue-name pdb-atom-reader)))
          (progn
            #+(or)(warn "The residue name changed from ~a to ~a but the residue number ~a didnt change"
                  (current-residue-name pdb-atom-reader)
                  residue-name
                  residue-number)
            t)
          nil)))

;;; Create a special readtable that treats quote #\' as a normal char
(defparameter *quote-readtable* (copy-readtable))
(set-syntax-from-char #\' #\A *quote-readtable* *readtable*)

(defun parse-fixed-double (num-str)
  (let ((neg (or (position #\- num-str)))
        (dot (position #\. num-str)))
    (multiple-value-bind (whole whole-end)
        (parse-integer num-str :start (if neg (1+ neg) 0) :end dot)
      (multiple-value-bind (frac frac-end)
          (parse-integer num-str :start (1+ whole-end))
        (let ((div (float (expt 10 (- frac-end (1+ whole-end))))))
      (let ((positive (+ (float whole) (/ frac div))))
        (if neg (- positive) positive)))))))

(defun parse-atom-name (line)
  (let ((name-string (string-trim '(#\space) (subseq line 12 16))))
    (intern name-string :keyword)))

(defun parse-line (line pdb-atom-reader read-atom-pos &optional big-z check-big-z)
  "* Arguments
- line :: String
- pdb-atom-reader :: pdb-atom-reader
- read-atom-pos :: boolean
- big-z :: boolean
- check-big-z :: boolean
* Description
Parse the PDB line in line and populate fields in pdb-atom-reader with the data.
If we need to read the atom position pass T for read-atom-pos. If it's a big-z PDB
file (the Z coordinate needs an extra floating point digit) then pass big-z = T.
If you want to check if it's a big-z file then pass check-big-z = T."
  (let ((core:*read-hook* nil))
    (labels ((my-read-from-string (str)
               ;; bypass eclector because it's too slow
               (with-input-from-string (sin str)
                 (core:fast-read sin)))
             (read-one-char (line start)
               (let ((ch (elt line start)))
                 (cond
                   ((char= ch #\space)
                    nil)
                   ((alpha-char-p ch)
                    (intern (string ch) :keyword))
                   ((digit-char-p ch)
                    (- (char-int ch) (char-int #\0)))
                   (t nil)))))
      (let* ((*package* (find-package :keyword))
             (head (string-right-trim '(#\space) (subseq line 0 (min 6 (length line))))))
        (cond 
          ((or (string= head "ATOM") (string= head "HETATM"))
           ;; If we are checking if it's a big-z pdb file then
           ;; test if the 52nd char is #\. - if it is - it's a big-z file
           (when check-big-z
             (when (char= (elt line 52) #\.)
               (format t "~a~%" line)
               (format t "(elt line 52) = ~a~%" (elt line 52))
               (setf big-z t)))
           (values
            (list* (if (string= head "ATOM") :atom :hetatm)
                   (parse-integer line :start 6 :end 12) ; atom-serial
                   (let ((*readtable* *quote-readtable*))
                     (parse-atom-name line))
                   (read-one-char line 16) ; alt-loc
                   (intern (string-trim '(#\space) (subseq line 17 20)) :keyword) ; residue-name
                   (read-one-char line 21)                ; chainid
                   (parse-integer line :start 22 :end 26) ; res-seq
                   (read-one-char line 26)                ; chainid
                   (if read-atom-pos
                       #+(or)(list 
                              (my-read-from-string line t nil :start 30 :end 38) ; x
                              (my-read-from-string line t nil :start 38 :end 46) ; y
                              (my-read-from-string line t nil :start 46 :end 54)) ; z
                       (list
                        (parse-fixed-double (subseq line 30 38))
                        (parse-fixed-double (subseq line 38 46))
                        (parse-fixed-double (subseq line 46 (if big-z 55 54))))
                       nil))
            big-z))
          ((string= head "TER") (list :ter))
          ((string= head "MTRIX" :start1 0 :end1 5)
           (list (intern head :keyword)
                 (my-read-from-string line t nil :start 7 :end 10)
                 (my-read-from-string line t nil :start 10 :end 20)
                 (my-read-from-string line t nil :start 20 :end 30)
                 (my-read-from-string line t nil :start 30 :end 40)
                 (my-read-from-string line t nil :start 45 :end 55)
                 (my-read-from-string line t nil :start 59 :end 60)))
          ((string= head "CONECT")
           (list* :conect
                  (parse-integer line :start 7 :end 11)
                  (loop for (start . end) in '((12 . 16) (17 . 21 ) (22 . 26 ) (27 . 31 ) (32 . 36 ) (37 . 41 ) (42 . 46 ) (47 . 51 ) (52 . 56 ) (57 . 61 ))
                        until (> start (length line))
                        collect (parse-integer line :start start :end end :junk-allowed t))))
          ((string= head "SSBOND")
           (list :ssbond
                 (read-one-char line 15)
                 (parse-integer (string-trim '(#\space #\tab) (subseq line 17 21)))
                 (read-one-char line 29)
                 (parse-integer (string-trim '(#\space #\tab) (subseq line 31 35)))))
          (t nil))))))



(defun pop-sequence-pdb-residue (pdb-atom-reader &optional (errorp t) error-val)
  "* Arguments
- pdb-atom-reader :: A pdb-reader
- errorp :: bool
- errorp :: T
* Description
Pop and return the next residue in the sequence.
If we run out of residues and errorp is NIL then return error-val
otherwise signal an error."
  (unless (< (sequence-index pdb-atom-reader) (length (elt (sequences pdb-atom-reader)
                                                           (sequences-index pdb-atom-reader))))
    (incf (sequences-index pdb-atom-reader))
    (setf (sequence-index pdb-atom-reader) 0)
    (return-from pop-sequence-pdb-residue error-val))
  (prog1
      (elt (elt (sequences pdb-atom-reader) (sequences-index pdb-atom-reader))
           (sequence-index pdb-atom-reader))
    (incf (sequence-index pdb-atom-reader))))

(defun finish-previous-sequence (pdb &optional (assign-tail t))
  "* Arguments
- pdb :: a pdb-scanner
- assign-tail :: bool
* Description
Finish off the current sequence and prepare to start another.
If assign-tail is NIL then don't assign the last read residue as
a tail - this is used when partial sequences are loaded to avoid
treating :main residues as :tail residues when that wont make sense and
create more problems."
  (when (current-reverse-sequence pdb)
    (when assign-tail
      (let ((tail-residue (car (current-reverse-sequence pdb))))
        (setf (context tail-residue) :tail)
        (try-to-assign-topology tail-residue (scanner pdb))))
    (push (nreverse (current-reverse-sequence pdb)) (reversed-sequences pdb))
    (setf (current-reverse-sequence pdb) nil)))

(defun try-to-assign-topology (res scanner)
  (let* ((topology (lookup-topology-using-pdb-name-and-context (residue-name res) (context res))))
    (if topology
        (progn
          (assert (not (symbolp topology)))
          (setf (topology res) topology))
        (let ((key (scanner-unknown-residue-key res)))
          (format t "Could not immediately identify topology for ~a~%" key)
          (unless (gethash key (unknown-residues scanner))
            (format t "Could not identify topology for ~a~%" key)
            (setf (gethash key (unknown-residues scanner)) t))))))

(defgeneric scan (system field scanner line))

(defun pdb-scanner-read-line (pdb eof-errorp eof)
  "* Arguments
- pdb : A pdb-scanner
* Description
Scan through the PDB file and build a high level description
of the contents that will be used later to read the PDB file.
The records that are read from the PDB file are:
ATOM - used to build a list of residues and whether they are :first :main :last.
TER  - Used to identify the ends of chains.
MTRIX- Used to build a list of matrices."
  (let* ((line (read-line (stream pdb) eof-errorp eof))
         (scanner (scanner pdb))
         (check-big-z t))
    (unless (eq line eof)
      (multiple-value-bind (line-data big-z)
          (parse-line line pdb nil nil check-big-z)
        (when line-data
          #+(or)(format t "Read:  ~a~%" line-data)
          (case (car line-data)
            ((:atom :hetatm)
             (if (and check-big-z (slot-boundp scanner 'big-z))
                 (unless (eq (big-z scanner) big-z)
                   (error "The big-z status of the PDB file changed"))
                 (progn
                   (setf (big-z scanner) big-z
                         check-big-z nil) ; only check big-z the first ATOM record
                   (when big-z (warn "This is a BIG-Z PDB file (The Z-coordinate contains an extra digit)"))))
             (destructuring-bind (head atom-serial atom-name alt-loc residue-name chain-id res-seq i-code)
                 line-data
               ;; Deal with the current line of data
               (when (new-residue-p residue-name res-seq i-code pdb)
                 #+(or)(progn
                         (format t "Starting new residue residue-name: ~a name-of-current-residue: ~a info: ~a~%" residue-name (current-residue-name pdb) (list res-seq i-code pdb))
                         (when (null (current-reverse-sequence pdb))
                           (format t "Starting new molecule~%")))
                 (setf (current-residue-number pdb) res-seq
                       (current-residue-name pdb) residue-name
                       (current-i-code pdb) i-code)
                 (let* ((context-guess (if (null (current-reverse-sequence pdb))
                                           :head
                                           :main))
                        (new-residue (make-instance 'pdb-residue
                                                    :chain-id chain-id
                                                    :res-seq res-seq
                                                    :i-code i-code
                                                    :name residue-name
                                                    :residue-name residue-name
                                                    :atom-names nil
                                                    :context context-guess
                                                    :atom-serial-first atom-serial)))
                   (setf (current-residue pdb) new-residue)
                   (unless (gethash residue-name (seen-residues (scanner pdb)))
                     ;;(format t "Creating hash-table for ~a~%" residue-name)
                     (setf (gethash residue-name (seen-residues (scanner pdb))) (make-hash-table)))
                   #+(or)(try-to-assign-topology new-residue (scanner pdb))
                   (when (previous-residue pdb)
                     (setf (atom-serial-last (previous-residue pdb)) (previous-atom-serial pdb)))
                   (push new-residue (current-reverse-sequence pdb))
                   ;; Set things up for the next new residue
                   (setf (previous-residue pdb) new-residue)))
               (let ((ht (gethash residue-name (seen-residues (scanner pdb)))))
                 (unless ht (error "Could not find residue-name ~a in the seen-residues: ~a" residue-name (alexandria:hash-table-keys (seen-residues (scanner pdb)))))
                 (setf (gethash atom-name ht) t))
               (push atom-name (atom-names (current-residue pdb)))
               ;; Now set things up for the next atom record
               (setf (previous-atom-serial pdb) atom-serial)))
            (:ter
             (finish-previous-sequence pdb))
            (:ssbond
             (destructuring-bind (head chain-id1 res-seq1 chain-id2 res-seq2)
                 line-data
               (push (make-instance 'disulphide
                                    :chain-id1 chain-id1
                                    :res-seq1 res-seq1
                                    :chain-id2 chain-id2
                                    :res-seq2 res-seq2) (disulphides scanner))))
            (:conect
             (destructuring-bind (head from &rest to)
                 line-data
               (push (make-instance 'connect
                                    :from from
                                    :to to)
                     (connects scanner))))
            (:mtrix1
             (destructuring-bind (head serial x y z t)
                 line-data
               (let ((matrix (geom:make-matrix-identity)))
                 (geom:at-row-col-put matrix 0 0 x)
                 (geom:at-row-col-put matrix 0 1 y)
                 (geom:at-row-col-put matrix 0 2 z)
                 (geom:at-row-col-put matrix 0 3 t)
                 (push matrix (matrices scanner)))))
            (:mtrix2
             (destructuring-bind (head serial x y z t)
                 line-data
               (let ((matrix (car (matrices scanner))))
                 (geom:at-row-col-put matrix 1 0 x)
                 (geom:at-row-col-put matrix 1 1 y)
                 (geom:at-row-col-put matrix 1 2 z)
                 (geom:at-row-col-put matrix 1 3 t))))
            (:mtrix3
             (destructuring-bind (head serial x y z t)
                 line-data
               (let ((matrix (car (matrices scanner))))
                 (geom:at-row-col-put matrix 2 0 x)
                 (geom:at-row-col-put matrix 2 1 y)
                 (geom:at-row-col-put matrix 2 2 z)
                 (geom:at-row-col-put matrix 2 3 t))))
            (otherwise nil)))))
    line))

(defgeneric adjust-residue-name (residue-name pdb-residue scanner system))

(defmethod adjust-residue-name (residue-name pdb-residue scanner system)
  )

(defmethod adjust-residue-name ((residue-name (eql :HIS)) pdb-residue scanner system)
  (let ((has-hd1 (member :hd1 (atom-names pdb-residue)))
        (has-he2 (member :he2 (atom-names pdb-residue))))
    (cond
      ((and has-hd1 has-he2)
       (setf (residue-name pdb-residue) :HIP))
      (has-hd1
       (setf (residue-name pdb-residue) :HID))
      (has-he2
       (setf (residue-name pdb-residue) :HIE)))))

(defun adjust-residue-names (scanner system)
  (loop for residues in (sequences scanner)
        do (loop for residue in residues
                 do (adjust-residue-name (residue-name residue) residue scanner system))))

(defun split-solvent (scanner)
  (let (solvents solutes)
    (loop for group in (sequences scanner)
          if (every (lambda (residue) (eq (name residue) :HOH)) group)
            do (setf solvents group)
          else do (push group solutes))
    (when solvents
      (loop for solvent in solvents
            do (setf (context solvent) :main)
            do (push (list solvent) solutes)))
    (setf (sequences scanner) (nreverse solutes))))

(defun assign-topologys (scanner)
  (loop for group in (sequences scanner)
        do (loop for residue in group
                 do (try-to-assign-topology residue scanner))))


(defun atom-names-match-topology (atom-names topology)
  (let* ((constitution (chem:get-constitution topology))
         (constitution-atoms (chem:get-constitution-atoms constitution))
         (constitution-atoms-as-list (chem:constitution-atoms-as-list constitution-atoms)))
    (let ((topology-atom-names (loop for ca in constitution-atoms-as-list
                                     collect (chem:atom-name ca))))
      (if (= (length atom-names) (length topology-atom-names))
          (loop for nm in topology-atom-names
                unless (member nm atom-names)
                  do (return-from atom-names-match-topology
                       (values nil topology-atom-names)))
          (values nil topology-atom-names))
      (values t topology-atom-names))))


(define-condition pdb-read-error (error)
  ((messages :initarg :messages :initform nil :reader messages))
  (:report
   (lambda (condition stream)
     (format stream "~a" (messages condition)))))

(defun sorted-mismatch (names-a names-b)
  (let ((names-in-both nil)
        unique-a
        unique-b)
    (loop for name in names-a
          when (member name names-b)
            do (push name names-in-both))
    (loop for name in names-a
          when (not (member name names-in-both))
            do (push name unique-a))
    (loop for name in names-b
          when (not (member name names-in-both))
            do (push name unique-b))
    (let ((sorted-unique-a (sort unique-a #'string< :key #'string))
          (sorted-unique-b (sort unique-b #'string< :key #'string)))
      (values sorted-unique-a sorted-unique-b))))
          
        
(defun validate-scanner (scanner &key ignore-missing-topology)
  (let* (saw-problems
         (*print-pretty* nil)
         (messages 
           (with-output-to-string (sout)
             (loop for sequence in (sequences scanner)
                   do (loop for pdb-res in sequence
                            for res-name = (residue-name pdb-res)
                            for res-seq = (res-seq pdb-res)
                            for chain-id = (chain-id pdb-res)
                            for atom-names = (atom-names pdb-res)
                            for topology = (topology pdb-res)
                            do (if topology
                                   (multiple-value-bind (match-p topology-atom-names)
                                       (atom-names-match-topology atom-names topology)
                                     (unless match-p
                                       (multiple-value-bind (unique-pdb unique-amber)
                                           (sorted-mismatch atom-names topology-atom-names)
                                       (format sout "Residue ~a ~a ~a - atom names mismatch~%   PDB: ~a and ~a in AMBER residue ~a~%" res-name chain-id res-seq unique-pdb unique-amber (chem:get-name topology) ))
                                       (setf saw-problems t)))
                                   (when (not ignore-missing-topology)
                                     (format sout "Residue ~a ~a ~a - is unknown~%" res-name chain-id res-seq)
                                     (setf saw-problems t))))))))
    (when saw-problems
      (error 'pdb-read-error :messages messages))
    ))
                                  


                                  
(defun scan-pdb-stream (fin &key progress system ignore-missing-topology)
  (let ((pdb-scanner (make-instance 'pdb-scanner
                                    :stream fin
                                    :scanner (make-instance 'scanner)))
        (bar (if progress
                 (cando:make-progress-bar :message "Scanned" :total (file-length fin) :divisions 100 :on progress)
                 nil)))
    (restart-case
        (progn
          (loop for x = (pdb-scanner-read-line pdb-scanner nil :eof)
                do (when progress (cando:progress-advance bar (file-position fin)))
                until (eq x :eof))
          (finish-previous-sequence pdb-scanner)
          (when progress (cando:progress-done bar)))
      (read-partial-sequence ()
        :report "Construct a partial structure"
        ;; Drop the current residue and proceed
        (format t "Terminating sequence read")
        (when (current-reverse-sequence pdb-scanner)
          (pop (current-reverse-sequence pdb-scanner)))
        (finish-previous-sequence pdb-scanner nil)))
    (let* ((sequences (nreverse (reversed-sequences pdb-scanner)))
           (scanner (scanner pdb-scanner)))
      (setf (sequences scanner) sequences)
      (finish-scanner scanner)
      (adjust-residue-names scanner system)
      (assign-topologys scanner)
      (split-solvent scanner)
      (let ((result-scanner (build-sequence scanner system)))
        (validate-scanner result-scanner :ignore-missing-topology ignore-missing-topology)
        result-scanner))))


(defun scan-pdb (filename &key progress system ignore-missing-topology)
  "* Arguments
- filename : A pathname
* Description
Scan the PDB file and use the ATOM records to build a list of residue sequences and matrices.
* Return
values residue-sequences matrices"
  (with-open-file (fin filename :direction :input)
    (scan-pdb-stream fin :progress progress
                         :system system
                         :ignore-missing-topology ignore-missing-topology)
    ))

(defun warn-of-unknown-topology (res scanner)
  (let* ((name (name res))
         (context (context res))
         (key (cons name context)))
    (unless (gethash key (unknown-residues scanner))
      (setf (gethash key (unknown-residues scanner)) t)
      (warn "No topology was found for ~a~%" key))))

(defun finish-scanner (scanner)
  (declare (optimize (debug 3)))
  ;; Reverse the matrices
  (setf (matrices scanner) (nreverse (matrices scanner)))
  ;; Form the disulfide bonds
  (mapc (lambda (dis)
          (let ((res1 (find-pdb-residue scanner (res-seq1 dis)
                                        (chain-id1 dis)))
                (res2 (find-pdb-residue scanner (res-seq2 dis)
                                        (chain-id2 dis))))
            (let* ((top1 (lookup-topology-using-pdb-name-and-context :CYX (context res1)))
                   (top2 (lookup-topology-using-pdb-name-and-context :CYX (context res2))))
              ;; If CYX topologies are available - use those
              (when top1 (setf (topology res1) top1))
              (when top2 (setf (topology res2) top2))
              (unless (topology res1)
                (warn "Could not identify a proper CYX residue for the first half of disulphide bond ~a" dis))
              (unless (topology res2)
                (warn "Could not identify a proper CYX residue for the second half of disulphide bond ~a" dis)))))
        (disulphides scanner))
  scanner)

(defun calculate-residue-sequence-number (res-seq i-code)
  "If i-code is not nil and as a string it is a single digit,
then multiply res-seq by 10 and add the i-code digit to it."
  (declare (integer res-seq))
  (cond
    ((null i-code)
     res-seq)
    ((integerp i-code)
     (+ (* res-seq 10) i-code))
    ((digit-char-p (elt (string i-code) 0))
     (+ (* res-seq 10) (parse-integer (string i-code))))
    (t res-seq)))

(defun read-and-process-line (reader eof-errorp eof big-z)
  (declare (optimize (debug 3)))
  "* Arguments
- reader : pdb-atom-reader
- eof-errorp : boolean
- eof : T
- big-z :: boolean
* Description 
Read the next line from the PDB file and process it, 
filling in the information in the pdb-atom-reader.
Pass big-z parse-line to tell it how to process the z-coordinate." 
  (let* ((line (read-line (stream reader) eof-errorp eof)))
    (if (eq line eof)
        eof
        (let ((line-data (parse-line line reader t big-z)))
          (when line-data
            ;;          (format t "Reading line: ~a~%" line-data)
            #+(or)(format *debug-io* "pdb.lisp read-and-process-line line-data -> ~s~%" line-data)
            (case (car line-data)
              ((:atom :hetatm)
               (destructuring-bind (head atom-serial atom-name alt-loc residue-name chain-id res-seq i-code x y z)
                   line-data
                 (when (new-residue-p residue-name res-seq i-code reader)
                   (setf (current-residue-number reader) res-seq
                         (current-i-code reader) i-code)
                   (let ((pdb-residue (pop-sequence-pdb-residue reader nil nil)))
                     (unless pdb-residue
                       (setf (molecule reader) nil)
;;;                       (format t "Starting a new molecule~%")
                       (setf (current-residue reader) nil)
                       (setf (current-topology reader) nil)
                       (setf pdb-residue (pop-sequence-pdb-residue reader nil nil))
;;;                       (format t "popped first residue: ~a~%" pdb-residue)
                       (ensure-molecule reader chain-id))
                     (setf (previous-topology reader) (current-topology reader)
                           (previous-residue reader) (current-residue reader)
                           (current-residue-name reader) residue-name
                           (current-topology reader) (topology pdb-residue))
                     (unless (eq (name pdb-residue) residue-name)
                       (error "There is a mismatch between the pdb-residue ~a and the expected residue ~a for ~a"
                              pdb-residue residue-name line-data))
                     (let ((prev-top (previous-topology reader))
                           (cur-top (current-topology reader)))
                       (if cur-top
                           ;; There is a topology - use it
                           (let ((cur-res (chem:build-residue-single-name cur-top)))
                             #+(or)(format *debug-io* "built residue with name ~a using topology ~a~%" cur-res cur-top)
                             (chem:set-id cur-res (calculate-residue-sequence-number res-seq i-code))
                             (setf (current-residue reader) cur-res)
                             (let ((prev-res (previous-residue reader)))
                               (chem:add-matter (ensure-molecule reader chain-id) cur-res)
                               (when prev-res
                                 (chem:connect-residues prev-top
                                                        prev-res
                                                        :+default
                                                        cur-top
                                                        cur-res
                                                        :-default))))
                           ;; There is no topology, create an empty residue
                           (let ((cur-res (chem:make-residue residue-name)))
                             (setf (current-residue reader) cur-res)
                             (chem:add-matter (ensure-molecule reader chain-id) cur-res) cur-res)))))
                 (let ((atom (or (chem:content-with-name-or-nil (current-residue reader) atom-name)
                                 (chem:content-with-name-or-nil (current-residue reader) (gethash atom-name *pdb-atom-map*)))))
                   (cond
                     (atom  ; there is an atom - fill it
                      ;; If there is a serial-to-atom hash-table
                      ;; then update it
                      (when (serial-to-atom reader)
                        (let ((in-connect (gethash atom-serial (serial-to-atom reader))))
                          (when in-connect
                            (setf (gethash atom-serial (serial-to-atom reader)) atom))))
                      (chem:set-id atom atom-serial)
                      (chem:set-position atom (geom:vec x y z))
                      (chem:setf-needs-build atom nil))
                     #+(or)((current-topology reader) ; there is a topology
                      (warn "Loaded atom ~a but amber form ~a does not recognize it"
                            (list atom-serial atom-name residue-name chain-id res-seq i-code)
                            (current-topology reader)))
                     (t                ; Add a new atom to the residue
                      (let ((atom (chem:make-atom atom-name (chem:element-from-atom-name-string (string atom-name)))))
                        (chem:set-id atom atom-serial)
                        (chem:add-matter (current-residue reader) atom)
                        (chem:set-position atom (geom:vec x y z))
                        (chem:setf-needs-build atom nil)
                        (when (serial-to-atom reader)
                          (setf (gethash atom-serial (serial-to-atom reader)) atom))))))))
              (:ter nil)
              (otherwise nil)))
          #+(or)(when (finish-molecule reader)
                  (setf (molecule reader) nil
                        (current-residue reader) nil
                        (current-topology reader) nil
                        (current-residue-number reader) nil
                        (current-i-code reader) nil))
          line))))

(defun connect-atoms-hash-table (scanner)
  (when (connects scanner)
    (let ((ht (make-hash-table :test #'eql)))
      (loop for connect in (connects scanner)
         do (setf (gethash (from connect) ht) :from)
         do (mapc (lambda (to) (setf (gethash to ht) :to))
                  (to connect)))
      ht)))


(defgeneric build-sequence (scanner system))

(defmethod build-sequence (scanner system)
  scanner)

(defparameter *serial-to-atoms* nil)

(defun molecules-from-sequences (sequences)
  (loop for sequence in sequences
        for index from 0
        for sequence-name = (cond
                              ((= (length sequence) 1)
                               (intern (format nil "~a_~d" (string (name (first sequence))) index) :keyword))
                              (t (intern (concatenate 'string "mol" (format nil "_~d" index)) :keyword)))
        collect (let ((mol (chem:make-molecule sequence-name)))
                  #+(or)(format *debug-io* "Creating molecule ~a for sequence ~a~%" mol sequence)
                  mol)))


(defun load-pdb-stream (fin &key filename scanner progress system ignore-missing-topology)
  (let* ((scanner (if scanner
                      scanner
                      (let ((new-scanner (build-sequence (scan-pdb-stream fin :progress progress
                                                                          :ignore-missing-topology ignore-missing-topology) system)))
                        (file-position fin 0)
                        new-scanner)))
         (serial-to-atoms (connect-atoms-hash-table scanner)))
    (setq *serial-to-atoms* serial-to-atoms)
    (let ((pdb-atom-reader (make-instance 'pdb-atom-reader :stream fin
                                                           :sequences (mapcar (lambda (seq)
                                                                                (copy-list seq))
                                                                              (sequences scanner))
                                                           :molecules (molecules-from-sequences (sequences scanner))
                                                           :connect-atoms (connects scanner)
                                                           :serial-to-atom serial-to-atoms)))
      (let ((bar (if progress
                     (cando:make-progress-bar :message "Load pdb" :total (file-length fin) :divisions 100 :on progress)
                     nil))
            (atom-table (make-hash-table :test #'eql)))
        (restart-case
            (progn
              (setf (molecule pdb-atom-reader) nil)
              (loop for x = (read-and-process-line pdb-atom-reader nil :eof (big-z scanner))
                    do (when bar (cando:progress-advance bar (file-position fin)))
                    until (eq x :eof))
              (when progress (format t "Loaded pdb~%")))
          (ran-out-of-sequence ()
            :report "Ran out of sequence while filling residues"
            (format t "Continuing with partial sequence~%")))
        (when bar (cando:progress-done bar)))
      (loop for connect in (connects scanner)
            for from-atom = (gethash (from connect) serial-to-atoms)
            for to-atoms = (mapcar (lambda (c) (gethash c serial-to-atoms)) (to connect))
            do (mapc (lambda (to-atom)
                       (when (and (typep from-atom 'chem:atom)
                                  (typep to-atom 'chem:atom)
                                  (not (chem:is-bonded-to from-atom to-atom)))
                         (chem:bond-to from-atom to-atom :single-bond)))
                     to-atoms))
      (let ((unbuilt-heavy-atoms 0)
            (aggregate (aggregate pdb-atom-reader)))
        (loop for molecule in (molecules pdb-atom-reader)
              do (chem:add-matter aggregate molecule))
        (chem:map-atoms
         nil
         (lambda (a)
           (when (and (not (eq (chem:get-element a) :H))
                      (chem:needs-build a))
             (incf unbuilt-heavy-atoms)))
         aggregate)
        (when (> unbuilt-heavy-atoms 0)
          (when progress
            (format t "There are ~a unbuilt heavy atoms - building them " unbuilt-heavy-atoms))
          (cando:simple-build-unbuilt-atoms aggregate))
        (when progress
          (format t "Building missing hydrogens~%"))
        (let ((built (cando:build-unbuilt-hydrogens aggregate)))
          (when progress
            (format t "Built ~d missing hydrogens~%" built)))
;;;            (cando:maybe-join-molecules-in-aggregate aggregate)
;;;            (cando:maybe-split-molecules-in-aggregate aggregate)
        (setf aggregate (classify-molecules aggregate system))
        (let ((name-only (if filename
                             (pathname-name (pathname filename))
                             "unknown")))
          (chem:set-name aggregate (intern name-only *package*)))
        (values aggregate scanner pdb-atom-reader)))))


(defun load-pdb (filename &key scanner progress system ignore-missing-topology)
  "    variable = loadPdb filename
      STRING                       _filename_

Load a Protein Databank format file with the file name _filename_.
The sequence numbers of the RESIDUEs will be determined from the
order of residues within the PDB file ATOM records.  For each
residue in the PDB file, LEaP searches the variables currently
defined for variable names that match the residue name.  If a
match is found, then the contents of the variable are copied into
the UNIT created for the PDB structure.  If no PDB `TER' card
separates the current residue from the previous one, a bond is
created between the connect1 ATOM of the previous residue and the
connect0 atom of the new one.  As atoms are read from the ATOM
records, their coordinates are written into the correspondingly
named ATOMs within the residue being built.  If the entire residue
is read and it is found that ATOM coordinates are missing, then
external coordinates are built from the internal coordinates that
were defined in the matching UNIT (residue) variable.  This allows
LEaP to build coordinates for hydrogens and lone pairs which are not
specified in PDB files.
"
  (let ((filename (leap.core:ensure-path filename)))
    (with-open-file (fin filename :direction :input)
      (load-pdb-stream fin :filename filename
                           :scanner scanner
                           :progress progress
                           :system system
                           :ignore-missing-topology ignore-missing-topology)
      )))

(defgeneric classify-molecules (aggregate system))

(defmethod classify-molecules (aggregate system)
  (let ((mol-id 0))
    (cando:do-molecules (molecule aggregate)
      (chem:set-id molecule (incf mol-id))
      (cond
        ((member (chem:get-name molecule) (list :hoh :wat))
         (chem:setf-molecule-type molecule :solvent))))
    aggregate))

#|
(defun read-line (pdb-atom-reader)
  
(defun read-atom (pdb-atom-reader)
  (

(get-macro-character #\A)
(let ((*readtable* (copy-readtable)))
  (set-syntax-from-char #\' #\A)
  (read-from-string "AB'C"))

(read-from-string "A'1")
(defun parse-atom (line)
  (
|#
