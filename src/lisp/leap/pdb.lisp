(in-package :leap.pdb)

(defclass line-filter ()
  ((chain-ids :initform :all :initarg :chain-ids :accessor chain-ids)
   (ignore-residues :initform (make-hash-table) :initarg :ignore-residues :accessor ignore-residues)
   (rename-residues :initform nil :initarg :rename-residues :accessor rename-residues)
   (ignore-atoms :initform (make-hash-table) :initarg :ignore-atoms :accessor ignore-atoms)
   (rename-atoms :initform nil :initarg :rename-atoms :accessor rename-atoms)))

(defvar *map-pdb-names-to-topology-names* (make-hash-table :test #'equal)
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

(defvar *pdb-atom-map* (make-hash-table :test #'eq)
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
  ((hetatmp :initform nil :initarg :hetatmp :accessor hetatmp)
   (start-lineno :initarg :start-lineno :accessor start-lineno)
   (sequences-index :initarg :sequences-index :accessor sequences-index)
   (sequence-index :initarg :sequence-index :accessor sequence-index)
   (chain-id :initarg :chain-id :accessor chain-id)
   (res-seq :initarg :res-seq :accessor res-seq)
   (i-code :initarg :i-code :accessor i-code)
   (topology :initform nil :accessor topology)
   (name :initarg :name :accessor name)
   (original-residue-name :initarg :original-residue-name :accessor original-residue-name)
   (residue-name :initarg :residue-name :accessor residue-name)
   (atom-names :initarg :atom-names :accessor atom-names)
   (context :initform nil :initarg :context :accessor context)
   (atom-serial-first :initform nil :initarg :atom-serial-first :accessor atom-serial-first)
   (atom-serial-last :initform nil :accessor atom-serial-last))
  (:documentation
     "Keep track of residue and context (:head :main :tail)"))


(defmethod print-object ((obj pdb-residue) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ln_~a|~a:~a chid_~a ~s ~a :atoms ~a"
            (name obj)
            (start-lineno obj)
            (sequences-index obj)
            (sequence-index obj)
            (chain-id obj)
            (context obj)
            (residue-name obj)
            (atom-names obj))))

(defclass disulphide ()
  ((chain-id1 :initarg :chain-id1 :reader chain-id1)
   (res-seq1 :initarg :res-seq1 :reader res-seq1)
   (chain-id2 :initarg :chain-id2 :reader chain-id2)
   (res-seq2 :initarg :res-seq2 :reader res-seq2)))

(defmethod print-object ((dis disulphide) stream)
  (print-unreadable-object (dis stream :type t)
    (format stream ":chain-id1 ~s :res-seq1 ~s :chain-id2 ~s :res-seq2 ~s"
            (chain-id1 dis)
            (res-seq1 dis)
            (chain-id2 dis)
            (res-seq2 dis))))

(defclass connect ()
  ((from :initarg :from :reader from)
   (to :initarg :to :reader to)))

(defmethod print-object ((obj connect) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (from obj) (to obj))))

(defun find-pdb-residue (scanner res-seq chain-id)
  (declare (optimize (debug 3)))
  (block found
    (loop for seq in (sequences scanner)
       do (loop for res in seq
             when (and (= res-seq (res-seq res))
                       (eq chain-id (chain-id res)))
             do (return-from found res)))))

(defclass pdb-reader ()
  ((current-residue-number :initform nil :initarg :current-residue-number :accessor current-residue-number)
   (current-residue :initform nil :initarg :current-residue :accessor current-residue)
   (current-i-code :initform nil :initarg :current-i-code :accessor current-i-code)))

(defclass scanner ()
  ((big-z :accessor big-z) ; If the pdb file has an extra big z coordinate
   (seen-residues :initform (make-hash-table) :accessor seen-residues)
   (sequences :initform nil :accessor sequences)
   (lineno-to-sequence-index :initform (make-hash-table :test 'eql) :accessor lineno-to-sequence-index)
   (matrices :initform nil  :accessor matrices)
   (disulphides :initform nil :initarg :disulphides :accessor disulphides)
   (unknown-residues :initform (make-hash-table :test #'equal) :accessor unknown-residues)
   (connects :initform nil :accessor connects)
   (problems :initform nil :initarg :problems :accessor problems)
   (line-filter :initform (make-instance 'line-filter) :initarg :filter :accessor line-filter)))

(defun scanner-unknown-residue-key (res)
  (cons (name res) (context res)))

(defclass pdb-scanner (pdb-reader scanner)
  ((filename :initarg :filename :accessor filename)
   (previous-atom-serial :initform nil :accessor previous-atom-serial)
   (previous-residue :initform nil :accessor previous-residue)
   (previous-hetatmp :initform nil :accessor previous-hetatmp)
   (previous-chain-id :initform nil :accessor previous-chain-id)
   (sequences-index :initform 0 :accessor sequences-index)
   (sequence-index :initform 0 :accessor sequence-index)
   (current-residue-name :initform nil :accessor current-residue-name)
   (current-hetatmp :initform nil :accessor current-hetatmp)
   (current-chain-id :initform nil :accessor current-chain-id)
   (current-reverse-sequence :initform nil :accessor current-reverse-sequence
                             :documentation "Accumulate the current sequence in reverse order")
   (reversed-sequences :initform nil :accessor reversed-sequences
                       :documentation "A list of sequences (proper order) in pushed in reverse order"))
  (:documentation " Keep track of sequence and matrix info while scanning a PDB file"))

(defun shallow-copy (object)
  (let* ((class (class-of object))
         (copy (allocate-instance class)))
    (dolist (slotd (clos:class-slots class) copy)
      (when (clos:slot-boundp-using-class class object slotd)
        (setf (clos:slot-value-using-class class copy slotd)
              (clos:slot-value-using-class class object slotd))))))


(defun shallow-copy-pdb-scanner (original)
  (let ((dup (shallow-copy original)))
    (let ((new-sequences (loop for sequence in (sequences dup)
                               collect (copy-seq sequence))))
      (setf (sequences dup) new-sequences))
    dup))

(defclass pdb-atom-reader (pdb-reader)
  ((sequences-index :initform 0 :accessor sequences-index)
   (sequence-index :initform 0 :accessor sequence-index)
   (sequences :initarg :sequences :accessor sequences)
   (previous-residue :initarg :previous-residue :accessor previous-residue)
   (previous-topology :initarg :previous-topology :accessor previous-topology)
   (previous-chain-id :initform nil :initarg :previous-chain-id :accessor previous-chain-id)
   (current-chain-id :initform nil :initarg :current-chain-id :accessor current-chain-id)
   (aggregate :initform (chem:make-aggregate) :reader aggregate)
   (molecules :initform nil :initarg :molecules :accessor molecules)
   (finish-molecule :initform nil :accessor finish-molecule)
   (molecule :initform nil :initarg :molecule :accessor molecule)
   (current-residue :initform nil :initarg :current-residue :accessor current-residue)
   (current-topology :initform nil :initarg :current-topology :accessor current-topology)
   (connect-atoms :initarg :connect-atoms :reader connect-atoms)
   (current-residue-name :initform nil :initarg :current-residue-name :accessor current-residue-name)
   (serial-to-atom :initarg :serial-to-atom :reader serial-to-atom)
   (pdb-scanner :initarg :pdb-scanner :accessor pdb-scanner)
   ))

(defun ensure-molecule (reader chain-id)
  (declare (ignore chain-id))
  (unless (molecule reader)
;;;    (format t "Getting molecule at ~a~%" (sequences-index reader))
    (setf (molecule reader) (elt (molecules reader) (sequences-index reader)))
    #+(or)(format *debug-io* "Added molecule ~a to aggregate~%" (molecule reader)))
  (molecule reader))

(defun new-residue-p (residue-name residue-sequence-number i-code chain-id pdb-atom-reader)
  (if (or (not (eq residue-sequence-number (current-residue-number pdb-atom-reader)))
          (not (eq i-code (current-i-code pdb-atom-reader)))
          (not (eq chain-id (previous-chain-id pdb-atom-reader))))
      t
      (if (not (eq residue-name (current-residue-name pdb-atom-reader)))
          (progn
            #+(or)(warn "The residue name changed from ~a to ~a but the residue number ~a didnt change"
                  (current-residue-name pdb-atom-reader)
                  residue-name
                  residue-sequence-number)
            t)
          nil)))

;;; Create a special readtable that treats quote #\' as a normal char
(defvar *quote-readtable* (copy-readtable))
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

(defclass pdb-line () ())

(defclass skip-line (pdb-line) ())

(defclass atom-line (pdb-line)
  ((lineno :initarg :lineno :accessor lineno)
   (atom-serial :initarg :atom-serial :accessor atom-serial)
   (atom-name :initarg :atom-name :accessor atom-name)
   (alt-loc :initarg :alt-loc :accessor alt-log)
   (residue-name :initarg :residue-name :accessor residue-name)
   (chain-id :initarg :chain-id :accessor chain-id)
   (residue-sequence-number :initarg :residue-sequence-number :accessor residue-sequence-number)
   (i-code :initarg :i-code :accessor i-code)
   (x :initarg :x :accessor x)
   (y :initarg :y :accessor y)
   (z :initarg :z :accessor z)
   ))

(defmethod print-object ((obj atom-line) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~4a: ~4a ~4a ~3a ~1a ~4a"
            (lineno obj)
            (atom-serial obj)
            (atom-name obj)
            (residue-name obj)
            (chain-id obj)
            (residue-sequence-number obj))))

(defclass hetatm-line (atom-line) ())

(defclass ter-line (pdb-line) ())

(defclass ssbond-line (pdb-line)
  ((chain-id1 :initarg :chain-id1 :accessor chain-id1)
   (res-seq1 :initarg :res-seq1 :accessor res-seq1)
   (chain-id2 :initarg :chain-id2 :accessor chain-id2)
   (res-seq2 :initarg :res-seq2 :accessor res-seq2)
   ))

(defclass conect-line (pdb-line)
  ((from :initarg :from :accessor from)
   (to :initarg :to :accessor to)))

(defclass mtrix-line (pdb-line)
  ((serial :initarg :serial :accessor serial)
   (x :initarg :x :accessor x)
   (y :initarg :y :accessor y)
   (z :initarg :z :accessor z)
   (tt :initarg :tt :accessor tt)))

(defclass mtrix1-line (mtrix-line) ())
(defclass mtrix2-line (mtrix-line) ())
(defclass mtrix3-line (mtrix-line) ())

(defun parse-line (line pdb-atom-reader read-atom-pos lineno &optional big-z check-big-z)
  "* Arguments
- line :: String
- pdb-atom-reader :: pdb-atom-reader
- read-atom-pos :: boolean
- lineno :: integer
- big-z :: boolean
- check-big-z :: boolean
* Description
Parse the PDB line in line and populate fields in pdb-atom-reader with the data.
If we need to read the atom position pass T for read-atom-pos. If it's a big-z PDB
file (the Z coordinate needs an extra floating point digit) then pass big-z = T.
If you want to check if it's a big-z file then pass check-big-z = T."
  (declare (ignore pdb-atom-reader))
  (labels ((read-one-char (line start)
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
      (let ((line
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
                  (make-instance (if (string= head "ATOM") 'atom-line 'hetatm-line)
                                 :lineno lineno
                                 :atom-serial (parse-integer line :start 6 :end 12) ; atom-serial
                                 :atom-name (let ((*readtable* *quote-readtable*))
                                              (parse-atom-name line))
                                 :alt-loc (read-one-char line 16) ; alt-loc
                                 :residue-name (intern (string-trim '(#\space) (subseq line 17 20)) :keyword) ; residue-name
                                 :chain-id (read-one-char line 21) ; chainid
                                 :residue-sequence-number (parse-integer line :start 22 :end 26) ; res-seq
                                 :i-code (read-one-char line 26) ; chainid/i-code?
                                 :x (when read-atom-pos
                                      (parse-fixed-double (subseq line 30 38)))
                                 :y (when read-atom-pos
                                      (parse-fixed-double (subseq line 38 46)))
                                 :z (when read-atom-pos
                                      (parse-fixed-double (subseq line 46 (if big-z 55 54)))))
                  big-z))
                ((string= head "TER") (make-instance 'ter-line))
                ((string= head "MTRIX" :start1 0 :end1 5)
                 (make-instance (case (elt head 5)
                                  (#\1 'mtrx1-line)
                                  (#\2 'mtrx2-line)
                                  (#\3 'mtrx3-line)
                                  (otherwise (error "Illegal MTRIX opcode ~a" head)))
                                :serial (read-from-string line t nil :start 7 :end 10)
                                :x (read-from-string line t nil :start 10 :end 20)
                                :y (read-from-string line t nil :start 20 :end 30)
                                :z (read-from-string line t nil :start 30 :end 40)
                                :tt (read-from-string line t nil :start 45 :end 55)
                                #|(read-from-string line t nil :start 59 :end 60)|# ))
                ((string= head "CONECT")
                 (make-instance 'conect-line
                                :from (parse-integer line :start 6 :end 11)
                                :to (loop for (start . end) in '((11 . 16) (16 . 21 ) (21 . 26 ) (26 . 31 ) (31 . 36 ) (36 . 41 ) (41 . 46 ) (46 . 51 ) (51 . 56 ) (56 . 61 ))
                                          until (> start (length line))
                                          collect (parse-integer line :start start :end end :junk-allowed t))))
                ((string= head "SSBOND")
                 (make-instance 'ssbond-line
                                :chain-id1 (read-one-char line 15)
                                :res-seq1 (parse-integer (string-trim '(#\space #\tab) (subseq line 17 21)))
                                :chain-id2 (read-one-char line 29)
                                :res-seq2 (parse-integer (string-trim '(#\space #\tab) (subseq line 31 35)))))
                (t nil))))
        line
        ))))


(defun lookup-pdb-residue (pdb-atom-reader sequences-index sequence-index)
  "* Arguments
- pdb-atom-reader :: A pdb-reader
- sequences-index :: Integer - which sequence
- sequence-index  :: Integer - index of pdb-residue in sequenc
   * Return
(values pdb-residue first-residue-of-molecule-p)
* Description
Lookup and return the pdb-residue in the sequence of sequences.
If this residue is the first residue of a molecule return first-residue-of-molecule-p = T
"
  (with-slots (pdb-scanner) pdb-atom-reader
    (let ((sequence (elt (sequences pdb-scanner) sequences-index)))
      (values (elt sequence sequence-index) (= sequence-index 0)))))


#+(or)
(defun pop-sequence-pdb-residue (pdb-atom-reader)
  "* Arguments
- pdb-atom-reader :: A pdb-reader
   * Return
(values pdb-residue first-residue-of-molecule-p)
* Description
Pop and return the next residue in the sequence.
If this residue is the first residue of a molecule return first-residue-of-molecule-p = T
"
  ;; If we have run out of pdb-residues in the current sequence
  ;;  then advance to the next sequence and return an error-val
  (let ((first-residue-of-molecule-p
          (if (< (sequence-index pdb-atom-reader) (length (elt (sequences pdb-atom-reader)
                                                               (sequences-index pdb-atom-reader))))
              nil                ; not the first residue of a molecule
              (progn
                (incf (sequences-index pdb-atom-reader))
                (setf (sequence-index pdb-atom-reader) 0)
                (format t "Advancing to next molecule ~%")
                t))))        ; this is the first residue of a molecule
    (let* ((sequence (elt (sequences pdb-atom-reader) (sequences-index pdb-atom-reader)))
           (pdb-residue (elt sequence (sequence-index pdb-atom-reader))))
      (incf (sequence-index pdb-atom-reader))
      (when (chem:verbose 2) (format t "Popped pdb-residue: ~a~%" pdb-residue))
      (values pdb-residue first-residue-of-molecule-p))))

(defun finish-previous-sequence (pdb-scanner &optional (assign-tail t))
  "* Arguments
- pdb-scanner :: a pdb-scanner
- assign-tail :: bool
* Description
Finish off the current sequence and prepare to start another.
If assign-tail is NIL then don't assign the last read residue as
a tail - this is used when partial sequences are loaded to avoid
treating :main residues as :tail residues when that wont make sense and
create more problems."
  (when (current-reverse-sequence pdb-scanner)
    (when assign-tail
      (let ((tail-residue (car (current-reverse-sequence pdb-scanner))))
        (setf (context tail-residue) :tail)
        (try-to-assign-topology tail-residue pdb-scanner)))
    (let ((seq (nreverse (current-reverse-sequence pdb-scanner))))
      (push seq (reversed-sequences pdb-scanner)))
    (incf (sequences-index pdb-scanner))
    (setf (sequence-index pdb-scanner) 0)
    (setf (previous-residue pdb-scanner) nil)
    (setf (previous-hetatmp pdb-scanner) nil)
    (setf (current-reverse-sequence pdb-scanner) nil)))

(defun try-to-assign-topology (res scanner)
  (let* ((topology (lookup-topology-using-pdb-name-and-context (residue-name res) (context res))))
    (if topology
        (progn
          (assert (not (symbolp topology)))
          (setf (topology res) topology))
        (let ((key (scanner-unknown-residue-key res)))
          #+(or)(format t "Could not immediately identify topology for ~a~%" key)
          (unless (gethash key (unknown-residues scanner))
            #+(or)(format t "Could not identify topology for ~a~%" key)
            (setf (gethash key (unknown-residues scanner)) t))))))

(defgeneric scan (system field scanner line))

(defun apply-filter-to-scanner (line-filter scanner)
  (check-type scanner pdb-scanner)
  (with-slots (chain-ids ignore-residues rename-residues ignore-atoms rename-atoms)
      line-filter
    (loop named outer-loop
          for cur-sequence = (sequences scanner) then (cdr cur-sequence)
          for sequence = (car cur-sequence)
          for selected-chain-id = (and sequence
                                       (> (length sequence) 0)
                                       (or (eq chain-ids :all)
                                           (member (chain-id (first sequence)) chain-ids)))
          if (null selected-chain-id)
            do (setf (car cur-sequence) nil)
          else
            do (loop named inner-loop
                     for inner-cur = sequence then (cdr inner-cur)
                     for pdb-residue = (car inner-cur)
                     for residue-name = (residue-name pdb-residue)
                     for atom-names = (atom-names pdb-residue)
                     for do-ignore-residue = (and ignore-residues
                                                  (gethash residue-name ignore-residues))
                     for do-renamed-residue = (and rename-residues
                                                   (gethash residue-name rename-residues))
                     for do-ignore-atoms = (and ignore-atoms
                                                (gethash residue-name ignore-atoms))
                     for do-rename-atoms = (and rename-atoms
                                                (gethash residue-name rename-atoms))
                     do (if do-ignore-residue
                            (setf (car inner-cur) nil)
                            (when (or do-ignore-atoms do-renamed-residue do-rename-atoms)
                              (let ((new-pdb-residue (shallow-copy pdb-residue)))
                                (when do-renamed-residue
                                  (let ((topology (lookup-topology-using-pdb-name-and-context do-renamed-residue (context pdb-residue))))
                                    (setf (name new-pdb-residue) do-renamed-residue
                                          (residue-name new-pdb-residue) do-renamed-residue
                                          (topology new-pdb-residue) topology)))
                                (let ((atom-names (if do-ignore-atoms
                                                      (loop for an in atom-names
                                                            unless (gethash an do-ignore-atoms)
                                                              collect an)
                                                      atom-names)))
                                  (setf (atom-names new-pdb-residue) atom-names)
                                  (when do-rename-atoms
                                    (let ((new-atoms (loop for atom-name in atom-names
                                                           for new-atom-name = (or (gethash atom-name do-rename-atoms) atom-name)
                                                           collect new-atom-name)))
                                      (setf (atom-names new-pdb-residue) new-atoms)))
                                  (setf (car inner-cur) new-pdb-residue)))))
                     when (null (cdr inner-cur))
                       do (return-from inner-loop nil))
          unless (consp (cdr cur-sequence))
            do (return-from outer-loop nil))
    ;; Filter out (nil)
    (setf (sequences scanner) (loop for sequence in (sequences scanner)
                                    if (and sequence (or (> (length sequence) 1) (first sequence)))
                                      collect sequence
                                    else
                                      collect nil))
    scanner))

(defun apply-filter-to-line (line-filter original)
  (unless line-filter
    (error "There is no line-filter"))
  (typecase original
    (atom-line
     (let ((data original))
       ;; Apply chain-id filter
       (setf data (if (chain-ids line-filter)
                      (cond
                        ((eq (chain-ids line-filter) :all)
                         data)
                        ((consp (chain-ids line-filter))
                         (if (member (chain-id data) (chain-ids line-filter))
                             data
                             (return-from apply-filter-to-line nil)))
                        (t (return-from apply-filter-to-line nil))
                        )))
       ;; data will be defined or we returned already
       (when (gethash (residue-name data) (ignore-residues line-filter))
         (return-from apply-filter-to-line nil))
       (let ((ignore-atoms-ht (gethash (residue-name original) (ignore-atoms line-filter))))
         (when ignore-atoms-ht
           (when (gethash (atom-name data) ignore-atoms-ht)
             (return-from apply-filter-to-line nil))))
       (when (or (rename-residues line-filter)
                 (rename-atoms line-filter))
         (let ((dup (shallow-copy data)))
           (setf data dup))
         (when (rename-residues line-filter)
           (let ((to-name (gethash (residue-name original) (rename-residues line-filter))))
             (when to-name
               (setf (residue-name data) to-name))))
         (when (rename-atoms line-filter)
           (let ((residue-atom-renames (gethash (residue-name original) (rename-atoms line-filter))))
             (when residue-atom-renames
               (let ((to-atom (gethash (atom-name original) residue-atom-renames)))
                 (when to-atom
                   (setf (atom-name data) to-atom)))))))
       data))
    (otherwise original)))

(defun pdb-scanner-read-line (fin pdb-scanner eof-errorp eof lineno)
  "* Arguments
- pdb-scanner : A pdb-scanner
* Description
Scan through the PDB file and build a high level description
of the contents that will be used later to read the PDB file.
The records that are read from the PDB file are:
ATOM - used to build a list of residues and whether they are :first :main :last.
TER  - Used to identify the ends of chains.
MTRIX- Used to build a list of matrices."
  (let* ((line (read-line fin eof-errorp eof))
         (scanner pdb-scanner)
         (check-big-z t))
    (if (eq line eof)
        (if eof-errorp
            (error "Encountered end of file")
            eof)
        (progn
          (multiple-value-bind (raw-line-data big-z)
              (parse-line line pdb-scanner nil lineno nil check-big-z)
            (let ((line-data raw-line-data))
              (when line-data
                #+(or)(format t "Read:  ~a~%" line-data)
                (typecase line-data
                  (atom-line
                   (let ((hetatmp (typep line-data 'hetatm-line)))
                     (setf (previous-hetatmp pdb-scanner) (current-hetatmp pdb-scanner)
                           (current-hetatmp pdb-scanner) hetatmp)
                     (if (and check-big-z (slot-boundp scanner 'big-z))
                         (unless (eq (big-z scanner) big-z)
                           (error "The big-z status of the PDB file changed"))
                         (progn
                           (setf (big-z scanner) big-z
                                 check-big-z nil) ; only check big-z the first ATOM record
                           (when big-z (warn "This is a BIG-Z PDB file (The Z-coordinate contains an extra digit)"))))
                     (with-slots (atom-serial atom-name alt-loc residue-name chain-id residue-sequence-number i-code)
                         line-data
                       (declare (ignore alt-loc))
                       (setf (previous-chain-id pdb-scanner) (current-chain-id pdb-scanner)
                             (current-chain-id pdb-scanner) chain-id)
                       ;; Deal with the current line of data
                       (when (new-residue-p residue-name residue-sequence-number i-code chain-id pdb-scanner)
                         #+(or)(progn
                           (format t "Starting new residue residue-name: ~a name-of-current-residue: ~a info: ~a~%"
                                   residue-name
                                   (current-residue-name pdb-scanner)
                                   (list residue-sequence-number i-code pdb-scanner))
                           (when (null (current-reverse-sequence pdb-scanner))
                             (format t "Starting new molecule~%"))
                           (when (and hetatmp (previous-hetatmp pdb-scanner))
                             (format t "Should be starting a new molecule HETATM~%")))
                         ;; Between HETATM residues there are no bonds - they are in separate molecules
                         (when (and hetatmp (previous-hetatmp pdb-scanner))
                           (finish-previous-sequence pdb-scanner))
                         (setf (current-residue-number pdb-scanner) residue-sequence-number
                               (current-residue-name pdb-scanner) residue-name
                               (current-i-code pdb-scanner) i-code)
                         (let* ((context-guess (if (null (current-reverse-sequence pdb-scanner))
                                                   :head
                                                   :main))
                                (new-residue (make-instance 'pdb-residue
                                                            :hetatmp hetatmp
                                                            :start-lineno lineno
                                                            :sequences-index (sequences-index pdb-scanner)
                                                            :sequence-index (prog1
                                                                                (sequence-index pdb-scanner)
                                                                              (incf (sequence-index pdb-scanner)))
                                                            :chain-id chain-id
                                                            :res-seq residue-sequence-number
                                                            :i-code i-code
                                                            :name residue-name
                                                            :original-residue-name residue-name
                                                            :residue-name residue-name
                                                            :atom-names nil
                                                            :context context-guess
                                                            :atom-serial-first atom-serial)))
                           (setf (current-residue pdb-scanner) new-residue)
                           (unless (gethash residue-name (seen-residues pdb-scanner))
                             ;;(format t "Creating hash-table for ~a~%" residue-name)
                             (setf (gethash residue-name (seen-residues pdb-scanner)) (make-hash-table)))
                           #+(or)(try-to-assign-topology new-residue (scanner pdb))
                           (when (previous-residue pdb-scanner)
                             (setf (atom-serial-last (previous-residue pdb-scanner)) (previous-atom-serial pdb-scanner)))
                           (push new-residue (current-reverse-sequence pdb-scanner))
                           ;; Set things up for the next new residue
                           (setf (previous-residue pdb-scanner) new-residue)))
                       (setf (gethash lineno (lineno-to-sequence-index pdb-scanner))
                             (cons (sequences-index (current-residue pdb-scanner))
                                   (sequence-index (current-residue pdb-scanner))))
                       (let ((ht (gethash residue-name (seen-residues pdb-scanner))))
                         (unless ht (error "Could not find residue-name ~a in the seen-residues: ~a"
                                           residue-name
                                           (alexandria:hash-table-keys (seen-residues pdb-scanner))))
                         (setf (gethash atom-name ht) t))
                       (push atom-name (atom-names (current-residue pdb-scanner)))
                       ;; Now set things up for the next atom record
                       (setf (previous-atom-serial pdb-scanner) atom-serial))))
                  (ter-line
                   (finish-previous-sequence pdb-scanner))
                  (ssbond-line
                   (with-slots (chain-id1 res-seq1 chain-id2 res-seq2)
                       line-data
                     (push (make-instance 'disulphide
                                          :chain-id1 chain-id1
                                          :res-seq1 res-seq1
                                          :chain-id2 chain-id2
                                          :res-seq2 res-seq2) (disulphides scanner))))
                  (conect-line
                   (with-slots (from to)
                       line-data
                     (push (make-instance 'connect
                                          :from from
                                          :to to)
                           (connects scanner))))
                  (mtrix1-line
                   (with-slots (serial x y z tt)
                       line-data
                     (declare (ignore serial))
                     (let ((matrix (geom:make-matrix-identity)))
                       (setf (geom:at matrix 0 0) x
                             (geom:at matrix 0 1) y
                             (geom:at matrix 0 2) z
                             (geom:at matrix 0 3) tt)
                       (push matrix (matrices scanner)))))
                  (mtrix2-line
                   (with-slots (serial x y z tt)
                       line-data
                     (declare (ignore serial))
                     (let ((matrix (car (matrices scanner))))
                       (setf (geom:at matrix 1 0) x
                             (geom:at matrix 1 1) y
                             (geom:at matrix 1 2) z
                             (geom:at matrix 1 3) tt))))
                  (mtrix3-line
                   (with-slots (serial x y z tt)
                       line-data
                     (declare (ignore serial))
                     (let ((matrix (car (matrices scanner))))
                       (setf (geom:at matrix 2 0) x
                             (geom:at matrix 2 1) y
                             (geom:at matrix 2 2) z
                             (geom:at matrix 2 3) tt))))
                  (otherwise nil)))))))))

(defgeneric adjust-residue-name (residue-name pdb-residue scanner system))

(defmethod adjust-residue-name (residue-name pdb-residue scanner system)
  (declare (ignore residue-name pdb-residue scanner system))
  )

(defmethod adjust-residue-name ((residue-name (eql :HIS)) pdb-residue scanner system)
  (declare (ignore scanner system))
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

#+(or)
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
  (let* ((constitution (chem:topology/get-constitution topology))
         (constitution-atoms (chem:constitution/get-constitution-atoms constitution))
         (constitution-atoms-as-list (chem:constitution-atoms-as-list constitution-atoms)))
    (let ((topology-atom-names (loop for ca in constitution-atoms-as-list
                                     collect (chem:atom-name ca))))
      (if (= (length atom-names) (length topology-atom-names))
          (loop for nm in topology-atom-names
                unless (member nm atom-names)
                  do (return-from atom-names-match-topology
                       (values nil topology-atom-names)))
          (return-from atom-names-match-topology
            (values nil topology-atom-names)))
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

(defclass pdb-problem ()
  ((pdb-res :initarg :pdb-res :accessor pdb-res)))

(defclass atoms-mismatch-problem (pdb-problem)
  ((unique-pdb-atom-names :initarg :unique-pdb-atom-names :accessor unique-pdb-atom-names)
   (unique-amber-atom-names :initarg :unique-amber-atom-names :accessor unique-amber-atom-names)
   ))

(defmethod print-object ((object atoms-mismatch-problem) stream)
  (if *print-pretty*
    (let* ((pdb-res (pdb-res object))
           (res-name (residue-name pdb-res))
           (res-seq (res-seq pdb-res))
           (chain-id (chain-id pdb-res))
           ;;(atom-names (atom-names (atom-names pdb-res)))
           (topology (topology pdb-res))
           (unique-pdb (unique-pdb-atom-names object))
           (unique-amber (unique-amber-atom-names object)))
    (declare (ignore res-name res-seq chain-id topology unique-pdb unique-amber))
      #+(or)(format stream "Residue ~a ~a ~a - atom names mismatch~%  atom name PDB: ~a and AMBER: ~a in AMBER residue ~a~%" res-name chain-id res-seq unique-pdb unique-amber (chem:get-name topology) ))
    (print-unreadable-object (object stream :type t))))

(defclass missing-topology (pdb-problem)
  ())

(defgeneric fix-problem (problem))

(defmethod fix-problem ((problem t))
  nil)

(defmethod fix-problem ((problem atoms-mismatch-problem))
  t
  )

(defun validate-scanner (scanner &key ignore-missing-topology)
  (let (problems)
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
                                (let ((problem (make-instance 'atoms-mismatch-problem
                                                              :pdb-res pdb-res
                                                              :unique-pdb-atom-names unique-pdb
                                                              :unique-amber-atom-names unique-amber)))
                                  (push problem problems)))))
                          (when (not ignore-missing-topology)
                            (let ((problem (make-instance 'missing-topology
                                                          :pdb-res pdb-res)))
                              (push problem problems))))))
    problems))

(defun scan-pdb-stream (fin pdb-scanner &key progress system ignore-missing-topology)
  (format t "In scan-pdb-stream~%")
  (let ((bar (if progress
                 (cando:make-progress-bar :message "Scanned" :total (file-length fin) :divisions 100 :on progress)
                 nil)))
    (restart-case
        (progn
          (loop for lineno from 1
                for x = (pdb-scanner-read-line fin pdb-scanner nil :eof lineno)
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
    (let ((sequences (reverse (reversed-sequences pdb-scanner))))
      (setf (sequences pdb-scanner) sequences)
      (finish-scanner pdb-scanner)
      (adjust-residue-names pdb-scanner system)
      (assign-topologys pdb-scanner)
      #+(or)(split-solvent pdb-scanner)
      (let* ((result-scanner (build-sequence pdb-scanner system))
             (problems (validate-scanner result-scanner :ignore-missing-topology ignore-missing-topology))
             (remaining-problems (loop for problem in problems
                                       unless (fix-problem problem)
                                         collect problem)))
        (when remaining-problems
          (format t "There is ~a problem~:P in the PDB file~%" (length remaining-problems)))
        (setf (problems result-scanner) remaining-problems)
        result-scanner)))
  pdb-scanner)


(defun scanPdb (filename &key progress system ignore-missing-topology)
  "* Arguments
- filename : A pathname
* Description
Scan the PDB file and use the ATOM records to build a list of residue sequences and matrices.
* Return
values residue-sequences matrices"
  (with-open-file (fin filename :direction :input)
    (let ((pdb-scanner (make-instance 'pdb-scanner
                                      :filename filename)))
      (scan-pdb-stream fin pdb-scanner
                       :progress progress
                       :system system
                       :ignore-missing-topology ignore-missing-topology)
      )))

(defun scanSelectChainIds (scanner list-or-all)
  (unless (or (eq list-or-all :all)
              (listp list-or-all))
    (error "You must specify :all or a list of chain ids"))
  (setf (chain-ids (line-filter scanner)) list-or-all))

(defun scanIgnoreResidues (scanner list)
  (let ((ht (make-hash-table)))
    (loop for entry in list
          unless (keywordp entry)
            do (leap.core:leap-error "You must provide a residue name - you passed ~a" entry)
          do (setf (gethash entry ht) t))
    (setf (ignore-residues (line-filter scanner)) ht)))

(defun scanRenameResidues (scanner list)
  (let ((ht (make-hash-table)))
    (loop for entry in list
          for from-res = (first entry)
          for to-res = (second entry)
          unless (and (= (length entry) 2)
                      (keywordp from-res)
                      (keywordp to-res))
            do (error "Illegal renameResidues entry: ~a" entry)
          do (setf (gethash from-res ht) to-res))
    (setf (rename-residues (line-filter scanner)) ht)))

(defun scanIgnoreAtoms (scanner list)
  ;; Parse the args into structured input
  (unless (listp list)
    (leap.core:leap-error "Illegal form for ignoreAtoms list ~s" list))
  (let ((ignore-ht (make-hash-table)))
    (loop for entry in list
          for residue-name = (first entry)
          for ignores = (second entry)
          for residue-atom-ht = (if residue-name
                                    (let ((ht (make-hash-table)))
                                      (setf (gethash residue-name ignore-ht) ht)
                                      ht)
                                    nil)
          unless (and (keywordp residue-name) ; residue-name
                      (consp (cdr entry))     ; atom ignore list
                      )
            do (leap.core:leap-error "Illegal form for ignoreAtoms entry ~s" entry)
          when residue-atom-ht
            do (loop for ignore in ignores
                     unless (keywordp ignore)
                       do (leap.core:leap-error "Illegal atom ignore ~s - provide atom names" ignore)
                     do (setf (gethash ignore residue-atom-ht) t)))
    (setf (ignore-atoms (line-filter scanner)) ignore-ht)))

(defun scanRenameAtoms (scanner list)
  ;; Parse the args into structured input
  (unless (listp list)
    (leap.core:leap-error "Illegal form for renameAtoms list ~s" list))
  (let ((rename-ht (make-hash-table)))
    (loop for entry in list
          for residue-name = (first entry)
          for renames = (cdr entry)
          for residue-atom-ht = (if residue-name
                                    (let ((ht (make-hash-table)))
                                      (setf (gethash residue-name rename-ht) ht)
                                      ht)
                                    nil)
          unless (and (keywordp residue-name) ; residue-name
                      (consp (cdr entry))     ; atom rename list
                      )
            do (leap.core:leap-error "Illegal form for renameAtoms entry ~s" entry)
          when residue-atom-ht
            do (loop for rename in renames
                     for atom-from = (first rename)
                     for atom-to = (second rename)
                     unless (= (length rename) 2)
                       do (leap.core:leap-error "Illegal atom rename ~s - provide a list of ( from-atom-name to-atom-name )" rename)
                     do (setf (gethash atom-from residue-atom-ht) atom-to)))
    (setf (rename-atoms (line-filter scanner)) rename-ht)))

(defun scanAmberCheck (pdb-scanner)
  "* Arguments
* Description
Check the scan of the PDB to see if its filters are sufficient to turn it into an AMBER PDB.
* Return
T if the it's an AMBER PDB and NIL if not."
  ;; Gather all the chainids
  (let* ((filter-scanner (make-filter-scanner pdb-scanner))
         (chain-ids-ht (make-hash-table))
         pdb-residues
         (amber-pdb t))
    (loop for sequence in (sequences filter-scanner)
          do (loop for pdb-residue in sequence
                   do (when pdb-residue
                        (setf (gethash (chain-id pdb-residue) chain-ids-ht) t)
                        (when (or (eq (chain-ids (line-filter pdb-scanner)) :all)
                                  (member (chain-id pdb-residue) (chain-ids (line-filter pdb-scanner))))
                          (push pdb-residue pdb-residues)))))
    (format t "Chain-ids: ~s~%" (chain-ids (line-filter pdb-scanner)))
    ;; Look for unrecognized residue names
    (let (unknown-residue-names)
      (loop for pres in pdb-residues
            for top = (lookup-topology-using-pdb-name-and-context (residue-name pres) (context pres))
            unless top
              do (progn
                   (setf amber-pdb nil)
                   (pushnew (cons (residue-name pres) (context pres)) unknown-residue-names :test 'equal)))
      (loop for residue-name-context in unknown-residue-names
            for residue-name = (car residue-name-context)
            for context = (cdr residue-name-context)
            do (format t "Unknown residue/context ~s/~s~%" residue-name context)))
    (when amber-pdb
      (let ((unknown-atom-names nil))
        (loop for pres in pdb-residues
              for original-residue-name = (original-residue-name pres)
              for residue-name = (residue-name pres)
              for context = (context pres)
              for top = (lookup-topology-using-pdb-name-and-context residue-name context)
              for atom-names = (atom-names pres)
              do (multiple-value-bind (match topology-atom-names)
                     (atom-names-match-topology atom-names top)
                   (unless match
                     (let ((extra-atom-names (loop for name in atom-names
                                                   unless (member name topology-atom-names)
                                                     collect name)))
                       (when extra-atom-names
                         (setf amber-pdb nil)
                         (pushnew (cons (cons original-residue-name context) extra-atom-names) unknown-atom-names :test 'equal))))))
        (loop for unk in unknown-atom-names
              for residue-name = (car (car unk))
              for context = (cdr (car unk))
              for names = (cdr unk)
              do (format t "Unknown atom names: ~s/~s -> ~s~%" residue-name context names))))
    (if amber-pdb
        (format t "This is an AMBER compatible PDB file after filters are applied.~%")
        (format t "This is NOT yet an AMBER compatible PDB file - keep adding filters.~%"))
    amber-pdb))

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
          (format t "Scanning disulphide: ~s~%" dis)
          (let ((res1 (find-pdb-residue scanner (res-seq1 dis)
                                        (chain-id1 dis)))
                (res2 (find-pdb-residue scanner (res-seq2 dis)
                                        (chain-id2 dis))))
            (format t "res1 = ~s~%" res1)
            (format t "res2 = ~s~%" res2)
            (let* ((top1 (lookup-topology-using-pdb-name-and-context :CYX (context res1)))
                   (top2 (lookup-topology-using-pdb-name-and-context :CYX (context res2))))
              (format t "top1 = ~s~%" top1)
              (format t "top2 = ~s~%" top2)
              ;; If CYX topologies are available - use those
              (when top1 (setf (topology res1) top1
                               (residue-name res1) :CYX))
              (when top2 (setf (topology res2) top2
                               (residue-name res2) :CYX))
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

(defun read-and-process-line (fin reader eof-errorp eof big-z lineno)
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
  (let ((line (read-line fin eof-errorp eof))
        pdb-residue
        sequences-index
        sequence-index
        first-residue-of-molecule-p)
    (if (eq line eof)
        (if eof-errorp
            (error "read-and-process-line encountered eof")
            eof)
        (let* ((raw-line-data (parse-line line reader t lineno big-z))
               (line-data (apply-filter-to-line (line-filter (pdb-scanner reader)) raw-line-data)))
          (when line-data
            ;;          (format t "Reading line: ~a~%" line-data)
            (when (chem:verbose 2) (format *debug-io* "pdb.lisp read-and-process-line~% pdb-residue -> ~a~% line-data -> ~s~%" pdb-residue line-data))
            (typecase line-data
              (atom-line
               (with-slots (atom-line-lineno atom-serial atom-name alt-loc residue-name chain-id residue-sequence-number i-code x y z)
                   line-data
                 (declare (ignore alt-loc))
                 (setf (previous-chain-id reader) (current-chain-id reader)
                       (current-chain-id reader) chain-id)
                 (when (chem:verbose 2) (format t "At top of atom residue-name: ~a residue-sequence-number: ~a i-code: ~a reader: ~a~%"
                                                residue-name residue-sequence-number i-code reader))
                 (when (new-residue-p residue-name residue-sequence-number i-code chain-id reader)
                   (setf (current-residue-number reader) residue-sequence-number
                         (current-i-code reader) i-code)
                   (let ((sequence-indices (gethash lineno (lineno-to-sequence-index (pdb-scanner reader)))))
                     (setf sequences-index (car sequence-indices)
                           sequence-index (cdr sequence-indices)))
                   (multiple-value-setq (pdb-residue first-residue-of-molecule-p)
                     (lookup-pdb-residue reader sequences-index sequence-index))
                   (when first-residue-of-molecule-p
                     (setf (molecule reader) nil)
                     (when (chem:verbose 2) (format t "Starting a new molecule~%"))
                     (setf (current-residue reader) nil)
                     (setf (current-topology reader) nil)
                     (when (chem:verbose 2) (format t "first residue of molecule: ~a~%" pdb-residue))
                     (ensure-molecule reader chain-id))
                   #+(or)(break "Check the reader: ~a~%" reader)
                   (setf (previous-topology reader) (current-topology reader)
                         (previous-residue reader) (current-residue reader)
                         (current-residue-name reader) residue-name
                         (current-topology reader) (topology pdb-residue))
                   (unless (eq (name pdb-residue) residue-name)
                     (error "There is a mismatch between~%   ~a and the expected residue:~% ~a for ~a"
                            pdb-residue residue-name line-data))
                   (let ((prev-top (previous-topology reader))
                         (cur-top (current-topology reader)))
                     (if cur-top
                         ;; There is a topology - use it
                         (let ((cur-res (chem:build-residue-single-name cur-top)))
                           (chem:set-pdb-name cur-res (residue-name pdb-residue))
                           #+(or)(format *debug-io* "built residue with name ~a using topology ~a~%" cur-res cur-top)
                           (chem:set-id cur-res (calculate-residue-sequence-number residue-sequence-number i-code))
                           (setf (current-residue reader) cur-res)
                           (let ((prev-res (previous-residue reader)))
                             (chem:add-matter (ensure-molecule reader chain-id) cur-res)
                             (when prev-res
                               (unless prev-top
                                 (error "Could not find topology for ~a" prev-res))
                               (unless cur-top
                                 (error "Could not find topology for ~a cur-top" cur-res))
                               (chem:connect-residues prev-top
                                                      prev-res
                                                      :+default
                                                      cur-top
                                                      cur-res
                                                      :-default))))
                         ;; There is no topology, create an empty residue
                         (let ((cur-res (chem:make-residue residue-name)))
                           (chem:set-pdb-name cur-res (residue-name pdb-residue))
                           (setf (current-residue reader) cur-res)
                           (chem:add-matter (ensure-molecule reader chain-id) cur-res) cur-res))))
                 (let ((atom (or (chem:content-with-name-or-nil (current-residue reader) atom-name)
                                 (chem:content-with-name-or-nil (current-residue reader) (gethash atom-name *pdb-atom-map*))
                                 )))
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
                                  (list atom-serial atom-name residue-name chain-id residue-sequence-number i-code)
                                  (current-topology reader)))
                     (t                ; Add a new atom to the residue
                      (let ((atom (chem:make-atom atom-name (chem:element-from-atom-name-string (string atom-name)))))
                        (chem:set-id atom atom-serial)
                        (chem:add-matter (current-residue reader) atom)
                        (chem:set-position atom (geom:vec x y z))
                        (chem:setf-needs-build atom nil)
                        (when (serial-to-atom reader)
                          (setf (gethash atom-serial (serial-to-atom reader)) atom))))))))
              (ter-line nil)
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
  (declare (ignore system))
  scanner)

(defvar *serial-to-atoms* nil)

(defun molecules-from-sequences (sequences)
  (loop for sequence in sequences
        for index from 0
        for make-molecule = (or (and (= (length sequence) 1) (first sequence))
                                (> (length sequence) 1))
        for sequence-name = (when make-molecule
                              (cond
                                ((and (= (length sequence) 1) (first sequence))
                                 (intern (format nil "~a_~d" (string (name (first sequence))) index) :keyword))
                                (t (intern (concatenate 'string "mol" (format nil "_~d" index)) :keyword))))
        when make-molecule
          collect (let ((mol (chem:make-molecule sequence-name)))
                    #+(or)(format *debug-io* "Creating molecule ~a for sequence ~a~%" mol sequence)
                    mol)))

(defun make-filter-scanner (pdb-scanner)
  (let ((dup (shallow-copy-pdb-scanner pdb-scanner))) ; Use a copy of the scanner
    (when (line-filter pdb-scanner)
      (apply-filter-to-scanner (line-filter pdb-scanner) dup))
    dup))

(defun load-pdb-stream (fin &key filename pdb-scanner progress system ignore-missing-topology)
  (let* ((pdb-scanner (if pdb-scanner
                          (make-filter-scanner pdb-scanner)
                          (let* ((new-scanner (build-sequence (let ((pdb-scanner (make-instance 'pdb-scanner
                                                                                                :filename filename)))
                                                                (scan-pdb-stream fin pdb-scanner
                                                                                 :progress progress
                                                                                 :ignore-missing-topology ignore-missing-topology))
                                                              system)))
                            (file-position fin 0)
                            new-scanner)))
         (serial-to-atoms (connect-atoms-hash-table pdb-scanner)))
    (setq *serial-to-atoms* serial-to-atoms)
    (let ((pdb-atom-reader (make-instance 'pdb-atom-reader :sequences (mapcar (lambda (seq)
                                                                                (copy-list seq))
                                                                              (sequences pdb-scanner))
                                                           :pdb-scanner pdb-scanner
                                                           :molecules (molecules-from-sequences (sequences pdb-scanner))
                                                           :connect-atoms (connects pdb-scanner)
                                                           :serial-to-atom serial-to-atoms)))
      (let ((bar (if progress
                     (cando:make-progress-bar :message "Load pdb" :total (file-length fin) :divisions 100 :on progress)
                     nil)))
        (restart-case
            (progn
              (setf (molecule pdb-atom-reader) nil)
              (loop for lineno from 1
                    for x = (read-and-process-line fin pdb-atom-reader nil :eof (big-z pdb-scanner) lineno)
<<<<<<< HEAD
                    do (when bar
=======
                    do (when (and (= 0 (floor lineno 100))
                                  bar)
>>>>>>> f7dfa02c (Improve leap's alignAxes)
                         (cando:progress-advance bar (file-position fin)))
                    until (eq x :eof))
              (when progress (format t "Loaded pdb~%")))
          (ran-out-of-sequence ()
            :report "Ran out of sequence while filling residues"
            (format t "Continuing with partial sequence~%")))
        (when bar (cando:progress-done bar)))
      (loop for connect in (connects pdb-scanner)
            for from-atom = (gethash (from connect) serial-to-atoms)
            for to-atoms = (mapcar (lambda (c) (gethash c serial-to-atoms)) (to connect))
            do (when (and from-atom (typep from-atom 'chem:atom))
                 (mapc (lambda (to-atom)
                         (when (and to-atom
                                    (typep to-atom 'chem:atom)
                                    (not (chem:is-bonded-to from-atom to-atom)))
                           (when (>= (chem:number-of-bonds from-atom) 4)
                             (error "CONECT ~a between ~a and ~a will fail" connect from-atom to-atom))
                           (when (>= (chem:number-of-bonds to-atom) 4)
                             (error "CONECT ~a will lead to too many bonds for to atom ~a" connect to-atom))
                           (format t "CONECT ~a to ~a~%" from-atom to-atom)
                           (chem:bond-to from-atom to-atom :single-bond)))
                       to-atoms)))
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
          (cando:simple-build-unbuilt-atoms aggregate)
          (when progress
            (format t "Built ~d heavy atom~:P~%" unbuilt-heavy-atoms)))
        (let ((built (cando:build-unbuilt-hydrogens aggregate)))
          (when progress
            (format t "Built ~d missing hydrogen~:P~%" built)))
;;;            (cando:maybe-join-molecules-in-aggregate aggregate)
;;;            (cando:maybe-split-molecules-in-aggregate aggregate)
        (setf aggregate (classify-molecules aggregate system))
        (let ((name-only (if filename
                             (pathname-name (pathname filename))
                             "unknown")))
          (chem:set-name aggregate (intern name-only *package*)))
        (values aggregate pdb-scanner pdb-atom-reader)))))


(defun loadPdb (filename-or-scanner &key (progress t) system ignore-missing-topology)
  "    variable = loadPdb filename
      STRING|SCANNER                       _filename_

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
  (let ((source (leap.core:lookup-variable filename-or-scanner)))
    (format t "source: ~s~%" source)
    (cond
      ((stringp source)
       (let ((filename (leap.core:ensure-path source)))
         (with-open-file (fin filename :direction :input)
           (load-pdb-stream fin :filename filename
                                :progress progress
                                :system system
                                :ignore-missing-topology ignore-missing-topology))))
      ((typep source 'pdb-scanner)
       (let* ((pdb-scanner source)
              (filename (filename pdb-scanner)))
         (with-open-file (fin filename :direction :input)
           (load-pdb-stream fin :filename filename
                                :pdb-scanner pdb-scanner
                                :progress progress
                                :system system
                                :ignore-missing-topology ignore-missing-topology))))
      (t (format t "You must provide a filename or a scanner~%")))))

(defgeneric classify-molecules (aggregate system))

(defmethod classify-molecules (aggregate system)
  (declare (ignore system))
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
