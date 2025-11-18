;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :sdf)

(define-condition sdf-parse-error (error)
  ((message
     :accessor sdf-parse-error-message
     :initarg :message
     :initform nil)
   (code
     :accessor sdf-parse-error-code
     :initarg :code
     :initform nil))
  (:report
    (lambda (condition stream)
      (format stream "~a~@[ encountered while parsing \"~a\"~].~%"
              (sdf-parse-error-message condition)
              (sdf-parse-error-code condition)))))

(esrap:defrule skippable
    (+ (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(esrap:defrule skippable?
    (* (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(parser.common-rules:defrule/s number
    (+ (esrap:character-ranges (#\0 #\9)))
  (:lambda (x)
    (parse-integer (esrap:text x))))

(parser.common-rules:defrule/s label
    (+ (esrap:character-ranges (#\A #\Z) (#\a #\z) #\?))
  (:text t))

(defparameter *elements*
  '("He" "Li" "Be" "Ne" "Na" "Mg" "Al" "Si" "Cl" "Ar" "Ca" "Sc" "Ti" "Cr" "Mn" "Fe" "Co"
    "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh"
    "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "Xe" "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu"
    "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Hf" "Ta" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl"
    "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac" "Th" "Pa" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es"
    "Fm" "Md" "No" "Lr"
    "H" "W" "U" "B" "C" "N" "O" "F" "P" "S" "K" "V" "I" "Y"))

(esrap:defrule atom-line
    (and parser.common-rules:float-literal
         parser.common-rules:whitespace/not-newline
         parser.common-rules:float-literal
         parser.common-rules:whitespace/not-newline
         parser.common-rules:float-literal
         parser.common-rules:whitespace/not-newline
         label
         (+ (and parser.common-rules:whitespace/not-newline parser.common-rules:integer-literal/decimal))
         parser.common-rules:whitespace*
         )
  (:destructure (x s1 y s2 z s3 label nrest space)
                (declare (ignore s0 s1 s2 s3 space))
                (list x y z label (mapcar #'second nrest))))

(esrap:defrule bond-line
    (and digits3 digits3 digits3 digits3 (* digits3) parser.common-rules:whitespace*)
  (:lambda (x)
    x))

(esrap:defrule m-chg-line
  (and "M" parser.common-rules:whitespace/not-newline
       "CHG" parser.common-rules:whitespace/not-newline
       number parser.common-rules:whitespace/not-newline (* (and number/s number/?s)))
  (:destructure (m s1 chg s2 num s3 pairs)
                (declare (ignore m s1 chg s2 s3))
                (list :chg num pairs)))

(esrap:defrule m-end-line
  (and "M" parser.common-rules:whitespace/not-newline
       "END")
  (:destructure (m s1 end)
                (declare (ignore m s1 end))
                (list :end)))

(esrap:defrule m-line
    (or m-chg-line m-end-line)
  (:lambda (data)
    data))

(esrap:defrule digit-or-space
    (or #\space #\0 #\1 #\2 #\3 #\4 #\5 #\6 #\7 #\8 #\9)
  (:lambda (data)
    (if (string= data " ")
        0
        (parse-integer data :junk-allowed t))))

(esrap:defrule digits3
    (and digit-or-space digit-or-space digit-or-space)
  (:destructure (d0 d1 d2)
                (+ (* 100 d0) (* 10 d1) d2)))

(esrap:defrule count-line
    (and digits3 digits3 digits3 digits3 digits3 digits3 (* #\space) (+ number/s) "V2000"))

#|
(esrap:parse 'count-line "8  8  0  0  0  0  0  0  0  0999 V2000")
(esrap:parse 'atom-line "  0.0000    4.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0")
|#



(defun build-atom-vector (atoms bonds)
  (let ((atom-vector (make-array (length atoms)))
        (map-numbers (make-hash-table)))
    (loop for atom-info in atoms
          for index from 0
          do (let* ((xpos (elt atom-info 0))
                    (ypos (elt atom-info 1))
                    (zpos (elt atom-info 2))
                    (label-string (elt atom-info 3))
                    (label-string (string label-string))
                    (match-element-string (find label-string *elements* :test #'string=))
                    (element-string (if match-element-string match-element-string "C"))
                    (element (intern element-string :keyword))
                    (label (if match-element-string nil (intern label-string :keyword)))
                    (atom-name (intern (format nil "~a~a" element-string index) :keyword))
                    (atom (chem:make-atom atom-name element)))
               (chem:set-position atom (geom:vec xpos ypos zpos))
               (when label (chem:set-property atom :label label))
               (setf (elt atom-vector index) atom)))
    (loop for bond-info in bonds
          for i1 = (1- (elt bond-info 0))
          for i2 = (1- (elt bond-info 1))
          for order-index = (elt bond-info 2)
          for order = (elt #(:no-bond :single-bond :double-bond :triple-bond) order-index)
          for atom1 = (elt atom-vector i1)
          for atom2 = (elt atom-vector i2)
          do (chem:bond-to atom1 atom2 order))
    atom-vector))

(defstruct m-line-state charges)

(defun parse-m-line (line m-line-state)
  (let ((m-line (esrap:parse 'm-line line)))
    (cond
      ((eq (car m-line) :chg)
       (loop for pair in (caddr m-line)
             for atom-index = (car pair)
             for charge = (cadr pair)
             do (push (cons atom-index charge) (m-line-state-charges m-line-state)))
       :chg)
      ((eq (car m-line) :end)
       :end)
      (t (error 'sdf-parse-error :message "Unknown m-line" :code line)))))

(defun parse-groups (sin eof-error-p eof)
  (flet ((terminating-read-line (sin)
           (let ((line (read-line sin eof-error-p eof)))
             (when (eq line eof)
               (return-from parse-groups eof))
             line)))
    (let* ((m-line-state (make-m-line-state))
           (name (terminating-read-line sin))
           (blank1 (terminating-read-line sin))
           (blank2 (terminating-read-line sin))
           (count-line (terminating-read-line sin))
           (header (esrap:parse 'count-line count-line))
           (num-atoms (first header))
           (num-bonds (second header))
           (atoms (loop for atom-index below num-atoms
                        for line = (string-trim '(#\space) (terminating-read-line sin))
                        collect (esrap:parse 'atom-line line)))
           (bonds (loop for bond-index below num-bonds
                        for line = (terminating-read-line sin)
                        collect (esrap:parse 'bond-line line)))
           (atom-vector (build-atom-vector atoms bonds)))
      ;; parse the M lines
      (loop for line = (string-trim '(#\space) (terminating-read-line sin))
            for result = (parse-m-line line m-line-state)
            until (eq result :end))
      (loop for index-charge in (m-line-state-charges m-line-state)
            for index = (1- (car index-charge))
            for charge = (cdr index-charge)
            do (chem:set-charge (elt atom-vector index) charge))
      (let ((atom-to-index (make-hash-table)))
        (loop for index below (length atom-vector)
              for atom = (aref atom-vector index)
              do (setf (gethash atom atom-to-index) index))
        (let* ((group (coerce atom-vector 'list))
               (groups (list group)))
          (values groups name))))))

(defun parse-mdl-molecule (fin eof-error-p eof)
  (multiple-value-bind (atom-groups name)
      (parse-groups fin eof-error-p eof)
    (cond
      ((eq atom-groups eof)
       (if eof-error-p
           (error 'end-of-file :stream fin)
           (return-from parse-mdl-molecule eof)))
      ((= (length atom-groups) 1)
       (let* ((atoms (first atom-groups))
              (molecule (chem:make-molecule nil))
              (residue (chem:make-residue nil)))
         (chem:add-matter molecule residue)
         (loop for group in atom-groups
               do (loop for atom in group
                        do (chem:add-matter residue atom)))
         (values molecule name)))
      (t (error 'sdf-parse-error :message "Deal with multiple molecules in sdf")))))

(defun parse-data-item (line fin eof-error-p eof)
  (let* ((data-name-start (position #\< line))
         (data-name-end (position #\> line :start data-name-start))
         (data-name (subseq line (1+ data-name-start) data-name-end))
         (data-val (read fin eof-error-p eof)))
    (values (intern data-name :keyword) data-val)))

(defun parse-sdf-section (fin eof-error-p eof)
  (multiple-value-bind (molecule name)
      (parse-mdl-molecule fin eof-error-p eof)
    (if (eq molecule eof)
        eof
        (progn
          (loop for line = (read-line fin eof-error-p eof)
                when (and (> (length line) 0) (char= (elt line 0) #\>))
                  do (multiple-value-bind (data-name data-val)
                         (parse-data-item line fin eof-error-p eof)
                       (chem:matter/set-property molecule data-name data-val))
                until (string= line "$$$$"))
          (values molecule name)))))

(defun ensure-names-unique (names-molecules)
  (loop
    (let ((sorted-by-name (sort names-molecules #'string< :key #'car))
          (collision nil)
          (counter 1))
      (loop with (prev-name . prev-molecule ) = (first sorted-by-name)
            for cur = (rest sorted-by-name) then (cdr cur)
            for (name . molecule) = (car cur)
            while cur
            do (when (string= name prev-name)
                 (let ((new-name (format nil "~a_~a" name (incf counter))))
                   #+(or)(format t "prev-name: ~a  name: ~a  new-name: ~a~%" prev-name name new-name)
                   (rplaca (car cur) new-name)))
               (setf prev-name name
                     prev-molecule molecule))
      (unless collision (return-from ensure-names-unique sorted-by-name)))))

(defun parse-sdf-file (fin)
  (let ((names-molecules
          (let (molecules names)
            (loop named load-molecules
                  do (multiple-value-bind (molecule name)
                         (parse-sdf-section fin nil :eof)
                       (if (eq molecule :eof)
                           (return-from load-molecules molecules)
                           (push (cons name molecule) molecules)))))))
    (let ((sorted-by-name (ensure-names-unique names-molecules)))
      (loop for (name . molecule) in sorted-by-name
            for residue = (chem:content-at molecule 0)
            for name-sym = (intern name :keyword)
            do (chem:set-name molecule name-sym)
               (chem:set-name residue name-sym))
      (mapcar #'cdr sorted-by-name))))

(defun load-sdf-as-list-of-molecules (fname)
  (with-open-file (fin fname :direction :input)
    (parse-sdf-file fin)))

(defun load-sdf-as-aggregate (fname)
  (let ((molecules (load-sdf-as-list-of-molecules fname))
        (aggregate (chem:make-aggregate)))
    (loop for mol in molecules
          do (chem:add-matter aggregate mol))
    aggregate))

(defun load (fname)
  "Load all of the molecules in the SDF file as a single aggregate"
  (load-sdf-as-aggregate fname))
  
