;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :sdf)

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
         label/s
         number/s number/s number/s number/s number/s number
         )
  (:destructure (x s1 y s2 z s3 label n0 n1 n2 n3 n4 n5)
                (declare (ignore s0 s1 s2 s3))
                (list x y z label n0 n1 n2 n3 n4 n5)))

(esrap:defrule bond-line
    (and number/s number/s number/s number/s number/s number)
  (:lambda (x)
    x))

(esrap:defrule header-line
    (and number/s number/s number/s number/s number/s number/s "999 V2000"))

#|
(esrap:parse 'header-line "8  8  0  0  0  0  0  0  0  0999 V2000")
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

(defun parse-groups (sin eof-error-p eof)
  (flet ((terminating-read-line (sin)
           (let ((line (read-line sin eof-error-p eof)))
             (when (eq line eof)
               (return-from parse-groups eof))
             line)))
    (let* ((name (terminating-read-line sin))
           (blank1 (terminating-read-line sin))
           (blank2 (terminating-read-line sin))
           (header-line (string-trim '(#\space) (terminating-read-line sin)))
           (header (esrap:parse 'header-line header-line))
           (num-atoms (first header))
           (num-bonds (second header))
           (atoms (loop for atom-index below num-atoms
                        for line = (string-trim '(#\space) (terminating-read-line sin))
                        collect (esrap:parse 'atom-line line)))
           (bonds (loop for bond-index below num-bonds
                        for line = (string-trim '(#\space) (terminating-read-line sin))
                        collect (esrap:parse 'bond-line line))) 
           (atom-vector (build-atom-vector atoms bonds)))
      (let ((atom-to-index (make-hash-table))
            groups)
        (loop for index below (length atom-vector)
              for atom = (aref atom-vector index)
              do (setf (gethash atom atom-to-index) index))
        (loop for index below (length atom-vector)
              for atom = (aref atom-vector index)
              when atom
                do (let* ((spanning-loop (chem:make-spanning-loop atom))
                          (group (loop for span-atom in (chem:all-atoms spanning-loop)
                                       for index = (gethash span-atom atom-to-index)
                                       do (setf (aref atom-vector index) nil)
                                       collect span-atom)))
                     (push group groups)))
        (values groups name)))))

(defun parse-mdl-molecule (fin eof-error-p eof)
  (multiple-value-bind (atom-groups name)
      (parse-groups fin eof-error-p eof)
    (cond
      ((eq atom-groups eof)
       (if eof-error-p
           (error "End of file hit in parse-mdl-molecule")
           (return-from parse-mdl-molecule eof)))
      ((= (length atom-groups) 1)
       (let* ((atoms (first atom-groups))
              (name-symbol (intern name :keyword))
              (molecule (chem:make-molecule name-symbol))
              (residue (chem:make-residue name-symbol)))
         (chem:add-matter molecule residue)
         (loop for group in atom-groups
               do (loop for atom in group
                        do (chem:add-matter residue atom)))
         molecule))
      (t (error "Deal with multiple molecules in sdf")))))

(defun parse-sdf-section (fin eof-error-p eof)
  (let ((molecule (parse-mdl-molecule fin eof-error-p eof)))
    (format t "Read molecule: ~a~%" molecule)
    (if (eq molecule eof)
        eof
        (progn
          (loop for line = (read-line fin eof-error-p eof)
                until (string= line "$$$$"))
          molecule))))

(defun parse-sdf-file (fin)
  (let ((molecules nil))
    (loop for molecule = (parse-sdf-section fin nil :eof)
          if (eq molecule :eof)
            do (return-from parse-sdf-file molecules)
          else
            do (push molecule molecules))))

(defun load-sdf-as-list-of-molecules (fname)
  (with-open-file (fin fname :direction :input)
    (parse-sdf-file fin)))

(defun load-sdf-as-aggregate (fname)
  (let ((molecules (load-sdf-as-list-of-molecules fname))
        (aggregate (chem:make-aggregate)))
    (loop for mol in molecules
          do (chem:add-matter aggregate mol))
    aggregate))
    
