;;; file: smiles.lisp
;;;
;;; Modified by Christian Schafmeister (2022)
;;;
;;; Copyright (c) 2008-2014 Cyrus Harmon (ch-lisp@bobobeach.com)
;;; All rights reserved.
;;;
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;;
;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.
;;;
;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.
;;;
;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


(in-package :smiles)

(defconstant +top-left+ 0)
(defconstant +bottom-left+ 1)
(defconstant +top-right+ 2)
(defconstant +bottom-right+ 3)

(defclass atom ()
  ()
  )

(defclass bond ()
  ())

(defclass atom-command ()
  ((atm :initarg :atm :accessor atm)
   (command :initarg :command :accessor command)))

;;; We need a way of explicitly storing configurations around double
;;; bonds and configurations at chiral centers. Hrm...
(defclass molecule ()
  ((name :initarg :name :accessor name :initform nil)
   (atom-name-hash :accessor atom-name-hash
                   :initform (make-hash-table :test 'equal))
   (chem-molecule :accessor chem-molecule :initarg :chem-molecule)
   (chem-residue :accessor chem-residue :initarg :chem-residue)
   (commands :accessor commands :initform nil)
   (chiral-arrangements :accessor chiral-arrangements :initform nil)
   (double-bond-arrangements :accessor double-bond-arrangements :initform nil))
  (:documentation "A class for representing molecules, that is a graph
  of atoms connected by bonds. It is usually the case that a molecule
  will be connected, but it is possible that certain
  operations (adding atoms, removing bonds, e.g.) will leave the
  molecule in a state where it is actually multiple disconnected
  components."))

(defun make-molecule (name)
  (let* ((cmol (chem:make-molecule name))
         (cres (chem:make-residue name))
         (mm (make-instance 'molecule
                            :chem-molecule cmol
                            :chem-residue cres)))
    (chem:add-matter cmol cres)
    mm))

(defun get-element (element)
  (intern element :keyword)
  )

(defparameter *aromatic-elements* '( (:|c| . :C)
                                    (:|n| . :N)
                                    (:|o| . :O)
                                    (:|p| . :P)
                                    (:|as| . :|As|)
                                    (:|se| . :|Se|)))

(defun add-atom (mol name orig-element)
  (let* ((element (or (cdr (assoc orig-element *aromatic-elements*))
                      orig-element))
         (atm (chem:make-atom name element)))
    (chem:add-matter (chem-residue mol) atm)
    atm))

(defclass spatial-arrangement () ())

(defclass vinyl-bond (spatial-arrangement)
  ((left-atom :initarg :left-atom :accessor left-atom)
   (right-atom :initarg :right-atom :accessor right-atom)
   (bond-direction :initarg :bond-direction :accessor bond-direction)))

(defmethod print-object ((geom vinyl-bond) stream)
  (print-unreadable-object (geom stream :type t)
    (format stream ":left-atom ~a :right-atom ~a :bond-direction ~a"
            (left-atom geom)
            (right-atom geom)
            (bond-direction geom))))

(defclass double-bond (spatial-arrangement)
  ((top-left-atom :accessor top-left-atom :initform nil)
   (bottom-left-atom :accessor bottom-left-atom :initform nil)
   (top-right-atom :accessor top-right-atom :initform nil)
   (bottom-right-atom :accessor bottom-right-atom :initform nil)
   (bond :initarg :bond :accessor bond)))

(defmethod print-object ((geom double-bond) stream)
  (print-unreadable-object (geom stream :type t)
    (format stream "~a/~a ~a=~a ~a/~a"
            (and (top-left-atom geom) (chem:get-name (top-left-atom geom)))
            (and (bottom-left-atom geom) (chem:get-name (bottom-left-atom geom)))
            (chem:get-name (chem:get-atom1 (bond geom)))
            (chem:get-name (chem:get-atom2 (bond geom)))
            (and (top-right-atom geom) (chem:get-name (top-right-atom geom)))
            (and (bottom-right-atom geom) (chem:get-name (bottom-right-atom geom))))))

(defclass tetrahedral-center (spatial-arrangement)
  ((chiral-atom :initarg :chiral-atom :accessor chiral-atom)
   (neighbors :accessor neighbors
              :initform (make-array 4 :initial-element nil :fill-pointer 0))
   ;; possible orientations are :clockwise and :anticlockwise
   (orientation :initarg :orientation :accessor orientation)))

(defun add-hydrogens (molecule &optional exclude-list)
  "Adds hydrogens to atoms in a molecule such that each atom has the
lowest normal valence consistent with the number of pre-existing bonds
for that atom. Returns two values, molecule and a list of the newly
added hydrogen ATOMs."
  (cando:do-atoms (atom (chem-molecule molecule))
    (chem:fill-in-implicit-hydrogens-with-residue atom (chem-residue molecule))))

(defparameter *bond-orders*
  '((:single-bond . 1)
    (:up-bond . 1.1)
    (:down-bond . 1.2)
    (:aromatic-bond . 1.5)
    (:double-bond . 2)
    (:triple-bond . 3)))

(defun get-bond-order-number (keyword-or-number)
  (etypecase keyword-or-number
    (keyword (cdr (assoc keyword-or-number *bond-orders*)))
    (number (cdr (rassoc keyword-or-number *bond-orders*)))))

(defun get-bond-order-keyword (keyword-or-number)
  (etypecase keyword-or-number
    (keyword (car (assoc keyword-or-number *bond-orders*)))
    (number (car (rassoc keyword-or-number *bond-orders*)))))

(defgeneric add-bond (molecule atom-identifier-1 atom-identifier-2
                      &key type order direction))

(defgeneric find-bond (molecule atom-identifier-1 atom-identifier-2))

(defgeneric remove-bond (molecule atom-or-bond &rest args))

(defgeneric find-bonds-containing (molecule atom))

(defmethod add-bond ((molecule molecule) atom-1 atom-2
                     &key (type :single-bond) (order 1) direction)
  (when direction
    (let ((bond (make-instance 'vinyl-bond
                               :left-atom atom-1
                               :right-atom atom-2
                               :bond-direction direction)))
      (push bond (double-bond-arrangements molecule))))
  (chem:bond-to atom-1 atom-2 type))

(defmethod remove-bond ((molecule molecule) (atom-identifier-1 atom) &rest args)
  (destructuring-bind (atom-identifier-2)
      args
    (error "(graph:remove-edge-between-nodes molecule atom-identifier-1 atom-identifier-2)")))

(defmethod remove-bond ((molecule molecule) (atom-identifier-1 string) &rest args)
  (destructuring-bind (atom-identifier-2)
      args
    (error "(graph:remove-edge-between-nodes molecule 
                                     (get-atom molecule atom-identifier-1)
                                     (get-atom molecule atom-identifier-2))")))

(defmethod remove-bond ((molecule molecule) (bond bond) &rest args)
  (declare (ignore args))
  (error "(graph:remove-edge molecule bond)"))

(defmethod find-bonds-containing ((molecule molecule) (atom atom))
  (error "(graph:find-edges-containing molecule atom)"))

(defmethod find-bond ((molecule molecule) (atom1 atom) (atom2 atom))
  (error "(graph:edgep molecule atom1 atom2)"))

#+(or)
(defun double-bond-configuration (molecule bond)
  (let ((atom1 (graph:node1 bond))
        (atom2 (graph:node2 bond)))
    (let ((atom1-bonds
            (remove bond (graph:find-edges-containing molecule atom1)))
          (atom2-bonds
            (remove bond (graph:find-edges-containing molecule atom2))))
      (cons atom1-bonds atom2-bonds))))

(defparameter *aromatic-atoms* '("c" "n" "o" "p" "s" "as" "se"))
(defparameter *organic-atoms* '("B" "C" "N" "O" "P" "S" "F" "Cl" "Br" "I"))

(defun aromatic-smiles-atom (string)
  (find string *aromatic-atoms* :test 'equal))

(defun in-double-bond (atm)
  (loop named bonds
        for bond in (chem:bonds-as-list atm)
        for atom1 = (chem:get-atom1 bond)
        for atom2 = (chem:get-atom2 bond)
        for order = (chem:get-order-from-atom bond atom1)
        when (eq order :double-bond)
          do (return-from in-double-bond bond))
  nil)

(defun remaining-atoms (atm vinyl1 vinyl2)
  "Collect the atoms bonded to atm that are not in a double bond with atm
and are not vinyl1 or vinyl2 (which may be nil)"
  (loop for bond in (chem:bonds-as-list atm)
        for order = (chem:get-order-from-atom bond atm)
        for other-atom = (chem:bond/get-other-atom bond atm)
        when (and (not (eq order :double-bond))
                  (not (or (eq other-atom vinyl1)
                           (eq other-atom vinyl2))))
            collect other-atom))

(defun parse-smiles-string (string &key name (add-implicit-hydrogens t))
  "Parses a SMILES description of a molecule and returns an instance
of the chem:MOLECULE class with the appropriate atoms and bonds."
  (let ((mol (make-molecule name))
        (ring-openings (make-hash-table))
        (element-count (make-hash-table))
        first-atom
        explicit-atoms
        aromatic)
    (labels ((read-number (stream)
               (parse-integer
                (coerce (loop for digit = (read-char stream)
                              for next = (peek-char nil stream nil)
                              collect digit
                              while (and next (digit-char-p next)))
                        'string)))
             (read-charge (stream charge-char)
               (let ((next-char (peek-char nil stream)))
                 (or
                  (let ((digit (digit-char-p next-char)))
                    (when digit (read-number stream)))
                  (let ((charge 1))
                    (when (and (not (char-equal next-char #\]))
                               (eql next-char charge-char))
                      (loop for char = (read-char stream)
                            do (incf charge)
                            while (eql charge-char (peek-char nil stream))))
                    charge))))
             (read-bracket-expression (stream)
               (let (isotope-number)
                 (let ((char (peek-char nil stream)))
                   (when (digit-char-p char)
                     (setf isotope-number (read-number stream)))
                   (let* ((element-string 
                            (coerce
                             (cons (read-char stream)
                                   (let ((next (peek-char nil stream)))
                                     (when (and (alpha-char-p next)
                                                (lower-case-p next))
                                       (cons (read-char stream) nil))))
                             'string))
                          (atom (add-molecule-atom
                                 (get-element element-string))))
                     (when (aromatic-smiles-atom element-string)
                       (setf aromatic t))
                     (when isotope-number
                       (let ((isotope (get-isotope (element atom) isotope-number)))
                         (setf (isotope atom) isotope)))
                     (loop for char = (peek-char nil stream)
                           do
                              (cond ((eq char #\])
                                     (read-char stream)
                                     (loop-finish))
                                    ((char-equal char #\H)
                                     (read-char stream)
                                     (let ((count
                                             (if (digit-char-p
                                                  (peek-char nil stream)) 
                                                 (read-number stream) 
                                                 1)))
                                       (dotimes (i count)
                                         (let ((h (add-molecule-atom (get-element "H"))))
                                           (add-bond mol atom h)))))
                                    ((char-equal char #\-)
                                     (read-char stream)
                                     (let ((charge (- (read-charge stream #\-))))
                                       (setf (charge atom) charge)))
                                    ((char-equal char #\+)
                                     (read-char stream)
                                     (let ((charge (read-charge stream #\+)))
                                       (setf (charge atom) charge)))
                                    ;; FIXME! Add support for @ and @@ here!!!
                                    ((char-equal char #\@)
                                     (read-char stream)
                                     (if (char-equal (peek-char nil stream) #\@)
                                         ;; it's clockwise
                                         (progn
                                           (read-char stream)
                                           (push (make-instance 'tetrahedral-center
                                                                :chiral-atom atom
                                                                :orientation :clockwise)
                                                 (chiral-arrangements mol)))
                                         ;; it's anticlockwise
                                         (push (make-instance 'tetrahedral-center
                                                              :chiral-atom atom
                                                              :orientation :anticlockwise)
                                               (chiral-arrangements mol))))
                                    ((char-equal char #\<)
                                     (read-char stream)
                                     (let* ((command (read-delimited-list #\> stream))
                                            (cmd (make-instance 'atom-command
                                                                :atm atom
                                                                :command command)))
                                       (push cmd (commands mol))))
                                    (t (let ((char (read-char stream)))
                                         (format t "unknown SMILES character! ~s~%" char)))))
                     (list (cons :explicit-atom atom))))))
             (add-molecule-atom (element)
               (let ((count (setf (gethash element element-count)
                                  (1+ (or (gethash element element-count) 0)))))
                 (let ((atom (add-atom mol 
                                       (intern (format nil "~A~A" element count) :keyword)
                                       element)))
                   ;; we need to squirrel away the first atom because
                   ;; the chirality w.r.t. implicit H atoms is
                   ;; different for the first atom.
                   (unless first-atom
                     (setf first-atom atom))
                   atom)))
             (read-branch (stream &optional source)
               (list (cons :branch
                           (loop for token = (read-smiles-stream stream source)
                                 while token))))
             (read-smiles-tokens (stream &optional source)
               (let ((char (read-char stream nil nil)))
                 (cond
                   ((null char) nil) 
                   ((eql char #\[)
                    (read-bracket-expression stream))
                   ((eql char #\()
                    (read-branch stream source))
                   ((eql char #\)) nil)
                   ((eql char #\-) (list (cons :bond :single-bond)))
                   ((eql char #\=) (list (cons :bond :double-bond)))
                   ((eql char #\#) (list (cons :bond :triple-bond)))
                   ((eql char #\:) (list (cons :bond :aromatic-bond)))
                   ((eql char #\/) (list (cons :bond :up-bond)))
                   ((eql char #\\) (list (cons :bond :down-bond)))
                   ((eql char #\.) (list (cons :disconnected nil)))
                   ((or (digit-char-p char)
                        (eql char #\%))
                    (let ((number (or (digit-char-p char)
                                      (read-number stream))))
                      (if number
                          (let* ((lookup (gethash number ring-openings)))
                            (if lookup
                                (progn
                                  (setf (gethash number ring-openings) nil)
                                  (list (cons :ring lookup)))
                                (list (cons :ring number))))
                          (error 'smiles-error
                                 :description "Couldn't read number!"))))

                   ;; Explicitly handle Cl and Br as they are the only
                   ;; multi-character atoms allowed in the SMILES
                   ;; organic subset. Other atoms have to be specified
                   ;; in brackets.
                   ((and (eql char #\C)
                         (eql (peek-char nil stream nil) #\l))
                    (read-char stream)
                    (list
                     (cons :atom (add-molecule-atom (get-element "Cl")))))
                   ((and (eql char #\B)
                         (eql (peek-char nil stream nil) #\r))
                    (read-char stream)
                    (list
                     (cons :atom (add-molecule-atom (get-element "Br")))))

                   (char
                    (unless (member (string char) *organic-atoms*
                                    :test #'string-equal) 
                      (error 'smiles-error
                             :description "Unexpected element symbol"))
                    (when (aromatic-smiles-atom (coerce (list char) 'string))
                      (setf aromatic t))
                    (list (cons :atom (add-molecule-atom
                                       (get-element
                                        ;; check for *, and use
                                        ;; element 0 for wildcard.
                                        (if (char-equal char #\*)
                                            0
                                            (string char))))))))))
             (read-smiles-stream (stream &optional last)
               (loop for tokens = (read-smiles-tokens stream last)
                     with bond-type = :single-bond
                     while tokens
                     do (loop for token in tokens
                              do (case (car token)
                                   (:disconnected (setf last nil))
                                   (:bond (setf bond-type (cdr token)))
                                   (:ring
                                    (let ((ring (cdr token)))
                                      (if (numberp ring)
                                          (setf (gethash ring ring-openings) last)
                                          (progn
                                            (loop named tets
                                                  for tet-center in (chiral-arrangements mol)
                                                  when (eq (chiral-atom tet-center) ring)
                                                    do (progn
                                                         (vector-push last (neighbors tet-center))
                                                         (return-from tets nil)))
                                            (add-bond mol ring last
                                                      :type (get-bond-order-keyword bond-type)
                                                      :order (get-bond-order-number bond-type))))))
                                   ((:atom :explicit-atom)
                                    (let ((atom (cdr token))
                                          bond-direction)
                                      (let ((last-arrangement
                                              (find last (chiral-arrangements mol) :key 'chiral-atom)))
                                        (when last-arrangement
                                          (vector-push atom (neighbors last-arrangement))))
                                      (let ((arrangement
                                              (find atom (chiral-arrangements mol) :key 'chiral-atom)))
                                        (when arrangement
                                          (vector-push last (neighbors arrangement))))
                                      (cond ((eql bond-type :up-bond)
                                             (setf bond-type :single-bond)
                                             (setf bond-direction :up-bond))
                                            ((eql bond-type :down-bond)
                                             (setf bond-type :single-bond)
                                             (setf bond-direction :down-bond)))
                                      (when last
                                        (apply #'add-bond mol last atom
                                               :type (get-bond-order-keyword bond-type)
                                               :order (get-bond-order-number bond-type)
                                               (when bond-direction
                                                 `(:direction ,bond-direction))))
                                      (setf bond-type (if aromatic
                                                          :aromatic-bond
                                                          :single-bond))
                                      (setf last atom)
                                      (when (eql (car token) :explicit-atom)
                                        (push atom explicit-atoms)))))))))
      (with-input-from-string (stream string)
        (read-smiles-stream stream)))
    ;; Now we need to do some post-processing.

    ;; 2. Find aromatic rings and replace Kekule structures with
    ;; explicitly aromatic rings. Unfortunately, we need to do this
    ;; BEFORE we add the hydrogens, I think, 
    
    #+nil 
    (multiple-value-bind (bonds cycles cycles-removed-mol)
        (graph:break-cycles mol)
      (print bonds)
      (print cycles)
      
      )
    ;; 1. Add implicit hydrogens
    (when add-implicit-hydrogens
      (let ((implicit-h-atoms (nth-value 1 (add-hydrogens mol explicit-atoms))))))

    ;; 1.5 Add implicit hydrogens to tetrahedral centers
    (loop for tet in (chiral-arrangements mol)
          for chiral-atom = (chiral-atom tet)
          for neighbors = (coerce (neighbors tet) 'list)
          for orientation = (orientation tet)
          for implicit-hydrogens = (loop named imp-hyd
                                        for bonded-neighbor in (chem:atom/bonded-atoms-as-list chiral-atom)
                                        for index from 0
                                        unless (member bonded-neighbor neighbors)
                                          collect bonded-neighbor)
          do (cond
               ((null implicit-hydrogens)
                ;; Do nothing
                )
               ((and (= (length implicit-hydrogens) 2)
                     (= (length (neighbors tet)) 2))
                (vector-push (first implicit-hydrogens) (neighbors tet))
                (vector-push (second implicit-hydrogens) (neighbors tet)))
               ((and (= (length implicit-hydrogens) 1)
                     (= (length (neighbors tet)) 3))
                (vector-push (first implicit-hydrogens) (neighbors tet)))
               ((and (null (elt (neighbors tet) 0))
                     (= (length implicit-hydrogens) 1))
                (setf (elt (neighbors tet) 0) (first implicit-hydrogens)))
               (t (error "I don't know where in ~s to put ~s" (neighbors tet) implicit-hydrogens)))
          do (progn
               (let ((neighbors (coerce (neighbors tet) 'list)))
                 (when (member nil neighbors)
                   (error "One of the neighbors of ~a is NIL: ~a" chiral-atom neighbors))
                 (chem:atom/reorder-bonds chiral-atom neighbors)
                 (chem:set-stereochemistry-type chiral-atom :chiral)
                 (ecase orientation
                   (:anticlockwise
                    (chem:set-configuration chiral-atom :right-handed))
                   (:clockwise
                    (chem:set-configuration chiral-atom :left-handed))))))

    ;; 3. Add double-bond-configurations for explicit spatial
    ;; arrangements around double bonds, based on the information
    ;; placed in the bonds when we were parsing them above.

    (let ((explicit-double-bonds (make-hash-table))
          )
      (loop for vinyl-bond in (double-bond-arrangements mol)
            for direction = (bond-direction vinyl-bond)
            for left-atom = (left-atom vinyl-bond)
            for right-atom = (right-atom vinyl-bond)
            for left-double-bond = (in-double-bond left-atom)
            for right-double-bond = (in-double-bond right-atom)
            ;; There is a double-bond to the left of the vinyl bond
            when left-double-bond
              do (let ((double-bond-geom (gethash left-double-bond explicit-double-bonds)))
                   (unless double-bond-geom
                     (setf double-bond-geom (make-instance 'double-bond
                                                           :bond left-double-bond)))
                   (setf (gethash left-double-bond explicit-double-bonds) double-bond-geom)
                   (ecase direction
                     (:up-bond
                      (setf (top-right-atom double-bond-geom) right-atom))
                     (:down-bond
                      (setf (bottom-right-atom double-bond-geom) right-atom))))
                 ;; There is a double-bond to the right of the vinyl bond
            when right-double-bond
              do (let ((double-bond-geom (gethash right-double-bond explicit-double-bonds)))
                   (unless double-bond-geom
                     (setf double-bond-geom (make-instance 'double-bond
                                                           :bond right-double-bond)))
                   (setf (gethash right-double-bond explicit-double-bonds) double-bond-geom)
                   (ecase direction
                     (:up-bond
                      (setf (bottom-left-atom double-bond-geom) left-atom))
                     (:down-bond
                      (setf (top-left-atom double-bond-geom) left-atom)))))
      (loop for dbond in (alexandria:hash-table-values explicit-double-bonds)
            for bond = (bond dbond)
            for left-atom = (chem:get-atom1 bond)
            for right-atom = (chem:get-atom2 bond)
            for left-remaining-atoms = (remaining-atoms left-atom (top-left-atom dbond) (bottom-left-atom dbond))
            for right-remaining-atoms = (remaining-atoms right-atom (top-right-atom dbond) (bottom-right-atom dbond))
            when left-remaining-atoms
              do (let ((original-remaining-atoms left-remaining-atoms))
                   (when (null (top-left-atom dbond))
                     (setf (top-left-atom dbond) (car left-remaining-atoms)
                           left-remaining-atoms (cdr left-remaining-atoms)))
                   (when (null (bottom-left-atom dbond))
                     (setf (bottom-left-atom dbond) (car left-remaining-atoms)
                           left-remaining-atoms (cdr left-remaining-atoms)))
                   (when left-remaining-atoms
                     (error "Inconsistent double-bond-geom ~a with left-remaining-atoms ~a" dbond original-remaining-atoms)))
            when right-remaining-atoms
              do (let ((original-remaining-atoms right-remaining-atoms))
                   (when (null (top-right-atom dbond))
                     (setf (top-right-atom dbond) (car right-remaining-atoms)
                           right-remaining-atoms (cdr right-remaining-atoms)))
                   (when (null (bottom-right-atom dbond))
                     (setf (bottom-right-atom dbond) (car right-remaining-atoms)
                           right-remaining-atoms (cdr right-remaining-atoms)))
                   (when right-remaining-atoms
                     (error "Inconsistent double-bond-geom ~a with right-remaining-atoms ~a" dbond original-remaining-atoms)))
            do (let ((atom-order (list (chem:get-name (top-left-atom dbond))
                                       (chem:get-name (bottom-left-atom dbond))
                                       (chem:get-name left-atom)
                                       (chem:get-name right-atom)
                                       (chem:get-name (top-right-atom dbond))
                                       (chem:get-name (bottom-right-atom dbond)))))
                 (chem:bond/set-property bond :double-bond-geometry atom-order)))
      )
    (values mol)))
