;;;
;;;    File: atom-tree.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
(in-package :atom-tree)

(defconstant +default-bond-angstroms+ 1.5)
(defconstant +default-angle-radians+ (* 109.5 0.0174533))
(defconstant +default-dihedral-radians+ (* 180.0 0.0174533))
(defclass atom-node ()
  ((target-atom :initarg target-atom :accessor target-atom)
   (target-pos :initform (geom:vec 0.0 0.0 0.0)
               :initarg :target-pos :accessor target-pos)
   (bond-angstroms :initform +default-bond-angstroms+
                   :initarg :bond-angstroms :accessor bond-angstroms)
   (bond-node :initarg :bond-node :accessor bond-node)
   (angle-radians :initform +default-angle-radians+
                  :initarg :angle-radians :accessor angle-radians)
   (angle-node :initarg :angle-node :accessor angle-node)
   (dihedral-radians :initarg :dihedral-radians :accessor dihedral-radians)
   (dihedral-node :initarg :dihedral-node :accessor dihedral-node)))

(defun build-origin (atom-node)
  (geom:set-all3 (target-pos atom-node) 0.0 0.0 0.0)
  (target-pos atom-node))

(defun build-from-bond (atom-node)
  (let* ((bond-pos (target-pos (bond-node atom-node))))
    (geom:v-set-using-bond (target-pos atom-node) (bond-angstroms atom-node) bond-pos)
    (target-pos atom-node)))

(defun build-from-bond-angle (atom-node)
  (let* ((bond-pos (target-pos (bond-node atom-node)))
         (angle-pos (target-pos (angle-node atom-node))))
    (geom:v-set-using-bond-angle (target-pos atom-node)
                                 (bond-angstroms atom-node) bond-pos
                                 (angle-radians atom-node) angle-pos)
    (target-pos atom-node)))

(defun build-from-bond-angle-dihedral (atom-node)
  (let* ((bond-pos (target-pos (bond-node atom-node)))
         (angle-pos (target-pos (angle-node atom-node)))
         (dihedral-pos (target-pos (dihedral-node atom-node))))
    (geom:v-set-using-bond-angle-dihedral (target-pos atom-node)
                                        (bond-angstroms atom-node) bond-pos
                                        (angle-radians atom-node) angle-pos
                                        (dihedral-radians atom-node) dihedral-pos)
    (target-pos atom-node)))

(defun build-atom-pos (atom-node)
  (cond ((dihedral-node atom-node) (build-from-bond-angle-dihedral atom-node))
        ((angle-node atom-node) (build-from-bond-angle atom-node))
        ((bond-node atom-node) (build-from-bond atom-node))
        (t (build-origin atom-node))))


(defclass atom-tree ()
  ((nodes :initarg :nodes :accessor nodes)
   (atoms-to-nodes :initform (make-hash-table) :accessor atoms-to-nodes)
   (coupled-dihedrals :initform (make-hash-table :test #'equal) :accessor coupled-dihedrals)))

(defun build-atom-tree-for-molecule (mol root-atom)
  (let ((spanning-loop (chem:make-spanning-loop root-atom)))
    (let* ((spanning-atoms (loop for x = (chem:next spanning-loop (constantly t))
                              while x
                              collect x))
           (tree (make-instance 'atom-tree)))
      (format t "spanning-atoms: ~a~%" spanning-atoms)
      (let ((tree-nodes (loop for a in spanning-atoms
                           for bond-atom = (chem:get-back-span spanning-loop a)
                           for bond-node = (and bond-atom (gethash bond-atom (atoms-to-nodes tree)))
                           for angle-atom = (and bond-atom (chem:get-back-span spanning-loop bond-atom))
                           for angle-node = (and angle-atom (gethash angle-atom (atoms-to-nodes tree)))
                           for dihedral-atom = (and angle-atom (chem:get-back-span spanning-loop angle-atom))
                           for dihedral-node = (and dihedral-atom (gethash dihedral-atom (atoms-to-nodes tree)))
                           for atom-node = (make-instance 'atom-node
                                                        :bond-node bond-node
                                                        :angle-node angle-node
                                                        :dihedral-node dihedral-node)
                           do (setf (target-atom atom-node) a)
                           do (setf (gethash a (atoms-to-nodes tree)) atom-node)
                           collect atom-node)))
        (setf (nodes tree) tree-nodes)
        tree))))

(defun build-atom-tree (mol &key root-property root-atom)
  (and root-property root-atom (error "Specify either root-property or root-atom but not both"))
  (when root-property
    (setq root-atom (car (select:atoms-with-property mol root-property))))
  (build-atom-tree-for-molecule mol root-atom))


(defun build-geometry-for-atom-tree (atom-tree)
  (loop for node in (nodes atom-tree)
     do (build-atom-pos node)))

(defun write-geometry-to-atoms (atom-tree)
  (loop for node in (nodes atom-tree)
       do (chem:set-position (target-atom node) (target-pos node))))
      

(defun find-coupled-dihedrals (atom-tree)
  (mapc (lambda (node)
          (let* ((bond-node (bond-node node))
                 (angle-node (angle-node node))
                 (dihedral-node (dihedral-node node))
                 (key (list bond-node angle-node dihedral-node))
                 (coupled (gethash key (coupled-dihedrals atom-tree))))
            (push node coupled)
            (setf (gethash key (coupled-dihedrals atom-tree)) coupled)))
        (nodes atom-tree)))

(defun assign-coupled-dihedral-angles (atom-tree)
  (maphash (lambda (k nodes)
             (case (length nodes)
               (3 (mapc (lambda (node angle)
                          (setf (dihedral-radians node) angle))
                        nodes (list (* 60.0 0.0174533)
                                    (* 180.0 0.0174533)
                                    (* 240.0 0.0174533))))
               (2 (mapc (lambda (node angle) (setf (dihedral-radians node) angle))
                        nodes (list (* 0.0 0.0174533)
                                    (* 180.0 0.0174533))))
               (1 (mapc (lambda (node angle) (setf (dihedral-radians node) angle))
                        nodes (list (* 0.0 0.0174533))))
               (t (error "Illegal number of atoms ~a share the dihedral ~a - ~a"
                         (length nodes)
                         (bond-node (car nodes))
                         (angle-node (car nodes))))))
           (coupled-dihedrals atom-tree)))

(and 1 2)

(defun extract-internal-coordinates-for-active-set (atom-tree active-set)
  (loop for node in (nodes atom-tree)
     for bond-atom = (and (bond-node node) (target-atom (bond-node node)))
     for angle-atom = (and (angle-node node) (target-atom (angle-node node)))
     for dihedral-atom = (and (dihedral-node node) (target-atom (dihedral-node node)))
     do (if (and bond-atom
                 angle-atom
                 dihedral-atom
                 (gethash bond-atom active-set)
                 (gethash angle-atom active-set)
                 (gethash dihedral-atom active-set))
            (let* ((ta (target-atom node))
                   (bond-angstroms (geom:vdistance (chem:get-position ta)
                                                          (chem:get-position bond-atom)))
                   (angle-radians (geom:vangle (chem:get-position ta)
                                                      (chem:get-position bond-atom)
                                                      (chem:get-position angle-atom)))
                   (dihedral-radians (geom:vdihedral (chem:get-position ta) ;
                                                    (chem:get-position bond-atom)
                                                    (chem:get-position angle-atom)
                                                    (chem:get-position dihedral-atom))))
              (setf (bond-angstroms node) bond-angstroms)
              (setf (angle-radians node) angle-radians)
              (setf (dihedral-radians node) dihedral-radians)))))
