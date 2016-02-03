(in-package :atom-tree)

(defconstant +default-bond-angstroms+ 1.5)
(defconstant +default-angle-radians+ (* 109.5 0.0174533))
(defconstant +default-dihedral-radians+ (* 180.0 0.0174533))
(defclass atom-node ()
  ((target-atom :initarg target-atom :accessor target-atom)
   (target-pos :initform (geom:make-v3 0.0 0.0 0.0)
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
    (geom:set-using-bond (target-pos atom-node) (bond-angstroms atom-node) bond-pos)
    (target-pos atom-node)))

(defun build-from-bond-angle (atom-node)
  (let* ((bond-pos (target-pos (bond-node atom-node)))
         (angle-pos (target-pos (angle-node atom-node))))
    (geom:set-using-bond-angle (target-pos atom-node)
                                 (bond-angstroms atom-node) bond-pos
                                 (angle-radians atom-node) angle-pos)
    (target-pos atom-node)))

(defun build-from-bond-angle-dihedral (atom-node)
  (let* ((bond-pos (target-pos (bond-node atom-node)))
         (angle-pos (target-pos (angle-node atom-node)))
         (dihedral-pos (target-pos (dihedral-node atom-node))))
    (geom:set-using-bond-angle-dihedral (target-pos atom-node)
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
  (let ((spanning-loop (core:make-cxx-object 'chem:spanning-loop)))
    (chem:set-top spanning-loop root-atom)
    (let* ((spanning-atoms (loop for x = (chem:next spanning-loop (constantly t))
                              while x
                              collect x))
           (tree (make-instance 'atom-tree)))
      (format t "spanning-atoms: ~a~%" spanning-atoms)
      (let ((tree-nodes (loop for a in spanning-atoms
                           for bond-atom = (chem:get-back-span a)
                           for bond-node = (and bond-atom (gethash bond-atom (atoms-to-nodes tree)))
                           for angle-atom = (and bond-atom (chem:get-back-span bond-atom))
                           for angle-node = (and angle-atom (gethash angle-atom (atoms-to-nodes tree)))
                           for dihedral-atom = (and angle-atom (chem:get-back-span angle-atom))
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
                   (bond-angstroms (geom:v3-distance (chem:get-position ta)
                                                          (chem:get-position bond-atom)))
                   (angle-radians (geom:v3-angle (chem:get-position ta)
                                                      (chem:get-position bond-atom)
                                                      (chem:get-position angle-atom)))
                   (dihedral-radians (geom:dihedral (chem:get-position ta) ;
                                                    (chem:get-position bond-atom)
                                                    (chem:get-position angle-atom)
                                                    (chem:get-position dihedral-atom))))
              (setf (bond-angstroms node) bond-angstroms)
              (setf (angle-radians node) angle-radians)
              (setf (dihedral-radians node) dihedral-radians)))))
       
