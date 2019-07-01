;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :fep)

(defun make-fep ()
  (make-instance 'fep:fep-calculation))

(defun version ()
  *version*)

(defun core-group-p (group)
  (loop for atom in group
        for label = (chem:matter-get-property-or-default atom :label nil)
        when (and label (find #\? (string label)))
          do (return-from core-group-p t))
  nil)

(defun group-name (group)
  (loop for atom in group
        for label = (chem:matter-get-property-or-default atom :label nil)
        when label
          do (return-from group-name label)))

(defgeneric parse-residue-groups (sketch))

(defmethod parse-residue-groups ((sketch string))
  (multiple-value-bind (groups map-names-numbers)
      (jsme:parse-groups sketch)
    (let (core
          (side-chains (make-hash-table)))
      (loop for group in groups
            if (core-group-p group)
              do (setf core group)
            else do (setf (gethash (group-name group) side-chains) group))
      (values core side-chains map-names-numbers))))

(defmethod parse-residue-groups ((sketch chem:chem-draw))
  (let* ((agg (chem:as-aggregate sketch))
         core-molecule
         (core-atom-count 0)
         side-chain-molecules)
    (let ((num-bonds 0)
          (sum-lengths 0.0))
      (chem:map-bonds nil (lambda (a1 a2 order)
                            (let* ((v1 (chem:get-position a1))
                                   (v2 (chem:get-position a2))
                                   (delta (geom:v- v1 v2))
                                   (len (geom:vlength delta)))
                              (setf sum-lengths (+ sum-lengths len)
                                    num-bonds (1+ num-bonds))))
                      agg)
      (let* ((average-length (/ sum-lengths num-bonds))
             (scale-factor (/ 1.5 average-length)))
        (format t "Average bond length: ~a~%" average-length)
        ;; Scale bonds to 1.5 Angstroms
        (chem:map-atoms nil (lambda (a)
                          (let ((pos (chem:get-position a)))
                            (chem:set-position a (geom:v* pos scale-factor))))
                        agg)))
    (chem:map-molecules nil (lambda (m)
                              (let ((num-atoms (chem:number-of-atoms m)))
                                (when (> num-atoms core-atom-count)
                                  (setf core-molecule m
                                        core-atom-count num-atoms))))
                        agg)
    (chem:map-molecules nil (lambda (m)
                              (unless (eq m core-molecule)
                                (push m side-chain-molecules)))
                        agg)
    ;; Handle the core-group
    (let ((core-group (chem:map-atoms 'list #'identity core-molecule))
          (map-names-numbers (make-hash-table)))
      (mapcar (lambda (a)
                (if (eq (chem:get-element a) :ne)
                    (let ((bonds (chem:bonds-as-list a)))
                      (unless (= (length bonds) 1)
                        (warn "There are the wrong number of bonds from ~a bonds: ~a" a bonds))
                      (let* ((bond (first bonds))
                             (other-atom (chem:get-other-atom bond a))
                             (attach-sym (intern (string-trim "?" (string (chem:get-name a))) :keyword)))
                        (chem:set-property other-atom :attach attach-sym))
                      (chem:set-property a :label (chem:get-name a)))
                    (let* ((id (chem:matter-get-property-or-default a :id nil)))
                      (when id
                        (setf (gethash (chem:get-name a) map-names-numbers) id)))))
              core-group)
      ;; Handle the side-chains
      (let ((side-chain-groups (make-hash-table)))
        (loop for side-chain-molecule in side-chain-molecules
              do (let* (side-chain-name
                        (side-chain-atoms (chem:map-atoms
                                           'list
                                           (lambda (a)
                                             (when (eq (chem:get-element a) :ne)
                                               (setf side-chain-name (chem:get-name a))
                                               (let* ((other-name (intern (subseq (string side-chain-name) 0 2) :keyword))
                                                      (bonds (chem:bonds-as-list a)))
                                                 (unless (= (length bonds) 1)
                                                   (warn "There are the wrong number of bonds ~a bonds: ~a" a bonds))
                                                 (let* ((bond (first bonds))
                                                        (other-atom (chem:get-other-atom bond a)))
                                                   (chem:set-name other-atom other-name)
                                                   (chem:set-property a :label (chem:get-name a)))))
                                             a)
                                           side-chain-molecule)))
                   (setf (gethash side-chain-name side-chain-groups) side-chain-atoms)))
        (values core-group side-chain-groups map-names-numbers)))))



(defun labeled-atom (atom)
  (first (chem:bonded-atoms-as-list atom)))

(defclass bond-vector ()
  ((start-vector :initarg :start-vector :accessor start-vector)
   (direction :initarg :direction :accessor direction)))

(defun build-core-residue (atoms)
  (let ((residue (chem:make-residue :cor))
        (out-atom-vectors (make-hash-table)))
    (loop for atom in atoms
          for label = (chem:matter-get-property-or-default atom :label nil)
          if label
            do (let* ((labeled-atom (labeled-atom atom))
                      (short-label (intern (subseq (string label) 0 2) :keyword))
                      (start (chem:get-position atom))
                      (dir (geom:v- (chem:get-position labeled-atom) start))
                      (bond-vector (make-instance 'bond-vector :start-vector start :direction dir)))
                 (chem:set-property labeled-atom :attach short-label)
                 (setf (gethash short-label out-atom-vectors) bond-vector)
                 (chem:remove-bond-to atom labeled-atom))
          else do (chem:add-matter residue atom))
    (values residue out-atom-vectors)))

(defun short-name (name)
  (intern (subseq (string name) 0 2) :keyword))

(defun atom-with-property (residue property property-value)
  (chem:map-atoms nil
                  (lambda (a)
                    (let ((pv (chem:matter-get-property-or-default a property nil)))
                      (when (eq pv property-value)
                        (return-from atom-with-property a))))
                  residue)
  nil)

(defun build-side-chain-residue (atoms side-chain-name core-residue core-vector)
  (let* ((residue (chem:make-residue side-chain-name))
         (short-side-chain-name (short-name side-chain-name))
         in-atom-vector attach-atom
         (core-atom (atom-with-property core-residue :attach short-side-chain-name)))
    (loop for atom in atoms
          for label = (chem:matter-get-property-or-default atom :label nil)
          if label
            do (let* ((labeled-atom (labeled-atom atom))
                      (short-label (intern (subseq (string label) 0 2) :keyword))
                      (start-vec (chem:get-position atom))
                      (dir (geom:v- (chem:get-position labeled-atom) start-vec)))
                 (setf attach-atom labeled-atom)
                 (chem:set-name labeled-atom short-label)
                 (setf in-atom-vector (make-instance 'bond-vector :start-vector start-vec :direction dir))
                 (chem:remove-bond-to atom labeled-atom))
          else do (chem:add-matter residue atom))
    (unless in-atom-vector
      (error "There is no in-atom-vector for atoms: ~a~%" atoms))
    (let* ((cross (geom:vcross (geom:vnormalized (direction in-atom-vector)) (geom:vnormalized (direction core-vector))))
           (angle (* 1.0 (asin (geom:vlength cross))))
           (residue-to-origin (geom:make-m4-translate (geom:v* (start-vector in-atom-vector) -1.0)))
           (residue-to-core (geom:make-m4-translate (chem:get-position core-atom)))
           (rotate (geom:make-m4-rotate-z angle)))
      (progn
        (chem:apply-transform-to-atoms residue residue-to-origin)
        (chem:apply-transform-to-atoms residue rotate)
        (chem:apply-transform-to-atoms residue residue-to-core)))
    (values residue in-atom-vector)))

(defvar *core-group* nil)
(defvar *side-chain-groups* nil)
(defvar *map-names-numbers* nil)

(defun build-topologys (sketch &key verbose)
  (let ((side-chains (make-hash-table))
        (topology-to-residue (make-hash-table)))
    (multiple-value-bind (core-group side-chain-groups map-names-numbers)
        (parse-residue-groups sketch)
      (setf *core-group* core-group
            *side-chain-groups* side-chain-groups
            *map-names-numbers* map-names-numbers)
      (format t "In build-topologys core-group: ~a~%" core-group)
      (format t "side-chain-groups: ~a~%" side-chain-groups)
      (format t "map-names-numbers: ~a~%" map-names-numbers)
      (multiple-value-bind (core-residue core-vectors)
          (build-core-residue core-group)
        (let ((core-topology (cando:make-simple-topology-from-residue core-residue)))
          (setf (gethash core-topology topology-to-residue) core-residue)
          (maphash (lambda (name vec)
                     (let* ((plug-name (intern (format nil "+~a" (string name)) :keyword))
                            (plug-atom (atom-with-property core-residue :attach name))
                            (plug-atom-name (chem:get-name plug-atom))
                            #+(or)(_ (format t "Making out plug plug-name: ~a    name: ~a~%" plug-name plug-atom-name))
                            (out-plug (chem:make-out-plug plug-name nil nil plug-atom-name :single-bond)))
                       (chem:add-plug core-topology (chem:get-name out-plug) out-plug)))
                   core-vectors)
          (loop for side-chain-name being the hash-keys in side-chain-groups using (hash-value side-chain-group)
                with side-chain-residue
                with side-chain-vector
                do (multiple-value-setq (side-chain-residue side-chain-vector)
                     (build-side-chain-residue side-chain-group
                                               side-chain-name
                                               core-residue
                                               (gethash (short-name side-chain-name) core-vectors)))
                do (let* ((side-chain-topology (cando:make-simple-topology-from-residue side-chain-residue))
                          (short-name (intern (subseq (string (chem:get-name side-chain-residue)) 0 2) :keyword))
                          (plug-name (intern (format nil "-~a" (string short-name)) :keyword))
                          (in-plug (chem:make-in-plug plug-name nil short-name :single-bond)))
                     (setf (gethash side-chain-topology topology-to-residue) side-chain-residue)
                     (chem:add-plug side-chain-topology plug-name in-plug)
                     (push (cons (chem:get-name side-chain-residue) side-chain-topology) (gethash short-name side-chains))))
          (values core-topology side-chains topology-to-residue map-names-numbers))))))



                                   
#+(or)
(defmethod build-topologys-chemdraw (chemdraw)
  (let ((side-chains (make-hash-table))
        (topology-to-residue (make-hash-table))
        (agg (chem:as-aggregate chemdraw)))
    (break "agg = ~a" agg)
    (multiple-value-bind (core-group side-chain-groups map-names-numbers)
        (parse-chemdraw chemdraw)
      (multiple-value-bind (core-residue core-vectors)
          (build-core-residue core-group)
        (let ((core-topology (cando:make-simple-topology-from-residue core-residue)))
          (setf (gethash core-topology topology-to-residue) core-residue)
          (maphash (lambda (name vec)
                     (let* ((plug-name (intern (format nil "+~a" (string name)) :keyword))
                            (plug-atom (atom-with-property core-residue :attach name))
                            (plug-atom-name (chem:get-name plug-atom))
                            #+(or)(_ (format t "Making out plug plug-name: ~a    name: ~a~%" plug-name plug-atom-name))
                            (out-plug (chem:make-out-plug plug-name nil nil plug-atom-name :single-bond)))
                       (chem:add-plug core-topology (chem:get-name out-plug) out-plug)))
                   core-vectors)
          (loop for side-chain-name being the hash-keys in side-chain-groups using (hash-value side-chain-group)
                with side-chain-residue
                with side-chain-vector
                do (multiple-value-setq (side-chain-residue side-chain-vector)
                     (build-side-chain-residue side-chain-group
                                               side-chain-name
                                               core-residue
                                               (gethash (short-name side-chain-name) core-vectors)))
                do (let* ((side-chain-topology (cando:make-simple-topology-from-residue side-chain-residue))
                          (short-name (intern (subseq (string (chem:get-name side-chain-residue)) 0 2) :keyword))
                          (plug-name (intern (format nil "-~a" (string short-name)) :keyword))
                          (in-plug (chem:make-in-plug plug-name nil short-name :single-bond)))
                     (setf (gethash side-chain-topology topology-to-residue) side-chain-residue)
                     (chem:add-plug side-chain-topology plug-name in-plug)
                     (push (cons (chem:get-name side-chain-residue) side-chain-topology) (gethash short-name side-chains))))
          (values core-topology side-chains topology-to-residue map-names-numbers))))))

(defun cartesian (l) 
  (cond ((null l) nil) 
        ((null (cdr l)) 
         (mapcar #'list (car l))) 
        (t (mapcan #'(lambda (x) (mapcar #'(lambda (y) (cons y x)) (car l))) 
                   (cartesian (cdr l))))))

(defvar *fep-time* 0)

(defun incf-fep-time ()
  (incf *fep-time*))

(defun print-object-readably-with-slots (obj stream)
  (format stream "#$(~s " (class-name (class-of obj)))
  (loop for slot in (clos:class-slots (class-of obj))
        for slot-name = (clos:slot-definition-name slot)
        for initargs = (clos:slot-definition-initargs slot)
        if (and (car initargs) (slot-boundp obj slot-name))
          do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
  (format stream ") "))

(defclass job-graph ()
  ((nodes :initform nil :initarg :nodes :accessor nodes)
   (morphs :initform nil :initarg :morphs :accessor morphs)))

(defmethod print-object ((obj job-graph) stream)
  (if *print-readably*
      (progn
        (format t "Printing job-graph~%")
        (print-object-readably-with-slots obj stream))
      (print-unreadable-object (obj stream)
        (format stream "~a" (class-name (class-of obj))))))

(defclass fep-morph ()
  ((source :initarg :source :accessor source)
   (target :initarg :target :accessor target)
   (stages :initarg :stages :initform 1 :accessor stages)
   (windows :initarg :windows :initform 11 :accessor windows)))

(defmethod print-object ((obj fep-morph) stream)
  (if *print-readably*
      (print-object-readably-with-slots obj stream)
      (print-unreadable-object (obj stream)
        (format stream "~a" (class-name (class-of obj))))))

(defun morph-string (fep-morph)
  (let ((source-name (name (source fep-morph)))
        (target-name (name (target fep-morph))))
    (format nil "~a-~a" source-name target-name)))

(defclass fep-structure ()
  ((name :initarg :name :accessor name)
   (drawing :initarg :drawing :accessor drawing)
   (core-residue :initarg :core-residue :accessor core-residue)
   (core-residue-name :initarg :core-residue-name :accessor core-residue-name)
   (side-chain-residue-names :initarg :side-chain-residue-names :accessor side-chain-residue-names)
   (core-atoms :initarg :core-atoms :accessor core-atoms)
   (side-chain-atoms :initarg :side-chain-atoms :accessor side-chain-atoms)
   (molecule :initarg :molecule :accessor molecule)
   (atom-order :initarg :atom-order :accessor atom-order)
   (net-charge :initform 0.0  :initarg :net-charge :accessor net-charge)
   (am1-charges :initarg :am1-charges :accessor am1-charges)
   (am1-bcc-charges :initarg :am1-bcc-charges :accessor am1-bcc-charges)))


(defmethod print-object ((obj fep-structure) stream)
  (if *print-readably*
      (print-object-readably-with-slots obj stream)
      (print-unreadable-object (obj stream)
        (format stream "~a ~a" (class-name (class-of obj)) (string (name obj))))))

(defun mol2-safe-residue-name (calculation name)
  name)

(defun mol2-safe-atom-name (calculation name)
  name)

(defclass calculation ()
  ((receptors :initform nil :initarg :receptors :accessor receptors)
   (ligands :initarg :ligands :accessor ligands)
   (mask-method :initform :default :initarg :mask-method :accessor mask-method)
   (core-topology :initarg :core-topology :accessor core-topology)
   (side-topologys :initarg :side-topologys :accessor side-topologys)
   (jobs :initarg :jobs :accessor jobs)
   (ti-lambdas :initarg :ti-lambdas :initform 11 :accessor ti-lambdas)
   (top-directory :initform (make-pathname :directory (list :relative "jobs"))
                  :initarg :top-directory :accessor top-directory )
   (stage :initform 0 :initarg :stage :accessor stage)
   (residue-name-to-pdb :initform (make-hash-table) :initarg :residue-name-to-pdb :accessor residue-name-to-pdb)
   ))

(defmethod print-object ((obj calculation) stream)
  (if *print-readably*
      (progn
        (format stream "#.(make-instance 'fep::calculation ")
        (loop for slot in (clos:class-slots (find-class 'fep::fep-calculation))
              for slot-name = (clos:slot-definition-name slot)
              for initargs = (clos:slot-definition-initargs slot)
              if (and (car initargs) (slot-boundp obj slot-name))
                do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
        (format stream ") "))
      (print-unreadable-object (obj stream)
        (format stream "feps"))))

(defun find-morph-with-name (name calculation)
  (loop for morph in (morphs (jobs calculation))
        when (string= (morph-string morph) (string name))
          do (return-from find-morph-with-name morph))
  (error "Could not find morph with name ~s in calculation with morphs: ~s" name (mapcar (lambda (m) (name m)) (morphs (jobs calculation)))))

(defparameter *top-guard* nil)

    
(defmacro with-top-directory ((calculation) &body body)
  `(progn
     ,@body))

(defmethod make-load-form ((object calculation) &optional env)
  (make-load-form-saving-slots object))

(defun unique-name (atom counters)
  (let* ((atom-name (chem:get-name atom))
         (element-name (string (chem:get-element atom)))
         (atom-name-string (string atom-name)))
    (let ((new-name (format nil "~a~a" element-name (incf (gethash element-name counters 0)))))
      (intern new-name :keyword))))

(defun combine-into-single-residue (molecule core-residue-name)
  (let* ((residues (chem:map-residues 'list #'identity molecule))
         (core-residue (unique-residue-with-name molecule core-residue-name))
         (side-residues (remove core-residue residues))
         (_ (when (= (length side-residues) (length residues))
              (error "The core residue ~s was not removed from ~s" core-residue residues)))
         (unique-atom-name-strings (chem:map-atoms 'list (lambda (latom) (string (chem:get-name latom))) core-residue))
         (element-name-counters (let ((counters (make-hash-table :test #'equal)))
                                  (loop for name in unique-atom-name-strings
                                        for first-digit-position = (position-if #'digit-char-p name :start 1)
                                        if first-digit-position
                                          do (let* ((element (subseq name 0 first-digit-position))
                                                    (num (parse-integer name :start first-digit-position))
                                                    (counters-num (gethash element counters 0)))
                                               (when (< counters-num num)
                                                 (setf (gethash element counters) (1+ num)))))
                                  counters)))
    ;; Ensure that all atom names in side-residues are unique
    (loop for residue in side-residues
          do (cando:do-atoms (atom1 residue)
               (let ((unique-atom-name (unique-name atom1 element-name-counters)))
                 (chem:set-name atom1 unique-atom-name))))
    ;; Now merge all of the atoms from side-residues into the core residue
    (loop for residue in side-residues
          for residue-atoms = (chem:map-atoms 'list #'identity residue)
          do (loop for atom2 in residue-atoms
                   do (chem:remove-atom-dont-delete-bonds residue atom2)
                      (chem:add-atom core-residue atom2)))
    ;; Now remove all side-residues from the molecule
    (loop for residue in side-residues
          do (chem:remove-residue molecule residue))))
    
(defun structure-name (tops)
  (let ((names (mapcar (lambda (x) (subseq (string (car x)) 1 3)) tops)))
    (intern (apply #'concatenate 'string (sort names #'string<)) :keyword)))

(defun set-coordinates-for-residue (res other-res)
  (chem:map-atoms
   nil
   (lambda (a)
     (let* ((other-atom (chem:atom-with-name other-res (chem:get-name a)))
            (pos (chem:get-position other-atom)))
       (chem:set-position a pos)))
   res))

(defun set-stereocenters (mol)
  (let ((stereocenters 
          (sort (cando:gather-stereocenters mol) #'string< :key #'chem:get-name)))
    (cando:set-stereoisomer-func stereocenters (constantly :S))))

(defun build-feps (core-topology side-topologys topology-to-residue map-names-numbers)
  (let ((top-groups nil))
    (maphash (lambda (name tops)
               (push tops top-groups))
             side-topologys)
    (let ((cartesian (cartesian top-groups)))
      (loop for cart in cartesian
            for name = (structure-name cart)
            for molecule = (chem:make-molecule name)
            for core-residue = (let ((res (chem:build-residue-single-name core-topology)))
                                 (set-coordinates-for-residue res (gethash core-topology topology-to-residue))
                                 res)
            collect (let (side-chain-residue-names)
                      (chem:add-matter molecule core-residue)
                      (loop for (name . top) in cart
                            for side-residue = (let ((res (chem:build-residue-single-name top)))
                                                 (set-coordinates-for-residue res (gethash top topology-to-residue))
                                                 res)
                            for in-plug = (chem:get-in-plug top)
                            for in-plug-name = (chem:get-name in-plug)
                            for out-plug-name = (chem:other-plug-name in-plug-name)
                            for side-atoms = (chem:map-atoms 'list #'identity side-residue)
                            do (progn
                                 (push (chem:get-name side-residue) side-chain-residue-names)
                                 (chem:add-matter molecule side-residue)
                                 (maphash (lambda (name number)
                                            (let ((atom (chem:atom-with-name core-residue name)))
                                              #+(or)(format t "Atom: ~a   number: ~a~%" atom number)
                                              (chem:set-property atom :map number)))
                                          map-names-numbers)
                                 (chem:connect-residues core-topology
                                                        core-residue
                                                        out-plug-name
                                                        top
                                                        side-residue
                                                        in-plug-name)))
                      (set-stereocenters molecule)
                      (make-instance 'fep-structure :name name :drawing molecule
                                                    :core-residue core-residue
                                                    :core-residue-name (chem:get-name core-residue)
                                                    :side-chain-residue-names side-chain-residue-names))))))

(defun setup-ligands (feps sketch &key verbose)
  (let ((ligands (multiple-value-bind (core-topology side-chains topology-to-residue map-names-numbers)
                     (build-topologys sketch :verbose verbose)
                   (setf (core-topology feps) core-topology)
                   (setf (side-topologys feps) side-chains)
                   (build-feps core-topology side-chains topology-to-residue map-names-numbers))))
    (loop for ligand in ligands
          for mol = (chem:matter-copy (drawing ligand))
          for side-chain-atoms = nil
          do (chem:fill-in-implicit-hydrogens mol)
             (ensure-unique-hydrogen-names mol)
             (cando:build-unbuilt-hydrogens mol)
             (setf (molecule ligand) mol)
             (let ((core-atoms (chem:map-atoms 'list #'identity
                                               (unique-residue-with-name mol (core-residue-name ligand)))))
               (setf (core-atoms ligand) core-atoms))
             (loop for side-chain-name in (side-chain-residue-names ligand)
                   for side-chain-residue = (unique-residue-with-name mol side-chain-name)
                   do (cando:do-atoms (a side-chain-residue)
                        (push a side-chain-atoms)))
             (setf (side-chain-atoms ligand) side-chain-atoms)
             (combine-into-single-residue mol (core-residue-name ligand)))
    (setf (ligands feps) ligands)))


(defun layout-ligands (calculation &key (xdelta 15.0) (ydelta 15.0) (accessor 'drawing))
  (let* ((ligands (copy-list (ligands calculation)))
         (numx (floor (sqrt (length ligands))))
         (numy (1+ numx))
         (agg (chem:make-aggregate nil)))
    (loop for ix from 0 below numx
          for x = (* (- ix (floor (/ numx 2))) xdelta)
          do (loop for iy from 0 below numy
                   for y = (* (- iy (floor (/ numy 2))) ydelta)
                   for transform = (geom:make-m4-translate (geom:vec x y 0.0))
                   do (let ((fep (pop ligands)))
                        (when fep
                          (let ((molecule (chem:matter-copy (funcall accessor fep))))
                            (chem:apply-transform-to-atoms molecule transform)
                            (chem:add-matter agg molecule))))))
    agg))

(defun layout-2 (name1 name2 feps &key (accessor 'drawing))
  (let ((mol1 (chem:matter-copy (funcall accessor (find name1 feps :key 'name))))
        (mol2 (chem:matter-copy (funcall accessor (find name2 feps :key 'name))))
        (agg (chem:make-aggregate 'two))
        (translate (geom:make-m4-translate (geom:vec 20.0 0.0 0.0))))
    (chem:add-matter agg mol2)
    (chem:apply-transform-to-atoms mol1 translate)
    (chem:add-matter agg mol1)
    agg))

(defun ensure-unique-hydrogen-names (molecule)
  (let ((hydrogen-counter 0))
    (cando:do-atoms (atom1 molecule)
      (when (eq (chem:get-element atom1) :H)
        (let ((hydrogen-name (intern (format nil "H~a" (incf hydrogen-counter)) :keyword)))
          (chem:set-name atom1 hydrogen-name))))))

(defun build-ligands (calculations)
  (let ((feps (ligands calculations)))
    (loop for fep in feps
          for mol = (molecule fep)
          do (format t "Building ~a~%" (name fep))
             (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) mol)
             (setf (molecule fep) mol)
          do (format t "build-structures fep: ~a   map-atoms: ~a~%" fep (sorted-map-atoms fep)))))

(defun minimize-ligands (calculation)
  (let ((feps (ligands calculation)))
    (loop for fep in feps
          for mol = (molecule fep)
          for agg = (chem:make-aggregate nil)
          do (chem:add-matter agg mol)
          do (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) agg)
          do (energy:minimize agg :cg-tolerance 0.0001 :max-tn-steps 100 :tn-tolerance 0.00001))))

(defun sorted-map-atoms (fep)
  (let ((molecule (molecule fep))
        (map-atoms nil))
    (chem:map-atoms nil (lambda (a)
                          (when (chem:matter-get-property-or-default a :map nil)
                            (push (cons (chem:matter-get-property a :map) a) map-atoms)))
                    molecule)
    (sort map-atoms #'< :key #'car)))


(defun anchor-to-pose (moveable-atoms fixed-pose-atoms &key stereochemical-restraints)
  (mapc (lambda (fixed-pose-atom moveable-atom)
          (let* ((pos (chem:get-position fixed-pose-atom))
                 (anchor-restraint (core:make-cxx-object
                                    'chem:restraint-anchor
                                    :atom moveable-atom
                                    :position pos
                                    :weight 10000.0)))
            (chem:clear-restraints moveable-atom)
            (chem:add-restraint moveable-atom anchor-restraint)
            ;; Assign stereochemical restraint if there is one that matches
            (let* ((map-id (chem:matter-get-property-or-default moveable-atom :map nil))
                   (stereochemistry-assignment (assoc map-id stereochemical-restraints)))
              (when stereochemistry-assignment
                (chem:set-stereochemistry-type moveable-atom :chiral)
                (format t "chem:get-stereochemistry-type ~a -> ~a~%" moveable-atom (chem:get-stereochemistry-type moveable-atom))
                (chem:set-configuration moveable-atom (cdr stereochemistry-assignment))
                (format t "Set stereochemistry of ~a to ~a~%" moveable-atom (cdr stereochemistry-assignment))))
            (format t "Anchored ~a to ~a restraints: ~a props: ~s~%" moveable-atom pos (chem:all-restraints moveable-atom) (chem:properties moveable-atom))))
        fixed-pose-atoms moveable-atoms))

(defun pose-ligands-using-pattern (calculation pattern docked-molecule &key stereochemical-restraints)
  (let* ((fixed-atoms-map (pattern-atoms pattern docked-molecule))
         (fixed-order (let (fo)
                        (maphash (lambda (index atom)
                                   (push index fo))
                                 fixed-atoms-map)
                        fo))
         (fixed-atoms (mapcar (lambda (index) (gethash index fixed-atoms-map)) fixed-order))
         (feps (ligands calculation)))
    (loop for fep in feps
          for molecule = (molecule fep)
          for moveable-atoms-map = (pattern-atoms pattern molecule)
          for moveable-atoms = (mapcar (lambda (index) (gethash index moveable-atoms-map)) fixed-order)
          do (format t "fep ~a~%" fep)
             (format t " moveable-atoms -> ~a~%" fep moveable-atoms)
             (format t "    fixed-atoms -> ~a~%" fixed-atoms)
          do (chem:superpose-one molecule moveable-atoms fixed-atoms)
          do (anchor-to-pose moveable-atoms fixed-atoms :stereochemical-restraints stereochemical-restraints)))
  (minimize-ligands calculation))



(defun pose-ligands (calculation fixed-atoms &key stereochemical-restraints)
  (let ((feps (ligands calculation)))
    (loop for fep in feps
          for molecule = (molecule fep)
          for moveable-atoms = (mapcar #'cdr (sorted-map-atoms fep))
          do (format t "fep ~a   moveable-atoms -> ~a~%" fep moveable-atoms)
          do (chem:superpose-one molecule moveable-atoms fixed-atoms)
          do (anchor-to-pose moveable-atoms fixed-atoms :stereochemical-restraints stereochemical-restraints)))
  (minimize-ligands calculation))



(defun build-initial-jobs (calculation &key (connections 3) stages windows)
  (let ((jobs (make-instance 'job-graph)))
    (let* ((unsorted-feps (copy-list (ligands calculation)))
           (sorted-feps (sort unsorted-feps #'< :key (lambda (x) (chem:number-of-atoms (drawing x))))))
      (loop for fep in sorted-feps
            for added = (nodes jobs)
            do (push fep (nodes jobs))
            do (loop for other in (subseq added 0 (min (length added) connections))
                     for stages-lambda = (append (when stages
                                                   (list :stages stages))
                                                 (when windows
                                                   (list :windows windows)))
                     ;; Force a direction for the calculation AA->AB and not AB->AA
                     do (if (string< (string (name fep)) (string (name other)))
                            (progn
                              (format t "Morph from source: ~s to target: ~s~%" fep other)
                              (push (apply #'make-instance 'fep-morph :source fep :target other stages-lambda) (morphs jobs)))
                            (progn
                              (format t "Morph from source: ~s to target: ~s~%" other fep)
                              (push (apply #'make-instance 'fep-morph :source other :target fep stages-lambda) (morphs jobs)))))))
    (setf (jobs calculation) jobs)))

#+(or)
(defun am1-file-name (fep type)
  (make-pathname :directory '(:relative "am1bcc") :name (format nil "am1-~a" (string (name fep))) :type type :defaults *default-pathname-defaults*))

(defun am1-calculation-complete (fep)
  "Return T if the AM1 charges for this FEP exist, are up to date and are available.
Otherwise return NIL."
  (let ((in-filename (am1-file-name fep "in"))
        (out-filename (am1-file-name fep "out"))
        (atom-order (atom-order fep)))
    (cond
      ((and (probe-file out-filename)
            (< (file-write-date in-filename) (file-write-date out-filename))
            (charges:read-am1-charges out-filename atom-order)))
      (t nil))))

(defun setup-am1-calculations (jupyter-job calculation &key (maxcyc 9999))
  (with-top-directory (calculation)
    (let ((ligands (ligands calculation)))
      (ensure-directories-exist *default-pathname-defaults*)
      (format t "Creating am1 scripts in ~a~%" (truename *default-pathname-defaults*))
      (let* (outputs
             (result (loop for ligand in ligands
                           for in-file = (make-instance 'sqm-input-file :name (name ligand))
                           for order-file = (make-instance 'sqm-atom-order-file :name (name ligand))
                           for atom-order = (let (atom-order)
                                              (write-file-if-it-has-changed
                                               (node-pathname in-file)
                                               (with-output-to-string (sout) 
                                                 (setf atom-order (charges:write-sqm-calculation sout (molecule ligand) :maxcyc maxcyc))))
                                              atom-order)
                           for job = (make-ligand-sqm-step ligand :sqm-input-file in-file)
                           do (setf (atom-order ligand) atom-order)
                           do (push (make-instance 'argument :option :output :node in-file) outputs)
                           collect job)))
        (format t "About to set outputs result -> ~s~%" result)
        (setf (outputs jupyter-job) outputs)
        (connect-graph jupyter-job)
        (mapc #'connect-graph result)
        result))))

(defun check-am1-calculations (calculation)
  (warn "This doesn't work right now - we need source tracking to identify the form that is the problem")
  #+(or)
  (let* ((status-entries (loop for fep in (ligands calculation)
                               for name = (string (name fep))
                               for status = (am1-calculation-complete fep)
                               collect (list name status)))
         (sorted-entries (sort status-entries #'string< :key #'car)))
    (error "Fix checking the am1-calculations")
    #+(or)(cl-jupyter-user:html (cl-markup:markup
                           (:table :border 3
                                   :cellpadding 4
                                   (loop for i below (length sorted-entries) by 5
                                         collect (cl-markup:markup
                                                  (:tr :align "right"
                                                       (loop for j from i below (+ i 5)
                                                             for entry = (elt sorted-entries j)
                                                             for name = (car entry)
                                                             for status = (cadr entry)
                                                             when entry
                                                               collect (cl-markup:markup
                                                                        (:td :bgcolor (if status "white" "red")
                                                                             (format nil "~a ~a" name (if status "Done" "Incomplete")))))))))))))



(defun read-am1-charges (calculation)
  (with-top-directory (calculation)
    (let ((count 0))
      (loop for ligand in (ligands calculation)
            for sqm-out = (make-instance 'sqm-output-file :name (name ligand))
            do (let ((am1-charge (charges::read-am1-charges (node-pathname sqm-out) (fep::atom-order ligand))))
                 (setf (am1-charges ligand) am1-charge))
            do (format t "file ~a done~%" (name ligand))
            do (incf count)))))

(defun balance-charges (atom-charge-hash-table net-charge)
  (let ((charge-sum 0.0)
        (num-charges 0))
    (maphash (lambda (atm charge)
               (format t "Charge: ~s ~f~%" atm charge)
               (incf charge-sum charge)
               (incf num-charges))
             atom-charge-hash-table)
    (format t "Total charge: ~f  number of charges: ~d  expected net-charge: ~f~%" charge-sum num-charges net-charge)
    (let ((charge-part (/ (- charge-sum net-charge) num-charges))
          (new-charge-total 0.0))
      (format t "Charge correction: ~f~%" charge-part)
      (maphash (lambda (atm charge)
                 (let ((new-charge (- charge charge-part)))
                   (format t "Updating charge for atom ~s from ~f to ~f~%" atm charge new-charge)
                   (setf (gethash atm atom-charge-hash-table) new-charge)
                   (incf new-charge-total new-charge)))
               atom-charge-hash-table)
      (format t "After adding net-charge the total charge is ~f~%" new-charge-total)
      (when (> (abs new-charge-total) 1.0e-5)
        (warn "After balance-charges the new-charge-total of the molecule is ~g~%" new-charge-total)))))

(defun calculate-am1-bcc-charges (calculation)
  (loop for fep in (ligands calculation)
        for am1-charges = (am1-charges fep)
        for bcc-corrections = (charges::calculate-bcc-corrections (fep::molecule fep))
        for am1-bcc-charges = (charges::combine-am1-bcc-charges (am1-charges fep) bcc-corrections)
        do (balance-charges am1-bcc-charges (net-charge fep))
           (setf (am1-bcc-charges fep) am1-bcc-charges)
           (maphash (lambda (atm charge)
                      (chem:set-charge atm charge))
                    am1-bcc-charges)))


(defun complex (calculation ligand-name &optional (receptor-index 0))
  (let ((receptor (elt (receptors calculation) receptor-index))
        (ligand (find ligand-name (ligands calculation) :test #'eq :key #'name)))
    (let ((agg (chem:matter-copy receptor)))
      (cando:remove-molecules agg 'cando:ligand)
      (chem:add-molecule agg (fep:molecule ligand))
      agg)))


(defun relax-receptor-around-ligand (receptor-aggregate ligand-molecule &key (weight 100.0))
  "Relax the receptor around the ligand."
  (cando:do-atoms (atom receptor-aggregate)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor)))
  (chem:add-matter receptor-aggregate ligand-molecule)
  (cando:do-atoms (atom ligand-molecule)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor))))

(defun pattern-atoms (smarts structure)
  "Given a smarts pattern and a structure, return a hash-table of tags to atoms"
  (cando:do-atoms (atom structure)
    (let (match)
      (when (setf match (chem:matches smarts atom))
        (return-from pattern-atoms (chem:tags-as-hashtable match))))))

(defclass ti-mask ()
  ((source :initarg :source :accessor source)
   (source-timask-residue-index :initarg :source-timask-residue-index :accessor source-timask-residue-index)
   (source-scmask-atom-names :initarg :source-scmask-atom-names :accessor source-scmask-atom-names)
   ;;   (source-timask :initarg :source-timask :accessor source-timask)
   ;;   (source-scmask :initarg :source-scmask :accessor source-scmask)
   (target :initarg :target :accessor target)
   (target-timask-residue-index :initarg :target-timask-residue-index :accessor target-timask-residue-index)
   (target-scmask-atom-names :initarg :target-scmask-atom-names :accessor target-scmask-atom-names)
   ;;   (target-timask :initarg :target-timask :accessor target-timask)
   ;;   (target-scmask :initarg :target-scmask :accessor target-scmask)
   ))


(defgeneric calculate-masks (morph kind)
  (:documentation "Use the source and target to calculate the masks. The kind argument allows the 
user to define new ways to calculate masks"))


(defun unique-residue-with-name (molecule name)
  (let (residues)
    (cando:do-residues (res molecule)
      (when (eq (chem:get-name res) name)
        (push res residues)))
    (if (= (length residues) 1)
        (first residues)
        (error "Only one residue in ~s can have the name ~s - found ~s" molecule name residues))))
    

(defmethod calculate-masks (morph (method (eql ':default)))
  "Really simple mask - the core and mutated residues are timask and the mutated residues are the scmask.
METHOD controls how the masks are calculated"
  (let ((source (source morph))
        (target (target morph)))
    (let ((source-core-residue (unique-residue-with-name (molecule source) (core-residue-name source)))
          (target-core-residue (unique-residue-with-name (molecule target) (core-residue-name target)))
          (source-softcore-atoms (side-chain-atoms source))
          (target-softcore-atoms (side-chain-atoms target))
          (combined-aggregate (cando:combine (molecule source) (molecule target))))
      (let ((sequenced-residues (chem:map-residues 'vector #'identity combined-aggregate))) ; vector of residues in order they will appear in topology file
        (let ((source-timask-residue-index (1+ (position source-core-residue sequenced-residues)))
              (source-scmask-atom-names (loop for atom in source-softcore-atoms
                                              collect (chem:get-name atom)))
              (target-timask-residue-index (1+ (position target-core-residue sequenced-residues)))
              (target-scmask-atom-names (loop for atom in target-softcore-atoms
                                              collect (chem:get-name atom))))
          (make-instance 'ti-mask
                         :source source
			 :source-timask-residue-index source-timask-residue-index
			 :source-scmask-atom-names source-scmask-atom-names
                         ;; :source-timask (format nil ":~d" source-timask-residue-index)
                         ;; :source-scmask (format nil ":~d@~{~d~^,~}" source-timask-residue-index source-scmask-atom-names)
                         :target target
			 :target-timask-residue-index target-timask-residue-index
			 :target-scmask-atom-names target-scmask-atom-names
                         ;; :target-timask (format nil ":~d" target-timask-residue-index)
                         ;; :target-scmask (format nil ":~d@~{~d~^,~}" target-timask-residue-index target-scmask-atom-names)
			 ))))))

(defun mask-substitutions (mask &optional stage)
  (unless (typep stage '(member nil :vdw-bonded :decharge :recharge))
    (error "stage ~s must be one of nil :vdw-bonded :decharge :recharge" stage))
  (list
   (cons :%TIMASK1% (format nil ":~d" (source-timask-residue-index mask)))
   (cons :%SCMASK1% (format nil ":~d@~{~d~^,~}" (source-timask-residue-index mask) (source-scmask-atom-names mask)))
   (cons :%TIMASK2% (format nil ":~d" (target-timask-residue-index mask)))
   (cons :%SCMASK2% (format nil ":~d@~{~d~^,~}" (target-timask-residue-index mask) (target-scmask-atom-names mask)))))

(defun crgmask-substitutions (mask stage)
  (unless (typep stage '(member nil :vdw-bonded :decharge :recharge))
    (error "stage ~s must be one of nil :vdw-bonded :decharge :recharge" stage))
  (cond
    ((eq stage :decharge)
     (List (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}" (target-timask-residue-index mask) (source-scmask-atom-names mask)))
           (cons :%IFSC% "0")))
    ((eq stage :vdw-bonded)
     (list (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}|:~d@~{~d~^,~}"
				    (source-timask-residue-index mask) (source-scmask-atom-names mask)
				    (target-timask-residue-index mask) (target-scmask-atom-names mask)))
           (cons :%IFSC% "1")))
    ((eq stage :recharge)
     (list (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}" (source-timask-residue-index mask) (target-scmask-atom-names mask)))
           (cons :%IFSC% "0")))
    (t nil)))



(defun average-core-atom-positions (source target)
  (loop for source-atom in (core-atoms source)
        for target-atom = (or (find (chem:get-name source-atom)
                                    (core-atoms target)
                                    :test (lambda (src-atom-name tgt-atom) (eq src-atom-name
                                                                               (chem:get-name tgt-atom))))
                              (error "Could not find a matching target atom with the name ~s" (chem:get-name source-atom)))
        for average-pos = (geom:v* (geom:v+ (chem:get-position source-atom)
                                            (chem:get-position target-atom))
                                   0.5)
        for source-delta = (geom:vlength (geom:v- (chem:get-position source-atom)
                                                  average-pos))
        for target-delta = (geom:vlength (geom:v- (chem:get-position target-atom)
                                                  average-pos))
        do (when (or (> source-delta 0.2) (> target-delta 0.2))
             (warn "The source atom ~s and target atom ~s are more than 0.2 A away from the average ~s~%" source-atom target-atom average-pos))
           (chem:set-position source-atom average-pos)
           (chem:set-position target-atom average-pos)))


(defun load-chem-draw-fep (filename)
  (handler-bind ((warning #'muffle-warning))
    (with-open-file (fin (open filename :direction :input))
      (chem:make-chem-draw fin nil nil))))


(defun scale-molecule (mol)
  (let ((total-length 0.0)
        (num-bonds 0))
    (loop for ligand in (ligands calculation)
          for mol = (drawing ligand)
          do (chem:map-bonds nil (lambda (a1 a2 bo)
                                   (let ((length (geom:vlength (geom:v- (chem:get-position a1)
                                                                        (chem:get-position a2)))))
                                     (when (> length 0.1)
                                       (incf num-bonds)
                                       (incf total-length length))))
                             mol))
    (let ((scale (/ 1.5 (/ total-length num-bonds))))
      (loop for ligand in (ligands calculation)
            for mol = (drawing ligand)
            do (cando:do-atoms (atm mol)
                 (chem:set-position atm (geom:v* (chem:get-position atm) scale)))))))
