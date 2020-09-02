;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :tirun)

(defparameter *sketch* "")

(defclass assembly ()
  ((core-group :initarg :core-group :accessor core-group)
   (core-attachment-indicators :initarg :core-attachment-indicators :accessor core-attachment-indicators)
   (core-attachment-atoms :initarg :core-attachment-atoms :accessor core-attachment-atoms)
   (core-residue :initarg :core-residue :accessor core-residue)
   (core-topology :initarg :core-topology :accessor core-topology)
   (core-out-vectors :initarg :core-out-vectors :accessor core-out-vectors)
   (side-chain-groups :initarg :side-chain-groups :accessor side-chain-groups)
   (side-chain-attachment-indicators :initarg :side-chain-attachment-indicators :accessor side-chain-attachment-indicators)
   (side-chain-attachment-names :initarg :side-chain-attachment-names :accessor side-chain-attachment-names)
   (topology-to-residue :initarg :topology-to-residue :accessor topology-to-residue)
   (side-chain-topologys :initform (make-hash-table) :initarg :side-chain-topologys :accessor side-chain-topologys)
   (map-atoms-numbers :initarg :map-atoms-numbers :accessor map-atoms-numbers)
   ))

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
         parser.common-rules:float-literal
         parser.common-rules:float-literal
         parser.common-rules:whitespace
         label/s
         number/s number/s number/s number/s number/s number/s number/s number/s number/s number/s number/s number
         )
  (:lambda (x)
    x))

(esrap:defrule bond-line
    (and number/s number/s number/s number/s number/s number/s number)
  (:lambda (x)
    x))

(esrap:defrule header-line
    (and number/s number/s number/s number/s number/s number/s number/s number/s number/s "0999 V2000"))

#|
(esrap:parse 'header-line "8  8  0  0  0  0  0  0  0  0999 V2000")
(esrap:parse 'atom-line "  0.0000    4.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0")
|#


(defun build-atom-vector (atoms bonds)
  (let ((atom-vector (make-array (length atoms)))
        (map-numbers (make-hash-table)))
    (loop for atom-info in atoms
          for index from 0
          for xpos = (elt atom-info 0)
          for ypos = (elt atom-info 1)
          for zpos = (elt atom-info 2)
          for raw-map-number = (elt atom-info 14)
          for map-number = (if (= raw-map-number 0) nil raw-map-number)
          for label-string = (string (elt atom-info 4))
          for match-element-string = (progn
                                       (find label-string *elements* :test #'string=))
          for element-string = (if match-element-string match-element-string "C")
          for element = (intern element-string :keyword)
          for label = (if match-element-string nil (intern label-string :keyword))
          for atom-name = (intern (format nil "~a~a" element-string index) :keyword)
          for atom = (chem:make-atom atom-name element)
          do (chem:set-position atom (geom:vec xpos ypos zpos))
          do (when label (chem:set-property atom :label label))
          do (when map-number (setf (gethash atom-name map-numbers) map-number))
          do (setf (elt atom-vector index) atom))
    (loop for bond-info in bonds
          for i1 = (1- (elt bond-info 0))
          for i2 = (1- (elt bond-info 1))
          for order-index = (elt bond-info 2)
          for order = (elt '(:no-bond :single-bond :double-bond :triple-bond) order-index)
          for atom1 = (elt atom-vector i1)
          for atom2 = (elt atom-vector i2)
          do (chem:bond-to atom1 atom2 order))
    (values atom-vector map-numbers)))

(defun parse-groups (mol-string)
  (with-input-from-string (sin mol-string)
    (let* ((smiles (read-line sin))
           (label (read-line sin))
           (blank (read-line sin))
           (header-line (string-trim '(#\space) (read-line sin)))
           (header (esrap:parse 'header-line header-line))
           (num-atoms (first header))
           (num-bonds (second header))
           (atoms (loop for atom-index below num-atoms
                        for line = (string-trim '(#\space) (read-line sin))
                        collect (esrap:parse 'atom-line line)))
           (bonds (loop for bond-index below num-bonds
                        for line = (string-trim '(#\space) (read-line sin))
                        collect (esrap:parse 'bond-line line))))
      (multiple-value-bind (atom-vector map-names-numbers)
          (build-atom-vector atoms bonds)
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
          (values groups map-names-numbers))))))



(defun build-attachment-residue (atoms attachment-indicators &key name direction)
  (let ((residue (chem:make-residue name))
        (connection-vectors (make-hash-table))
        (attachment-atoms (make-hash-table))
        (indicator-atoms (make-hash-table)))
    (maphash (lambda (label atom)
               (let* ((attachment-atom (let ((bonded-atoms (chem:bonded-atoms-as-list atom)))
                                         (unless (= (length bonded-atoms) 1)
                                           (error "There are too many bonded atoms to ~a" atom))
                                         (first bonded-atoms)))
                      (short-label (intern (subseq (string label) 0 2) :keyword))
                      (start (chem:get-position atom))
                      (dir (geom:v- (chem:get-position attachment-atom) start))
                      (bond-vector (if (eq direction :out)
                                       (let* ((start (chem:get-position attachment-atom))
                                              (end (chem:get-position atom))
                                              (dir (geom:v- end start)))
                                         (make-instance 'bond-vector :start-vector start :direction dir))
                                       (let* ((start (chem:get-position atom))
                                              (end (chem:get-position attachment-atom))
                                              (dir (geom:v- end start)))
                                         (make-instance 'bond-vector :start-vector start :direction dir)))))
                 (setf (gethash short-label attachment-atoms) attachment-atom)
                 (setf (gethash short-label connection-vectors) bond-vector)
                 (chem:remove-bond-to atom attachment-atom)
                 (setf (gethash atom indicator-atoms) t)))
             attachment-indicators)
    (loop for atom in atoms
          when (null (gethash atom indicator-atoms))
            do (chem:add-matter residue atom))
    (values residue connection-vectors attachment-atoms)))


(defgeneric parse-residue-groups (sketch))

(defmethod parse-residue-groups ((sketch chem:chem-draw))
  (let* ((agg (chem:as-aggregate sketch)))
    (parse-residue-groups agg)))

(defmethod parse-residue-groups ((agg chem:aggregate))
  (let* (core-molecule
         side-chain-molecules
         (core-atom-count 0)
         (map-atoms-numbers (make-hash-table)))
    (let ((num-bonds 0)
          (sum-lengths 0.0)
          (index 0))
      (chem:map-atoms nil (lambda (a)
                            (setf (gethash a map-atoms-numbers) (incf index)))
                      agg)
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
    ;; Find the core-molecule - it has the property :id :core
    (chem:map-molecules nil (lambda (m)
                              (let ((num-atoms (chem:number-of-atoms m)))
                                (when (eq (chem:matter-get-property-or-default m :id nil) :core)
                                  (setf core-molecule m))))
                        agg)
    ;; Everything else is a side-chain molecules
    (chem:map-molecules nil (lambda (m)
                              (unless (eq m core-molecule)
                                (push m side-chain-molecules)))
                        agg)
    ;; Handle the core-group
    (let ((core-group (chem:map-atoms 'list #'identity core-molecule))
          (core-attachment-indicators (make-hash-table))
          (map-names-numbers (make-hash-table)))
      (mapcar (lambda (a)
                (if (member (chem:get-element a) (list :du :ne))
                    (let ((bonds (chem:bonds-as-list a)))
                      (unless (= (length bonds) 1)
                        (warn "There are the wrong number of bonds from ~a bonds: ~a" a bonds))
                      (setf (gethash (chem:get-name a) core-attachment-indicators) a))))
              core-group)
      (multiple-value-bind (core-residue out-vectors core-attachment-atoms)
          (build-attachment-residue core-group core-attachment-indicators :name :cor :direction :out)
        ;; Handle the side-chains
        (let ((side-chain-groups (make-hash-table))
              (side-chain-attachment-indicators (make-hash-table)))
          (loop for side-chain-molecule in side-chain-molecules
                do (let* (side-chain-name
                          (side-chain-atoms (chem:map-atoms
                                             'list
                                             (lambda (a)
                                               (when (member (chem:get-element a) (list :du :ne))
                                                 (setf side-chain-name (chem:get-name a))
                                                 (format t "side-chain-name: ~s~%" side-chain-name)
                                                 (setf (gethash side-chain-name side-chain-attachment-indicators) a))
                                               a)
                                             side-chain-molecule)))
                     (setf (gethash side-chain-name side-chain-groups) side-chain-atoms)))
          ;; Sort out the side-chains now
          (make-instance 'assembly
                         :map-atoms-numbers map-atoms-numbers
                         :core-group core-group
                         :core-attachment-atoms core-attachment-atoms
                         :core-attachment-indicators core-attachment-indicators
                         :core-residue core-residue
                         :core-out-vectors out-vectors
                         :side-chain-groups side-chain-groups
                         :side-chain-attachment-indicators side-chain-attachment-indicators
                         ))))))

(defun group-name (group)
  (loop for atom in group
        for label = (chem:matter-get-property-or-default atom :label nil)
        when label
          do (return-from group-name label)))

(defvar *core-group* nil)
(defvar *side-chain-groups* nil)
(defvar *map-names-numbers* nil)

(defun build-assembly (sketch &key verbose)
  (let* ((topology-to-residue (make-hash-table))
         (assembly (parse-residue-groups sketch))
         (core-topology (cando:make-simple-topology-from-residue (core-residue assembly))))
    (setf (core-topology assembly) core-topology)
    (maphash (lambda (name vec)
               (let* ((plug-name (intern (format nil "+~a" (string name)) :keyword))
                      (plug-atom (gethash name (core-attachment-atoms assembly)))
                      (plug-atom-name (chem:get-name plug-atom))
                      #+(or)(_ (format t "Making out plug plug-name: ~a    name: ~a~%" plug-name plug-atom-name))
                      (out-plug (chem:make-out-plug plug-name nil nil plug-atom-name :single-bond)))
                 (chem:add-plug core-topology (chem:get-name out-plug) out-plug)))
             (core-out-vectors assembly))
    (let ((side-chain-attachment-names (make-hash-table)))
      (loop for side-chain-name being the hash-keys in (side-chain-groups assembly) using (hash-value side-chain-group)
            with side-chain-residue
            with side-chain-vector
            with attachment-atom
            do (multiple-value-setq (side-chain-residue side-chain-vector attachment-atom)
                 (build-side-chain-residue side-chain-attachment-names
                                           side-chain-name
                                           side-chain-group
                                           (side-chain-attachment-indicators assembly)
                                           (gethash (short-name side-chain-name) (core-out-vectors assembly))
                                           (core-attachment-atoms assembly)))
            do (let* ((side-chain-topology (cando:make-simple-topology-from-residue side-chain-residue))
                      (short-name (short-name (chem:get-name side-chain-residue)))
                      (plug-name (intern (format nil "-~a" (string short-name)) :keyword))
                      (in-plug (chem:make-in-plug plug-name nil (chem:get-name attachment-atom) :single-bond)))
                 (setf (gethash side-chain-topology topology-to-residue) side-chain-residue)
                 (chem:add-plug side-chain-topology plug-name in-plug)
                 (push side-chain-topology (gethash short-name (side-chain-topologys assembly)))))
      (setf (topology-to-residue assembly) topology-to-residue)
      (setf (side-chain-attachment-names assembly) side-chain-attachment-names))
    assembly))
                                   
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


(defun labeled-atom (atom)
  (first (chem:bonded-atoms-as-list atom)))

(defclass bond-vector ()
  ((start-vector :initarg :start-vector :accessor start-vector)
   (direction :initarg :direction :accessor direction)))

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

(defun build-side-chain-residue (side-chain-names side-chain-name atoms side-chain-attachment-indicators core-vector core-attachment-atoms)
  (let* ((residue (chem:make-residue side-chain-name))
         (short-side-chain-name (short-name side-chain-name))
         in-atom-vector
         attach-atom
         (core-atom (gethash short-side-chain-name core-attachment-atoms)))
    (push side-chain-name (gethash short-side-chain-name side-chain-names))
    (loop for atom in atoms
          for label = (chem:get-name atom)
          if (eq label side-chain-name)
            do (let* ((attachment-atom (let ((bonded-neighbors (chem:bonded-atoms-as-list atom)))
                                         (unless (= (length bonded-neighbors) 1)
                                           (error "There are too many bonded neighbors to ~a" atom))
                                         (first bonded-neighbors)))
                      (start-vec (chem:get-position atom))
                      (dir (geom:v- (chem:get-position attachment-atom) start-vec)))
                 (setf in-atom-vector (make-instance 'bond-vector :start-vector start-vec :direction dir))
                 (chem:remove-bond-to atom attachment-atom)
                 (setf attach-atom attachment-atom))
          else
            do (chem:add-matter residue atom))
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
    (values residue in-atom-vector attach-atom)))

(defun build-tiruns (assembly)
  ;;; core-topology side-topologys topology-to-residue map-names-numbers)
  (let ((top-groups nil))
    (maphash (lambda (name tops)
               (push tops top-groups))
             (side-chain-topologys assembly))
    (let ((cartesian (cartesian-product top-groups)))
      (loop for cart in cartesian
            for name = (structure-name cart)
            for molecule = (chem:make-molecule name)
            for core-residue = (let ((res (chem:build-residue-single-name (core-topology assembly))))
                                 (set-coordinates-for-residue res (core-residue assembly))
                                 res)
            collect (let (side-chain-residue-names)
                      (chem:add-matter molecule core-residue)
                      (loop for top in cart
                            for name = (chem:get-name top)
                            for side-residue = (let ((res (chem:build-residue-single-name top)))
                                                 (set-coordinates-for-residue res (gethash top (topology-to-residue assembly)))
                                                 res)
                            for in-plug = (chem:get-in-plug top)
                            for in-plug-name = (chem:get-name in-plug)
                            for out-plug-name = (chem:other-plug-name in-plug-name)
                            for side-atoms = (chem:map-atoms 'list #'identity side-residue)
                            do (progn
                                 (push (chem:get-name side-residue) side-chain-residue-names)
                                 (chem:add-matter molecule side-residue)
                                 #+(or)
                                 (maphash (lambda (name number)
                                            (let ((atom (chem:atom-with-name core-residue name)))
                                              #+(or)(format t "Atom: ~a   number: ~a~%" atom number)
                                              (chem:set-property atom :map number)))
                                          (map-names-numbers assembly))
                                 (chem:connect-residues (core-topology assembly)
                                                        core-residue
                                                        out-plug-name
                                                        top
                                                        side-residue
                                                        in-plug-name)))
                      (set-stereocenters molecule)
                      (make-instance 'tirun-structure :name name :drawing molecule
                                                    :core-residue core-residue
                                                    :core-residue-name (chem:get-name core-residue)
                                                    :side-chain-residue-names side-chain-residue-names))))))

(defun setup-ligands (tiruns sketch &key verbose)
  (let* ((assembly (build-assembly sketch :verbose verbose)))
    (setf (core-topology tiruns) (core-topology assembly))
    (setf (side-topologys tiruns) (side-chain-topologys assembly))
    (let ((ligands (build-tiruns assembly)))
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
      (setf (ligands tiruns) ligands))))

(defun ensure-unique-hydrogen-names (molecule)
  (let ((hydrogen-counter 0))
    (cando:do-atoms (atom1 molecule)
      (when (eq (chem:get-element atom1) :H)
        (let ((hydrogen-name (intern (format nil "H~a" (incf hydrogen-counter)) :keyword)))
          (chem:set-name atom1 hydrogen-name))))))

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
  (let ((names (mapcar (lambda (x) (subseq (string (chem:get-name x)) 1 3)) tops)))
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
                   do (let ((tirun (pop ligands)))
                        (when tirun
                          (let ((molecule (chem:matter-copy (funcall accessor tirun))))
                            (chem:apply-transform-to-atoms molecule transform)
                            (chem:add-matter agg molecule))))))
    agg))

(defun layout-2 (name1 name2 tiruns &key (accessor 'drawing))
  (let ((mol1 (chem:matter-copy (funcall accessor (find name1 tiruns :key 'name))))
        (mol2 (chem:matter-copy (funcall accessor (find name2 tiruns :key 'name))))
        (agg (chem:make-aggregate 'two))
        (translate (geom:make-m4-translate (geom:vec 20.0 0.0 0.0))))
    (chem:add-matter agg mol2)
    (chem:apply-transform-to-atoms mol1 translate)
    (chem:add-matter agg mol1)
    agg))

(defun build-ligands (calculations)
  (let ((tiruns (ligands calculations)))
    (loop for tirun in tiruns
          for mol = (molecule tirun)
          do (format t "Building ~a~%" (name tirun))
             (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) mol)
             (setf (molecule tirun) mol)
          do (format t "build-structures tirun: ~a   map-atoms: ~a~%" tirun (sorted-map-atoms tirun)))))

(defun minimize-ligands (calculation)
  (let ((tiruns (ligands calculation)))
    (loop for tirun in tiruns
          for mol = (molecule tirun)
          for agg = (chem:make-aggregate nil)
          do (chem:add-matter agg mol)
          do (funcall (find-symbol "ASSIGN-ATOM-TYPES" :leap) agg)
          do (energy:minimize agg :cg-tolerance 0.0001 :max-tn-steps 100 :tn-tolerance 0.00001))))

(defun sorted-map-atoms (tirun assembly)
  (let ((molecule (molecule tirun))
        (map-atoms nil))
    (chem:map-atoms nil (lambda (a)
                          (let ((num (gethash a (map-atoms-numbers assembly))))
                            (when num
                              (push (cons num a) map-atoms))))
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
         (tiruns (ligands calculation)))
    (loop for tirun in tiruns
          for molecule = (molecule tirun)
          for moveable-atoms-map = (pattern-atoms pattern molecule)
          for moveable-atoms = (mapcar (lambda (index) (gethash index moveable-atoms-map)) fixed-order)
          do (format t "tirun ~a~%" tirun)
             (format t " ~a moveable-atoms -> ~a~%" tirun moveable-atoms)
             (format t "    fixed-atoms -> ~a~%" fixed-atoms)
          do (chem:superpose-one molecule moveable-atoms fixed-atoms)
          do (anchor-to-pose moveable-atoms fixed-atoms :stereochemical-restraints stereochemical-restraints)))
  (minimize-ligands calculation))



(defun pose-ligands (calculation fixed-atoms &key stereochemical-restraints)
  (let ((tiruns (ligands calculation)))
    (loop for tirun in tiruns
          for molecule = (molecule tirun)
          for moveable-atoms = (mapcar #'cdr (sorted-map-atoms tirun))
          do (format t "tirun ~a   moveable-atoms -> ~a~%" tirun moveable-atoms)
          do (chem:superpose-one molecule moveable-atoms fixed-atoms)
          do (anchor-to-pose moveable-atoms fixed-atoms :stereochemical-restraints stereochemical-restraints)))
  (minimize-ligands calculation))



(defun build-job-nodes (calculation)
  (let ((jobs (make-instance 'job-graph)))
    (let* ((unsorted-tiruns (copy-list (ligands calculation)))
           (sorted-tiruns (sort unsorted-tiruns #'< :key (lambda (x) (chem:number-of-atoms (drawing x))))))
      (loop for tirun in sorted-tiruns
            for added = (nodes jobs)
            do (push tirun (nodes jobs))))
    (setf (jobs calculation) jobs)))

(defgeneric connect-job-nodes (calculation type &key &allow-other-keys))

(defun add-job-edge (source-node target-node calculation)
  (let ((stages-lambda (append (when (ti-stages calculation)
                                 (list :stages (ti-stages calculation)))
                               (when (ti-lambdas calculation)
                                 (list :lambda-values
                                       (if *testing*
                                           *testing-lambdas*
                                           (funcall (ti-lambdas calculation) *system*)))))))
    (push (apply #'make-instance
                 'tirun-morph
                 :source source-node
                 :target target-node
                 :morph-mask (calculate-masks source-node
                                              target-node
                                              (mask-method calculation))
                 stages-lambda)
          (morphs (jobs calculation)))))
  
(defmethod connect-job-nodes (calculation (type (eql :simple)) &key (connections 3) stages windows)
  (let* ((jobs (jobs calculation))
         (nodes (nodes jobs))
         (unsorted-tiruns (copy-list nodes))
         (sorted-tiruns (sort unsorted-tiruns #'< :key (lambda (x) (chem:number-of-atoms (drawing x)))))
         added)
    (loop for tirun in sorted-tiruns
          do (push tirun added)
          do (loop for other in (subseq (cdr added) 0 (min (length (cdr added)) connections))
                   do (add-job-edge tirun other calculation)))))
