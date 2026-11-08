(defpackage #:topology-adjustment-tests
  (:use #:cl)
  (:export #:run-tests))

(in-package #:topology-adjustment-tests)

(defun signals-error-p (thunk)
  (handler-case (progn (funcall thunk) nil) (error () t)))

(defun near (a b)
  (assert (< (abs (- a b)) 1d-5) () "Expected ~s, got ~s" b a))

(defun make-fixture (mode tetrahedral-p)
  "A small joint tree with N bonded to B, H, and a C in another residue.
In :RING mode C is a separate tree branch. Otherwise it is N's child.
The tetrahedral fixture adds H2. No force field or rotamer library is needed."
  (let* ((aggregate (chem:make-aggregate))
         (molecule (chem:make-molecule :test))
         (residue (chem:make-residue :local))
         (remote (chem:make-residue :remote))
         (names (if tetrahedral-p '(:A :B :N :H :H2) '(:A :B :N :H)))
         (atoms (mapcar (lambda (name)
                          (let ((atom (chem:make-atom name
                                        (case name ((:H :H2) :H) (:N :N) (t :C)))))
                            (chem:add-atom residue atom) atom)) names))
         (carbon (chem:make-atom :C :C))
         (atom-table (chem:make-atom-table))
         (local-monomer (make-instance 'topology:monomer :monomers '(:local)))
         (remote-monomer (make-instance 'topology:monomer :monomers '(:remote)))
         (local-topology (make-instance 'topology:topology :name :local))
         (remote-topology (make-instance 'topology:topology :name :remote))
         (plug-name (if (eq mode :incoming) :-amide :+amide))
         (remote-plug-name (case mode (:ring :+ring-amide) (:incoming :+amide) (t :-amide))))
    (chem:add-atom remote carbon)
    (chem:add-matter molecule residue)
    (chem:add-matter molecule remote)
    (chem:add-molecule aggregate molecule)
    (chem:bond-to (first atoms) (second atoms) :single-bond)
    (chem:bond-to (second atoms) (third atoms) :single-bond)
    (dolist (atom (cdddr atoms)) (chem:bond-to (third atoms) atom :single-bond))
    (chem:bond-to (third atoms) carbon :single-bond)
    (when (eq mode :ring) (chem:bond-to (first atoms) carbon :single-bond))
    (chem:construct-from-molecule atom-table molecule nil nil (make-hash-table :test #'eq))
    (setf (gethash plug-name (topology:plugs local-topology))
          (topology::make-out-plug plug-name :N :single-bond)
          (gethash remote-plug-name (topology:plugs remote-topology))
          (topology::make-out-plug remote-plug-name :C :single-bond))
    (topology::register-coupling-with-monomers
     (case mode
       (:ring (make-instance 'topology:ring-coupling
                :monomer1 local-monomer :plug1 plug-name
                :monomer2 remote-monomer :plug2 remote-plug-name))
       (:incoming (make-instance 'topology:directional-coupling
                    :source-monomer remote-monomer :source-plug-name remote-plug-name
                    :target-monomer local-monomer :target-plug-name plug-name))
       (t (make-instance 'topology:directional-coupling
            :source-monomer local-monomer :source-plug-name plug-name
            :target-monomer remote-monomer :target-plug-name remote-plug-name))))
    (let* ((joints (coerce (loop for name in names for i from 0
                                collect (if (< i 3)
                                            (kin:make-xyz-joint (list 0 0 i) name atom-table)
                                            (kin:make-bonded-joint (list 0 0 i) name atom-table)))
                           'vector))
           (other (if (eq mode :ring)
                      (kin:make-xyz-joint '(0 1 0) :C atom-table)
                      (kin:make-bonded-joint '(0 1 0) :C atom-table)))
           (root (aref joints 0))
           (center (aref joints 2))
           (target (aref joints 3))
           (positions (make-hash-table :test #'eq))
           (assembler
             (make-instance 'topology::assembler-base :aggregate aggregate
               :monomer-positions positions
               :ataggregate
               (make-instance 'topology:ataggregate :atmolecules
                 (vector (make-instance 'topology::atmolecule :atresidues
                           (vector (make-instance 'topology::atresidue
                                     :topology local-topology :joints joints)
                                   (make-instance 'topology::atresidue
                                     :topology remote-topology :joints (vector other))))))))
           (internals (topology:make-internals-for-assembler assembler))
           (coords (topology:make-coordinates-for-assembler assembler)))
      (setf (gethash local-monomer positions)
            (make-instance 'topology::monomer-position :molecule-index 0 :residue-index 0)
            (gethash remote-monomer positions)
            (make-instance 'topology::monomer-position :molecule-index 0 :residue-index 1))
      (kin:joint/add-child root (aref joints 1))
      (kin:joint/add-child (aref joints 1) center)
      (loop for i from 3 below (length joints) do (kin:joint/add-child center (aref joints i)))
      (kin:joint/add-child (if (eq mode :ring) root center) other)
      (loop for joint across joints
            for point in '((0d0 1d0 0d0) (0d0 0d0 0d0) (1d0 0d0 0d0)
                           (1d0 0d0 1d0) (1d0 0d0 -1d0))
            do (loop for value in point for i from (kin:joint/position-index-x3 joint)
                     do (setf (aref coords i) value)))
      (loop for value in (list 1.5d0 (/ (sqrt 3d0) 2d0) 0d0)
            for i from (kin:joint/position-index-x3 other)
            do (setf (aref coords i) value))
      (kin:update-internal-coords root internals coords)
      (values assembler internals coords root target other plug-name))))

(defun test-placement (mode tetrahedral-p)
  (multiple-value-bind (assembler internals coords root target other plug-name)
      (make-fixture mode tetrahedral-p)
    (let* ((spec (if tetrahedral-p
                     `(:tetrahedral-hydrogens :other-hydrogen :H2
                       :references (:B (:plug ,plug-name)) :orientation :positive)
                     `(:planar :references (:B (:plug ,plug-name)))))
           (adjustment (topology::instantiate-joint-adjustment spec target))
           (center (kin:parent target))
           (ic (kin:joint/position-index-x3 center))
           (ih (kin:joint/position-index-x3 target)))
      (topology:initialize-adjustment adjustment assembler)
      (assert (eq other (second (topology::reference-joints adjustment))))
      (when (eq mode :ring) (assert (not (eq center (kin:parent other)))))
      (topology:internal-adjust adjustment assembler internals)
      (kin:update-xyz-coords root internals coords)
      (near (geom:calculate-distance-array ic ih coords) 1d0)
      (if tetrahedral-p
          (progn
            ;; B x C points down; positive orientation puts H on that side.
            (assert (< (aref coords (+ ih 2)) (aref coords (+ ic 2))))
            (near (geom:calculate-angle-array
                   ih ic (kin:joint/position-index-x3
                          (second (topology::hydrogen-pair adjustment))) coords)
                  (acos (/ -1d0 3d0))))
          (progn
            (near (aref coords (+ ih 2)) (aref coords (+ ic 2)))
            (dolist (reference (topology::reference-joints adjustment))
              (near (geom:calculate-angle-array ih ic (kin:joint/position-index-x3 reference) coords)
                    (* 2/3 pi)))
            ;; Reversing reference order gives exactly the same placement.
            (setf (topology::reference-joints adjustment)
                  (reverse (topology::reference-joints adjustment)))
            (let ((before (copy-seq coords)))
              (topology:internal-adjust adjustment assembler internals)
              (kin:update-xyz-coords root internals coords)
              (map nil #'near coords before))
            (dolist (references `((:B :H) (:B :B) (:B (:plug :+missing))))
              (assert (signals-error-p
                       (lambda ()
                         (topology:initialize-adjustment
                          (topology::instantiate-joint-adjustment
                           `(:planar :references ,references) target) assembler))))))))))

(defun run-tests ()
  (topology::test-new-adjustments)
  (assert (eq :external (nth-value 1 (find-symbol "WALK-ATAGGREGATE-JOINTS" :topology))))
  (dolist (mode '(:outgoing :incoming :ring))
    (test-placement mode nil)
    (test-placement mode t))
  (format t "~&ADJUSTMENT-PLACEMENT-TESTS-PASSED~%")
  t)
