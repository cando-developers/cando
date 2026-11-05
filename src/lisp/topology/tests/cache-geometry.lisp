(in-package #:topology-shape-key-tests)

(defun test-cp-amplitude-boundaries ()
  ;; Regression for the C++ amplitude extraction bug near +/-pi: computing
  ;; amplitude by dividing sine components loses precision at the half turn.
  (dolist (phase (list 0d0 1d-16 (- 1d-16) pi (- pi) (- pi 1d-12) (+ (- pi) 1d-12)))
    (let ((zjs (chem:cpring-5-cremer-pople-to-zjs 0.2d0 phase)))
      (multiple-value-bind (q phi) (chem:cpring-5-zjs-to-cremer-pople zjs)
        (assert (< (abs (- q 0.2d0)) 1d-12) ()
                "CP amplitude at phase ~s: expected 0.2, measured ~s" phase q)
        (assert (< (abs (atan (sin (- phi phase)) (cos (- phi phase)))) 1d-12)))))
  t)

(defun test-ring-cache-reconstruction (joints)
  (dolist (size '(5 6))
    (let* ((cp-size (if (= size 5) 2 3))
           (driver (make-instance (if (= size 5) 'topology:ring5-driver 'topology:ring6-driver)
                                 :ring-joints (subseq joints 0 size) :exocyclic-joints nil
                                 :cremer-pople-range
                                 (make-instance 'topology::range :begin 0 :end cp-size)
                                 :lengths-range
                                 (make-instance 'topology::range :begin cp-size :end (+ cp-size size))
                                 :angles-range
                                 (make-instance 'topology::range :begin (+ cp-size size)
                                                                :end (+ cp-size (* 2 size)))))
           (declaration (first (topology:parse-cache-declarations '((:ring :pucker :anchor)))))
           (component (topology::resolve-cache-component declaration :test-monomer
                                                         (second joints) (list driver)))
           (policy (make-instance 'ten-degree-test-foldamer))
           (cpinternals (make-array (+ cp-size (* 2 size)) :element-type 'double-float))
           (externals (make-array (+ 3 (reduce #'max joints :key #'kin:joint/position-index-x3))
                                  :element-type 'double-float :initial-element 0d0)))
      ;; A nonplanar regular polygon supplies finite bond geometry. The selector
      ;; is deliberately not the first ring joint.
      (loop for joint in (topology::ring-joints driver)
            for index from 0
            for phase = (* 2d0 pi (/ index size))
            for offset = (kin:joint/position-index-x3 joint)
            do (setf (aref externals offset) (* 1.4d0 (cos phase))
                     (aref externals (+ offset 1)) (* 1.4d0 (sin phase))
                     (aref externals (+ offset 2)) (* 0.15d0 (cos (+ 0.4d0 (* 2 phase))))))
      (topology::driver-extract-cpinternals-from-coordinates driver cpinternals externals)
      (let* ((cache (topology:measure-binned-shape-key-cache policy (list component) externals))
             (bins (topology:binned-shape-key-component cache :pucker))
             (supporting-geometry (subseq cpinternals cp-size)))
        (topology:write-ring-cache-component policy component bins cpinternals)
        (assert (equalp supporting-geometry (subseq cpinternals cp-size)))
        (topology:build-only-ring-externals-from-cpinternals nil driver nil externals cpinternals)
        (let ((reconstructed (topology:measure-binned-shape-key-cache
                              policy (list component) externals)))
          (assert (equal bins (topology:binned-shape-key-component reconstructed :pucker)) ()
                  "Ring ~d reconstruction: expected bins ~s, measured ~s; parameters ~s"
                  size bins (topology:binned-shape-key-component reconstructed :pucker)
                  (topology:measure-cache-component component externals)))))))

(defun test-xyz-cache-dihedral (residue atom-table)
  "Amber's Cartesian joints can provide torsions without being internal drivers."
  (let* ((joints (loop for index below 4
                       collect (kin:make-xyz-joint
                                (list 0 0 index)
                                (chem:get-name (chem:content-at residue index)) atom-table)))
         (declaration (first (topology:parse-cache-declarations '((:dihedral :phi :C)))))
         (externals (make-array
                     (+ 3 (reduce #'max joints :key #'kin:joint/position-index-x3))
                     :element-type 'double-float :initial-element 0d0)))
    (loop for (parent child) on joints while child do (kin:joint/add-child parent child))
    (dolist (joint (subseq joints 0 3))
      (assert (signals-error-p
               (lambda () (topology::resolve-cache-component declaration :amber joint nil)))))
    (loop for joint in joints
          for point in '((0d0 1d0 0d0) (0d0 0d0 0d0) (1d0 0d0 0d0) (1d0 0d0 1d0))
          do (loop for value in point
                   for offset from (kin:joint/position-index-x3 joint)
                   do (setf (aref externals offset) value)))
    (let* ((component (topology::resolve-cache-component declaration :amber (fourth joints) nil))
           (cache (topology:measure-binned-shape-key-cache
                   (make-instance 'ten-degree-test-foldamer) (list component) externals)))
      (assert (member (first (topology:binned-shape-key-component cache :phi)) '(9 27))))))

(defun test-cache-geometry ()
  (let* ((aggregate (chem:load-mol2
                     (asdf:system-relative-pathname
                      "topology" #p"../regression-tests/data/hexapeptide.mol2")))
         (molecule (chem:content-at aggregate 0))
         (residue (chem:content-at molecule 0))
         (atom-table (chem:make-atom-table)))
    (chem:construct-from-molecule atom-table molecule nil nil (make-hash-table :test #'eq))
    (test-xyz-cache-dihedral residue atom-table)
    (let* ((joints (loop for index below 6
                         for name = (chem:get-name (chem:content-at residue index))
                         collect (if (zerop index)
                                     (kin:make-jump-joint (list 0 0 index) name atom-table)
                                     (kin:make-bonded-joint (list 0 0 index) name atom-table))))
           (ring5 (make-instance 'topology:ring5-driver :ring-joints (subseq joints 0 5)))
           (ring6 (make-instance 'topology:ring6-driver :ring-joints joints))
           (declarations (topology:parse-cache-declarations
                          '((:ring :pucker :anchor) (:dihedral :bb1 :angle))))
           (ring-declaration (first declarations))
           (angle-declaration (second declarations)))
      (loop for (parent child) on joints while child do (kin:joint/add-child parent child))
      (test-ring-cache-reconstruction joints)
      ;; Choosing a different unique anchor must not change the CP coordinate order.
      (dolist (ring (list ring5 ring6))
        (dolist (anchor (topology::ring-joints ring))
          (let ((component (topology::resolve-cache-component
                            ring-declaration :test-monomer anchor (list ring))))
            (assert (eq ring (topology:resolved-cache-driver component)))
            (assert (eq (first joints) (first (topology::ring-joints
                                              (topology:resolved-cache-driver component))))))))
      ;; Shared anchors are ambiguous; a unique anchor selects the correct ring
      ;; regardless of driver ordering or the cache field's name.
      (assert (signals-error-p
               (lambda () (topology::resolve-cache-component
                           ring-declaration :test-monomer (first joints) (list ring5 ring6)))))
      (dolist (drivers (list (list ring5 ring6) (list ring6 ring5)))
        (let ((component (topology::resolve-cache-component
                          ring-declaration :test-monomer (sixth joints) drivers)))
          (assert (eq :ring6 (topology:resolved-cache-kind component)))
          (assert (eq ring6 (topology:resolved-cache-driver component)))))
      (dolist (selector (list nil (sixth joints)))
        (assert (signals-error-p
                 (lambda () (topology::resolve-cache-component
                             ring-declaration :test-monomer selector (list ring5))))))
      ;; Reject roots and shallow joints before reaching the coordinate extractor.
      (dolist (joint (subseq joints 0 3))
        (assert (signals-error-p
                 (lambda () (topology::resolve-cache-component
                             angle-declaration :test-monomer joint nil)))))
      (let* ((component (topology::resolve-cache-component
                         angle-declaration :test-monomer (fourth joints) nil))
             (externals (make-array
                         (+ 3 (reduce #'max joints :key #'kin:joint/position-index-x3))
                         :element-type 'double-float :initial-element 0.0d0)))
        ;; Four points with a right-angle torsion; either sign convention bins
        ;; to a quarter turn. Test measurement through the public cache API.
        (loop for joint in joints
              for point in '((0d0 1d0 0d0) (0d0 0d0 0d0)
                             (1d0 0d0 0d0) (1d0 0d0 1d0))
              do (loop for value in point
                       for offset from (kin:joint/position-index-x3 joint)
                       do (setf (aref externals offset) value)))
        (let* ((cache (topology:measure-binned-shape-key-cache
                       (make-instance 'ten-degree-test-foldamer) (list component) externals))
               (bins (topology:binned-shape-key-component cache :bb1)))
          (assert (member (first bins) '(9 27)))
          (assert (equal '(:bb1) (topology:binned-shape-key-names cache)))))))
  ;; Decoding must touch only CP parameters, preserving the bond geometry
  ;; that reconstruction obtains from the supporting backbone rotamer.
  (dolist (kind '(:ring5 :ring6))
    (let* ((size (ecase kind (:ring5 2) (:ring6 3)))
           (driver (make-instance (ecase kind (:ring5 'topology:ring5-driver)
                                              (:ring6 'topology:ring6-driver))
                                  :cremer-pople-range
                                  (make-instance 'topology::range :begin 2 :end (+ 2 size))))
           (component (make-instance 'topology:resolved-ring-cache-component
                                     :kind kind :driver driver))
           (bins (ecase kind (:ring5 '(2 35)) (:ring6 '(3 9 35))))
           (cpinternals (make-array 8 :element-type 'double-float :initial-element 17d0))
           (policy (make-instance 'ten-degree-test-foldamer)))
      (topology:write-ring-cache-component policy component bins cpinternals)
      (assert (equal bins (topology:bin-shape-component
                          policy kind (coerce (subseq cpinternals 2 (+ 2 size)) 'list))))
      (assert (every (lambda (value) (= 17d0 value)) (subseq cpinternals 0 2)))
      (assert (every (lambda (value) (= 17d0 value)) (subseq cpinternals (+ 2 size))))))
  (format t "CACHE-GEOMETRY-TESTS-PASSED~%")
  t)
