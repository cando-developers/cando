(in-package #:clasp-tests)


(test-true sign-dihedral (> (geom:calculate-dihedral (geom:vec -1.0 1.0 0.0) (geom:vec -1.0 0.0 0.0) (geom:vec 0.0 0.0 0.0) (geom:vec 0.0 0.0 1.0)) 0.0))
(test-true sign-inv-dihedral (< (geom:vz (geom:point-from-stub-and-internal-coordinates (geom:make-matrix t) 1.0 (/ pi 2.0) (/ pi 2.0))) 0.0))
(test-true sign-other-dihedral (> (geom:vz (geom:build-using-bond-angle-dihedral
                                            1.0 (geom:vec 0.0 0.0 0.0)
                                            (/ pi 2.0) (geom:vec -1.0 0.0 0.0)
                                            (/ pi 2.0) (geom:vec -1.0 1.0 0.0)))
                                  0.0))
(test-true sign-plane-vector-angle1 (< (geom:plane-vector-angle  1.0  1.0) 0.0))
(test-true sign-plane-vector-angle2 (< (geom:plane-vector-angle -1.0  1.0) 0.0))
(test-true sign-plane-vector-angle3 (> (geom:plane-vector-angle  1.0 -1.0) 0.0))
(test-true sign-plane-vector-angle4 (> (geom:plane-vector-angle -1.0 -1.0) 0.0))
