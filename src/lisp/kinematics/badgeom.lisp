(in-package :kin.work)

(defun bad-carbon-geometry-p (aggregate)
  ;; Check that sp3 carbons have good geometry
  (labels
      ((calculate-angle-degrees (va vx vb)
         (/ (geom:calculate-angle va vx vb) 0.0174533))
       (bad-sp3-angle-p (angle)
         (or (< angle 90.0) (> angle 120.0)))
       (bad-sp2-angle-p (angle)
         (or (< angle 100.0) (> angle 134.0)))
       (bad-angles-p (aggregate)
         (cando:do-atoms (atm aggregate)
           (let ((vx (chem:get-position atm))
                 (bonds (chem:bonds-as-list atm)))
             (cond
               ((= (length bonds) 4)
                (let* ((aa (chem:get-other-atom (first bonds) atm))
                       (ab (chem:get-other-atom (second bonds) atm))
                       (ac (chem:get-other-atom (third bonds) atm))
                       (ad (chem:get-other-atom (fourth bonds) atm))
                       (va (chem:get-position aa))
                       (vb (chem:get-position ab))
                       (vc (chem:get-position ac))
                       (vd (chem:get-position ad))
                       (axb (calculate-angle-degrees va vx vb))
                       (axc (calculate-angle-degrees va vx vc))
                       (axd (calculate-angle-degrees va vx vd))
                       (bxc (calculate-angle-degrees vb vx vc))
                       (bxd (calculate-angle-degrees vb vx vd))
                       (cxd (calculate-angle-degrees vc vx vd)))
                  #+(or)(progn
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ab) axb)
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ac) axc)
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ad) axd)
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ac) bxc)
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ad) bxd)
                          (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ac) (chem:get-name atm) (chem:get-name ad) cxd))
                  (when (or (bad-sp3-angle-p axb)
                            (bad-sp3-angle-p axc)
                            (bad-sp3-angle-p axd)
                            (bad-sp3-angle-p bxc)
                            (bad-sp3-angle-p bxd)
                            (bad-sp3-angle-p cxd))
                    (return-from bad-angles-p t))
                  ))
               ((= (length bonds) 3)
                (let* ((aa (chem:get-other-atom (first bonds) atm))
                       (ab (chem:get-other-atom (second bonds) atm))
                       (ac (chem:get-other-atom (third bonds) atm))
                       (va (chem:get-position aa))
                       (vb (chem:get-position ab))
                       (vc (chem:get-position ac))
                       (axb (calculate-angle-degrees va vx vb))
                       (axc (calculate-angle-degrees va vx vc))
                       (bxc (calculate-angle-degrees vb vx vc)))
                  #+(or)(progn
                    (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ab) axb)
                    (format t "~a ~a ~a : ~5,3f~%" (chem:get-name aa) (chem:get-name atm) (chem:get-name ac) axc)
                    (format t "~a ~a ~a : ~5,3f~%" (chem:get-name ab) (chem:get-name atm) (chem:get-name ac) bxc))
                  (when (or (bad-sp2-angle-p axb)
                            (bad-sp2-angle-p axc)
                            (bad-sp2-angle-p bxc))
                    (return-from bad-angles-p t))
                  ))               )))))
    (or (bad-angles-p aggregate)
        )))

(defparameter *dkp* (chem:make-chem-info-graph (chem:compile-smarts "[C&H1:0]1-[N:1](-[C:2])-[C:3](=[O:4])-[C&H0:5]-[N:6](-[C&H2:7])-[C:8]1(=[O:9])")))

(defun bad-dkp-geometry-p (agg)
  (let* ((mol (cando:mol agg 0))
         (mol-graph (chem:make-molecule-graph-from-molecule mol))
         (matches (chem:boost-graph-vf2 *dkp* mol-graph)))
    (loop for match in matches
          for c1 = (elt match 2)
          for n1 = (elt match 1)
          for co1 = (elt match 3)
          for oc1 = (elt match 4)
          for c1n = (chem:get-name c1)
          for n1n = (chem:get-name n1)
          for co1n = (chem:get-name co1)
          for oc1n = (chem:get-name oc1)
          for dihed1 = (abs (/ (chem:calculate-dihedral-for-atoms c1 n1 co1 oc1) 0.0174533))
          for c2 = (elt match 7)
          for n2 = (elt match 6)
          for co2 = (elt match 8)
          for oc2 = (elt match 9)
          for c2n = (chem:get-name c2)
          for n2n = (chem:get-name n2)
          for co2n = (chem:get-name co2)
          for oc2n = (chem:get-name oc2)
          for dihed2 = (abs (/ (chem:calculate-dihedral-for-atoms c2 n2 co2 oc2) 0.0174533))
          when (or (> dihed1 45.0)
                   (> dihed2 45.0))
            do (return-from bad-dkp-geometry-p t)))
  nil)

(defparameter *planar-aromatic* (chem:make-chem-info-graph (chem:compile-smarts "[*:0]1:[C&H1:1](~[#1:2]):[C&H1:3](-[#1:4]):[*:5]:[*]:[*]1")))

(defun bad-aromatic-geometry-p (agg)
  (unless chem:*current-rings* (error "The chem:*current-rings* dynamic variable must be defined - use (chem:identify-rings matter)"))
  (aromaticity:with-aromaticity-information (agg :am1bcc)
    (let* ((mol (cando:mol agg 0))
           (mol-graph (chem:make-molecule-graph-from-molecule mol))
           (matches (chem:boost-graph-vf2 *planar-aromatic* mol-graph)))
      (loop for match in matches
            for h1 = (elt match 2)
            for c1 = (elt match 1)
            for c2 = (elt match 3)
            for h2 = (elt match 4)
            for h1n = (chem:get-name h1)
            for c1n = (chem:get-name c1)
            for h2n = (chem:get-name h2)
            for c2n = (chem:get-name c2)
            for ph1 = (chem:get-position h1)
            for pc1 = (chem:get-position c1)
            for pc2 = (chem:get-position c2)
            for ph2 = (chem:get-position h2)
            for torsion = (geom:calculate-dihedral ph1 pc1 pc2 ph2)
            when (> (abs (/ torsion 0.0174533)) 30.0)
              do (return t)
            finally (return nil)))))


(defun bad-geometry-p (agg)
  "Test for different kinds of common bad geometry 
- return true if there was any bad geometry and nil if there wasn't"
  (or (bad-carbon-geometry-p agg)
      (let ((chem:*current-rings* (chem:identify-rings agg)))
        (bad-aromatic-geometry-p agg)
        (bad-dkp-geometry-p agg))))
