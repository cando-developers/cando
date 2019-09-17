
(defun canonicalize-lipids (molecules top-atom-names bottom-atom-names)
  (let (aligned-molecules
        (top (geom:vec 0.0 0.0 0.0))
        (bottom (geom:vec 0.0 0.0 0.0)))
    (loop for cmol in molecules
          for mol = (chem:matter-copy cmol)
          for cog = (chem:geometric-center mol)
          for m4 = (geom:make-m4-translate (geom:v* cog -1.0))
          do (chem:apply-transform-to-atoms mol m4)
          do (push mol aligned-molecules))
    (flet ((average-atom-pos (mol atom-names)
             (loop named top-calc
                   with ttop = (geom:vec 0.0 0.0 0.0)
                   for name in atom-names
                   for atom = (chem:first-atom-with-name mol name)
                   do (setf ttop (geom:v+ ttop (chem:get-position atom)))
                   finally (return-from top-calc (geom:v* ttop (/ 1.0 (length atom-names)))))))
      (loop with zdir = (geom:vec 0.0 0.0 1.0)
            for mol in aligned-molecules
            for top = (average-atom-pos mol top-atom-names)
            for bot = (average-atom-pos mol bottom-atom-names)
            for delta = (geom:vnormalized (geom:v- top bot))
            for dot = (geom:vdot delta zdir)
            for angle = (acos dot)
            for cross = (geom:vcross delta zdir)
            for cross-length = (geom:vlength cross)
            for cross-dir = (if (or (< cross-length 0.035)  ; 2 deg
                                    (> cross-length 3.107)) ; 178deg
                                (geom:vec 1.0 0.0 0.0)
                                (geom:vnormalized cross))
            for transform = (geom:make-m4-rotate-axis (- angle) cross)
            do (chem:apply-transform-to-atoms mol transform)))
    aligned-molecules))


(defun delta (lipids zside)
  (let ((lipid-vec (coerce lipids 'vector)))
    (loop for i1 from 0 below (1- (length lipid-vec))
          for mol1 = (elt lipid-vec i1)
          for i1-center = (chem:geometric-center mol1)
          when (= (float-sign zside) (float-sign (geom:vz i1-center)))
            do (loop for i2 from (1+ i1) below (length lipid-vec)
                     for mol2 = (elt lipid-vec i2)
                     for i2-center = (chem:geometric-center mol2)
                     when (= (float-sign zside) (float-sign (geom:vz i2-center)))
                       do (let ((delta (geom:v- i2-center i1-center)))
                            (when (and (> (geom:vx delta) 0.0)
                                       #+(or)(> (geom:vy delta) 0.0)
                                       (< (geom:vlength delta) 20.0))
                              (format t "~,2f~a~,2f~%" (geom:vx delta) #\tab (geom:vy delta))))))))
                         
                     
