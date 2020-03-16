(in-package :sdf)

(defun write-sdf-to-file (aggregate filename)
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (write-sdf-stream aggregate fout)))

(defun write-sdf-stream (aggregate stream)
  (format stream "CHEMBL1093087~%
                    3D~%
Structure written by MMmdl.~%")
  (let ((atom-to-index (make-hash-table)))
    (let atoms (chem:map-atoms 'vector
                               (lambda (atm) atm)
                               aggregate)
    (bonds (chem:map-bonds 'vector
                           (lambda (atm1 atm2 bond-order)
                             (list atm1 atm2 bond-order))
                           aggregate))
    ((coordinates (make-array (* 3 (length atoms)) :element-type 'double-float)))
  (loop for index below (length atoms)
        for atm = (elt atoms index)  ;;get element from vector at position
        do (setf (gethash atm atom-to-index) index) ;;key=atm value=index
                 for pos = (chem:get-position atm)
                 for coord-index = (* 3 index)
                 do (setf (elt coordinates coord-index) (geom:vx pos)
                   (elt coordinates (+ 1 coord-index)) (geom:vy pos)
                   (elt coordinates (+ 2 coord-index)) (geom:vz pos))
                    (format stream "~a"
                            (coordinates))))))
                        
