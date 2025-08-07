(in-package :cando)


(defun align-atoms (&key moveable-atoms fixed-atoms moveable-matter)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and apply the transformation to MOVEABLE-MATTER (an aggregate, molecule, or residue)."
  (chem:align-matter moveable-matter moveable-atoms fixed-atoms))

(defun align-transform (&key moveable-atoms fixed-atoms)
  "Align the atoms in MOVEABLE-ATOMS-SEQUENCE to the FIXED-ATOMS-SEQUENCE and return the transform"
  (chem:align-transform moveable-atoms fixed-atoms))


(defun heavy-atoms-index3 (aggregate &key (test t))
  "Return a vector of index3 for heavy-atoms in AGGREGATE that pass the test"
  (let ((heavy-atoms-index3 (make-array 256 :element-type 'ext:byte32 :adjustable t :fill-pointer 0))
        (chosen nil)
        (index3 0))
    ;; Build the indexes3
    (chem:do-molecules (mol aggregate)
      (chem:do-residues (res mol)
        (chem:do-atoms (atm res)
          (unless (chem:element-is-hydrogen (chem:atom/get-element atm))
            (when (or (eq t test) (funcall test atm))
              (push atm chosen)
              (vector-push-extend index3 heavy-atoms-index3)
              (incf index3 3))))))
    (values (copy-seq heavy-atoms-index3) (nreverse chosen))))


(defun align-heavy-atoms (aggregate coordinates)
  "Return a vector of coordinates, one for each of COORDINATES.  Use the heavy atoms in AGGREGATE for the alignment"
  (let* ((heavy-atoms-index3 (heavy-atoms-index3 aggregate))
         (transforms (loop with fixed-coords = (aref coordinates 0)
                           for index from 0 below (length coordinates)
                           for moveable-coords = (aref coordinates index)
                           collect (chem:nvectors-align fixed-coords moveable-coords heavy-atoms-index3)))
         (results (loop with results = (make-array (length transforms))
                        for one-coords across coordinates
                        for transform in transforms
                        for index from 0
                        for transformed = (chem:nvector-transform one-coords transform)
                        do (setf (aref results index) transformed)
                        finally (return results))))
    results))


(defclass one-coords ()
  ((aggregate :initarg :aggregate :reader aggregate)
   (one-coords :initarg :coords :reader one-coords)))


(defun coordinates-with-rms-differences (coordinates index3 &key (rms-threshold 0.2) key)
  "Return a list of COORDINATES that are different from each other by at least RMS-THRESHOLD considering only the INDEX3 coordinates.
If KEY is passed then use KEY to get the coordinates from COORDINATES"
  (flet ((different-p (unique-coords one-coords &optional key)
           "Return T if the ONE-COORDS are more than RMS-THRESHOLD different when compared to everything in UNIQUE-COORDS"
           (if key
               (loop for one-unique-coords-struct in unique-coords
                     for one-unique-coords = (funcall key one-unique-coords-struct)
                     for transform = (chem:nvectors-align one-unique-coords one-coords index3)
                     for rms = (chem:nvector-root-mean-square-difference-transformed one-unique-coords one-coords transform index3)
                     when (< rms rms-threshold)
                       do (return-from different-p nil))
               (loop for one-unique-coords in unique-coords
                     for transform = (chem:nvectors-align one-unique-coords one-coords index3)
                     for rms = (chem:nvector-root-mean-square-difference-transformed one-unique-coords one-coords transform index3)
                     when (< rms rms-threshold)
                       do (return-from different-p nil)))
           t))
    (let ((unique-coords nil))
      (if key
          (loop for coords-struct in coordinates
                for one-coords = (funcall key coords-struct)
                when (different-p unique-coords one-coords key)
                  do (push coords-struct unique-coords))
          (loop for maybe-one-coords in coordinates
                when (different-p unique-coords one-coords)
                  do (push one-coords unique-coords)))
      unique-coords)))

(defun aggregates-with-rms-differences (aggregates &key (rms-threshold 0.2))
  "Return a list of aggregates that are different from each other by at least RMS-THRESHOLD"
  (let ((coords-to-aggregate (make-hash-table))
        (all-coords (make-array (length aggregates)))
        (index3 (heavy-atoms-index3 (first aggregates))))
    (loop for aggregate in aggregates
          for one-coords = (chem:matter/extract-coordinates aggregate)
          for index from 0
          do (setf (aref all-coords index) one-coords
                   (gethash one-coords coords-to-aggregate) aggregate))
    (let ((unique-coords (coordinates-with-rms-differences all-coords index3 :rms-threshold rms-threshold)))
      (loop for coords across unique-coords
            collect (gethash coords coords-to-aggregate)))))

