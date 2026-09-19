(in-package #:cando)

(defun geometry-residue-maps (matter)
  (let ((atom-residues (make-hash-table :test #'eq))
        (residue-indices (make-hash-table :test #'eq))
        (index 0))
    (chem:do-residues (residue matter)
      (setf (gethash residue residue-indices) (incf index))
      (chem:do-atoms (atom residue)
        (when (gethash atom atom-residues)
          (error "Atom appears in multiple residue entries: ~s" atom))
        (setf (gethash atom atom-residues) residue)))
    (values atom-residues residue-indices)))

(defun write-geometry-atom (stream atom atom-residues residue-indices)
  (let ((residue (gethash atom atom-residues)))
    (unless residue (error "Atom lies outside mapped residues: ~s" atom))
    (let ((number (chem:get-file-sequence-number residue)))
      (format stream "residue[~a] ~a:~a (traversal-index=~d)"
              (if (= number #xffffffff) "unknown" number)
              (chem:get-name residue) (chem:get-name atom)
              (gethash residue residue-indices)))))

(defun report-long-bonds (matter cutoff &key (stream *standard-output*))
  "Print bonds longer than CUTOFF Angstrom using atoms' current positions.
Each endpoint is identified by its file sequence number (CHEM:GET-FILE-SEQUENCE-NUMBER),
residue name and atom name, followed by its one-based traversal index.
The PDB reader stores the original residue sequence number without encoding
insertion codes. Sequence numbers may repeat; traversal indices are
unique within MATTER. Return the number of bonds reported.
Uses the deduplicating bond loop and does not modify MATTER."
  (unless (and (realp cutoff) (< 0 cutoff most-positive-double-float))
    (error "Bond length cutoff must be finite and positive: ~s" cutoff))
  (multiple-value-bind (atom-residues residue-indices) (geometry-residue-maps matter)
   (let ((count 0))
    (chem:map-bonds
     nil
     (lambda (a b order bond)
       (declare (ignore order bond))
       (let ((ra (gethash a atom-residues))
             (rb (gethash b atom-residues))
             (distance (geom:calculate-distance (chem:get-position a)
                                                 (chem:get-position b))))
         (unless (and ra rb)
           (error "Bond endpoint lies outside mapped residues: ~s -- ~s" a b))
         (unless (and (realp distance) (<= 0 distance most-positive-double-float))
           (error "Nonfinite bond length: ~s -- ~s" a b))
         (when (> distance cutoff)
           (incf count)
           (format stream "~&~,4f A  " distance)
           (write-geometry-atom stream a atom-residues residue-indices)
           (write-string " -- " stream)
           (write-geometry-atom stream b atom-residues residue-indices)
           (terpri stream))))
     matter)
    (finish-output stream)
    count)))

(defun report-vdw-clashes (matter &key (ratio-cutoff 0.6d0)
                                     (include-14 t) (stream *standard-output*))
  "Report pairs closer than RATIO-CUTOFF times their summed elemental VDW radii.
Uses current atom positions and CHEM:VDW-RADIUS-FOR-ELEMENT, not force-field
Lennard-Jones parameters. This is a geometric diagnostic, not an energy estimate.
Exclude 1-2 and 1-3 pairs; include 1-4 pairs unless INCLUDE-14 is NIL.
Print each pair once, worst distance/radius ratio first. Report distance,
radius sum, overlap, ratio, and both residue/atom identities. Return clash count.
Includes hydrogens and intermolecular contacts. Uses the deduplicating bond
loop to build connectivity. Does not modify MATTER. Pair scanning is O(N^2)."
  (unless (and (realp ratio-cutoff) (< 0 ratio-cutoff 1))
    (error "VDW clash ratio cutoff must be between zero and one: ~s" ratio-cutoff))
  (multiple-value-bind (atom-residues residue-indices) (geometry-residue-maps matter)
    (let* ((atoms (chem:map-atoms 'vector #'identity matter))
           (n (length atoms))
           (indices (make-hash-table :test #'eq))
           (neighbors (make-array n :initial-element nil))
           (positions (make-array n))
           (radii (make-array n))
           (clashes nil))
      (loop for atom across atoms for i from 0
            for radius = (chem:vdw-radius-for-element (chem:get-element atom))
            do (unless (and (realp radius) (< 0 radius most-positive-double-float))
                 (error "Invalid elemental VDW radius for ~s: ~s" atom radius))
               (setf (gethash atom indices) i
                     (aref positions i) (chem:get-position atom)
                     (aref radii i) radius))
      (chem:map-bonds
       nil (lambda (a b order bond)
             (declare (ignore order bond))
             (let ((i (gethash a indices)) (j (gethash b indices)))
               (unless (and i j) (error "Bond endpoint outside MATTER: ~s -- ~s" a b))
               (pushnew j (aref neighbors i))
               (pushnew i (aref neighbors j))))
       matter)
      (dotimes (i n)
        ;; Breadth-first distances avoid mistakes in rings or duplicated bonds.
        (let ((excluded (make-hash-table :test #'eql)) (frontier (list i)))
          (setf (gethash i excluded) t)
          (dotimes (depth (if include-14 2 3))
            (let ((next nil))
              (dolist (j frontier)
                (dolist (k (aref neighbors j))
                  (unless (gethash k excluded)
                    (setf (gethash k excluded) t)
                    (push k next))))
              (setf frontier next)))
          (loop for j from (1+ i) below n
                unless (gethash j excluded)
                  do (let* ((distance (geom:calculate-distance (aref positions i) (aref positions j)))
                            (radius-sum (+ (aref radii i) (aref radii j))))
                       (unless (and (realp distance) (<= 0 distance most-positive-double-float))
                         (error "Nonfinite distance: ~s -- ~s" (aref atoms i) (aref atoms j)))
                       (when (< distance (* ratio-cutoff radius-sum))
                         (push (list (/ distance radius-sum) distance radius-sum i j) clashes))))))
      (setf clashes (stable-sort clashes #'< :key #'first))
      (dolist (clash clashes)
        (destructuring-bind (ratio distance radius-sum i j) clash
          (format stream "~&VDW clash: ~,4f A  radii-sum=~,4f A overlap=~,4f A ratio=~,3f  "
                  distance radius-sum (- radius-sum distance) ratio)
          (write-geometry-atom stream (aref atoms i) atom-residues residue-indices)
          (write-string " -- " stream)
          (write-geometry-atom stream (aref atoms j) atom-residues residue-indices)
          (terpri stream)))
      (finish-output stream)
      (length clashes))))

(defun report-geometry (matter &key (bond-cutoff 2.0d0) (vdw-ratio-cutoff 0.6d0)
                                  (include-14 t) (stream *standard-output*))
  "Print long-bond and geometric VDW-clash reports using current atom positions.
Return a plist containing :LONG-BONDS and :VDW-CLASHES counts. Cutoffs are in
Angstrom for bonds and a dimensionless radius-sum fraction for clashes."
  (format stream "~&=== Long bonds (> ~,3f A) ===~%" bond-cutoff)
  (let ((bonds (report-long-bonds matter bond-cutoff :stream stream)))
    (format stream "~&=== VDW clashes (distance/radii-sum < ~,3f; 1-4 ~a) ===~%"
            vdw-ratio-cutoff (if include-14 "included" "excluded"))
    (let ((clashes (report-vdw-clashes matter :ratio-cutoff vdw-ratio-cutoff
                                     :include-14 include-14 :stream stream)))
      (format stream "~&Geometry summary: ~d long bonds, ~d VDW clashes.~%" bonds clashes)
      (finish-output stream)
      (list :long-bonds bonds :vdw-clashes clashes))))
