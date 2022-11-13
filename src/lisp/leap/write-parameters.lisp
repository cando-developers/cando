(in-package :leap.write-parameters)


(defun merge-force-fields (force-field-name)
  (let ((combined-force-field (chem:find-force-field force-field-name))
        (merged-force-field (chem:force-field/make)))
    (loop for partial-force-field in (chem:combined-force-field/force-fields-as-list combined-force-field)
          do (chem:force-field/force-field-merge merged-force-field partial-force-field))
    merged-force-field))


(defun ensure-string (obj)
  (if obj
      (string obj)
      ""))

(defun rad-to-deg (rad)
  (/ rad 0.0174533))

(defun write-atom-parameters (ffnonbond-db stream)
  (format stream "MASS~%")
  (let ((ffnonbond-term-vector (chem:ffnonbond-db/term-vector ffnonbond-db)))
    (loop for index below (length ffnonbond-term-vector)
          for nonbond = (aref ffnonbond-term-vector index)
          for type = (chem:ffnonbond/get-type nonbond)
          for mass = (chem:ffnonbond/get-mass nonbond)
          for polarizability = (chem:ffnonbond/get-polarizability nonbond)
          for comment = (chem:ffnonbond/get-type-comment nonbond)
          do (sys:bformat stream "%-2s %-8.3f   %8.3f   %s%N"
                          (string type)
                          mass
                          polarizability
                          (if comment comment "")))
    ))

(defun write-stretch-parameters (ffstretch-db stream)
  (format stream "~%BOND~%")
  (let ((ffstretch-parms (chem:ffparameter-base-db/parameters ffstretch-db)))
    (maphash (lambda (key ffstretch)
               (declare (ignore key))
               (multiple-value-bind (type1 type2)
                   (chem:ffstretch/get-types ffstretch)
                 (sys:bformat stream "%-2s-%-2s%8.2f%8.3f    %s%N"
                              (string type1)
                              (string type2)
                              (chem:ffstretch/get-kb-kcal-per-angstrom-squared ffstretch)
                              (chem:ffstretch/get-r0-angstrom ffstretch)
                              (ensure-string (chem:ffparameter/get-comment ffstretch)))
                 ))
             ffstretch-parms)))

(defun write-angle-parameters (ffangle-db stream)
  (format stream "~%ANGLE~%")
  (let ((ffangle-parms (chem:ffparameter-base-db/parameters ffangle-db)))
    (maphash (lambda (key ffangle)
               (declare (ignore key))
               (multiple-value-bind (type1 type2 type3)
                   (chem:ffangle/get-types ffangle)
                 (sys:bformat stream "%-2s-%-2s-%-2s%9.3f%12.3f    %s%N"
                              (string type1)
                              (string type2)
                              (string type3)
                              (chem:ffangle/get-k2-kcal-per-radian-squared ffangle)
                              (rad-to-deg (chem:ffangle/get-angle-radian ffangle))
                              (ensure-string (chem:ffparameter/get-comment ffangle))
                              )))
             ffangle-parms)))


(defun write-ptor-periodicity (ptor periodicity max-periodicity stream)
  (multiple-value-bind (type1 type2 type3 type4)
      (chem:ffptor/get-types ptor)
    (let ((pk (* (chem:ffptor/get-idivf ptor periodicity)
                 (chem:ffptor/get-v-kcal ptor periodicity)))) ; See leap::parse-ptor-line
                                        ; From https://ambermd.org/FileFormats.php#topology
                                        ; The factor by which the torsional barrier is divided.
                                        ; Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                                        ; details. Basically, the actual torsional potential is
                                        ; (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
      (sys:bformat stream "%-2s-%-2s-%-2s-%-2s%4d%9.3f%14.3f%15.3f   %s%N"
                   (if type1 (string type1) "X")
                   (string type2)
                   (string type3)
                   (if type1 (string type4) "X")
                   (chem:ffptor/get-idivf ptor periodicity)
                   pk
                   (chem:ffptor/get-phase-degrees ptor periodicity)
                   (float (if (= periodicity max-periodicity)
                              periodicity ; periodicity > 0 indicates this is the last one
                              (- periodicity)) 0s0) ; periodicity < 0 indicates there are more coming
                   (if (= periodicity max-periodicity)
                       (ensure-string (chem:ffparameter/get-comment ptor))
                       "")))))

(defun write-ptor-parameters (ffptor-db stream)
  (format stream "~%DIHE~%")
  (let* ((ffptor-parms (chem:ffparameter-base-db/parameters ffptor-db))
         (ffptor-general-terms nil)
         (ffptor-specific-terms nil))
    ;; Separate out the general and specific terms
    (loop for ptor being each hash-values of ffptor-parms
          do (multiple-value-bind (t1 t2 t3 t4)
                 (chem:ffptor/get-types ptor)
               (declare (ignore t2 t3))
               (cond
                 ((and t1 t4)
                  (push ptor ffptor-specific-terms))
                 ((null (or t1 t4))
                  (push ptor ffptor-general-terms))
                 (t (error "ptor ~a is neither general nor specific" ptor)))))
    ;; Write out the general terms
    (loop for ptor in ffptor-general-terms
          do (loop with max-periodicity = (chem:ffptor/max-periodicity ptor)
                   for periodicity from 1 to max-periodicity
                   when (chem:ffptor/has-periodicity ptor periodicity)
                     do (write-ptor-periodicity ptor periodicity max-periodicity stream)))
    ;; Write out the specific terms
    (loop for ptor in ffptor-specific-terms
          do (loop with max-periodicity = (chem:ffptor/max-periodicity ptor)
                   for periodicity from 1 to max-periodicity
                   when (chem:ffptor/has-periodicity ptor periodicity)
                     do (write-ptor-periodicity ptor periodicity max-periodicity stream)))))

(defun write-itor-periodicity (itor periodicity max-periodicity stream)
  (multiple-value-bind (type1 type2 type3 type4)
      (chem:ffitor/get-types itor)
    (let ((pk (chem:ffitor/get-v-kcal itor periodicity)))
                                        ; See leap::parse-itor-line
                                        ; From https://ambermd.org/FileFormats.php#topology
                                        ; The factor by which the torsional barrier is divided.
                                        ; Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                                        ; details. Basically, the actual torsional potential is
                                        ; (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
      (sys:bformat stream "%-2s-%-2s-%-2s-%-2s%9.3f%14.3f%15.3f   %s%N"
                   (if type1 (string type1) "X")
                   (if type2 (string type2) "X")
                   (string type3)
                   (if type1 (string type4) "X")
                   pk
                   (chem:ffitor/get-phase-degrees itor periodicity)
                   (float (if (= periodicity max-periodicity)
                              periodicity ; periodicity > 0 indicates this is the last one
                              (- periodicity)) 0s0) ; periodicity < 0 indicates there are more coming
                   (if (= periodicity max-periodicity)
                       (ensure-string (chem:ffparameter/get-comment itor))
                       "")))))

(defun write-itor-parameters (ffitor-db stream)
  (format stream "~%IMPROPER~%")
  (let* ((ffitor-parms (chem:ffparameter-base-db/parameters ffitor-db))
         (ffitor-general-terms nil)
         (ffitor-specific-terms nil))
    ;; Separate out the general and specific terms
    (loop for itor being each hash-values of ffitor-parms
          do (multiple-value-bind (t1 t2 t3 t4)
                 (chem:ffitor/get-types itor)
               (declare (ignore t3))
               (cond
                 ((and t1 t2 t4)
                  (push itor ffitor-specific-terms))
                 (t
                  (push itor ffitor-general-terms)))))
    ;; Write out the general terms
    (loop for itor in ffitor-general-terms
          do (loop with max-periodicity = (chem:ffitor/max-periodicity itor)
                   for periodicity from 1 to max-periodicity
                   when (chem:ffitor/has-periodicity itor periodicity)
                     do (write-itor-periodicity itor periodicity max-periodicity stream)))
    ;; Write out the specific terms
    (loop for itor in ffitor-specific-terms
          do (loop with max-periodicity = (chem:ffitor/max-periodicity itor)
                   for periodicity from 1 to max-periodicity
                   when (chem:ffitor/has-periodicity itor periodicity)
                     do (write-itor-periodicity itor periodicity max-periodicity stream)))))

(defun write-nonbond-parameters (ffnonbond-db stream)
  (let ((ffnonbond-term-vector (chem:ffnonbond-db/term-vector ffnonbond-db)))
    (loop for index below (length ffnonbond-term-vector)
          for nonbond = (aref ffnonbond-term-vector index)
          for type = (chem:ffnonbond/get-type nonbond)
          for radius = (chem:ffnonbond/get-radius-angstroms nonbond)
          for epsilon = (chem:ffnonbond/get-epsilon-kcal nonbond)
          do (sys:bformat stream "  %-2s%16.4f%8.4f%N"
                          (string type)
                          radius
                          epsilon ))))

(defun write-force-field-parameters (force-field &optional (stream t))
  (write-atom-parameters (chem:get-nonbond-db force-field) stream)
  (sys:bformat stream "%N")
  (write-stretch-parameters (chem:get-stretch-db force-field) stream)
  (sys:bformat stream "%N")
  (write-angle-parameters (chem:get-angle-db force-field) stream)
  (sys:bformat stream "%N")
  (write-ptor-parameters (chem:get-ptor-db force-field) stream)
  (sys:bformat stream "%N")
  (write-itor-parameters (chem:get-itor-db force-field) stream)
  (sys:bformat stream "%N")
  (write-nonbond-parameters (chem:get-nonbond-db force-field) stream)
  (sys:bformat stream "%N")
  (sys:bformat stream "END%N")
  )

