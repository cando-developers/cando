

(in-package :leap)


;;(defparameter *ff* nil)
;;(export '*ff*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Define the amber-system class for leap specific force field configuration
;;;

(defclass amber-system () ())
(defparameter *amber-system* (make-instance 'amber-system))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Provide force field readers
;;;

(defmethod chem:parse-ptor-db (stream-in (system amber-system))
  (let ((ffptordb (core:make-cxx-object 'chem:ffptor-db))
        ptor)
    (loop for line = (read-line stream-in nil :eof)
       until (eq line :eof)
       until (= (length line) 0)
       do (multiple-value-bind (atom1 atom2 atom3 atom4 barrier phase periodicity scee scnb comment)
              (parse-ptor-line line)
            (let ((iperiod (floor (+ (abs periodicity) 0.5))))
              (unless (<= 1 iperiod 6)
                (error "Bad multiplicity ~s" periodicity))
              (unless ptor
                (setf ptor (core:make-cxx-object 'chem:ffptor :type1 atom1 :type2 atom2 :type3 atom3 :type4 atom4 :level :parameterized)))
              (chem:set-v-k-cal ptor iperiod barrier)
              (chem:set-phase-rad ptor iperiod phase)
              (when (and (> scee 0.0) (> scnb 0.0))
                (chem:set-scee ptor scee)
                (chem:set-scnb ptor scnb))
              (when (> (length comment) 0)
                (chem:set-ptor-comment ptor comment)))
            (when (> periodicity 0.0)
              (chem:add-ffptor ffptordb ptor)
              (setf ptor nil))))
    ffptordb))

(defun parse-ptor-line (line)
  (let* ((atom1 (intern (string-trim " " (subseq line 0 2)) :keyword))
        (atom2 (intern (string-trim " " (subseq line 3 5)) :keyword))
        (atom3 (intern (string-trim " " (subseq line 6 8)) :keyword))
        (atom4 (intern (string-trim " " (subseq line 9 11)) :keyword))
        (divider (parse-integer (subseq line 11 15)))
        (barrier-str (string-trim " " (subseq line 15 30)))
        (barrier (progn
                   (fortran::read-double-float-from-string barrier-str)))
        (phase (fortran::read-double-float-from-string (string-trim " " (subseq line 30 45))))
        (periodicity (fortran::read-double-float-from-string (string-trim " " (subseq line 45 (min 60 (length line))))))
        (comment (if (> (length line) 60) (subseq line 60) ""))
        (scee -1.0)
        (scnb -1.0))
    ;; Parse "SCEE=1.0 SCNB=1.0 2-methl..."
    ;;        .123456789.123456789
    (if (and (> (length comment) 5) (string-equal comment "SCEE=" :start1 0 :end1 5))
        (let ((scee-end (position #\space comment :start 5)))
          (setf scee (fortran::read-double-float-from-string (subseq comment 5 scee-end)))
          (let ((scnb-start (position-if #'alpha-char-p comment :start scee-end)))
            (if (string-equal comment "SCNB=" :start1 scnb-start :end1 14)
                (setf scnb (fortran::read-double-float-from-string (subseq comment 14 17)))))))
    (let ((barrier (/ barrier divider)))
      (values atom1 atom2 atom3 atom4 barrier phase periodicity scee scnb comment))))


;;;   Oh-Cg-Cg-Oh   1          -1.00           0.00          -1.  SCEE=1.0 SCNB=1.0 2,3-propanediol, 2-methyl 2,3-butanediol, etc KNK
;;; --> (values :|Oh| :|Cg|) :|Cg| :|Oh| -1.0
;;(cons :ph1 0.0) (cons :scee 1.0) (cons :scnb 1.0) (cons :comment "2,3-propanediol, 2-methyl 2,3-butanediol, etc KNK"))








#||
(defmethod parse-ptor-db-line (stream-in (system cando-system))
  ;; A function that reads one line and returns a FFPtor_sp object
  (let ((ffptor-db (core:make-cxx-object 'chem:ffptor-db)))
    ;;While not blank line
    
    ;; Parse a bunch of stuff and build an alist for one ptor
    (loop for line = (read-line stream-in)
       unless (= (length line) 0)
         do (let (parameters)
    ))
    (let ((ffptor (apply core:make-cxx-object 'chem:ffptor parameters)))
      (chem:add-ffptor ffptor-db ffptor))
    ffptor-db))
|#

    
