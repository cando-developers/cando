

(in-package :leap)


;;(defparameter *ff* nil)
;;(export '*ff*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Define the amber-system class for leap specific force field configuration
;;;

(defclass amber-system () ())
(defparameter *amber-system* (make-instance 'amber-system))

(define-condition bad-multiplicity ()
  ((line :initarg :line :accessor bad-multiplicity-line)
   (multiplicity :initarg :multiplicity :accessor bad-multiplicity-multiplicity)))

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
                 (if (not (<= 1 iperiod 6))
                     (error "bad multiplicity in line: ~a" line) #+(or)(signal  'bad-multiplicity :line line :multiplicity periodicity)
                     (progn
                       (unless ptor
                         (setf ptor (core:make-cxx-object 'chem:ffptor :type1 atom1 :type2 atom2 :type3 atom3 :type4 atom4 :level :parameterized)))
                       (chem:set-v-k-cal ptor iperiod barrier)
                       (chem:set-phase-rad ptor iperiod phase)
                       (when (and (> scee 0.0) (> scnb 0.0))
                         (chem:set-scee ptor scee)
                         (chem:set-scnb ptor scnb))
                       (when (> (length comment) 0)
                         (chem:set-ptor-comment ptor comment))
                       (when (> periodicity 0.0)
                         (chem:add-ffptor ffptordb ptor)
                         (setf ptor nil)))))))
    ffptordb))

(defconstant +degrees-to-radians+ 0.0174533)

(defun parse-ptor-line (line)
  (flet ((not-white-space-p (c) (char>= c #\!))
         (white-space-p (c) (char< c #\!)))
    (let* ((atom1-name (string-trim " " (subseq line 0 2)))
           (atom1 (if (string= atom1-name "X") NIL (intern atom1-name :keyword)))
           (atom2 (intern (string-trim " " (subseq line 3 5)) :keyword))
           (atom3 (intern (string-trim " " (subseq line 6 8)) :keyword))
           (atom4-name (string-trim " " (subseq line 9 11)))
           (atom4 (if (string= atom4-name "X") NIL (intern atom4-name :keyword)))
           (divider-start (position-if #'digit-char-p line :start 11))
           (divider-end (position #\space line :start divider-start))
           (divider (round (fortran::parse-double-float (subseq line divider-start divider-end))))
           (barrier-start (position-if #'not-white-space-p line :start divider-end))
           (barrier-end (position-if #'white-space-p line :start barrier-start))
           (barrier (fortran::parse-double-float (subseq line barrier-start barrier-end)))
           (phase-start (position-if #'not-white-space-p line :start barrier-end))
           (phase-end (position-if #'white-space-p line :start phase-start))
           (phase-deg (fortran::parse-double-float (subseq line phase-start phase-end)))
           (phase (* phase-deg +degrees-to-radians+))
           (periodicity-start (position-if #'not-white-space-p line :start phase-end))
           (periodicity-end (position-if #'white-space-p line :start periodicity-start))
           (periodicity (fortran::parse-double-float (subseq line periodicity-start periodicity-end)))
           (comment-start (if periodicity-end (position-if #'not-white-space-p line :start periodicity-end) nil))
           (comment (if comment-start (subseq line comment-start) ""))
           (scee -1.0)
           (scnb -1.0))
      ;; Parse "SCEE=1.0 SCNB=1.0 2-methl..."
      ;;        .123456789.123456789
      (if (and (> (length comment) 5) (string-equal comment "SCEE=" :start1 0 :end1 5))
          (let ((scee-end (position #\space comment :start 5)))
            (setf scee (fortran::parse-double-float (subseq comment 5 scee-end)))
            (let ((scnb-start (position-if #'alpha-char-p comment :start scee-end)))
              (if (string-equal comment "SCNB=" :start1 scnb-start :end1 14)
                  (setf scnb (fortran::parse-double-float (subseq comment 14 17)))))))
      (let ((barrier (/ barrier divider)))
        (values atom1 atom2 atom3 atom4 barrier phase periodicity scee scnb comment)))))



;;; ------------------------------------------------------------
;;;
;;; Initialize sybyl type rules
;;;

(defvar chem:*sybyl-type-assignment-rules* nil
  "This symbol is defined and used in mol2.cc")

(defun chem:initialize-sybyl-type-rules ()
  (setf chem:*sybyl-type-assignment-rules*
        (with-open-file (fin (merge-pathnames "ATOMTYPE_SYBYL.DEF"
                                              (translate-logical-pathname #P"source-dir:extensions;cando;src;data;force-field;"))
                             :direction :input)
          (leap.antechamber-type-definition-parser:read-antechamber-type-rules fin))))


(eval-when (:execute :load-toplevel)
  ;;; Initialize the sybyl-type-rules
  (chem:initialize-sybyl-type-rules))


(defmethod chem:compute-merged-nonbond-force-field-for-aggregate (aggregate)
  (let* ((aggregate-force-field-name (chem:force-field-name aggregate))
         (aggregate-force-field (leap.core:nonbond-force-field-component aggregate-force-field-name))
         (force-field-names (chem:map-molecules
                             'list
                             (lambda (molecule)
                               (chem:force-field-name molecule))
                             aggregate))
         (unique-force-field-names (remove-duplicates force-field-names))
         (force-fields-and-name-alist (loop for name in unique-force-field-names
                                            for nonbond-force-field = (leap.core:nonbond-force-field-component name)
                                            collect (cons nonbond-force-field name))))
    (chem:combine-nonbond-force-fields aggregate-force-field force-fields-and-name-alist)))
