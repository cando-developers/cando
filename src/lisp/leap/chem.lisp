

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

(esrap:defrule nonb-type-name
    (and (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9))
         (* (or (esrap:character-ranges (#\! #\~)))) ; entire ascii printable range
         )
  (:text t))



(esrap:defrule bond-angle-dihe-type-name
    (and (esrap:character-ranges (#\a #\z) (#\A #\Z) (#\0 #\9))
         (* (or (esrap:character-ranges (#\! #\,))
                (esrap:character-ranges (#\. #\~)))) ; entire ascii printable range, exclude #\- (dash)
         )
  (:text t))

(esrap:defrule ptor-types
    (and bond-angle-dihe-type-name parser.common-rules:whitespace* "-" parser.common-rules:whitespace*
         bond-angle-dihe-type-name parser.common-rules:whitespace* "-" parser.common-rules:whitespace*
         bond-angle-dihe-type-name parser.common-rules:whitespace* "-" parser.common-rules:whitespace*
         bond-angle-dihe-type-name
         )
  (:destructure (t1 s1a d1 s1b t2 s2a d2 s2b t3 s3a d3 s3b t4)
    (declare (ignore d1 d2 d3 s1a s1b s2a s2b s3a s3b))
    (list t1 t2 t3 t4)))

(esrap:defrule ptor-rule
  (and ptor-types parser.common-rules:whitespace+
       (esrap:? (and parser.common-rules:integer-literal/decimal
                     parser.common-rules:whitespace+))
          parser.common-rules:float-literal
          parser.common-rules:whitespace+
          parser.common-rules:float-literal
          parser.common-rules:whitespace+
       (or parser.common-rules:float-literal
           parser.common-rules:integer-literal/decimal)
       parser.common-rules:whitespace*
       (esrap:? (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
                     (* esrap:character))))
  (:destructure (types blank0 opt-idivf pk blank1 phase blank2 pn blank3 comment)
                (declare (ignore blank0 blank1 blank2 blank3))
                (list types
                      (if opt-idivf
                          (first opt-idivf)
                          1)
                      pk
                      phase
                      pn
                      (esrap:text comment))))

(esrap:defrule scee
  (and (or "SCEE" "scee")
       parser.common-rules:whitespace*
       "="
       parser.common-rules:whitespace*
       parser.common-rules:float-literal)
  (:destructure (scee space1 eq1 space2 value)
                value))

(esrap:defrule scnb
  (and (or "SCNB" "scnb")
       parser.common-rules:whitespace*
       "="
       parser.common-rules:whitespace*
       parser.common-rules:float-literal)
  (:destructure (scnb space1 eq1 space2 value)
                value))

(esrap:defrule scee-scnb-comment
  (and scee
       parser.common-rules:whitespace+
       scnb
       (or (and parser.common-rules:whitespace+
                (and (esrap:character-ranges (#\a #\z) (#\A #\Z))
                     (* esrap:character)))
           parser.common-rules:whitespace*))
  (:destructure (scee blank1 scnb comment)
                (declare (ignore blank1))
                (list scee scnb (esrap:text (second comment)))))

(defun maybe-parse-scee-scnb (comment)
  (when (stringp comment)
    (when (or (search "scee" comment)
              (search "SCEE" comment))
      (handler-case (esrap:parse 'scee-scnb-comment comment)
        (esrap:parse-error (err)
          ;; ignore parse errors - it's a comment
          nil)))))


(defun parse-ptor-line (line)
  (let ((parse (esrap:parse 'ptor-rule line)))
    (destructuring-bind (type-names idivf barrier phase periodicity comment)
        parse
      (destructuring-bind (atom1 atom2 atom3 atom4)
          (mapcar (lambda (name) (intern name :keyword)) type-names)
        (when (eq atom1 :x)
          (setf atom1 nil))
        (when (eq atom4 :x)
          (setf atom4 nil))
        (let ((scee -1.0)
              (scnb -1.0))
          ;; Parse "SCEE=1.0 SCNB=1.0 2-methl..."
          ;;        .123456789.123456789
          (let ((maybe-scee-scnb (maybe-parse-scee-scnb comment)))
            (when maybe-scee-scnb
              (destructuring-bind (got-scee got-scnb comment)
                  maybe-scee-scnb
                (setf scee got-scee
                      scnb got-scnb))))
          (let ((barrier (/ barrier idivf)))
            (values atom1 atom2 atom3 atom4 barrier phase periodicity scee scnb comment)))))))



(esrap:defrule nonbond-line
    (and parser.common-rules:whitespace*
         nonb-type-name parser.common-rules:whitespace+
         parser.common-rules:float-literal
         parser.common-rules:whitespace+
         parser.common-rules:float-literal
         (esrap:? (and parser.common-rules:whitespace+
                       (* esrap:character))))
  (:destructure (blank0 type blank1 radius blank2 edep maybe-comment)
                (declare (ignore blank0 blank1 blank2 maybe-comment))
                (list (intern type :keyword) radius edep)))

(defun chem:parse-nonbond-line (line)
  "Return :end when we hit the END keyword"
  (when (stringp line)
    (let ((tline (string-trim " " line)))
      (cond
        ((string= tline "END")
         :end)
        ((> (length tline) 0)
         (esrap:parse 'nonbond-line (string-trim " " tline)))
        (t nil)))))

  
#+(or)
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
           (divider (fortran::parse-double-float (subseq line divider-start divider-end)))
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
  #-nosmarts(chem:initialize-sybyl-type-rules)
  #+nosmarts(format t "Not calling chem:initialize-sybyl-type-rules~%"))


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
