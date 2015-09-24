(require :asdf)
(require :cando)

(require :cando)

(core:fieldsp (core:make-builtin 'chem:in-plug))

(dolist (x *all-cxx-classes*)
  (when (eq (symbol-package x) #.(find-package :CHEM))
    (format t "Building ~a ..." x)
    (let ((obj (make-builtin x)))
      (when (core:fieldsp obj) (format t "---- FIELDS -----")))
    (format t "...done~%")))

(let ((x (make-builtin 'chem:chem-info-match))))

(progn
  (setf *print-pretty* nil)
  (setf *default-pathname-defaults* (pathname "/Users/meister/Development/clasp/projects/cando/work/forceField/")))

(defmacro enable-print-object (aclass)
  `(progn
     (defmethod cl:print-object ((object ,aclass) stream)
       (if *print-readably*
           (progn
             (core:print-builtin-object object stream))
           (call-next-method)))))

(defparameter +chem-classes+ '(chem:fftypes-db
                                chem:chem-info))
(progn
  (dolist (c +chem-classes+)
    (eval `(enable-print-object ,c))))

(defparameter parms (chem:make-read-amber-parameters))

(with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
  (chem:read-types parms fin))

#+(or)(let ((*print-readably* t)
            (*print-circle* t)
            (*print-pretty* nil))
        (prin1 (chem:get-type-rules parms)))

(defparameter *x* (chem:get-type-rules parms))
(defparameter *rule* (chem:fftypes-get-rule *x* 0))

(chem:fftypes-number-of-rules *x*)
(let ((*print-readably* t)
      (*print-circle* t)
      (*print-pretty* nil))
  (prin1 *rule*))



(setq *print-readably* t)
(setq *print-pretty* nil)
(print (core:make-builtin 'chem:atom :name 'chem-keyword::C1 :element 'chem-keyword::C :properties nil :contents #()))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Reader
;;;

(defun read-builtin (stream char n)
  (declare (ignore char))
  (let ((description (read stream t nil t)))
    (apply #'core:make-builtin (car description) (cdr description))))

(set-dispatch-macro-character #\# #\I #'read-builtin)
(get-dispatch-macro-character #\# #\I)


#I(chem:residue-list)

(let ((*print-readably* t))
  (print *x*))


(setf *test* #I(CHEM:ATOM-TEST CHEM-KEYWORD::INT 3 CHEM-KEYWORD::TEST CHEM-KEYWORD::|SAPRingSize| ))

(let ((*print-readably* t)) (print *test*))

(defparameter *x* (core:make-builtin 'chem:chem-info))
(chem:compile-smarts *x* "C1CCCCC=1")
(let ((*print-circle* t))
  (print *x*))

(setf *print-readably* t)


(chem:as-smarts *rule*)

(setf *print-circle* t)
(defun read-builtin (stream char n)
  (declare (ignore char))
  (let ((x (read stream t)))
    (format t "About to make-instance with |~s|~%" x)
    (apply #'make-instance x)))

(set-macro-character #\{ #'read-builtin)

(defclass foo () ((a :initarg :a :accessor foo-a)))

(defparameter *x* (make-instance 'foo :a 1))
(setf (foo-a *x*) *x*)

(defmethod print-object ((x foo) stream)
  (if *print-readably*
      (format stream "{ ~_~s ~2i ~_~s ~_~W }" (class-name (class-of x)) :a (foo-a x))
      (call-next-method)))

(let ((*print-readably* t)
      (*print-circle* t))
  (prin1 *x*))

(defparameter *z* (with-input-from-string (sin "#1={foo :a #1# }  1") (read sin t)))
(foo-a *z*)


(print "Hello")





(apropos "cxx")

(cxx-instance-p 3.0)
