
(in-package :chem)


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
