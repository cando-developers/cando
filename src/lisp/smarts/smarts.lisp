(in-package :smarts)

(defun make-chem-info (&key tests smarts)
  (let ((ci (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts ci smarts)
    (if tests
        (chem:define-tests ci tests))
    ci))

(defun make-chem-info (&key tests smarts)
  (let ((ci (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts ci smarts)
    (if tests
        (chem:define-tests ci tests))
    ci))


(defmethod architecture.builder-protocol:make-node ((builder (eql 'chem:chem-info))
                                                    (head (eql :atom))
                                                    &rest args
                                                    &key kind symbol bounds)
  (let ((sym (intern (string-upcase symbol) :keyword)))
    (ecase kind
      (:organic
       (chem:create-sapelement sym)
       #+(or)(core:make-cxx-object 'chem:atom-test :sym sym :test :sapelement)
       )
      (:aromatic
       (chem:create-sap-aromatic-element sym)
       #+(or)(core:make-cxx-object 'chem:atom-test :sym sym :test :saparomatic-element)))))

(defmethod architecture.builder-protocol:make-node ((builder (eql 'chem:chem-info))
                                                    (head (eql :chain))
                                                    &rest args)
  
  (format t ":chain make-node head: ~s args: ~s~%" head args))


#+(or)(defmethod architecture.builder-protocol:relate ((builder (eql 'chem:chem-info))
                                                 (head (eql :atom))
                                                 left
                                                 right)
  (format t ":atom relate head: ~s left: ~s right:~s~%" head left right))

#+(or)(defmethod architecture.builder-protocol:relate ((builder (eql 'chem:chem-info))
                                                 (head (eql :element))
                                                 left
                                                 right)
  (format t ":element relate head: ~s left: ~s right:~s~%" head left right))




(defgeneric build (head tree &rest args))

(defmethod build ((head (eql :atom)) tree &rest args)
  (cond
    ((consp tree)
     (format t "Tree: ~s~%" (car (car tree)))
     (apply 'build (car (car tree))))
    (t (let* ((kind (getf args :kind))
              (symbol-name (getf args :symbol))
              (symbol (intern (string-upcase symbol-name) :keyword))
              (bounds (getf args :bounds)))
         (ecase kind
           (:organic
            (core:make-cxx-object 'chem:atom-test :sym symbol :test :sapelement))
           (:aromatic
            (core:make-cxx-object 'chem:atom-test :sym symbol :test :saparomatic-element))
           )))))

(defmethod build ((head (eql :atom)) (tree cons) &rest args)
     (format t "Tree: ~s~%" (car (car tree)))
  (apply 'build (car (car tree))))

(defmethod build ((head (eql :atom)) (tree null) &rest args)
  (let* ((kind (getf args :kind))
         (symbol-name (getf args :symbol))
         (symbol (intern (string-upcase symbol-name) :keyword))
         (bounds (getf args :bounds)))
    (ecase kind
      (:organic
       (core:make-cxx-object 'chem:atom-test :sym symbol :test :sapelement))
      (:aromatic
       (core:make-cxx-object 'chem:atom-test :sym symbol :test :saparomatic-element))
      )))
