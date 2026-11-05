(in-package #:topology)

(defclass cache-component-declaration (cando.serialize:serializable)
  ((kind :initarg :kind :reader cache-component-kind)
   (name :initarg :name :reader cache-component-name)
   (atom-name :initarg :atom-name :reader cache-component-atom-name))
  (:documentation "An explicit (kind cache-name atom-selector) declaration."))

(defmethod initialize-instance :after ((entry cache-component-declaration) &key)
  (unless (and (member (cache-component-kind entry) '(:ring :dihedral))
               (cache-component-name entry)
               (symbolp (cache-component-name entry))
               (symbol-package (cache-component-name entry))
               (keywordp (cache-component-atom-name entry)))
    (error "Invalid cache component declaration")))

(defun parse-cache-declarations (expressions &key atom-exists-p)
  "Parse :CACHE entries, normalizing atom names as DEFINE-TOPOLOGY does.
Optional ATOM-EXISTS-P validates each atom against the topology's graph.
Ring uniqueness and joint ancestry are checked later against assembled geometry."
  (unless (binned-shape-key-proper-list-p expressions)
    (error ":cache must be a proper list"))
  (let ((seen (make-hash-table :test #'eq)))
    (loop for expression in expressions
          collect
          (progn
            (unless (and (binned-shape-key-proper-list-p expression)
                         (= 3 (length expression)))
              (error "Cache entry must be (kind name atom): ~s" expression))
            (destructuring-bind (kind name atom) expression
              (unless (and atom (symbolp atom))
                (error "Invalid cache atom selector ~s" atom))
              (when (gethash name seen)
                (error "Duplicate cache declaration name ~s" name))
              (setf (gethash name seen) t)
              (let ((atom-name (intern (symbol-name atom) :keyword)))
                (when (and atom-exists-p (not (funcall atom-exists-p atom-name)))
                  (error "Cache field ~s selects missing atom ~s" name atom-name))
                (make-instance 'cache-component-declaration
                               :kind kind :name name :atom-name atom-name)))))))

(defun topology-cache-declarations (topology)
  "Return the explicitly declared cache components for TOPOLOGY."
  (getf (residue-properties (constitution topology)) :cache))
