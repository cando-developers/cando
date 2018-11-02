
(in-package :cando)

(defparameter chem:*cando-database* (chem:make-cando-database :default)
  "This is where topology's are stored")

(defun register-topology (name topology)
  (chem:setf-find-topology chem:*cando-database* name topology))

(defun lookup-topology (name &optional errorp)
  (multiple-value-bind (topology foundp)
      (chem:find-topology name errorp)
    topology))
    

(defun register-entity (name entity)
  (error "Implement register-entity")
  )

(defun lookup-entity (name &optional (errorp t))
  (let ((maybe-topology (lookup-topology name nil)))
    (if maybe-topology
        maybe-topology
        (when errorp
          (error "Implement find-entity for ~a - it isn't a topology" name)))))
