
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
  (error "Implement find-entity")
  #+(or)(chem:find-entity name errorp))
