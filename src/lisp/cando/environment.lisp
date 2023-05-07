
(in-package :cando)

-(defun lookup-topology (name &optional errorp)
  (multiple-value-bind (topology foundp)
      (chem:find-topology #+(or)(intern (string-upcase (string name)) :keyword) name errorp)
    topology))

(defun walk-topologys (func)
  (chem:walk-topologys chem:*cando-database* func))

(defun list-topologys ()
  (chem:walk-topologys chem:*cando-database* (lambda (topology &rest args) (print topology))))

(defun register-entity (name entity)
  (error "Implement register-entity")
  )

(defun lookup-entity (name &optional (errorp t))
  (let ((maybe-topology (lookup-topology name nil)))
    (if maybe-topology
        maybe-topology
        (when errorp
          (error "Implement find-entity for ~a - it isn't a topology" name)))))
