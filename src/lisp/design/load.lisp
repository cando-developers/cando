(in-package :design.load)


(defun build-topology (name constitution plugs)
  (let* ((topology (chem:make-topology name constitution plugs))
         (constitution-atoms (chem:get-constitution-atoms constitution))
         (constitution-atom-list (chem:constitution-atoms-as-list constitution-atoms)))
    (loop for constitution-atom in constitution-atom-list
          for constitution-atom-index = (chem:constitution-atoms-index constitution-atoms (chem:atom-name constitution-atom))
          for topology-atom-info = (chem:make-topology-atom-info nil 0.0 constitution-atom-index)
          do (chem:topology-add-atom-info topology topology-atom-info))
    topology))

(defun load-cdxml (filename)
  (let ((pname (merge-pathnames filename)))
    (unless (probe-file pname)
      (error "File not found: ~a" pname))
    (let* ((chemdraw (cando:load-chem-draw pname))
           (prepare-topologys (design.joint-tree::extract-prepare-topologys chemdraw))
           (topologys (loop for prep-top in prepare-topologys
                            for joint-template = (design.joint-tree::build-internal-coordinate-joint-template-tree prep-top)
                            for name = (design.joint-tree:name prep-top)
                            for constitution = (design.joint-tree:constitution prep-top)
                            for plugs = (design.joint-tree:plugs prep-top)
                            for topology = (build-topology name constitution plugs)
                            do (chem:add-topology constitution topology)
                               (chem:set-property topology :joint-template joint-template)
                            collect topology)))
      (loop for topology in topologys
            do (leap.core:register-variable (chem:get-name topology) topology))
      (values topologys chemdraw prepare-topologys))))

