(in-package :kin)


(defun build-topology-from-prep-top (prep-top)
  (let* ((name (design.joint-tree:name prep-top))
         (constitution (design.joint-tree:constitution prep-top))
         (plugs (design.joint-tree:plugs prep-top))
         (stereo-information (design.joint-tree:stereo-information prep-top))
         (topology (chem:make-topology name constitution plugs))
         (stereo-information (chem:get-stereo-information constitution))
         (constitution-atoms (chem:|Constitution_O::getConstitutionAtoms| constitution))
;;;         (constitution-atom-list (chem:constitution-atoms-as-list constitution-atoms)))
         (stereoisomer-atoms-list (loop for stereoisomer in (chem:stereoisomers-as-list stereo-information)
                                        for stereoisomer-name = (chem:get-name stereoisomer)
                                        for stereoisomer-atoms = (chem:make-stereoisomer-atoms stereoisomer-name constitution-atoms)
                                        collect stereoisomer-atoms)))
    (chem:set-stereoisomer-atoms topology :absolute stereoisomer-atoms-list)
    topology))

(defun setup (design)
  (design.joint-tree:extract-prepare-topologys design))

(defvar *prep-tops* nil)

(defun load-cdxml (filename)
  (let ((pname (merge-pathnames filename)))
    (unless (probe-file pname)
      (error "File not found: ~a" pname))
    (let* ((chemdraw (cando:load-chem-draw pname))
           (codes (chem:chem-draw-code chemdraw)))
      ;; Evaluate code - this must define the design:configure function
      (loop for code in codes
            do (eval code))
      (let* ((design (design:configure chemdraw))
             (prepare-topologys (design.joint-tree:extract-prepare-topologys design)))
        (setf *prep-tops* prepare-topologys)
        (when chem:*verbose*
          (format t "prepare-topologys = ~a~%" prepare-topologys))
        (let ((topologys (loop for prep-top in prepare-topologys
                               for joint-template = (design.joint-tree:build-internal-coordinate-joint-template-tree prep-top)
                               for topology = (build-topology-from-prep-top prep-top)
                               for constitution = (design.joint-tree:constitution prep-top)
                               do (chem:add-topology constitution topology)
                                  (chem:set-property topology :joint-template joint-template)
                                  (cando:register-topology topology)
                               collect topology)))
          (loop for topology in topologys
                do (cando:register-topology topology))
          (format t "Returning from load-cdxml~%")
          (make-instance 'design:grammar
                         :topologys topologys
                         :cap-name-map (design:cap-name-map design)))))))

