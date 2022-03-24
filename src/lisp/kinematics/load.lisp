(in-package :kin)


(defun build-topology-from-prep-top (prep-top)
  (let* ((name (name prep-top))
         (constitution (constitution prep-top))
         (plugs (plugs prep-top))
         (stereo-information (stereo-information prep-top))
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
  (extract-prepare-topologys design))

(defvar *prep-tops* nil)

(defun load-cdxml (filename)
  (let ((pname (merge-pathnames filename)))
    (unless (probe-file pname)
      (error "File not found: ~a" pname))
    (let* ((chemdraw (cando:load-chem-draw pname))
           (codes (chem:chem-draw-code chemdraw)))
      ;; Evaluate code - this must define the configure function
      (loop for code in codes
            do (eval code))
      (let* ((design (configure chemdraw))
             (prepare-topologys (extract-prepare-topologys design)))
        (setf *prep-tops* prepare-topologys)
        (when chem:*verbose*
          (format t "prepare-topologys = ~a~%" prepare-topologys))
        (let ((topologys (loop for prep-top in prepare-topologys
                               for joint-template = (joint-template prep-top)
                               for topology = (build-topology-from-prep-top prep-top)
                               for constitution = (constitution prep-top)
                               do (chem:add-topology constitution topology)
                                  (chem:set-property topology :joint-template joint-template)
                                  (cando:register-topology topology)
                               collect topology)))
          (format t "Returning from load-cdxml~%")
          (make-instance 'grammar
                         :topologys topologys
                         :cap-name-map (cap-name-map design)))))))


(defun load-fragments* (filename)
  (unless (probe-file filename)
    (error "The file ~a was not found" filename))
  (let ((*package* (find-package :kin.fragment))
        (kin.fragment:*fragments* nil))
    (cl:load filename)
    (values kin.fragment:*fragments* kin.fragment:*cap-map*)))

(defun load-fragments (filename)
  (multiple-value-bind (fragments cap-map)
      (load-fragments* filename)
    (let* ((prep-tops (loop for fragment in fragments
                            for prep-top = (kin.fragment:extract-prepare-topology fragment)
                            collect prep-top))
           (topologys (loop for prep-top in prep-tops
                            for joint-template = (joint-template prep-top)
                            for topology = (build-topology-from-prep-top prep-top)
                            for constitution = (constitution prep-top)
                            do (chem:add-topology constitution topology)
                               (chem:set-property topology :joint-template joint-template)
                               (cando:register-topology topology)
                            collect topology)))
      (setf *prep-tops* prep-tops)
      (make-instance 'grammar
                     :topologys topologys
                     :cap-name-map cap-map))))


