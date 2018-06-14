(in-package :design.load)

(defun build-stereoisomer (name chiral-atoms stereoisomer-index)
  "Build a stereoisomer. Assemble a name for the stereiosomer using the name and
the stereoisomer-index.   The stereoisomer-index is like a bit-vector that indicates
the stereochemistry of this stereoisomer (0 is :S and 1 is :R).
The name is constructed by appending \"{chiral-atom-name/(S|R),...}\".
So if name is \"ALA\" and stereoisomer-index is 1 the name becomes ALA{CA/S}."
  (if chiral-atoms
      (let* ((configurations (loop for chiral-index below (length chiral-atoms)
                                   for config = (if (logbitp chiral-index stereoisomer-index) :S :R)
                                   collect config))
             (stereo-configurations (mapcar (lambda (atom config) (chem:make-stereo-configuration (chem:get-name atom) config))
                                            chiral-atoms configurations))
             (new-name-string (format nil "~A[~{~A~^,~}]"
                                      name
                                      (mapcar (lambda (atom config)
                                                (format nil "~A/~A"
                                                        (string (chem:get-name atom))
                                                        (string config))) chiral-atoms configurations)))
             (new-name (intern new-name-string :keyword)))
        #+(or)
        (progn
          (format t "name -> ~a~%" name)
          (format t "chiral-atoms -> ~a~%" chiral-atoms)
          (format t "stereoisomer-index -> ~a~%" stereoisomer-index)
          (format t "new-name -> ~a~%" new-name))
        (chem:make-stereoisomer new-name nil stereoisomer-index stereo-configurations))
      (progn
        #+(or)(format t "Single stereoisomer name -> ~a~%" name)
        (chem:make-stereoisomer name name 0 nil))))

(defun build-stereoinformation (prep-topology)
  "Build a residue from the constitution-atoms in **prep-topology**.
Then identify all stereocenters sorted by distance from the in-plug/origin-plug atom.
Build stereoinformation from all of this and return it."
  (with-accessors ((name design.joint-tree:name)
                   (constitution-atoms design.joint-tree:constitution-atoms)
                   (plugs design.joint-tree:plugs)) prep-topology
    (let* ((residue (chem:build-residue constitution-atoms))
           (in-plug (find-if 'chem:get-is-in plugs)))
;;;      (chem:define-stereochemical-configurations-for-all-atoms residue)
      (unless in-plug
        (error "There must be an in-plug in the residue - the only plugs available are ~s" plugs))
      (let* ((root-atom-name (chem:root-atom-name in-plug))
             (root-atom (chem:atom-with-name residue root-atom-name))
             (spanning-loop (chem:make-spanning-loop root-atom))
             (atoms (chem:all-atoms spanning-loop))
             (chiral-atoms (remove-if (lambda (x) (null (eq (chem:get-stereochemistry-type x) :chiral))) atoms))
             (number-of-chiral-atoms (length chiral-atoms))
             (number-of-stereoisomers (expt 2 number-of-chiral-atoms)))
        #+(or)
        (progn
          (format t "residue: ~s~%" residue)
          (format t "plugs: ~s ~%" plugs)
          (format t "in plugs: ~s~%" (find-if 'chem:get-is-in plugs))
          (format t "root atom: ~s~%" root-atom)
          (format t "chiral atoms: ~s~%" chiral-atoms)
          (terpri)
          (terpri))
        (let ((stereoisomers (loop for stereoisomer-index below number-of-stereoisomers
                                   for stereoisomer = (build-stereoisomer name chiral-atoms stereoisomer-index)
                                   collect stereoisomer)))
          (chem:make-stereoinformation stereoisomers nil))))))

(defun load-cdxml (filename)
  (let ((pname (merge-pathnames filename)))
    (unless (probe-file pname)
      (error "File not found: ~a" pname))
    (let* ((chemdraw (cando:load-chem-draw pname))
           (prepare-topologys (design.joint-tree::extract-prepare-topologys chemdraw))
           (topologys (loop for prep-top in prepare-topologys
                            for joint-template = (design.joint-tree::build-internal-coordinate-joint-template-tree prep-top)
                            for stereoinfo = (build-stereoinformation prep-top)
                            for constitution = (chem:make-constitution (design.joint-tree:name prep-top)
                                                                       "no-comment"
                                                                       (design.joint-tree:constitution-atoms prep-top)
                                                                       stereoinfo
                                                                       (design.joint-tree:plugs prep-top)
                                                                       nil)
                            for top = (chem:make-topology (design.joint-tree::name prep-top)
                                                          constitution
                                                          (design.joint-tree:plugs prep-top))
                            do (chem:set-property top :joint-template joint-template)
                            collect top)))
      (loop for topology in topologys
            do (leap.core:register-variable (chem:get-name topology) topology))
      (values topologys chemdraw prepare-topologys))))

