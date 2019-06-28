;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :fep)

(defclass fep-calculation (calculation)
  ((solvent-box :initform :tip3pbox :initarg :solvent-box :accessor solvent-box)
   (solvent-buffer :initform 12.0 :initarg :solvent-buffer :accessor solvent-buffer)
   (solvent-closeness :initform 0.75 :initarg :solvent-closeness :accessor solvent-closeness)
   (solvent-residue-name :initform :WAT :initarg :solvent-residue-name :accessor solvent-residue-name)
   (script-0-setup :initform 'default-script-0-setup :accessor script-0-setup)
   (script-1-leap :initform 'default-script-1-leap :accessor script-1-leap)))

(defmethod print-object ((obj fep-calculation) stream)
  (if *print-readably*
      (progn
        (format stream "#$(fep::fep-calculation ")
        (loop for slot in (clos:class-slots (find-class 'fep::fep-calculation))
              for slot-name = (clos:slot-definition-name slot)
              for initargs = (clos:slot-definition-initargs slot)
              if (and (car initargs) (slot-boundp obj slot-name))
                do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
        (format stream ") "))
      (print-unreadable-object (obj stream)
        (format stream "fep-calculation"))))

(defmacro powerloop ((&rest clauses) &rest final)
  (if (null clauses)
      `(progn ,@final)
    `(loop ,@(first clauses) do (powerloop (,@(rest clauses)) ,@final))))
 
(defun make-script-1-leap (calculation &key input-feps-file)
  (with-top-directory (calculation)
    (let (work-list morph-jobs)
      (powerloop
       ((for receptor in (receptors calculation)
             for jobs = (jobs calculation))
        (for morph in (morphs jobs)
             for source = (source morph)
             for target = (target morph)))
       (let (analysis-jobs)
         (loop
           for side in '(:ligand :complex)
           for side-name = (format nil "~a-vdw-bonded" (string-downcase side))
           for input-topology-file = (make-instance 'morph-side-topology-file
                                                    :morph morph :side side
                                                    :name side-name
                                                    :extension "parm7")
           for coord-node = (make-instance 'morph-side-coordinate-file
                                           :morph morph :side side
                                           :name side-name
                                           :extension "rst7")
           for mol2-node = (make-instance 'morph-side-mol2-file
                                         :morph morph :side side
                                         :name side-name)
           for script = (make-instance 'morph-side-script
                                       :morph morph
                                       :side side
                                       :script *solvate-addion-morph-side-script*
                                       :name "solvate-addion"
                                       :extension "lisp")
           for solvate-addion-job = (connect-graph
                                     (make-instance 'morph-side-cando-job
                                                    :morph morph
                                                    :side side
                                                    :script script
                                                    :inputs (arguments :%INPUT% input-feps-file)
                                                    :outputs (arguments :%COORDINATES% coord-node
                                                                        :%TOPOLOGY% input-topology-file
                                                                        :%MOL2% mol2-node)
                                                    :makefile-clause (standard-cando-makefile-clause script)))
           for morph-side-prepare-job = (make-morph-side-prepare morph side
                                                                 :input-topology-file input-topology-file
                                                                 :input-coordinate-file coord-node)
           for strip-job = (make-morph-side-strip morph side
                                                  :input-topology-file input-topology-file
                                                  :input-coordinate-file (output-file morph-side-prepare-job :-r))
           do (let (stage-jobs)
                (loop for stage in (case (stages morph)
                                     (1 '(:vdw-bonded))
                                     (3 '(:decharge :vdw-bonded :recharge))
                                     (otherwise (error "Illegal number of stages ~a - only 1 or 3 are allowed for morph ~s" (stages morph) morph)))
                      with source = (output-file strip-job :%SOURCE%)
                      with target = (output-file strip-job :%TARGET%)
                      with solvated = (output-file strip-job :%SOLVATED%)
                      for script-source = (ecase stage
                                            (:decharge *decharge*)
                                            (:vdw-bonded nil)
                                            (:recharge *recharge*))
                      for script = (unless (eq stage :vdw)
                                     (make-instance 'morph-side-stage-script-file :morph morph :side side :stage stage :script script-source
                                                                                  :name (string-downcase stage)
                                                                                  :extension "lisp"))
                      for inputs = (ecase stage
                                     (:decharge (arguments :%FEPS% input-feps-file
                                                           :%SOLVATED% solvated
                                                           :%SOURCE% source))
                                     (:vdw-bonded nil)
                                     (:recharge (arguments :%FEPS% input-feps-file
                                                           :%SOLVATED% solvated
                                                           :%TARGET% target)))
                      for outputs = (ecase stage
                                      (:decharge
                                       (arguments
                                        :%DECHARGE-MOL2% (make-instance 'morph-side-stage-mol2-file
                                                                     :morph morph
                                                                     :side side
                                                                     :stage stage
                                                                     :name "decharge")
                                        :%DECHARGE-TOPOLOGY% (make-instance 'morph-side-stage-topology-file
                                                                          :morph morph
                                                                          :side side
                                                                          :stage stage
                                                                          :name "decharge")
                                        :%DECHARGE-COORDINATES% (make-instance 'morph-side-stage-coordinates-file
                                                                             :morph morph
                                                                             :side side
                                                                             :stage stage
                                                                             :name "decharge")))
                                      (:vdw-bonded nil)
                                      (:recharge
                                       (arguments
                                        :%RECHARGE-MOL2% (make-instance 'morph-side-stage-mol2-file
                                                                     :morph morph
                                                                     :side side
                                                                     :stage stage
                                                                     :name "recharge")
                                        :%RECHARGE-TOPOLOGY% (make-instance 'morph-side-stage-topology-file
                                                                          :morph morph
                                                                          :side side
                                                                          :stage stage
                                                                          :name "recharge")
                                        :%RECHARGE-COORDINATES% (make-instance 'morph-side-stage-coordinates-file
                                                                             :morph morph
                                                                             :side side
                                                                             :stage stage
                                                                             :name "recharge"))))
                      for stage-job = (unless (eq stage :vdw-bonded)
                                        (connect-graph
                                         (make-instance 'morph-side-stage-cando-job
                                                        :morph morph :side side :stage stage :script script
                                                        :inputs inputs
                                                        :outputs outputs
                                                        :makefile-clause (standard-cando-makefile-clause script))))
                      do (let (decharge-jobs
                               vdw-jobs
                               recharge-jobs
                               (lambda-values (loop for window-index from 0 below (windows morph)
                                                    collect (/ (float window-index) (1- (windows morph))))))
                           (loop for lambda-value in lambda-values
                                 for lambda-label = (format nil "~5,3f" lambda-value)
                                 do (case stage
                                      (:decharge (let* ((input-topology-file (output-file stage-job :%DECHARGE-TOPOLOGY%))
                                                        (input-coordinate-file (output-file stage-job :%DECHARGE-COORDINATES%))
                                                        (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                     :input-topology-file input-topology-file
                                                                                     :input-coordinate-file input-coordinate-file))
                                                        (input-coordinate-file (output-file heat-job :-r)))
                                                   (push (make-ti-step morph side stage lambda-label lambda-values
                                                                       :input-topology-file input-topology-file
                                                                       :input-coordinate-file input-coordinate-file)
                                                         decharge-jobs)))
                                      (:vdw-bonded (let* ((input-topology-file input-topology-file)
                                                          (input-coordinate-file (output-file morph-side-prepare-job :-r))
                                                          (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                       :input-topology-file input-topology-file
                                                                                       :input-coordinate-file input-coordinate-file))
                                                          (input-coordinate-file (output-file heat-job :-r)))
                                                     (push (make-ti-step morph side stage lambda-label lambda-values
                                                                         :input-topology-file input-topology-file
                                                                         :input-coordinate-file input-coordinate-file)
                                                           vdw-jobs)))
                                      (:recharge (let* ((input-topology-file (output-file stage-job :%RECHARGE-TOPOLOGY%))
                                                        (input-coordinate-file (output-file stage-job :%RECHARGE-COORDINATES%))
                                                        (heat-job (make-heat-ti-step morph side stage lambda-label lambda-values
                                                                                     :input-topology-file input-topology-file
                                                                                     :input-coordinate-file input-coordinate-file))
                                                        (input-coordinate-file (output-file heat-job :-r)))
                                                   (push (make-ti-step morph side stage lambda-label lambda-values
                                                                       :input-topology-file input-topology-file
                                                                       :input-coordinate-file input-coordinate-file)
                                                         recharge-jobs)))))
                           (let ((script (make-instance 'python-script-file
                                                        :script *python-getdvdl*
                                                        :name "getdvdl"
                                                        :extension "py")))
                             (push (connect-graph
                                    (make-instance 'morph-side-stage-python-job
                                                   :morph morph :side side :stage stage :script script
                                                   :inputs (case stage
                                                              (:decharge (arguments :. (mapcar (lambda (job) (output-file job :-e)) decharge-jobs)))
                                                              (:vdw-bonded (arguments :. (mapcar (lambda (job) (output-file job :-e)) vdw-jobs)))
                                                              (:recharge (arguments :. (mapcar (lambda (job) (output-file job :-e)) recharge-jobs))))
                                                   :outputs (arguments :%STAGE-ANALYSIS% (make-instance 'morph-side-stage-file
                                                                                                      :morph morph :side side :stage stage
                                                                                                      :name "dvdl" :extension "dat"))
                                                   :makefile-clause (standard-makefile-clause "python2 getdvdl.py :%STAGE-ANALYSIS% :%INPUTS%")))
                                    stage-jobs))))
                ;; combine stage-jobs
                (let ((side-script (make-instance 'morph-side-script
                                                  :morph morph :side side
                                                  :script *combine-stages*
                                                  :name "side-script"
                                                  :extension "lisp")))
                  (push (connect-graph
                         (make-instance 'morph-side-cando-job
                                        :morph morph :side side
                                        :script side-script
                                        :inputs (arguments :. (mapcar (lambda (job) (output-file job :%STAGE-ANALYSIS%)) stage-jobs))
                                        :outputs (arguments :%SIDE-ANALYSIS% (make-instance 'morph-side-file
                                                                                          :morph morph :side side
                                                                                          :name "side-analysis" :extension "dat"))
                                        :makefile-clause (standard-cando-makefile-clause side-script :add-inputs t)))
                        analysis-jobs))))
         ;; combine analysis-jobs
         (let ((morph-script (make-instance 'morph-script
                                      :morph morph
                                      :script *combine-sides*
                                      :name "morph-script"
                                      :extension "lisp")))
           (push (connect-graph
                  (make-instance 'morph-cando-job
                                 :morph morph
                                 :script morph-script
                                 :inputs (arguments :. (mapcar (lambda (job) (output-file job :%SIDE-ANALYSIS%)) analysis-jobs))
                                 :outputs (arguments :%MORPH-ANALYSIS% (make-instance 'morph-file
                                                                                    :morph morph
                                                                                    :name "morph-analysis" :extension "dat"))
                                 :makefile-clause (standard-cando-makefile-clause morph-script :add-inputs t)))
                 morph-jobs))))
      morph-jobs)))

(defmethod generate-jobs (calculation)
  (let ((*default-pathname-defaults* (merge-pathnames (top-directory calculation) *default-pathname-defaults*)))
    (let* ((jupyter-job (make-instance 'jupyter-job))
           (am1-jobs (setup-am1-calculations jupyter-job calculation))
           (feps-precharge (make-instance 'feps-file :name "precharge")))
      (fep:save-feps calculation (node-pathname feps-precharge))
      (push (make-instance 'argument :option :feps-precharge :node feps-precharge) (outputs jupyter-job))
      (let* ((script (make-instance 'cando-script-file
                                    :name "charge"
                                    :script *cando-charge-script*))
             (feps-out (make-instance 'feps-file :name "postcharge")))
        (connect-graph
         (make-instance 'cando-job
                        :inputs (apply #'arguments
                                       :%FEPS% feps-precharge
                                       (loop for am1-job in am1-jobs
                                             for output = (output-file am1-job :-o)
                                             append (list :%INPUT% output)))
                        :outputs (arguments :%OUTPUT% feps-out)
                        :script script
                        :makefile-clause (standard-cando-makefile-clause script)))
        (let ((morph-jobs (make-script-1-leap calculation :input-feps-file feps-out)))
          ;; Do more preparation
          (generate-all-code calculation (list jupyter-job) (mapcar (lambda (job) (output-file job :%MORPH-ANALYSIS%)) morph-jobs)))
        jupyter-job))))
