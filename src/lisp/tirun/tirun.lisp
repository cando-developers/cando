;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :tirun)

(defparameter *system* nil "The default system is NIL")


(defun make-tirun ()
  (make-instance 'tirun:tirun-calculation))

(defun add-receptor(tirun receptor)
  (pushnew receptor (receptors tirun)))

(defun version ()
  *version*)

(defun cartesian-product (l) 
  (cond ((null l) nil) 
        ((null (cdr l)) 
         (mapcar #'list (car l))) 
        (t (mapcan #'(lambda (x) (mapcar #'(lambda (y) (cons y x)) (car l))) 
                   (cartesian-product (cdr l))))))

(defvar *tirun-time* 0)

(defun incf-tirun-time ()
  (incf *tirun-time*))

(defun print-object-readably-with-slots (obj stream)
  (format stream "#$(~s " (class-name (class-of obj)))
  (loop for slot in (clos:class-slots (class-of obj))
        for slot-name = (clos:slot-definition-name slot)
        for initargs = (clos:slot-definition-initargs slot)
        if (and (car initargs) (slot-boundp obj slot-name))
          do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
  (format stream ") "))

(defclass job-graph ()
  ((nodes :initform nil :initarg :nodes :accessor nodes)
   (morphs :initform nil :initarg :morphs :accessor morphs)))

(defmethod print-object ((obj job-graph) stream)
  (if *print-readably*
      (progn
        (print-object-readably-with-slots obj stream))
      (print-unreadable-object (obj stream)
        (format stream "~a" (class-name (class-of obj))))))

(defgeneric generate-lambda-values (system &key divisions))

(defmethod generate-lambda-values ((system null) &key (divisions 10))
  "The default method for generate-lambda-values - 0.0 ... 1.0 with equal sized divisions"
  (loop for lambda-int from 0 to divisions
        for lambda-val = (/ (float lambda-int) (float divisions))
        collect lambda-val))

        
(defclass tirun-morph ()
  ((amber-version :initform 18 :initarg :amber-version :accessor amber-version)
   (source :initarg :source :accessor source)
   (target :initarg :target :accessor target)
   (morph-mask :initarg :morph-mask :accessor morph-mask)
   (stages :initarg :stages :initform 1 :accessor stages)
   (lambda-values :type list :initarg :lambda-values :initform (generate-lambda-values *system*) :accessor lambda-values)))

(defmethod print-object ((obj tirun-morph) stream)
  (if *print-readably*
      (print-object-readably-with-slots obj stream)
      (print-unreadable-object (obj stream)
        (format stream "~a" (class-name (class-of obj))))))

(defun morph-string (tirun-morph)
  (let ((source-name (name (source tirun-morph)))
        (target-name (name (target tirun-morph))))
    (format nil "~a-~a" source-name target-name)))

(defclass simple-tirun-structure ()
  ((name :initarg :name :accessor name)
   (drawing :initarg :drawing :accessor drawing)
   (core-residue :initarg :core-residue :accessor core-residue)
   (core-residue-name :initarg :core-residue-name :accessor core-residue-name)
   (core-atoms :initarg :core-atoms :accessor core-atoms)
   (molecule :initarg :molecule :accessor molecule)
   (atom-order :initarg :atom-order :accessor atom-order)
   (net-charge :initform 0.0  :initarg :net-charge :accessor net-charge)
   (am1-charges :initarg :am1-charges :accessor am1-charges)
   (am1-bcc-charges :initarg :am1-bcc-charges :accessor am1-bcc-charges)))

(defmethod print-object ((object simple-tirun-structure) stream)
  (print-unreadable-object (object stream)
    (format stream "~a ~a" (class-name (class-of object)) (name object))))

(defclass tirun-structure (simple-tirun-structure)
  ((side-chain-residue-names :initarg :side-chain-residue-names :accessor side-chain-residue-names)
   (side-chain-atoms :initarg :side-chain-atoms :accessor side-chain-atoms)))


(defmethod print-object ((obj simple-tirun-structure) stream)
  (if *print-readably*
      (print-object-readably-with-slots obj stream)
      (print-unreadable-object (obj stream)
        (format stream "~a ~a" (class-name (class-of obj)) (string (name obj))))))

(defun mol2-safe-residue-name (calculation name)
  name)

(defun mol2-safe-atom-name (calculation name)
  name)

(defvar *default-calculation-settings*
  '(("Simulation Input Parameters"
     ((:%DECHARGE-RECHARGE-TI-IN.NSTLIM% "Decharge/recharge nstlim" (3000000 3000)) ; if *testing* must be >= 3000
      (:%PREPARE-MIN-IN.MAXCYC% "Number of steepest decent minimization steps" (50000))
      (:%VDW-TI-IN.NSTLIM% "Vdw nstlim" (3000000 3000)) ; if *testing* must be >= 3000
      (:%PREPARE-HEAT-IN.NSTLIM% "PREPARE-HEAT-IN.NSTLIM" (100000 3000))
      (:%PREPARE-PRESS-IN.NSTLIM% "PREPARE-PRESS-IN.NSTLIM" (2500000 3000))
      (:%DECHARGE-RECHARGE-HEAT-IN.NSTLIM% "DECHARGE-RECHARGE-HEAT-IN.NSTLIM" (100000 3000))
      (:%VDW-HEAT-IN.NSTLIM% "VDW-HEAT-IN.NSTLIM" (100000 3000))
      (:%PREPARE-PRESS-IN-DT% "NPT with SHAKE timestep (fs)" 0.002)
      (:%PREPARE-REST-WT% "Equilibration harmonic restraint weight (kcal/mol/A^2)" 5.00)
      (:%DT% "TI without SHAKE timestep (fs) [DT]" 0.001)
      (:%NTC% "SHAKE-on (2) SHAKE-off (1) [NTC]" 2)
      (:%NTF% "Force calculation-on (1) -off (2) [NTF]" 1)
      (:%TEMP0% "Constant temperature" 300.0 )
      (:%SCALPHA% "Softcore alpha in potential for 1 and 3 step TI" 0.5)
      (:%SCBETA% "Softcore beta in Coulomb eq. for 1 step TI" 12.0) 
      ))))

(defun default-calculation-settings (&key (default-settings *default-calculation-settings*) testing)
  (let (all-settings)
    (loop for batch in default-settings
          for name = (first batch)
          for settings = (second batch)
          do (loop for setting in settings
                   for setting-key = (first setting)
                   for setting-name = (second setting)
                   for setting-value = (third setting)
                   do (push (cons setting-key (if (consp setting-value)
                                                  (if testing
                                                      (second setting-value)
                                                      (first setting-value))
                                                  setting-value))
                            all-settings)))
    all-settings))


(defparameter *testing-lambdas* (loop for idx from 0 to 5 collect (/ (float idx) 5.0)))

(defclass calculation ()
  ((receptors :initform nil :initarg :receptors :accessor receptors)
   (ligands :initarg :ligands :accessor ligands)
   (mask-method :initform :default :initarg :mask-method :accessor mask-method)
   (core-topology :initform nil :initarg :core-topology :accessor core-topology)
   (side-topologys :initform nil :initarg :side-topologys :accessor side-topologys)
   (jobs :initform (make-instance 'job-graph) :initarg :jobs :accessor jobs)
   (ti-stages :initarg :ti-stages :initform 1 :accessor ti-stages)
   (ti-lambdas :initarg :ti-lambdas :initform 'generate-lambda-values :accessor ti-lambdas)
   (settings :initform (default-calculation-settings) :initarg :settings :accessor settings)
   (top-directory :initform (make-pathname :directory (list :relative "jobs"))
                  :initarg :top-directory :accessor top-directory )
   (stage :initform 0 :initarg :stage :accessor stage)
   (residue-name-to-pdb-alist :initform nil :initarg :residue-name-to-pdb-alist :accessor residue-name-to-pdb-alist)
   ))


(defmethod print-object ((obj calculation) stream)
  (if *print-readably*
      (progn
        (format stream "#.(make-instance 'tirun::calculation ")
        (loop for slot in (clos:class-slots (find-class 'calculation))
              for slot-name = (clos:slot-definition-name slot)
              for initargs = (clos:slot-definition-initargs slot)
              if (and (car initargs) (slot-boundp obj slot-name))
                do (format stream "~s ~s " (car initargs) (slot-value obj slot-name)))
        (format stream ") "))
      (print-unreadable-object (obj stream)
        (format stream "tiruns"))))

(defun find-morph-with-name (name calculation)
  (loop for morph in (morphs (jobs calculation))
        when (string= (morph-string morph) (string name))
          do (return-from find-morph-with-name morph))
  (error "Could not find morph with name ~s in calculation with morphs: ~s" name (mapcar (lambda (m) (name m)) (morphs (jobs calculation)))))

(defparameter *top-guard* nil)


(defparameter *write-files* nil)
(defmacro with-top-directory ((calculation) &body body)
  `(let ((*write-files* t))
     ,@body))

(defmethod make-load-form ((object calculation) &optional env)
  (make-load-form-saving-slots object))


(defun am1-file-name (tirun type)
  (make-pathname :directory '(:relative "am1bcc") :name (format nil "am1-~a" (string (name tirun))) :type type :defaults *default-pathname-defaults*))

(defun am1-calculation-complete (tirun)
  "Return T if the AM1 charges for this TIRUN exist, are up to date and are available.
Otherwise return NIL."
  (let ((in-filename (am1-file-name tirun "in"))
        (out-filename (am1-file-name tirun "out"))
        (atom-order (atom-order tirun)))
    (cond
      ((and (probe-file out-filename)
            (< (file-write-date in-filename) (file-write-date out-filename))
            (charges:read-am1-charges out-filename atom-order)))
      (t nil))))

(defun setup-am1-calculations (jupyter-job calculation &key (maxcyc 9999))
  (with-top-directory (calculation)
    (let* ((ligands (ligands calculation))
           (jobs (jobs calculation))
           (morphs (morphs jobs))
           (only-ligands-in-morphs (let (ligs)
                                     (loop for morph in morphs
                                           do (pushnew (source morph) ligs)
                                           do (pushnew (target morph) ligs))
                                     ligs)))
      (ensure-jobs-directories-exist *default-pathname-defaults*)
      (let* (outputs
             (result (loop for ligand in only-ligands-in-morphs
                           for in-file = (make-instance 'sqm-input-file :name (name ligand))
                           for order-file = (make-instance 'sqm-atom-order-file :name (name ligand))
                           for molecule-charge = (let ((charge 0.0))
                                                   (cando:do-atoms (atom (molecule ligand))
                                                     (incf charge (chem:get-charge atom)))
                                                   charge)
                           for atom-order = (let (atom-order)
                                              (write-file-if-it-has-changed
                                               (node-pathname in-file)
                                               (with-output-to-string (sout) 
                                                 (setf atom-order (charges:write-sqm-calculation sout (molecule ligand) :maxcyc maxcyc :qm-charge molecule-charge))))
                                              atom-order)
                           for job = (make-ligand-sqm-step ligand :sqm-input-file in-file)
                           do (setf (atom-order ligand) atom-order)
                           do (push (make-instance 'argument :option :output :node in-file) outputs)
                           collect job)))
        (setf (outputs jupyter-job) outputs)
        (connect-graph jupyter-job)
        (mapc #'connect-graph result)
        result))))

(defun check-am1-calculations (calculation)
  (warn "This doesn't work right now - we need source tracking to identify the form that is the problem")
  #+(or)
  (let* ((status-entries (loop for tirun in (ligands calculation)
                               for name = (string (name tirun))
                               for status = (am1-calculation-complete tirun)
                               collect (list name status)))
         (sorted-entries (sort status-entries #'string< :key #'car)))
    (error "Fix checking the am1-calculations")
    #+(or)(jupyter:html (cl-markup:markup
                           (:table :border 3
                                   :cellpadding 4
                                   (loop for i below (length sorted-entries) by 5
                                         collect (cl-markup:markup
                                                  (:tr :align "right"
                                                       (loop for j from i below (+ i 5)
                                                             for entry = (elt sorted-entries j)
                                                             for name = (car entry)
                                                             for status = (cadr entry)
                                                             when entry
                                                               collect (cl-markup:markup
                                                                        (:td :bgcolor (if status "white" "red")
                                                                             (format nil "~a ~a" name (if status "Done" "Incomplete")))))))))))))



(defun read-am1-charges (calculation)
  (with-top-directory (calculation)
	(let* ((count 0)
           (jobs (jobs calculation))
           (morphs (morphs jobs))
           (only-ligands-in-morphs (let (ligs)
                                     (loop for morph in morphs
                                           do (pushnew (source morph) ligs)
                                           do (pushnew (target morph) ligs))
                                     ligs)))
      (loop for ligand in only-ligands-in-morphs
            for sqm-out = (make-instance 'sqm-output-file :name (name ligand))
            do (let ((am1-charge (charges::read-am1-charges (node-pathname sqm-out) (tirun::atom-order ligand))))
                 (setf (am1-charges ligand) am1-charge))
            do (format t "file ~a done~%" (name ligand))
            do (incf count)))))

(defun balance-charges (atom-charge-hash-table net-charge)
  (let ((charge-sum 0.0)
        (num-charges 0))
    (maphash (lambda (atm charge)
               (incf charge-sum charge)
               (incf num-charges))
             atom-charge-hash-table)
    (let ((charge-part (/ (- charge-sum net-charge) num-charges))
          (new-charge-total 0.0))
      (maphash (lambda (atm charge)
                 (let ((new-charge (- charge charge-part)))
                   (setf (gethash atm atom-charge-hash-table) new-charge)
                   (incf new-charge-total new-charge)))
               atom-charge-hash-table)
      (when (> (abs new-charge-total) 1.0e-5)
        (warn "After balance-charges the new-charge-total of the molecule is ~g~%" new-charge-total)))))

(defun calculate-am1-bcc-charges (calculation)
   (let* ((jobs (jobs calculation))
	 (morphs (morphs jobs))
	 (only-ligands-in-morphs (let (ligs)
				   (loop for morph in morphs
					 do (pushnew (source morph) ligs)
					 do (pushnew (target morph) ligs))
			           ligs)))
  (loop for tirun in only-ligands-in-morphs
  	for am1-charges = (am1-charges tirun)
        for bcc-corrections = (charges::calculate-bcc-corrections (tirun::molecule tirun))
        for am1-bcc-charges = (charges::combine-am1-bcc-charges (am1-charges tirun) bcc-corrections)
        do (balance-charges am1-bcc-charges (net-charge tirun))
           (setf (am1-bcc-charges tirun) am1-bcc-charges)
           (maphash (lambda (atm charge)
                      (chem:set-charge atm charge))
                    am1-bcc-charges))))


(defun build-complex (calculation ligand-name &optional (receptor-index 0))
  (let ((receptor (elt (receptors calculation) receptor-index))
        (ligand (find ligand-name (ligands calculation) :test #'eq :key #'name)))
    (let ((agg (chem:matter-copy receptor)))
      (cando:remove-molecules agg 'cando:ligand)
      (chem:add-molecule agg (tirun:molecule ligand))
      agg)))


(defun relax-receptor-around-ligand (receptor-aggregate ligand-molecule &key (weight 100.0))
  "Relax the receptor around the ligand."
  (cando:do-atoms (atom receptor-aggregate)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor)))
  (chem:add-matter receptor-aggregate ligand-molecule)
  (cando:do-atoms (atom ligand-molecule)
    (let ((restraint-anchor (core:make-cxx-object
                             'chem:restraint-anchor
                             :atom atom
                             :position (chem:get-position atom)
                             :weight weight)))
      (chem:clear-restraints atom)
      (chem:add-restraint atom restraint-anchor))))

(defun pattern-atoms (smarts structure)
  "Given a smarts pattern and a structure, return a hash-table of tags to atoms"
  (cando:do-atoms (atom structure)
    (let (match)
      (when (setf match (chem:matches smarts atom))
        (return-from pattern-atoms (chem:tags-as-hashtable match))))))

(defclass ti-mask ()
  ((equivalent-atom-names :initarg :equivalent-atom-names :accessor equivalent-atom-names)
   (source :initarg :source :accessor source)
   (source-timask-residue-index :initarg :source-timask-residue-index :accessor source-timask-residue-index)
   (source-scmask-atom-names :initarg :source-scmask-atom-names :accessor source-scmask-atom-names)
   ;;   (source-timask :initarg :source-timask :accessor source-timask)
   ;;   (source-scmask :initarg :source-scmask :accessor source-scmask)
   (target :initarg :target :accessor target)
   (target-timask-residue-index :initarg :target-timask-residue-index :accessor target-timask-residue-index)
   (target-scmask-atom-names :initarg :target-scmask-atom-names :accessor target-scmask-atom-names)
   ;;   (target-timask :initarg :target-timask :accessor target-timask)
   ;;   (target-scmask :initarg :target-scmask :accessor target-scmask)
   ))

(defmethod print-object ((obj ti-mask) stream)
  (if *print-readably*
      (progn
        (print-object-readably-with-slots obj stream))
      (print-unreadable-object (obj stream)
        (format stream "~a" (class-name (class-of obj))))))


(defgeneric calculate-masks (source-structure target-structure kind)
  (:documentation "Use the source and target to calculate the masks. The kind argument allows the 
user to define new ways to calculate masks"))


(defun unique-residue-with-name (molecule name)
  (let (residues)
    (cando:do-residues (res molecule)
      (when (eq (chem:get-name res) name)
        (push res residues)))
    (if (= (length residues) 1)
        (first residues)
        (error "Only one residue in ~s can have the name ~s - found ~s" molecule name residues))))
    

(defmethod calculate-masks (source target (method (eql ':default)))
  "Really simple mask - the core and mutated residues are timask and the mutated residues are the scmask.
METHOD controls how the masks are calculated"
  (let ((source-core-residue (unique-residue-with-name (molecule source) (core-residue-name source)))
        (target-core-residue (unique-residue-with-name (molecule target) (core-residue-name target)))
        (source-softcore-atoms (side-chain-atoms source))
        (target-softcore-atoms (side-chain-atoms target))
        (combined-aggregate (cando:combine (molecule source) (molecule target))))
    (let ((sequenced-residues (chem:map-residues 'vector #'identity combined-aggregate))) ; vector of residues in order they will appear in topology file
      (let* ((source-timask-residue-index (1+ (position source-core-residue sequenced-residues)))
            (source-scmask-atom-names (loop for atom in source-softcore-atoms
                                            collect (chem:get-name atom)))
            (target-timask-residue-index (1+ (position target-core-residue sequenced-residues)))
            (target-scmask-atom-names (loop for atom in target-softcore-atoms
                                            collect (chem:get-name atom)))
            (source-core-atom-names (let (names)
                                      (cando:do-atoms (atom source-core-residue)
                                        (unless (member (chem:get-name atom) source-scmask-atom-names)
                                          (push (chem:get-name atom) names)))
                                      names))
            (target-core-atom-names (let (names)
                                      (cando:do-atoms (atom target-core-residue)
                                        (unless (member (chem:get-name atom) target-scmask-atom-names)
                                          (push (chem:get-name atom) names)))
                                      names)))
        (unless (and (= (length source-core-atom-names)
                        (length target-core-atom-names))
                     (every (lambda (name)
                              (member name target-core-atom-names))
                            source-core-atom-names))
          (error "The source-core-atom-names ~s must match the target-core-atom-names ~s - or you need to define equivalent pairs.  source-softcore-atoms ~a   target-softcore-atoms ~a"
                 source-core-atom-names
                 target-core-atom-names
                 source-softcore-atoms
                 target-softcore-atoms))
        (make-instance 'ti-mask
                       ;; core atoms by the default method must have the same names in
                       ;; source and target
                       :equivalent-atom-names (mapcar 'cons source-core-atom-names source-core-atom-names)
                       :source source
		       :source-timask-residue-index source-timask-residue-index
		       :source-scmask-atom-names source-scmask-atom-names
                       ;; :source-timask (format nil ":~d" source-timask-residue-index)
                       ;; :source-scmask (format nil ":~d@~{~d~^,~}" source-timask-residue-index source-scmask-atom-names)
                       :target target
		       :target-timask-residue-index target-timask-residue-index
		       :target-scmask-atom-names target-scmask-atom-names
                       ;; :target-timask (format nil ":~d" target-timask-residue-index)
                       ;; :target-scmask (format nil ":~d@~{~d~^,~}" target-timask-residue-index target-scmask-atom-names)
		       )))))

(defun mask-substitutions (mask &optional stage)
  (unless (typep stage '(member nil :vdw-bonded :decharge :recharge))
    (error "stage ~s must be one of nil :vdw-bonded :decharge :recharge" stage))
  (list
   (cons :%TIMASK1% (format nil ":~d" (source-timask-residue-index mask)))
   (cons :%SCMASK1% (format nil ":~d@~{~d~^,~}" (source-timask-residue-index mask) (source-scmask-atom-names mask)))
   (cons :%TIMASK2% (format nil ":~d" (target-timask-residue-index mask)))
   (cons :%SCMASK2% (format nil ":~d@~{~d~^,~}" (target-timask-residue-index mask) (target-scmask-atom-names mask)))))

(defun crgmask-substitutions (mask stage)
  (unless (typep stage '(member nil :vdw-bonded :decharge :recharge))
    (error "stage ~s must be one of nil :vdw-bonded :decharge :recharge" stage))
  (cond
    ((eq stage :decharge)
     (List (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}" (target-timask-residue-index mask) (source-scmask-atom-names mask)))
           (cons :%IFSC% "0")))
    ((eq stage :vdw-bonded)
     (list (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}|:~d@~{~d~^,~}"
				    (source-timask-residue-index mask) (source-scmask-atom-names mask)
				    (target-timask-residue-index mask) (target-scmask-atom-names mask)))
           (cons :%IFSC% "1")))
    ((eq stage :recharge)
     (list (cons :%CRGMASK% (format nil ":~d@~{~d~^,~}" (source-timask-residue-index mask) (target-scmask-atom-names mask)))
           (cons :%IFSC% "0")))
    (t nil)))



(defun average-core-atom-positions (source target equivalent-atom-names)
  (loop for equiv-pair in equivalent-atom-names
        for source-atom-name = (car equiv-pair)
        for target-atom-name = (cdr equiv-pair)
        for source-atom = (or (find source-atom-name (core-atoms source) :test (lambda (src-atom-name src-atom)
                                                                                 (eq src-atom-name
                                                                                     (chem:get-name src-atom))))
                              (error "Could not find source atom ~a in ~a" source-atom-name (core-atoms source)))
        for target-atom = (or (find target-atom-name (core-atoms target) :test (lambda (tgt-atom-name tgt-atom)
                                                                                 (eq tgt-atom-name
                                                                                     (chem:get-name tgt-atom))))
                              (error "Could not find target atom ~a in ~a" target-atom-name (core-atoms target)))
        for average-pos = (geom:v* (geom:v+ (chem:get-position source-atom)
                                            (chem:get-position target-atom))
                                   0.5)
        for source-delta = (geom:vlength (geom:v- (chem:get-position source-atom)
                                                  average-pos))
        for target-delta = (geom:vlength (geom:v- (chem:get-position target-atom)
                                                  average-pos))
        do (when (or (> source-delta 0.1) (> target-delta 0.1))
             (warn "The source atom ~s and target atom ~s are more than 0.1 A away from the average ~s~%" source-atom target-atom average-pos))
           (chem:set-position source-atom average-pos)
           (chem:set-position target-atom average-pos)))


(defun load-chem-draw-tirun (filename)
  (handler-bind ((warning #'muffle-warning))
    (with-open-file (fin (open filename :direction :input))
      (chem:make-chem-draw fin nil nil))))

#+(or)
(defun scale-molecule (mol)
  (let ((total-length 0.0)
        (num-bonds 0))
    (loop for ligand in (ligands calculation)
          for mol = (drawing ligand)
          do (chem:map-bonds nil (lambda (a1 a2 bo)
                                   (let ((length (geom:vlength (geom:v- (chem:get-position a1)
                                                                        (chem:get-position a2)))))
                                     (when (> length 0.1)
                                       (incf num-bonds)
                                       (incf total-length length))))
                             mol))
    (let ((scale (/ 1.5 (/ total-length num-bonds))))
      (loop for ligand in (ligands calculation)
            for mol = (drawing ligand)
            do (cando:do-atoms (atm mol)
                 (chem:set-position atm (geom:v* (chem:get-position atm) scale)))))))

(defun validate-atom-types (matter)
  (cando:do-atoms (atom matter)
    (let ((type (chem:get-type atom)))
      (when (or (null type) (eq :du type))
        (error "Illegal type ~a for ~a in ~a" type atom matter)))))

(defun check-calculation-atom-types (calculation)
  (loop for ligand in (ligands calculation)
        for molecule = (molecule ligand)
        do (leap:assign-atom-types molecule)
        do (validate-atom-types molecule)
           (format t "ligand ~a has valid atom types~%" molecule))
  (loop for receptor in (receptors calculation)
        do (validate-atom-types receptor)
           (format t "Receptor ~a has valid atom types~%" receptor)))

(defun ensure-jobs-directories-exist (pathname)
        (ensure-directories-exist pathname))
