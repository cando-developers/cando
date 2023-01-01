(in-package :foldamer)

(define-condition no-matching-context ()
  ((context :initarg :context :accessor context))
  (:report (lambda (condition stream)
             (format stream "no-matching-context ~a" (context condition)))))

(defclass node ()
  ((name :initarg :name :accessor name)))

(defmethod print-object ((obj node) stream)
  (print-unreadable-object (obj stream :type t :identity t)
    (format stream "~a" (name obj))))

(defclass edge ()
  ((from-node :initarg :from-node :accessor from-node)
   (to-node :initarg :to-node :accessor to-node)
   (name :initarg :name :accessor name)))

(defmethod print-object ((obj edge) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a" (from-node obj) (name obj) (to-node obj))))

(defclass dag ()
  ((nodes :initform nil :accessor nodes)
   (edges :initform nil :accessor edges)))

(defun ensure-monomer (node map)
  (let ((monomer (gethash node map)))
    (unless monomer
      (error "Could not find monomer for node ~a" node))
    monomer))

(defun oligomer-space-from-dag (foldamer dag focus-node topology-groups)
  (let ((node-to-monomer (make-hash-table))
        (oligomer-space (make-instance 'topology:oligomer-space
                                       :foldamer foldamer)))
    (loop for node in (nodes dag)
          for name = (name node)
          for names = (gethash name topology-groups)
          for monomer = (make-instance 'topology:monomer
                                       :monomers names
                                       :id name)
          do (setf (gethash node node-to-monomer) monomer)
          do (vector-push-extend monomer (topology:monomers oligomer-space)))
    (loop for edge in (edges dag)
          for source-node = (from-node edge)
          for target-node = (to-node edge)
          for source-monomer = (ensure-monomer source-node node-to-monomer)
          for target-monomer = (ensure-monomer target-node node-to-monomer)
          for coupling-name = (name edge)
          for source-plug-name = (topology:out-plug-name coupling-name)
          for target-plug-name = (topology:in-plug-name coupling-name)
          for directional-coupling = (make-instance 'topology:directional-coupling
                                                    :name coupling-name
                                                    :source-monomer source-monomer
                                                    :target-monomer target-monomer
                                                    :source-plug-name source-plug-name
                                                    :target-plug-name target-plug-name)
          do (setf (gethash source-plug-name (topology:couplings source-monomer))
                   directional-coupling
                   (gethash target-plug-name (topology:couplings target-monomer))
                   directional-coupling)
          do (vector-push-extend directional-coupling (topology:couplings oligomer-space)))
    (values oligomer-space (gethash focus-node node-to-monomer))))

(defun ensure-name (maybe-name)
  (cond
    ((and (consp maybe-name) (eq (car maybe-name) :cap))
     (cadr maybe-name))
    ((symbolp maybe-name)
     maybe-name)
    (t (error "Illegal ensure-name ~s" maybe-name))))

(defun parse-recursive (sub-tree prev-node dag)
  (cond
    ((null sub-tree))
    ((consp (car sub-tree))
     (parse-recursive (car sub-tree) prev-node dag)
     (parse-recursive (cdr sub-tree) prev-node dag))
    ((symbolp (car sub-tree))
     (let* ((plug-name (car sub-tree))
            (node (make-instance 'node :name (ensure-name (cadr sub-tree)))))
       (cond
         ((topology:in-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :name (topology:coupling-name plug-name)
                                     :from-node node
                                     :to-node prev-node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         ((topology:out-plug-name-p plug-name)
          (let ((edge (make-instance 'edge
                                     :name (topology:coupling-name plug-name)
                                     :from-node prev-node
                                     :to-node node)))
            (push node (nodes dag))
            (push edge (edges dag))))
         (t (error "Illegal plug-name ~s" plug-name)))
       (parse-recursive (cddr sub-tree) node dag)))
    (t (error "Illegal sub-tree ~s" sub-tree))))


(defun parse-for-oligomer-space (foldamer tree)
  (let ((dag (make-instance 'dag))
        (node (make-instance 'node :name (ensure-name (car tree)))))
    (push node (nodes dag))
    (parse-recursive (cdr tree) node dag)
    (oligomer-space-from-dag foldamer dag node topology::*topology-groups*)))

(defclass training-oligomer-space ()
  ((expression :initarg :expression :accessor expression)
   (oligomer-space :initarg :oligomer-space :accessor oligomer-space)
   (focus-monomer :initarg :focus-monomer :accessor focus-monomer)
   (monomer-context-matcher :initarg :monomer-context-matcher :accessor monomer-context-matcher)))

(cando:make-class-save-load
 training-oligomer-space
 :print-unreadably
 (lambda (obj stream)
   (print-unreadable-object (obj stream :type t)
     (format stream "~s" (expression obj)))))

(defclass foldamer ()
  ((topologys :initform nil :initarg :topologys :accessor topologys)
   (training-oligomer-spaces :initarg :training-oligomer-spaces :accessor training-oligomer-spaces)))

(cando:make-class-save-load
 foldamer)

(defparameter *foldamers* (make-hash-table))

(defun define-foldamer (name contexts)
  (let ((total-sequences 0)
        (topologys nil)
        (foldamer (make-instance 'foldamer)))
    (let ((training-oligomer-spaces
            (loop for context in contexts
                  collect (multiple-value-bind (oligomer-space focus-monomer)
                              (parse-for-oligomer-space foldamer context)
                            (format t "trainers ~a ~a~%" (topology:number-of-sequences oligomer-space) context)
                            (let ((topologys-in-oligomer-space (topology:topologys-in-oligomer-space oligomer-space)))
                              (loop for topology in topologys-in-oligomer-space
                                    do (pushnew topology topologys)))
                            (incf total-sequences (topology:number-of-sequences oligomer-space))
                            (let ((monomer-context-matcher (monomer-context:parse context)))
                              (make-instance 'training-oligomer-space
                                             :expression context
                                             :oligomer-space oligomer-space
                                             :focus-monomer focus-monomer
                                             :monomer-context-matcher monomer-context-matcher))))))
      (reinitialize-instance foldamer
                             :topologys topologys
                             :training-oligomer-spaces training-oligomer-spaces)
      (format t "Total trainers: ~a~%" total-sequences)
      (setf (gethash name *foldamers*)
            foldamer))))


(defun find-oligomer-for-monomer-context (foldamer monomer-context)
  (let ((training-spaces (training-oligomer-spaces foldamer)))
    (loop for training-space in training-spaces
          for oligomer-space = (oligomer-space training-space)
          for focus-monomer = (focus-monomer training-space)
          for monomer-context-matcher = (monomer-context-matcher training-space)
          for num-sequences = (topology:number-of-sequences oligomer-space)
          do (loop for num-seq below num-sequences
                   for oligomer = (topology:make-oligomer oligomer-space num-seq)
                   for match = (monomer-context:match monomer-context-matcher focus-monomer oligomer)
                   for trainer-context = (monomer-context:match-as-string match)
                   when (string= monomer-context trainer-context)
                     do (return-from find-oligomer-for-monomer-context (values oligomer focus-monomer))))))


(defun calculate-files (trainer-context &optional root-pathname)
  (let* ((data-dir (if root-pathname (merge-pathnames #P"data/" root-pathname)
                       #P"data/"))
         (output-dir (if root-pathname (merge-pathnames #P"output/" root-pathname)
                         #P"output/"))
         (input-file (make-pathname :name trainer-context :type "input" :defaults data-dir))
         (sdf-file (make-pathname :name trainer-context :type "sdf" :defaults output-dir))
         (internals-file (make-pathname :name trainer-context :type "internals" :defaults output-dir))
         (log-file (make-pathname :name trainer-context :type "log" :defaults output-dir))
         (svg-file (make-pathname :name trainer-context :type "svg" :defaults output-dir))
         (done-file (make-pathname :name trainer-context :type "done" :defaults output-dir)))
    (values input-file done-file sdf-file internals-file log-file svg-file)))

(defun generate-training-oligomers (foldamer path &key force-save)
  (ensure-directories-exist path)
  (let ((foldamer-path (merge-pathnames #P"foldamer.dat" path)))
    (cando:save-cando foldamer foldamer-path))
  (let ((training-spaces (training-oligomer-spaces foldamer))
        (index 0))
    (with-open-file (fmakefile (merge-pathnames "makefile" path) :direction :output :if-exists :supersede)
      (let ((all-done-files nil))
        (loop for training-space in training-spaces
              for oligomer-space = (oligomer-space training-space)
              for focus-monomer = (focus-monomer training-space)
              for monomer-context-matcher = (monomer-context-matcher training-space)
              for num-sequences = (topology:number-of-sequences oligomer-space)
              do (loop for num-seq below num-sequences
                       for oligomer = (topology:make-oligomer oligomer-space num-seq)
                       for match = (monomer-context:match monomer-context-matcher focus-monomer oligomer)
                       unless match
                         do (error "Bad match")
                       do (incf index)
                       do (let* ((trainer-name (monomer-context:match-as-string match)))
                            (multiple-value-bind (input-file done-file)
                                (calculate-files trainer-name path)
                              (multiple-value-bind (local-input-file local-done-file)
                                  (calculate-files trainer-name)
                                (ensure-directories-exist input-file)
                                (when (or force-save (null (probe-file input-file)))
                                  (with-open-file (fout input-file :direction :output :if-exists :supersede)
                                    (format fout "(ql:quickload :topology)~%")
                                    (format fout "(format t \"Building trainer in file ~~s~~%\" *load-pathname*)~%")
                                    (format fout "(defparameter agg (foldamer:build-trainer ~s))~%" trainer-name)
                                    (format fout "(sys:exit 0)~%")
                                    )
                                  (push local-done-file all-done-files)
                                  (format fmakefile "~a : ~a~%" local-done-file local-input-file)
                                  (format fmakefile "~a$(CLASP) -t c -f cclasp -l ~s~%~%" #\tab (namestring local-input-file))
                                  ))))))
        (format fmakefile "all: ~{~a ~}~%" all-done-files)
        ))))

(defun register-topologys (foldamer)
  (loop for topology in (topologys foldamer)
        do (cando:register-topology topology (topology:name topology))))

(defun build-trainer (trainer-context &key (steps 3) (load-pathname *load-pathname*))
  (format t "Loading foldamer.dat for building ~s~%" trainer-context)
  (let ((root-pathname (make-pathname :directory (butlast (pathname-directory load-pathname)))))
    (multiple-value-bind (input-file done-file sdf-file internals-file log-file svg-file)
        (calculate-files trainer-context root-pathname)
      (declare (ignore input-file))
      (ensure-directories-exist sdf-file)
      (with-open-file (fsdf sdf-file :direction :output
                                     :if-exists :append
                                     :if-does-not-exist :create)
        (with-open-file (flog log-file :direction :output
                                       :if-exists :append
                                       :if-does-not-exist :create)
          (ensure-directories-exist log-file)
          (let* ((fragment-pool (if (probe-file internals-file)
                                    (topology:load-fragment-pool internals-file)
                                    (make-instance 'topology:fragment-pool
                                                   :monomer-context trainer-context)))
                 (foldamer-dat-pathname (merge-pathnames #P"foldamer.dat" root-pathname))
                 (foldamer (cando:load-cando foldamer-dat-pathname)))
            (register-topologys foldamer)
            (multiple-value-bind (oligomer focus-monomer)
                (find-oligomer-for-monomer-context foldamer trainer-context)
              (format flog "Building trainer for monomer context ~a~%" (dump-local-monomer-context focus-monomer))
              (let* ((conf (topology:make-conformation oligomer :focus-monomer focus-monomer))
                     (agg (topology:aggregate conf))
                     (molecule (cando:mol agg 0))
                     (monomer-positions (topology:monomer-positions conf))
                     (focus-residue-index (gethash focus-monomer monomer-positions))
                     (ataggregate (topology:ataggregate conf))
                     (atmolecule (aref (topology:atmolecules ataggregate) 0))
                     (focus-atresidue (aref (topology:atresidues atmolecule) focus-residue-index))
                     (total-count (topology:total-count fragment-pool)))
                (with-open-file (fout svg-file :direction :output :if-exists :supersede)
                  (let* ((sketch2d (sketch2d:sketch2d molecule))
                         (svg (sketch2d:svg sketch2d)))
                    (cl-svg:stream-out fout (sketch2d:render-svg-scene svg))))
                (flet ((to-deg (rad)
                         (/ rad 0.0174533)))
                  (format flog "internals ~a~%" trainer-context)
                  (handler-case
                      (loop for count below steps
                            do (handler-case
                                   (handler-bind
                                       ((chem:minimizer-error (lambda (err)
                                                                (let ((save-filename (make-pathname :name (format nil "~a-~a" (pathname-name flog) count)
                                                                                                    :type "cando"
                                                                                                    :defaults flog)))
                                                                (warn "build-trainer - the minimizer reported: ~a - writing to ~a" err save-filename)
                                                                (invoke-restart 'save-and-skip-rest-of-minimization save-filename))))
                                     (cando:starting-geometry-with-restarts agg)))
                                 (smirnoff:missing-dihedral (err)
                                   (let ((save-filename (make-pathname :type "cando" :defaults flog)))
                                     (format flog "Missing dihedral ~a - saving molecule to ~s~%" err save-filename)
                                     (cando:save-cando (smirnoff:molecule err) save-filename))
                                   (signal err)))
                            do (when (not (topology:bad-geometry-p agg))
                                 (topology::copy-atom-positions-into-joints conf)
                                 (topology::update-joint-tree-internal-coordinates conf)
                                 (let* ((internals (loop for joint across (topology:joints focus-atresidue)
                                                         for name = (kin:joint/name joint)
                                                         collect (cond
                                                                   ((typep joint 'kin:jump-joint)
                                                                    (make-instance 'topology:jump-internal
                                                                                   :name name))
                                                                   ((typep joint 'kin:complex-bonded-joint)
                                                                    (let ((distance (kin:bonded-joint/get-distance joint))
                                                                          (angle (kin:bonded-joint/get-theta joint))
                                                                          (dihedral (kin:bonded-joint/get-phi joint)))
                                                                      (make-instance 'topology:complex-bonded-internal
                                                                                     :name name
                                                                                     :bond distance
                                                                                     :angle angle
                                                                                     :dihedral dihedral)))
                                                                   ((typep joint 'kin:bonded-joint)
                                                                    (let ((distance (kin:bonded-joint/get-distance joint))
                                                                          (angle (kin:bonded-joint/get-theta joint))
                                                                          (dihedral (kin:bonded-joint/get-phi joint)))
                                                                      (make-instance 'topology:bonded-internal
                                                                                     :name name
                                                                                     :bond distance
                                                                                     :angle angle
                                                                                     :dihedral dihedral)))
                                                                   (t (format flog "unknown-joint ~a ~a~%" name joint)))))
                                        (fragment-internals (make-instance 'topology:fragment-internals
                                                                           :index total-count
                                                                           :internals internals)))
                                   (let ((seen-index (topology:seen-fragment-internals fragment-pool fragment-internals)))
                                     (if (topology:good-fragment-internals fragment-internals)
                                         (if (not seen-index)
                                             (progn
                                               (sdf:write-sdf-stream agg fsdf)
                                               (format flog "Saving new conformation ~a~%" total-count)
                                               (topology:dump-fragment-internals fragment-internals flog)
                                               (push fragment-internals (topology:fragments fragment-pool)))
                                             (progn
                                               (format flog "Ignoring conformation ~a - seen before at ~a~%" total-count seen-index)
                                               (topology:dump-fragment-internals fragment-internals flog)
                                               ))
                                         (progn
                                           (format flog "Ignoring conformation ~a - bad geometry~%" total-count)
                                           (topology:dump-fragment-internals fragment-internals flog)
                                           )))))
                            do (incf total-count)
                            )
                    (error (err)
                      (format flog "There was an error ~a- terminating~%" err) ;
                      (return-from build-trainer nil) ;
                      ))
                  )
                (setf (topology:total-count fragment-pool) total-count)
                (topology:save-fragment-pool fragment-pool internals-file)
                (with-open-file (fout done-file :direction :output :if-exists :supersede)
                  (format fout "finished-steps ~a~%" total-count)
                  (let ((fragment-pool-length (length (topology:fragments fragment-pool))))
                    (format fout "hits ~a~%" fragment-pool-length)
                    (format fout "hit-to-total-ratio ~5,4f~%" (/ (float fragment-pool-length 1.0d0) (float total-count 1.0d0)))))))))))))


(defun prepare-to-build-trainer (&key (smirnoff #P"~/Development/openff-sage/inputs-and-results/optimizations/vdw-v1/forcefield/force-field.offxml" ))
  (leap:load-smirnoff-params smirnoff))

(defun assemble-fragment-pool-map (filename)
  (let* ((foldamer-filename (merge-pathnames #P"foldamer.dat" filename))
         (foldamer (cando:load-cando foldamer-filename))
         (fragment-pool-map (make-instance 'topology:fragment-pool-map))
         (training-oligomer-spaces (training-oligomer-spaces foldamer))
         (fragment-pool-map (make-instance 'topology:fragment-pool-map)))
    (loop for training-space in training-oligomer-spaces
          for oligomer-space = (oligomer-space training-space)
          for focus-monomer = (focus-monomer training-space)
          for monomer-context-matcher = (monomer-context-matcher training-space)
          for num-sequences = (topology:number-of-sequences oligomer-space)
          do (loop for num-seq below num-sequences
                   for oligomer = (topology:make-oligomer oligomer-space num-seq)
                   for match = (monomer-context:match monomer-context-matcher focus-monomer oligomer)
                   unless match
                     do (error "Bad match")
                   do (let* ((trainer-name (monomer-context:match-as-string match)))
                        (format t "Loading internals for ~a~%" trainer-name)
                        (finish-output)
                        (multiple-value-bind (input-file done-file sdf-file internals-file)
                            (calculate-files trainer-name filename)
                          (declare (ignore input-file done-file sdf-file))
                          (let ((fragment-pool (topology:load-fragment-pool internals-file)))
                            (setf (gethash trainer-name (topology:monomer-context-to-fragment-pool fragment-pool-map))
                                  fragment-pool))))))
    fragment-pool-map))

(defun recursive-dump-local-monomer-context (monomer prev-monomer depth)
  (when (>= depth 0)
    (let (rest)
      (maphash (lambda (key coupling)
                 (let ((other-monomer (topology:other-monomer coupling monomer)))
                   (unless (eq other-monomer prev-monomer)
                     (let ((other-tree (recursive-dump-local-monomer-context other-monomer monomer (1- depth))))
                       (when other-tree
                         (push (cons key other-tree) rest))))
                   ))
                 (topology:couplings monomer))
               (let ((sorted-rest (sort rest #'string< :key (lambda (val) (string (car val))))))
                 (list* (topology:monomers monomer) (mapcar (lambda (alist) (list (car alist) (cdr alist))) sorted-rest))))))

(defun dump-local-monomer-context (monomer)
  (recursive-dump-local-monomer-context monomer nil 2))

(defun foldamer-monomer-context (monomer oligomer-or-space foldamer)
  "Return the monomer-context string for the monomer in the oligomer that is made of the foldamer"
  (loop for training-oligomer-space in (training-oligomer-spaces foldamer)
        for monomer-context-matcher = (monomer-context-matcher training-oligomer-space)
        for match = (monomer-context:match monomer-context-matcher monomer oligomer-or-space)
        when match
          do (progn
               (return-from foldamer-monomer-context (monomer-context:match-as-string match))))
  (error 'no-matching-context
         :context (dump-local-monomer-context monomer)))


(defun check-foldamer-describes-oligomer-space (foldamer oligomer-space)
  (loop for monomer across (topology:monomers oligomer-space)
        for monomer-context = (foldamer-monomer-context monomer oligomer-space foldamer)
        when monomer-context
          do (format t "monomer-context ~a~%   is matched by ~a~%" (recursive-dump-local-monomer-context monomer nil 1) monomer-context)
        )
  (format t "foldamer completely describes oligomer-space~%"))
