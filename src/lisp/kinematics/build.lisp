(in-package :cando-user)

(defparameter kin.work::*build-conformation*
  (format nil "~{~s~%~}~%"
          (list
           `(ql:quickload :kinematics)
           `(format t "Starting build-conformation~%")
           `(kin.work:build-conformation)
           `(sys:quit)
           )))
  
(in-package :kin.work)

(defun separate-arguments (args)
  (let ((results (make-hash-table :test 'equal))
        (sep-pos (loop for arg in args
                       for pos from 0
                       when (string= arg "--")
                         do (return pos))))
    (loop for cur = (nthcdr (1+ sep-pos) args) then (cddr cur)
          for key = (car cur)
          for value = (cadr cur)
          until (null cur)
          do (setf (gethash key results) value))
    results))


(defun conformation-was-seen (results conf-internal cycle logout)
  (when (kin:entries results)
    (let ((pos (position conf-internal (kin:entries results)
                         :test 'kin:internals-match)))
      (if pos
          (progn
            #+(or)(format t "Cycle: ~a - conformation was seen in history at pos ~a of ~a existing conformations~%"
                    cycle pos (length (kin:entries results)))
            #+(or)(break "Check pos: ~a" pos) 
            t
            )
          (progn
            (format logout "Cycle: ~a - there is a new conformation with a history of ~a existing conformations~%" cycle (length (kin:entries results)))
            #+(or)(break "A new conformation was found")
            nil)
          ))))

(defparameter *maxcycles* 10000)
(defparameter *max-searches-since-last-new-conformation* 10)
(defun outer-search-conformations (args trainer results)
  (let* ((trainer-fn (kin:trainer-filename trainer))
         (log-filename (pathname (gethash "--log" args)))
         (allsdf-filename (pathname (gethash "--allsdf" args)))
         (uniquesdf-filename (pathname (gethash "--uniquesdf" args)))
         (sketch-fn (pathname (gethash "--sketch" args)))
         (conformation (kin:conformation trainer))
         (focus-residue-sequence-number (kin:focus-residue-sequence-number trainer))
         (superposer (core:make-cxx-object 'chem:superpose-engine))
         (aggregate (kin:aggregate conformation))
         (molecule (cando:mol aggregate 0))
         (focus-residue (cando:res molecule focus-residue-sequence-number))
         (superpose-atoms (chem:map-atoms 'list #'identity focus-residue))
         (cycle 0)
         (bad-geometry 0)
         (searches-since-last-new-conformation 0)
         (new-conformations 0)
         first-conformation-seen
         )
    (let* ((sketch (sketch2d:sketch2d (chem:content-at aggregate 0)))
           (svg (sketch2d:svg sketch))
           (as-string (sketch2d:render-svg-to-string svg)))
      (with-open-file (fout sketch-fn :direction :output)
        (princ as-string fout)))
    (format t "Opening log file: ~a~%" log-filename)
    (with-open-file (logout log-filename :direction :output)
      (with-open-file (allsdfout allsdf-filename :direction :output)
        (with-open-file (uniquesdfout uniquesdf-filename :direction :output)
          (tagbody
           top
             (incf cycle)
             (when (> cycle *maxcycles*)
               (format logout "Terminating because we hit *maxcycles*[~a]~%"
                       *maxcycles*)
               (go done))
             (format logout "Cycle: ~a/~a since-last-new: ~a/~a results: ~a ~a~%"
                     cycle
                     *maxcycles*
                     searches-since-last-new-conformation
                     *max-searches-since-last-new-conformation*
                     (length (kin:entries results))
                     trainer-fn)
             (finish-output logout)
             (cando:jostle aggregate)
             (energy:minimize aggregate :print-intermediate-results nil)
             (when (bad-geometry-p aggregate)
               (incf bad-geometry)
               (go top))
             ;; superpose atoms
             (if (null first-conformation-seen)
                 (let* ((center-of-residue (chem:geometric-center focus-residue))
                        (transform (geom:make-m4-translate (geom:v* center-of-residue -1.0))))
                   (chem:apply-transform-to-atoms aggregate transform)
                   (let ((fixed-coords (chem:make-simple-vector-coordinate-from-atom-list superpose-atoms)))
                     (chem:set-fixed-all-points superposer fixed-coords))
                   (setf first-conformation-seen t))
                 (let* ((moveable-coords (chem:make-simple-vector-coordinate-from-atom-list superpose-atoms)))
                   (chem:set-moveable-all-points superposer moveable-coords)
                   (let ((transform (chem:superpose superposer)))
                     (chem:apply-transform-to-atoms aggregate transform))))
             (sdf:write-sdf-stream aggregate allsdfout :name (format nil "all-~a" cycle))
             (kin:copy-atom-positions-into-joints conformation)
             (kin:update-joint-tree-internal-coordinates conformation)
             (let ((internals (kin:extract-rotamer trainer :index cycle)))
               (when (> searches-since-last-new-conformation *max-searches-since-last-new-conformation*)
                 (format logout "Terminating because we hit *max-searches-since-last-new-conformation*[~a] test~%"
                         *max-searches-since-last-new-conformation*)
                 (go done))
               (when (conformation-was-seen results internals cycle logout)
                 (incf searches-since-last-new-conformation)
                 (go top))
               (sdf:write-sdf-stream aggregate uniquesdfout :name (format nil "unique-~a" cycle))
               (incf new-conformations)
               (setf searches-since-last-new-conformation 0)
               (push internals (kin:entries results))
               )
             (go top)
           done
             ))
        results))))

(defun build-conformation-args (args)
  (let* ((hargs (separate-arguments args))
         (grammar-fn (gethash "--grammar" hargs))
         (trainer-fn (gethash "--trainer" hargs))
         (output-fn (gethash "--output" hargs))
         (grammar (cando:load-cando grammar-fn))
         (trainer (cando:load-cando trainer-fn))
         results)
    (loop for topology in (kin:topologys grammar)
          do (cando:register-topology topology))
    (format t "Loaded grammar ~a -> ~a~%" grammar-fn grammar)
    (format t "Loaded trainer ~a -> ~a~%" trainer-fn trainer)
    (leap:easy-gaff2)
    (setf results (make-instance 'kin:rotamer-list))
    (let* ((conformation (kin:conformation trainer))
           (aggregate (kin:aggregate conformation)))
      (leap:assign-atom-types aggregate)
      (outer-search-conformations hargs trainer results)
      (ensure-directories-exist output-fn)
      (cando:save-cando results output-fn)
      #+(or)(kin:save-rotamer-list results output-fn))))

(defparameter *number-of-cycles* 10)
(defun build-conformation ()
  (let ((args (loop for ai from 1 below (core:argc)
                    collect (core:argv ai))))
    (build-conformation-args args)))


(defun build-conformation-str (str)
  (let* ((args (core:split str " ")))
    (format t "Args: ~s~%" args)
    (build-conformation-args args)))



(defun assemble-rotamer-library (project)
  (let* ((grammar-node (kin.tasks:make-grammar-node))
         (grammar-pn (kin.tasks::project-node-pathname grammar-node project))
         (grammar (cando:load-cando grammar-pn))
         (rotamer-map (make-hash-table :test 'equal)))
    (loop for topology in (kin:topologys grammar)
          do (cando:register-topology topology))
    (loop for pathname in (kin:trainer-node-pathnames grammar)
          for dir = (second (pathname-directory pathname))
          for trainer-node = (kin.tasks::make-trainer-node dir)
          for trainer = (cando:load-cando (kin.tasks:project-node-pathname trainer-node project))
          for conf-node = (kin.tasks::make-conformation-node trainer-node)
          for conf = (cando:load-cando (kin.tasks:project-node-pathname conf-node project))
          for context = (kin:context trainer)
          do (format t "trainer-node: ~a~%" trainer-node)
          do (format t "context: ~a~%" context)
          do (format t "conf-node: ~a~%" conf-node)
          do (format t "conf: ~a~%" conf)
          do (setf (gethash context rotamer-map) conf))
    (make-instance 'kin:rotamer-library
                   :context-to-rotamer-list rotamer-map)))


