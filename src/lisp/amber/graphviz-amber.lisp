(in-package :graphviz-amber)

(defmethod draw-node (id (node amber:job) stream)
  (if (amber:script node)
      (format stream "~a [label = \"~a\", style=\"filled,solid\",fillcolor=\"yellow\"  ];~%" id (amber:filename (amber:script node)))
      (format stream "~a [label = \"~a\", style=\"filled,solid\",fillcolor=\"#FFE0E0\" ];~%" id (class-name (class-of node)))))

#+(or)
(defmethod draw-node (id (node amber:cpptraj-job-mixin) stream)
  (format stream "~a [label = \"~a\" ];~%" id (amber:name (amber:script node))))

(defgeneric label (node))

(defmethod label ((node amber:node-file))
  (format nil "~{~a~%~}" (core:split (namestring (amber:node-pathname node)) "/")))

(defmethod draw-node (id (node amber:argument) stream)
  (format stream "~a [label = \"~a\",shape=rectangle];~%" id (label (amber:node node))))

(defmethod draw-edge (source-id target-id label stream)
  (format stream "~a -> ~a [label = \"~a\"];~%" source-id target-id label))

(defun maybe-draw-node (key node stream id-map)
  (unless key
    (error "key is nil for node ~s~%" node))
  (let ((id (gethash key id-map)))
    (if id
        id
        (let ((tid (gensym)))
          (draw-node tid node stream)
          (setf (gethash key id-map) tid)
          tid))))

(defun safe-definers (argument)
  (if (eq (amber:option argument) :.)
      (let (definers)
        (mapc (lambda (one)
                (loop for definer in (amber:definers one)
                      do (pushnew definer definers :test #'eq)))
              (amber:node argument))
        definers)
      (amber:definers (amber:node argument))))

(defun safe-users (argument)
  (if (eq (amber:option argument) :.)
      (let (users)
        (mapc (lambda (one)
                (loop for definer in (amber:users one)
                      do (pushnew definer users :test #'eq)))
              (amber:node argument))
        users)
      (amber:users (amber:node argument))))

(defun safe-inputs (job)
  (loop for argument in (amber:inputs job)
        if (consp (amber:node argument))
          append (mapcar
                  (lambda (node-file) (make-instance 'amber::argument
                                                    :option :.part
                                                    :node node-file))
                  (amber:node argument))
        else
          append (list argument)))

(defun maybe-draw-file (key node stream id-map all-files)
  (unless key
    (error "key is nil for node ~s~%" node))
  (let ((definers (safe-definers node))
        (users (safe-users node)))
    (when (and (null all-files) definers users)
      (let ((id (gethash key id-map)))
        (if id
            id
            (let ((tid (gensym)))
              (draw-node tid node stream)
              (setf (gethash key id-map) tid)
              tid))))))

(defun draw-one-job (job stream id-map all-files)
  (let ((job-id (maybe-draw-node job job stream id-map)))
    (loop for input in (safe-inputs job)
          for drawn-input = (maybe-draw-file (amber:node input) input stream id-map all-files)
          do (when drawn-input (draw-edge drawn-input job-id (amber:option input) stream)))
    (loop for output in (amber:outputs job)
          for drawn-output = (maybe-draw-file (amber:node output) output stream id-map all-files)
          do (when drawn-output (draw-edge job-id drawn-output (amber:option output) stream)))))

(defun gather-jobs (job seen-jobs)
  (unless (gethash job seen-jobs)
    (setf (gethash job seen-jobs) t)
    (loop for output in (amber:outputs job)
          do (loop for user in (amber:users output)
                   do (gather-jobs user seen-jobs)))))

(defun draw-graph-stream (job stream &optional all-files)
  (format stream "digraph G {~%")
  (let ((unique-jobs (make-hash-table))
        (id-map (make-hash-table)))
    (gather-jobs job unique-jobs)
    (maphash (lambda (job dummy)
               (declare (ignore dummy))
               (draw-one-job job stream id-map all-files))
             unique-jobs))
  (format stream "}~%"))

#+(or)
(defun draw-graph-stream (job stream &optional all-files)
  (let ((unique-jobs (make-hash-table))
        (id-map (make-hash-table)))
    (gather-jobs job unique-jobs)
    (format stream "digraph G {~%")
    (maphash (lambda (job dummy)
               (declare (ignore dummy))
               (draw-one-job job stream id-map all-files))
             unique-jobs)
    (format stream "}~%")))
