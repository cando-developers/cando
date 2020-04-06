;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :ti-rundot)

(defmethod draw-node (id (node ti-run:job) stream)
  (if (ti-run:script node)
      (format stream "~a [label = \"~a.~a\", style=\"filled,solid\",fillcolor=\"yellow\" ];~%" id (ti-run:name (ti-run:script node)) (ti-run:extension (ti-run:script node)))
      (format stream "~a [label = \"~a\", style=\"filled,solid\",fillcolor=\"#FFE0E0\" ];~%" id (class-name (class-of node)))))

(defmethod draw-node (id (node ti-run:amber-job-mixin) stream)
  (format stream "~a [label = \"~a\" ];~%" id (ti-run:name (ti-run:script node))))

(defgeneric label (node))

(defmethod label ((node ti-run:node-file))
  (format nil "~{~a~%~}" (core:split (namestring (ti-run:node-pathname node)) "/")))

(defmethod draw-node (id (node ti-run:argument) stream)
  (format stream "~a [label = \"~a\",shape=rectangle];~%" id (label (ti-run:node node))))

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
  (if (eq (ti-run:option argument) :.)
      (let (definers)
        (mapc (lambda (one)
                (loop for definer in (ti-run:definers one)
                      do (pushnew definer definers :test #'eq)))
              (ti-run:node argument))
        definers)
      (ti-run:definers (ti-run:node argument))))

(defun safe-users (argument)
  (if (eq (ti-run:option argument) :.)
      (let (users)
        (mapc (lambda (one)
                (loop for definer in (ti-run:users one)
                      do (pushnew definer users :test #'eq)))
              (ti-run:node argument))
        users)
      (ti-run:users (ti-run:node argument))))

(defun safe-inputs (job)
  (loop for argument in (ti-run:inputs job)
        if (consp (ti-run:node argument))
          append (mapcar
                  (lambda (node-file) (make-instance 'ti-run::argument
                                                    :option :.part
                                                    :node node-file))
                  (ti-run:node argument))
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
          for drawn-input = (maybe-draw-file (ti-run:node input) input stream id-map all-files)
          do (when drawn-input (draw-edge drawn-input job-id (ti-run:option input) stream)))
    (loop for output in (ti-run:outputs job)
          for drawn-output = (maybe-draw-file (ti-run:node output) output stream id-map all-files)
          do (when drawn-output (draw-edge job-id drawn-output (ti-run:option output) stream)))))

(defun gather-jobs (job seen-jobs)
  (unless (gethash job seen-jobs)
    (setf (gethash job seen-jobs) t)
    (loop for output in (ti-run:outputs job)
          do (loop for user in (ti-run:users output)
                   do (gather-jobs user seen-jobs)))))

(defun draw-graph-stream (jobs stream &optional all-files)
  (let ((unique-jobs (make-hash-table))
        (id-map (make-hash-table)))
    (loop for job in jobs
          do (gather-jobs job unique-jobs))
    (format stream "digraph G {~%")
    (maphash (lambda (job dummy)
               (declare (ignore dummy))
               (draw-one-job job stream id-map all-files))
             unique-jobs)
    (format stream "}~%")))
