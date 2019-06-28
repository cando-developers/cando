;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :fepdot)

(defmethod draw-node (id (node fep:job) stream)
  (if (fep:script node)
      (format stream "~a [label = \"~a.~a\", style=\"filled,solid\",fillcolor=\"yellow\" ];~%" id (fep:name (fep:script node)) (fep:extension (fep:script node)))
      (format stream "~a [label = \"~a\", style=\"filled,solid\",fillcolor=\"#FFE0E0\" ];~%" id (class-name (class-of node)))))

(defmethod draw-node (id (node fep:amber-job-mixin) stream)
  (format stream "~a [label = \"~a\" ];~%" id (fep:name (fep:script node))))

(defmethod draw-node (id (node fep:cpptraj-job-mixin) stream)
  (format stream "~a [label = \"~a\" ];~%" id (fep:name (fep:script node))))

(defgeneric label (node))

(defmethod label ((node fep:node-file))
  (format nil "~{~a~%~}" (core:split (namestring (fep:node-pathname node)) "/")))

(defmethod draw-node (id (node fep:argument) stream)
  (format stream "~a [label = \"~a\",shape=rectangle];~%" id (label (fep:node node))))

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
  (if (eq (fep:option argument) :.)
      (let (definers)
        (mapc (lambda (one)
                (loop for definer in (fep:definers one)
                      do (pushnew definer definers :test #'eq)))
              (fep:node argument))
        definers)
      (fep:definers (fep:node argument))))

(defun safe-users (argument)
  (if (eq (fep:option argument) :.)
      (let (users)
        (mapc (lambda (one)
                (loop for definer in (fep:users one)
                      do (pushnew definer users :test #'eq)))
              (fep:node argument))
        users)
      (fep:users (fep:node argument))))

(defun safe-inputs (job)
  (loop for argument in (fep:inputs job)
        if (consp (fep:node argument))
          append (mapcar
                  (lambda (node-file) (make-instance 'fep::argument
                                                    :option :.part
                                                    :node node-file))
                  (fep:node argument))
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
          for drawn-input = (maybe-draw-file (fep:node input) input stream id-map all-files)
          do (when drawn-input (draw-edge drawn-input job-id (fep:option input) stream)))
    (loop for output in (fep:outputs job)
          for drawn-output = (maybe-draw-file (fep:node output) output stream id-map all-files)
          do (when drawn-output (draw-edge job-id drawn-output (fep:option output) stream)))))

(defun gather-jobs (job seen-jobs)
  (unless (gethash job seen-jobs)
    (setf (gethash job seen-jobs) t)
    (loop for output in (fep:outputs job)
          do (loop for user in (fep:users output)
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
