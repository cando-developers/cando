(in-package #:topology-shape-key-tests)

(defun piece-rotamer-summary (pairs)
  (map 'vector (lambda (pair)
                 (cons (car pair) (aref (topology:internals-values (cdr pair)) 0))) pairs))

(defun test-shape-key-consumers ()
  (let* ((database (make-instance 'topology:rotamers-database :foldamer-name :test))
         (rotamers (make-instance 'topology:sidechain-rotamers
                                  :monomer-context :test-context
                                  :rotamer-vector
                                  (map 'vector (lambda (value)
                                                 (make-instance 'topology:sidechain-rotamer
                                                                :internals-values (vector value)))
                                       '(1d0 2d0 3d0))))
         (map (topology:shape-key-to-index rotamers))
         (key (topology:make-binned-shape-key '((:ring 2 13) (:bb1 12))))
         (equal-key (topology:make-binned-shape-key '((:bb1 12) (:ring 2 13))))
         (other-key (topology:make-binned-shape-key '((:bb1 13) (:ring 2 13)))))
    (setf (gethash :test-context (topology:context-to-rotamers database)) rotamers
          (topology:binned-shape-key-map-value key map) #(0 2)
          (topology:binned-shape-key-map-value other-key map) #(1))
    (assert (= 2 (topology::max-rotamer-count rotamers)))
    (assert (equalp #((0 . 1d0) (2 . 3d0))
                    (piece-rotamer-summary
                     (topology::allowed-rotamers-for-piece :test-context equal-key database :test))))
    (assert (equalp #((1 . 2d0))
                    (piece-rotamer-summary
                     (topology::allowed-rotamers-for-piece :test-context other-key database :test))))
    (uiop:with-temporary-file (:pathname pathname)
      (cando.serialize:save-cando database pathname)
      (dolist (native '(t nil))
        (let ((restored (cando.serialize:load-cando pathname :use-native-reader native)))
          (assert (equalp #((0 . 1d0) (2 . 3d0))
                          (piece-rotamer-summary
                           (topology::allowed-rotamers-for-piece
                            :test-context equal-key restored :test))))))))
  (format t "SHAPE-KEY-CONSUMER-TESTS-PASSED~%")
  t)
