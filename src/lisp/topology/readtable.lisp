(in-package :topology)

(named-readtables:defreadtable topology:syntax
  (:dispatch-macro-char #\# #\l
                        #'(lambda (stream char n)
                            (declare (ignore n))
                            (let ((open-ch (read-char stream)))
                              (unless (char= open-ch #\()
                                (error "First char after #l reader macro must be ("))
                              (let ((plug-name (let ((pn (make-array 16 :adjustable t :fill-pointer 0)))
                                                 (loop for ch = (read-char stream)
                                                       when (char= ch #\.)
                                                         do (return pn)
                                                       do (vector-push-extend ch pn)))))
                                (let ((dot-ch (read-char stream)))
                                  (unless (char= dot-ch #\.)
                                    (error "Separator char in #l reader macro must be a dot (.) - got ~c" dot-ch)))
                                (let ((digit-ch (read-char stream))
                                      (plug-atom-id nil))))))))

