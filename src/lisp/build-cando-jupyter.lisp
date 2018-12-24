(require :asdf)
(progn
  (defun all-subdirs (dir)
    (let (dirs)
      (labels ((trav (d)
                 (dolist (d (uiop:subdirectories d))
                   (push d dirs)
                   (trav d))))
        (trav dir))
      dirs))

;;; Add the cando hostname
  (format t "Setting CANDO hostname~%")

;;; Add directories for ASDF to search for systems
  (let* ((topdir (translate-logical-pathname #P"source-dir:extensions;cando;src;lisp;"))
         (dirs (all-subdirs topdir)))
    (format t "Pushing topdir: ~a~%" topdir)
    (push topdir asdf:*central-registry*)
    (dolist (dir dirs)
      (format t "Pushing dir: ~a~%" dir)
      (push dir asdf:*central-registry*)))

  (progn
    (format t "Pushing dir: ~a~%" *default-pathname-defaults*)
    (push *default-pathname-defaults* asdf:*central-registry*))
  )
(load "quicklisp:setup.lisp")
(ql:quickload "build-cando-jupyter" :verbose t)
