
(require :asdf)
(load "quicklisp:setup.lisp")
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
  (progn
    (setf (logical-pathname-translations "cando")
          '(("**;*.*" "source-dir:extensions;cando;src;**;*.*"))))

;;; Add directories for ASDF to search for systems
  (let* ((topdir (translate-logical-pathname #P"cando:lisp;"))
         (dirs (all-subdirs topdir)))
    (push topdir asdf:*central-registry*)
    (dolist (dir dirs)
      (format t "Pushing dir: ~a~%" dir)
      (push dir asdf:*central-registry*)))

  (progn
    (format t "Pushing dir: ~a~%" *default-pathname-defaults*)
    (push *default-pathname-defaults* asdf:*central-registry*))
  )
(ql:quickload "build-cando" :verbose t)
