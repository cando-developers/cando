;;;
;;;    File: start-cando.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
#+(or)
(eval-when (:load-toplevel :execute :compile-toplevel)
  (format t "Turning on core:*echo-repl-read* in start-cando.lisp~%")
  (setq core::*echo-repl-read* t))

(in-package :common-lisp-user)

#+(or)
(format t "Starting start-cando.lisp script - arguments: ~s~%" core:*command-line-arguments*)

;;; Load the ASDF package manager
(progn
  (require :asdf))
#+(or)
(format t "Loaded asdf version ~s~%" (asdf/upgrade:asdf-version))

;;; Register eclector with asdf so we don't try and load it again
(asdf:register-immutable-system :eclector)

(defun all-subdirs (dir)
  (let (dirs)
    (labels ((trav (d)
               (dolist (d (uiop:subdirectories d))
                 (push d dirs)
                 (trav d))))
      (trav dir))
    dirs))

;;; Add the cando hostname
(progn
  (setf (logical-pathname-translations "cando")
        '(("**;*.*" "source-dir:extensions;cando;src;**;*.*"))))

;;; Add directories for ASDF to search for systems
(let* ((topdir (translate-logical-pathname #P"cando:lisp;"))
       (dirs (all-subdirs topdir)))
  (push topdir asdf:*central-registry*)
  (dolist (dir dirs)
    (push dir asdf:*central-registry*)))

(progn
  (push *default-pathname-defaults* asdf:*central-registry*))

;;;(make-package :cando)

(let ((amber-home
        (namestring (uiop:ensure-directory-pathname (or (ext:getenv "AMBERHOME") "/amber/")))))
  (setf (logical-pathname-translations "amber")
        (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*")))))

;;; Setup or startup the Cando system 
;;; If :setup-cando is in *features* then don't load the cando system
(progn
  #+(or)
  (handler-bind ((error (lambda (&rest args) (sys:safe-backtrace))))
    (load "quicklisp:setup.lisp"))
  (load "quicklisp:setup.lisp"))

(progn
  (funcall (find-symbol "QUICKLOAD" :ql) "cando-user" :silent (null (core:is-interactive-lisp)))
  ;; Ensure that all threads start in the :CANDO-USER package
  (core:symbol-global-value-set '*package* (find-package :cando-user))
  )

;; Startup lparallels
(cando:lazy-setup-lparallel-kernel)

