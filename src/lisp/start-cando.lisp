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
(in-package :common-lisp-user)

(format t "Running start-cando.lisp script~%")

;;; Load the ASDF package manager
(require :asdf)

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

(make-package :cando)

(let ((amber-home
       (namestring (uiop:ensure-directory-pathname (or (ext:getenv "AMBERHOME") "/amber/")))))
  (setf (logical-pathname-translations "amber")
        (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
  (format t "Setting *amber-home* -> ~a~%" amber-home))

(defvar cando::*root-directory*)
(progn
  (if (member :docker *features*)
    (setf cando::*root-directory* (pathname "/src/"))
    (setf cando::*root-directory* (uiop:ensure-directory-pathname (ext:getenv "HOME"))))
  (format t "Setting *root-directory* -> ~a~%" cando::*root-directory*))

;;; Setup or startup the Cando system 
;;; If :setup-cando is in *features* then don't load the cando system
(format t "Starting Cando~%")
(progn
  (if (member :setup-cando *features*)
      (progn
        (format t "Skipping (asdf:loadsystem :cando)~%"))
      (progn
        (format t "Loading Cando system.~%")
        (load "~/quicklisp/setup.lisp")
        (funcall (find-symbol "QUICKLOAD" :ql) "cando-user")
        ;; Ensure that all threads start in the :CANDO-USER package
        (core:symbol-global-value-set '*package* (find-package :cando-user))))
  (core:process-command-line-load-eval-sequence)
  (core::tpl))
