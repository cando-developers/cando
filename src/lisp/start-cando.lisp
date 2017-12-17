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
(eval-when (:load-toplevel :execute :compile-toplevel)
  (setq core::*echo-repl-read* t))

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

;;;(make-package :cando)

(let ((amber-home
        (namestring (uiop:ensure-directory-pathname (or (ext:getenv "AMBERHOME") "/amber/")))))
  (setf (logical-pathname-translations "amber")
        (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
  (format t "Setting *amber-home* -> ~a~%" amber-home))

;;; Setup to run slime if we are in a jupyter notebook
#+jupyter
(defparameter *started-swank* nil)
#+jupyter
(defun ext::start-swank ()
  ;; Bad!  This is hard-coded to work with docker
  (if *started-swank*
      (format t "Swank is already running~%")
      (progn
        (let ((swank-loader (or (probe-file (format nil "~a/~a" (or (ext:getenv "HOME") "/") "/slime/swank-loader.lisp"))
                                (probe-file (format nil "~a/~a" (or (ext:getenv "SLIME_HOME") "/") "/swank-loader.lisp"))
                                (error "Cannot find swank - set SLIME_HOME environment variable and try again"))))
          (format t "swank-loader -> ~a~%" swank-loader)
          (load swank-loader))
        (let ((swank-loader-init (find-symbol "INIT" "SWANK-LOADER")))
          (funcall swank-loader-init :delete nil :reload nil :load-contribs nil))
        (let ((swank-create-server (find-symbol "CREATE-SERVER" "SWANK")))
          (mp:process-run-function 'swank-main
                                   (lambda () (funcall swank-create-server
                                                       :port 4005
                                                       :interface "0.0.0.0")))
          (setf *started-swank* t)))))
#+jupyter
(export 'ext::start-swank :ext)

;;; Setup or startup the Cando system 
;;; If :setup-cando is in *features* then don't load the cando system
(format t "Starting Cando~%")
(format t "*features* -> ~a~%" *features*)

#-cando-jupyter
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
  #+(or)(core:process-command-line-load-eval-sequence)
  #+(or)(core::tpl))
