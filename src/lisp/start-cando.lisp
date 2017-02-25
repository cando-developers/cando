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
(progn
  (setf (logical-pathname-translations "cando")
        '(("**;*.*" "source-dir:extensions;cando;src;**;*.*"))))

;;; Add directories for ASDF to search for systems
(let* ((topdir (translate-logical-pathname #P"cando:lisp;cando;"))
       (dirs (all-subdirs topdir)))
  (push topdir asdf:*central-registry*)
  (dolist (dir dirs)
    (format t "Pushing dir: ~a~%" dir)
    (push dir asdf:*central-registry*)))

;;; Setup or startup the Cando system 
;;; If :setup-cando is in *features* then don't load the cando system
(progn
  (if (member :setup-cando *features*)
      (format t "Setting up cando without (asdf:loadsystem :cando)~%")
      (progn
        (format t "Starting Cando - loading cando system.~%")
        (asdf:load-system "cando")))
  (in-package :cando-user)
  (core:process-command-line-load-eval-sequence)
  (core::tpl))
