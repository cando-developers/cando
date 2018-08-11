(in-package :cando-user)

(defun show (matter)
  (cond
    ((typep matter 'chem:aggregate))
    ((typep matter 'chem:molecule)
     (let ((agg (chem:make-aggregate nil)))
       (chem:add-matter agg matter)
       (setf matter agg)))
    ((typep matter 'chem:residue)
     (let ((agg (chem:make-aggregate nil))
           (mol (chem:make-molecule nil)))
       (chem:add-matter agg mol)
       (chem:add-matter mol matter)
       (setf matter agg)))
    (t (error "You cannot show a ~a" matter)))
  (cando-jupyter:show-aggregate matter))

(defun repr (widget representation &optional (selection "all"))
  (funcall (find-symbol "ADD-REPRESENTATION" :nglv) widget representation :selection selection))

(defun cl-jupyter-kernel-start (&optional connection-file-name)
  (let ((amber-home
          (namestring (or (if (ext:getenv "AMBERHOME")
                              (probe-file (ext:getenv "AMBERHOME"))
                              "/usr/local/amber/")
                          (probe-file "/home/app/amber16-data/")
                          "/dev/null"))))
    (setf (logical-pathname-translations "amber")
          (list (list "**;*.*" (concatenate 'string amber-home "/**/*.*"))))
    (format t "Setting amber host pathname translation -> ~a~%" amber-home))
  (let ((*readtable* (copy-readtable)))
    (set-macro-character #\~ #'leap.core:leap-lookup-variable-reader-macro)
    (let ((cl-jup (find-symbol "KERNEL-START" "CL-JUPYTER")))
      (if cl-jup
          (if connection-file-name
              (funcall cl-jup connection-file-name)
              (funcall cl-jup))
          (error "cl-jupyter is not installed")))))


