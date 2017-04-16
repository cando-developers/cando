(add-path "/Users/meister/Development/amber/dat/leap/cmd/")
(add-path "/Users/meister/Development/amber/dat/leap/lib/")
(add-path "/Users/meister/Development/amber/dat/leap/dat/")
(add-path "/Users/meister/Development/amber/dat/leap/parm/")
(source "leaprc.gaff")

(leap.core:search-path #P"gaff.dat")
leap.core::*path*


(in-package :cando-user)
(progn
  (asdf:load-system "amber")
  (use-package :leap :cando-user)
  (load "~/Downloads/chromatinfiber-11million/setup.lisp"))

(amber:load-off "~/Development/amber/dat/leap/lib/solvents.lib")


(progn
  (defparameter *1mil-layout* (leap.pdb:scanpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :progress t))
  (defparameter *1mil* (leap.pdb:loadpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :layout *1mil-layout* :progress t)))

(cando:chimera *1mil*)


;;; ------------------------------------------------------------
;;;
;;; Test fortran code

;;; Compile the (defpackage #:fortran) in packages.lisp
;;; and then compile the fortran.lisp file

(apply #'open ("/tmp/test.for" :direction :output))

(fortran:with-fortran-output-file (fout "/tmp/test.for" :direction :output)
  (fortran:format 4 "%10d")
  (dotimes (i 100) (fortran:write i))
  (fortran:format 4 "%10.3f")
  (dotimes (i 50) (fortran:write (float (+ i 0.1)))))


(defun setup-amber-force-field ()
  (let ((parms (chem:make-read-amber-parameters)))
    (with-open-file (fin "~/Development/amber/dat/leap/parm/parm99.dat")
      (chem:read-parameters parms fin))
    (chem:get-force-field parms)))
(defparameter *amber* (setup-amber-force-field))


(progn
  (defparameter *ef* (core:make-cxx-object 'chem:energy-function))
  (defparameter *ff* (energy:setup-amber))
  (chem:generate-standard-energy-function-tables *ef* *1mil* *amber* :show-progress t)
  (chem:summarize-terms *ef*))

(defparameter *m* (chem:get-missing-parameters *ef*))
(loop for x in *m*
   when (typep x 'chem:ffstretch)
     collect x)

(chem:map-residues
 nil
 (lambda (r)
   (chem:map-atoms
    nil
    (lambda (a)
      (when (eq (chem:get-type a) :c8)
        (print (list (chem:get-name r) a ))))
    r))
   *1mil*)

(print (chem:get-number-of-atoms (chem:atom-table *ef*)))
(defun excluded-atom-list (ef)
  (let ((atom-table (chem:atom-table ef))
        (excluded-atoms-list (core:make-native-vector<int>))
        (excluded-atoms-num (core:make-native-vector<int>)))
    (loop for i below (chem:get-number-of-atoms atom-table)
       for num = (chem:push-back-excluded-atom-indices-and-sort atom-table excluded-atoms-list i)
       do (core:native-vector<int>-push-back excluded-atoms-num num))
    (values excluded-atoms-num excluded-atoms-list)))
(chem:number-of-terms (chem:get-stretch-component *ef*))
(chem:safe-amber-energy-stretch-term (chem:get-stretch-component *ef*) 99999999)


(defvar *ean*)
(defvar *eal*)
(time (multiple-value-setq (*ean* *eal*) (excluded-atom-list *ef* )))
(core:native-vector<int>-size *ean*)
(core:native-vector<int>-size *eal*)

(loop for x below 100
   do (format t "~a -> ~a~%" x (core:native-vector<int>-elt *eal* x)))

(defparameter *eal* (core:make-native-vector<int>))
(chem:push-back-excluded-atom-indices-and-sort (chem:atom-table *ef*) *eal* 1)
(loop for i below (core:native-vector<int>-size *eal*)
   collect (core:native-vector<int>-elt *eal* i))


(apropos "excluded")
(print (chem:excluded-atom-list (chem:atom-table *ef*) 0))

(chem:atom-
(chem:generate-restraint-energy-function-tables *ef* *1mil* *ff* :show-progress t)
(print "Done")



(chem:contents-as-list *1mil*)

(defparameter *1mil-layout* (leap.pdb:scanpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :progress t))
(defparameter *1mil* (leap.pdb:loadpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :layout *1mil-layout* :progress t))


;;; -------------
;;; Test the leap parser
;;;
(asdf:load-system "amber")

(leap.core:add-path #P"/Users/meister/Development/amber/dat/leap/prep/")

(find-symbol "ADD-PATH" :leap.core)


(progn
  (add-path #P"/Users/meister/Development/amber/dat/leap/prep/")
  (add-path #P"/Users/meister/Development/amber/dat/leap/lib/")
  (add-path #P"/Users/meister/Development/amber/dat/leap/parm/")
  (add-path #P"/Users/meister/Development/amber/dat/leap/cmd/"))

(source #P"leaprc.protein.fb15")

(defparameter pep (

(apply 'leap.commands::leap.log-file '(:|leap.log|))
(leap.commands::leap.log-file '(goof :|foo.test|))

(apropos "system-index")
(ql-dist::system-index)

(fdefinition 'leap.commands::leap.log-file)
(apropos "load-script")


(with-input-from-string (sin "
loadOff test.lib
# This is a comment

x = { \"A\"
      \"B\"
      \"C\" }
test = 3
;asdfas
x = 4")
  (loop for x = (leap.commands::leap-parse-one-line sin nil :eof)
     until (eq x :eof)
     when x
     do (print x)))


(with-input-from-string (sin "a
b
c
")
  (read-delimited-list #\r sin))

