(in-package :cando-user)
(progn
  (asdf:load-system "amber")
  (use-package :amber :cando-user)
  (load "~/Downloads/chromatinfiber-11million/setup.lisp"))

(amber:load-off "~/Development/amber/dat/leap/lib/solvents.lib")


(progn
  (defparameter *1mil-layout* (amber.pdb:scanpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :progress t))
  (defparameter *1mil* (amber.pdb:loadpdb "~/Downloads/chromatin-1mil/chromsmall.pdb" :layout *1mil-layout* :progress t)))

(cando:chimera *1mil*)

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

(defparameter *1mil-layout* (amber.pdb:scanpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :progress t))
(defparameter *1mil* (amber.pdb:loadpdb "~/Downloads/chromatin-1mil/chrom40.pdb" :layout *1mil-layout* :progress t))


;;; -------------
;;; Test the leap parser
;;;
(asdf:load-system "amber")
(progn
  (amber:add-path #P"/Users/meister/Development/amber/dat/leap/prep/")
  (amber:add-path #P"/Users/meister/Development/amber/dat/leap/lib/")
  (amber:add-path #P"/Users/meister/Development/amber/dat/leap/parm/")
  (amber:add-path #P"/Users/meister/Development/amber/dat/leap/cmd/"))

(amber.leap::source #P"leaprc.protein.fb15")
(apply 'amber.leap::leap.log-file '(:|leap.log|))
(amber.leap::leap.log-file '(goof :|foo.test|))

(apropos "system-index")
(ql-dist::system-index)

(fdefinition 'amber.leap::leap.log-file)
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
  (loop for x = (amber.leap::leap-parse-one-line sin nil :eof)
     until (eq x :eof)
     when x
     do (print x)))


(with-input-from-string (sin "a
b
c
")
  (read-delimited-list #\r sin))

