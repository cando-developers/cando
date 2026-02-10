#+tests(in-package #:clasp-tests)

(defparameter *rosetta-nonbond-test-force-field* :rosetta-nonbond-test)
(defparameter *rosetta-nonbond-test-type* :rnb)
(defparameter *rosetta-nonbond-rrep-scale* 0.6000000238418579d0)

(defun rosetta-nonbond--ensure-force-field ()
  (let* ((force-field (chem:force-field/make))
         (nonbond-db (chem:get-nonbond-db force-field))
         (ffnonbond (chem:make-ffnonbond *rosetta-nonbond-test-type*)))
    (chem:ffnonbond/set-radius-angstroms ffnonbond 1.5d0)
    (chem:ffnonbond/set-epsilon-kcal ffnonbond 0.2d0)
    (chem:ffnonbond-db-add nonbond-db ffnonbond)
    (when (chem:find-force-field *rosetta-nonbond-test-force-field* nil)
      (leap:clear-force-field *rosetta-nonbond-test-force-field*))
    (leap.core:add-force-field-or-modification
     force-field
     :force-field-name *rosetta-nonbond-test-force-field*
     :combined-force-field-class-name 'chem:combined-force-field)
    ffnonbond))

(defun rosetta-nonbond--make-two-atom-aggregate (force-field-name atom-type)
  (let* ((agg (chem:make-aggregate))
         (mol1 (chem:make-molecule :m1))
         (mol2 (chem:make-molecule :m2))
         (res1 (chem:make-residue :r1))
         (res2 (chem:make-residue :r2))
         (a1 (chem:make-atom :a1 :c))
         (a2 (chem:make-atom :a2 :c)))
    (chem:add-atom res1 a1)
    (chem:add-atom res2 a2)
    (chem:add-matter mol1 res1)
    (chem:add-matter mol2 res2)
    (chem:add-molecule agg mol1)
    (chem:add-molecule agg mol2)
    (chem:setf-force-field-name agg force-field-name)
    (chem:setf-force-field-name mol1 force-field-name)
    (chem:setf-force-field-name mol2 force-field-name)
    (chem:set-property a1 :given-atom-type atom-type)
    (chem:set-property a2 :given-atom-type atom-type)
    (values agg a1 a2)))

(defun rosetta-nonbond--term (sigma epsilon rswitch rcut)
  (let* ((rrep (* *rosetta-nonbond-rrep-scale* sigma))
         (sr-rep (/ sigma rrep))
         (sr2-rep (* sr-rep sr-rep))
         (sr6-rep (* sr2-rep sr2-rep sr2-rep))
         (sr12-rep (* sr6-rep sr6-rep))
         (rep-rrep (* epsilon (+ 1.0d0 sr12-rep (* -2.0d0 sr6-rep))))
         (drep-rrep (* epsilon (+ (* -12.0d0 sr12-rep) (* 12.0d0 sr6-rep)) (/ 1.0d0 rrep)))
         (mij drep-rrep)
         (bij (- rep-rrep (* mij rrep)))
         (sr-sw (/ sigma rswitch))
         (sr2-sw (* sr-sw sr-sw))
         (sr6-sw (* sr2-sw sr2-sw sr2-sw))
         (sr12-sw (* sr6-sw sr6-sw))
         (atr-sw (* epsilon (+ sr12-sw (* -2.0d0 sr6-sw))))
         (d-atr-sw (* epsilon (+ (* -12.0d0 sr12-sw) (* 12.0d0 sr6-sw)) (/ 1.0d0 rswitch)))
         (dr (- rcut rswitch))
         (aa-shift (/ (+ (* 2.0d0 atr-sw) (* dr d-atr-sw)) (* dr dr dr)))
         (bb-shift (/ (+ (* -2.0d0 dr d-atr-sw) (* -3.0d0 atr-sw)) (* dr dr)))
         (cc-shift d-atr-sw)
         (dd-shift atr-sw)
         (aa aa-shift)
         (bb (+ bb-shift (* -3.0d0 aa-shift rswitch)))
         (cc (+ cc-shift (* -2.0d0 bb-shift rswitch) (* 3.0d0 aa-shift rswitch rswitch)))
         (dd (+ dd-shift (* bb-shift rswitch rswitch)
                (* -1.0d0 aa-shift rswitch rswitch rswitch)
                (* -1.0d0 cc-shift rswitch))))
    (list :sigma sigma :epsilon epsilon :mij mij :bij bij :aa aa :bb bb :cc cc :dd dd)))

(defun rosetta-nonbond--term-from-ffnonbond (ffnonbond rswitch rcut)
  (let* ((radius (chem:ffnonbond/get-radius-angstroms ffnonbond))
         (epsilon (chem:ffnonbond/get-epsilon-kcal ffnonbond))
         (rstar (+ radius radius))
         (epsilonij (sqrt (* epsilon epsilon)))
         (rt6 (expt rstar 6))
         (rt12 (* rt6 rt6))
         (parm-a (* epsilonij rt12))
         (parm-c (* 2.0d0 epsilonij rt6))
         (sigma6 (/ (* 2.0d0 parm-a) parm-c))
         (sigma (expt sigma6 (/ 1.0d0 6.0d0)))
         (epsilon-term (/ (* parm-c parm-c) (* 4.0d0 parm-a))))
    (rosetta-nonbond--term sigma epsilon-term rswitch rcut)))

(defun rosetta-nonbond--expected-energy (rr term rep-weight rswitch rcut)
  (let* ((sigma (getf term :sigma))
         (epsilon (getf term :epsilon))
         (mij (getf term :mij))
         (bij (getf term :bij))
         (aa (getf term :aa))
         (bb (getf term :bb))
         (cc (getf term :cc))
         (dd (getf term :dd)))
    (cond
      ((<= rr (* *rosetta-nonbond-rrep-scale* sigma))
       (+ (- epsilon) (* rep-weight (+ bij (* mij rr)))))
      ((<= rr sigma)
       (let* ((sr (/ sigma rr))
              (sr2 (* sr sr))
              (sr6 (* sr2 sr2 sr2))
              (sr12 (* sr6 sr6))
              (fa-rep-epsilon (* epsilon (+ sr12 (* -2.0d0 sr6))))
              (fa-rep (+ epsilon fa-rep-epsilon))
              (fa-atr (- epsilon)))
         (+ fa-atr (* rep-weight fa-rep))))
      ((<= rr rswitch)
       (let* ((sr (/ sigma rr))
              (sr2 (* sr sr))
              (sr6 (* sr2 sr2 sr2))
              (sr12 (* sr6 sr6)))
         (* epsilon (+ sr12 (* -2.0d0 sr6)))))
      ((<= rr rcut)
       (+ dd (* cc rr) (* rr rr (+ bb (* aa rr)))))
      (t 0.0d0))))

(defun rosetta-nonbond--component (energy-function)
  (or (find-if (lambda (comp) (typep comp 'chem:energy-rosetta-nonbond))
               (chem:all-components energy-function))
      (error "No rosetta nonbond component found in energy-function")))

(defun rosetta-nonbond-generate-xy-pairs (&key (rep-weights (list 0.5d0 1.0d0 2.0d0))
                                               (start 0.1d0)
                                               (end 10.0d0)
                                               (step 0.1d0)
                                               (rswitch 4.5d0)
                                               (rcut 6.0d0)
                                               (energy-source :actual)
                                               stream)
  "Generate x,y pairs for the rosetta nonbond between two atoms.

Returns an alist of (rep-weight . ((distance energy) ...)).
When STREAM is non-nil, prints \"distance energy\" lines to STREAM and
separates each rep-weight series with a comment line."
  (let* ((ffnonbond (rosetta-nonbond--ensure-force-field))
         (term (rosetta-nonbond--term-from-ffnonbond ffnonbond rswitch rcut)))
    (loop for rep in rep-weights
          collect
          (multiple-value-bind (agg a1 a2)
              (rosetta-nonbond--make-two-atom-aggregate *rosetta-nonbond-test-force-field*
                                                        *rosetta-nonbond-test-type*)
            (let* ((setup (list :rosetta (list 'chem:energy-rosetta-nonbond :rep-weight rep)))
                   (energy-function (chem:make-energy-function :matter agg
                                                               :use-excluded-atoms nil
                                                               :assign-types nil
                                                               :setup setup))
                   (component (rosetta-nonbond--component energy-function))
                   (pos (chem:make-nvector (chem:get-nvector-size energy-function)))
                   (pairs nil)
                   (energy-rosetta-nonbond (chem:find-component energy-function 'chem:energy-rosetta-nonbond))
                   )
              (format t "energy-rosetta-nonbond ~s~%" (core:encode energy-rosetta-nonbond))
              (when stream
                (format stream "# rep-weight ~,6f~%" rep))
              (loop for dist from start to end by step
                    do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                       (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                       (chem:load-coordinates-into-vector energy-function pos)
                       (let* ((energy (ecase energy-source
                                        (:actual
                                         (chem:energy-component-evaluate-energy energy-function component pos))
                                        (:expected
                                         (rosetta-nonbond--expected-energy dist term rep rswitch rcut)))))
                         (push (list dist energy) pairs)
                         (when stream
                           (format stream "~,6f ~,12f~%" dist energy))))
              (setf pairs (nreverse pairs))
              (cons rep pairs))))))

(defun run-rosetta-nonbond-generate-xy-pairs ()
  (with-open-file (fout "/tmp/rosetta-nonbond.xy" :direction :output)
    (rosetta-nonbond-generate-xy-pairs
     :rep-weights (list 0.0)
     :start 2.5
     :end 4.0
     :step 0.01
     :stream fout)))

(core:set-simd-width 1)

(let* ((ffnonbond (rosetta-nonbond--ensure-force-field))
       (rswitch 4.5d0)
       (rcut 6.0d0)
       (term (rosetta-nonbond--term-from-ffnonbond ffnonbond rswitch rcut))
       (sigma (getf term :sigma))
       (distances (list (* 0.5d0 sigma)
                        (* 0.8d0 sigma)
                        (* 1.2d0 sigma)
                        (+ rswitch 0.3d0)
                        (+ rcut 0.3d0)))
       (rep-weights (list 0.5d0 1.0d0 2.0d0))
       (tol 1.0d-8)
       (all-ok t))
  (loop for rep in rep-weights
        do (multiple-value-bind (agg a1 a2)
               (rosetta-nonbond--make-two-atom-aggregate *rosetta-nonbond-test-force-field*
                                                         *rosetta-nonbond-test-type*)
             (let* ((setup (list :rosetta (list 'chem:energy-rosetta-nonbond :rep-weight rep)))
                    (energy-function (chem:make-energy-function :matter agg
                                                                :use-excluded-atoms nil
                                                                :assign-types nil
                                                                :setup setup))
                    (component (rosetta-nonbond--component energy-function))
                    (pos (chem:make-nvector (chem:get-nvector-size energy-function))))
               (loop for dist in distances
                     do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                        (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                        (chem:load-coordinates-into-vector energy-function pos)
                        (let* ((actual (chem:energy-component-evaluate-energy energy-function component pos))
                               (expected (rosetta-nonbond--expected-energy dist term rep rswitch rcut))
                               (delta (abs (- actual expected))))
                          (when (>= delta tol)
                            #-tests(format t "Problem for rep(~s) dist(~s) delta(~s) should be less than ~s~%" rep dist delta tol)
                            (setf all-ok nil)))))))
  #+tests
  (test-true rosetta-nonbond-energy-scan all-ok))
