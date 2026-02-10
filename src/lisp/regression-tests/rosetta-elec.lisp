#+tests(in-package #:clasp-tests)

(defparameter *rosetta-elec-test-force-field* :rosetta-elec-test)
(defparameter *rosetta-elec-test-type* :rnb)

(defparameter *rosetta-elec-default-params*
  (list :elec-weight 1.0d0
        :eps-core 6.0d0
        :eps-solvent 80.0d0
        :rmin 1.45d0
        :rlow 1.85d0
        :rhi 4.5d0
        :rcut 5.5d0
        :rpairlist 7.5d0))

(defun rosetta-elec--params (&key
                               (elec-weight (getf *rosetta-elec-default-params* :elec-weight))
                               (eps-core (getf *rosetta-elec-default-params* :eps-core))
                               (eps-solvent (getf *rosetta-elec-default-params* :eps-solvent))
                               (rmin (getf *rosetta-elec-default-params* :rmin))
                               (rlow (getf *rosetta-elec-default-params* :rlow))
                               (rhi (getf *rosetta-elec-default-params* :rhi))
                               (rcut (getf *rosetta-elec-default-params* :rcut))
                               (rpairlist (getf *rosetta-elec-default-params* :rpairlist)))
  (list :elec-weight elec-weight
        :eps-core eps-core
        :eps-solvent eps-solvent
        :rmin rmin
        :rlow rlow
        :rhi rhi
        :rcut rcut
        :rpairlist rpairlist))

(defun rosetta-elec--ensure-force-field ()
  (let* ((force-field (chem:force-field/make))
         (nonbond-db (chem:get-nonbond-db force-field))
         (ffnonbond (chem:make-ffnonbond *rosetta-elec-test-type*)))
    (chem:ffnonbond/set-radius-angstroms ffnonbond 1.5d0)
    (chem:ffnonbond/set-epsilon-kcal ffnonbond 0.2d0)
    (chem:ffnonbond-db-add nonbond-db ffnonbond)
    (when (chem:find-force-field *rosetta-elec-test-force-field* nil)
      (leap:clear-force-field *rosetta-elec-test-force-field*))
    (leap.core:add-force-field-or-modification
     force-field
     :force-field-name *rosetta-elec-test-force-field*
     :combined-force-field-class-name 'chem:combined-force-field)
    ffnonbond))

(defun rosetta-elec--make-two-atom-aggregate (force-field-name atom-type &key (charge1 1.0d0) (charge2 -1.0d0))
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
    (chem:set-charge a1 charge1)
    (chem:set-charge a2 charge2)
    (values agg a1 a2)))

(defun rosetta-elec--amber-charge-conversion ()
  (cond ((boundp 'chem:*amber-charge-conversion-18.2223*)
         (symbol-value 'chem:*amber-charge-conversion-18.2223*))
        (t 18.2223d0)))

(defun rosetta-elec--kqq (charge1 charge2)
  (let* ((conv (rosetta-elec--amber-charge-conversion))
         (scale (* conv conv)))
    (* scale charge1 charge2)))

(defun rosetta-elec--term (params kqq)
  (let* ((eps-core (getf params :eps-core))
         (eps-solvent (getf params :eps-solvent))
         (rmin (getf params :rmin))
         (rlow (getf params :rlow))
         (rhi (getf params :rhi))
         (rcut (getf params :rcut))
         (rx-rmin (* 0.25d0 rmin))
         (ex-rmin (exp (- rx-rmin)))
         (gx-rmin (* ex-rmin (+ 1.0d0 rx-rmin (* 0.5d0 rx-rmin rx-rmin))))
         (ep-rmin (+ eps-solvent (* gx-rmin (- eps-core eps-solvent))))
         (e-rmin (* kqq (/ 1.0d0 ep-rmin)))
         (rx-rlow (* 0.25d0 rlow))
         (ex-rlow (exp (- rx-rlow)))
         (gx-rlow (* ex-rlow (+ 1.0d0 rx-rlow (* 0.5d0 rx-rlow rx-rlow))))
         (ep-rlow (+ eps-solvent (* gx-rlow (- eps-core eps-solvent))))
         (e-rlow (* kqq (/ 1.0d0 ep-rlow)))
         (gx-p-rlow (* -0.5d0 ex-rlow rx-rlow rx-rlow))
         (ep-p-rlow (* 0.25d0 gx-p-rlow (- eps-core eps-solvent)))
         (de-rlow (* ep-p-rlow (/ 1.0d0 (* ep-rlow ep-rlow)) (- kqq)))
         (rx-rhi (* 0.25d0 rhi))
         (ex-rhi (exp (- rx-rhi)))
         (gx-rhi (* ex-rhi (+ 1.0d0 rx-rhi (* 0.5d0 rx-rhi rx-rhi))))
         (ep-rhi (+ eps-solvent (* gx-rhi (- eps-core eps-solvent))))
         (e-rhi (* kqq (/ 1.0d0 ep-rhi)))
         (gx-p-rhi (* -0.5d0 ex-rhi rx-rhi rx-rhi))
         (ep-p-rhi (* 0.25d0 gx-p-rhi (- eps-core eps-solvent)))
         (de-rhi (* ep-p-rhi (/ 1.0d0 (* ep-rhi ep-rhi)) (- kqq)))
         (dr-low (- rlow rmin))
         (dr-high (- rcut rhi))
         (aa-low-shifted (* (+ (* -2.0d0 e-rlow) (* 2.0d0 e-rmin) (* de-rlow dr-low))
                            (/ 1.0d0 (* dr-low dr-low dr-low))))
         (bb-low-shifted (* (+ (* -3.0d0 e-rmin) (* 3.0d0 e-rlow) (* -1.0d0 de-rlow dr-low))
                            (/ 1.0d0 (* dr-low dr-low))))
         (cc-low-shifted 0.0d0)
         (dd-low-shifted e-rmin)
         (aa-low aa-low-shifted)
         (bb-low (+ bb-low-shifted (* -3.0d0 aa-low-shifted rmin)))
         (cc-low (+ cc-low-shifted (* -2.0d0 bb-low-shifted rmin) (* 3.0d0 aa-low-shifted rmin rmin)))
         (dd-low (+ dd-low-shifted
                    (* bb-low-shifted rmin rmin)
                    (* -1.0d0 aa-low-shifted rmin rmin rmin)
                    (* -1.0d0 cc-low-shifted rmin)))
         (aa-high-shifted (* (+ (* 2.0d0 e-rhi) (* de-rhi dr-high))
                             (/ 1.0d0 (* dr-high dr-high dr-high))))
         (bb-high-shifted (* (+ (* -2.0d0 de-rhi dr-high) (* -3.0d0 e-rhi))
                             (/ 1.0d0 (* dr-high dr-high))))
         (cc-high-shifted de-rhi)
         (dd-high-shifted e-rhi)
         (aa-high aa-high-shifted)
         (bb-high (+ bb-high-shifted (* -3.0d0 aa-high-shifted rhi)))
         (cc-high (+ cc-high-shifted (* -2.0d0 bb-high-shifted rhi) (* 3.0d0 aa-high-shifted rhi rhi)))
         (dd-high (+ dd-high-shifted
                     (* bb-high-shifted rhi rhi)
                     (* -1.0d0 aa-high-shifted rhi rhi rhi)
                     (* -1.0d0 cc-high-shifted rhi))))
    (list :kqq kqq
          :e-rmin e-rmin
          :aa-low aa-low
          :bb-low bb-low
          :cc-low cc-low
          :dd-low dd-low
          :aa-high aa-high
          :bb-high bb-high
          :cc-high cc-high
          :dd-high dd-high)))

(defun rosetta-elec--expected-energy (rr term params)
  (let* ((elec-weight (getf params :elec-weight))
         (eps-core (getf params :eps-core))
         (eps-solvent (getf params :eps-solvent))
         (rmin (getf params :rmin))
         (rlow (getf params :rlow))
         (rhi (getf params :rhi))
         (rcut (getf params :rcut))
         (kqq (getf term :kqq))
         (e-rmin (getf term :e-rmin))
         (aa-low (getf term :aa-low))
         (bb-low (getf term :bb-low))
         (cc-low (getf term :cc-low))
         (dd-low (getf term :dd-low))
         (aa-high (getf term :aa-high))
         (bb-high (getf term :bb-high))
         (cc-high (getf term :cc-high))
         (dd-high (getf term :dd-high)))
    (cond
      ((<= rr rmin)
       (* elec-weight e-rmin))
      ((<= rr rlow)
       (* elec-weight (+ dd-low (* cc-low rr) (* rr rr (+ bb-low (* aa-low rr))))))
      ((<= rr rhi)
       (let* ((rx (* 0.25d0 rr))
              (ex (exp (- rx)))
              (gx (* ex (+ 1.0d0 rx (* 0.5d0 rx rx))))
              (ep (+ eps-solvent (* (- eps-core eps-solvent) gx)))
              (ee (* kqq (/ 1.0d0 ep))))
         (* elec-weight ee)))
      ((<= rr rcut)
       (* elec-weight (+ dd-high (* cc-high rr) (* rr rr (+ bb-high (* aa-high rr))))))
      (t 0.0d0))))

(defun rosetta-elec--component (energy-function)
  (or (find-if (lambda (comp) (typep comp 'chem:energy-rosetta-elec))
               (chem:all-components energy-function))
      (error "No rosetta elec component found in energy-function")))

(defun rosetta-elec--setup (params)
  (list :rosetta
        (list 'chem:energy-rosetta-elec
              :elec-weight (getf params :elec-weight)
              :eps-core (getf params :eps-core)
              :eps-solvent (getf params :eps-solvent)
              :rmin (getf params :rmin)
              :rlow (getf params :rlow)
              :rhi (getf params :rhi)
              :rcut (getf params :rcut)
              :rpairlist (getf params :rpairlist))))

(defun rosetta-elec-generate-xy-pairs (&key (elec-weights (list 0.5d0 1.0d0 2.0d0))
                                            (charge1 1.0d0)
                                            (charge2 -1.0d0)
                                            (start 0.1d0)
                                            (end 10.0d0)
                                            (step 0.1d0)
                                            (eps-core (getf *rosetta-elec-default-params* :eps-core))
                                            (eps-solvent (getf *rosetta-elec-default-params* :eps-solvent))
                                            (rmin (getf *rosetta-elec-default-params* :rmin))
                                            (rlow (getf *rosetta-elec-default-params* :rlow))
                                            (rhi (getf *rosetta-elec-default-params* :rhi))
                                            (rcut (getf *rosetta-elec-default-params* :rcut))
                                            (rpairlist (getf *rosetta-elec-default-params* :rpairlist))
                                            (energy-source :actual)
                                            stream)
  "Generate x,y pairs for the rosetta elec between two atoms.

Returns an alist of (elec-weight . ((distance energy) ...)).
When STREAM is non-nil, prints \"distance energy\" lines to STREAM and
separates each elec-weight series with a comment line."
  (rosetta-elec--ensure-force-field)
  (loop for weight in elec-weights
        collect
        (multiple-value-bind (agg a1 a2)
            (rosetta-elec--make-two-atom-aggregate *rosetta-elec-test-force-field*
                                                   *rosetta-elec-test-type*
                                                   :charge1 charge1
                                                   :charge2 charge2)
          (let* ((params (rosetta-elec--params
                          :elec-weight weight
                          :eps-core eps-core
                          :eps-solvent eps-solvent
                          :rmin rmin
                          :rlow rlow
                          :rhi rhi
                          :rcut rcut
                          :rpairlist rpairlist))
                 (setup (rosetta-elec--setup params))
                 (energy-function (chem:make-energy-function :matter agg
                                                             :use-excluded-atoms nil
                                                             :assign-types nil
                                                             :setup setup))
                 (component (rosetta-elec--component energy-function))
                 (pos (chem:make-nvector (chem:get-nvector-size energy-function)))
                 (kqq (rosetta-elec--kqq charge1 charge2))
                 (term (rosetta-elec--term params kqq))
                 (pairs nil))
            (when stream
              (format stream "# elec-weight ~,6f~%" weight))
            (loop for dist from start to end by step
                  do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                     (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                     (chem:load-coordinates-into-vector energy-function pos)
                     (let* ((energy (ecase energy-source
                                      (:actual
                                       (chem:energy-component-evaluate-energy energy-function component pos))
                                      (:expected
                                       (rosetta-elec--expected-energy dist term params)))))
                       (push (list dist energy) pairs)
                       (when stream
                         (format stream "~,6f ~,12f~%" dist energy))))
            (setf pairs (nreverse pairs))
            (cons weight pairs)))))

(defun run-rosetta-elec-generate-xy-pairs ()
  (with-open-file (fout "/tmp/rosetta-elec.xy" :direction :output)
    (rosetta-elec-generate-xy-pairs
     :elec-weights (list 1.25d0)
     :start 0.2d0
     :end 8.0d0
     :step 0.05d0
     :stream fout)))



(core:set-simd-width 1)

(let* ((ffnonbond (rosetta-elec--ensure-force-field))
       (charge1 1.0d0)
       (charge2 -1.0d0)
       (params (rosetta-elec--params))
       (rmin (getf params :rmin))
       (rlow (getf params :rlow))
       (rhi (getf params :rhi))
       (rcut (getf params :rcut))
       (distances (list (* 0.8d0 rmin)
                        (* 0.5d0 (+ rmin rlow))
                        (* 0.5d0 (+ rlow rhi))
                        (* 0.5d0 (+ rhi rcut))
                        (+ rcut 0.5d0)))
       (elec-weights (list 0.5d0 1.0d0 2.0d0))
       (tol 1.0d-8)
       (all-ok t))
  (declare (ignore ffnonbond))
  (loop for weight in elec-weights
        do (multiple-value-bind (agg a1 a2)
               (rosetta-elec--make-two-atom-aggregate *rosetta-elec-test-force-field*
                                                      *rosetta-elec-test-type*
                                                      :charge1 charge1
                                                      :charge2 charge2)
             (let* ((weight-params (rosetta-elec--params
                                    :elec-weight weight
                                    :eps-core (getf params :eps-core)
                                    :eps-solvent (getf params :eps-solvent)
                                    :rmin rmin
                                    :rlow rlow
                                    :rhi rhi
                                    :rcut rcut
                                    :rpairlist (getf params :rpairlist)))
                    (setup (rosetta-elec--setup weight-params))
                    (energy-function (chem:make-energy-function :matter agg
                                                                :use-excluded-atoms nil
                                                                :assign-types nil
                                                                :setup setup))
                    (component (rosetta-elec--component energy-function))
                    (pos (chem:make-nvector (chem:get-nvector-size energy-function)))
                    (kqq (rosetta-elec--kqq charge1 charge2))
                    (term (rosetta-elec--term weight-params kqq)))
               (loop for dist in distances
                     do (chem:set-position a1 (geom:vec 0.0d0 0.0d0 0.0d0))
                        (chem:set-position a2 (geom:vec dist 0.0d0 0.0d0))
                        (chem:load-coordinates-into-vector energy-function pos)
                        (let* ((actual (chem:energy-component-evaluate-energy energy-function component pos))
                               (expected (rosetta-elec--expected-energy dist term weight-params))
                               (delta (abs (- actual expected))))
                          (when (>= delta tol)
                            #-tests(format t "Problem for elec-weight(~s) dist(~s) delta(~s) should be less than ~s~%"
                                           weight dist delta tol)
                            (setf all-ok nil)))))))
  #+tests
  (test-true rosetta-elec-energy-scan all-ok))
