;;;
;;;    File: molecules.lisp
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
(in-package :cando)


(defun scramble-positions (matter &key (center (geom:vec 0.0 0.0 0.0)) (width 10.0)
                           &aux (half-width (/ width 2.0)))
  (flet ((scramble-atom (atom)
           (let ((pos (geom:vec
                       (+ (- (random width) half-width) (geom:get-x center))
                       (+ (- (random width) half-width) (geom:get-y center))
                       (+ (- (random width) half-width) (geom:get-z center)))))
             (chem:set-position atom pos))))
    (cond
      ((typep matter 'chem:matter)
       (chem:map-atoms
        nil
        (lambda (a) (scramble-atom a)) matter))
      ((hash-table-p matter)
       (maphash (lambda (atm value)
                  (declare (ignore value))
                  (scramble-atom atm)) matter))
      ((listp matter)
       (mapc (lambda (a) (scramble-atom a)) matter))
      (t (error "Add support to scramble positions for ~s" matter)))))
              
  
(defun dump-atoms (matter)
  (chem:map-atoms nil
                  (lambda (x)
                    (format t "~6a ~a ~a ~a ~s~%"
                            (chem:get-name x)
                            (chem:get-type x)
                            (chem:get-configuration x)
                            (chem:get-position x)
                            (chem:properties x)))
                  matter))

(defun dump-stereocenters (centers)
  (dotimes (i (length centers))
    (let ((a (elt centers i)))
      (format t "~6a ~a~%" (chem:get-name a) (chem:get-configuration a)))))

(defun set-stereoisomer (atom-vec index &key show)
  (dotimes (i (length atom-vec))
    (let ((config (if (logbitp i index) :S :R))
          (atom (elt atom-vec i)))
      (chem:set-configuration atom config)))
  (when show
    (let ((bits (ceiling (length atom-vec))))
      (format t "======== Stereoisomer ~3d/~v,'0b~%" index bits index)
      (dump-stereocenters atom-vec))))

(defun gather-stereocenters (matter)
  (let (chiral-atoms)
    (chem:map-atoms nil (lambda (atm)
                          (when (eq (chem:get-stereochemistry-type atm) :chiral)
                            (push atm chiral-atoms)))
                    matter)
    (make-array (length chiral-atoms) :initial-contents chiral-atoms)))

(defun set-stereoisomer-func (atom-vec func &key show)
  "* Arguments
- atom-vec :: A vector of chiral centers
- func :: A function that takes an atom and returns 'chem:S or 'chem:R
* Description
Set the stereochemistry of a collection of stereocenters using a function that returns :S or :R for each atom."
  (dotimes (i (length atom-vec))
    (let* ((atom (elt atom-vec i))
           (config (funcall func (elt atom-vec i))))
      (chem:set-configuration atom config)))
  (when show
    (dump-stereocenters atom-vec)))

(defun set-stereoisomer-mapping (matter atom-name-to-config &key show)
  "* Arguments
- matter : The aggregate/molecule/residue that contains the named atoms.
- atom-name-to-config :: An alist that maps atom names to :r or :s
* Description
Set the stereochemistry of a collection of stereocenters using a alist of atom names to stereochemistry.
Example:  (set-stereoisomer-mapping *agg* '((:C1 :R) (:C2 :S))"
  (loop for (name . config) in atom-name-to-config
     do (let ((atom (chem:first-atom-with-name matter name)))
          (when show (format t "Atom name: ~a  atom: ~a config: ~a~%" name atom config))
          (chem:set-configuration atom config))))

;;; Set a list of stereocenters using an integer
;;; A 0 bit is :S and 1 bit is :R
;;; The value 13 is #b1101  and it sets the configuration of the atoms
;;;    to (:R :S :R :R ).
;;;  The least significant bit is the first stereocenter and the most significant bit is the last.
(defun set-stereoisomer-using-number (list-of-centers num)
  (loop for atom across list-of-centers
        for tnum = num then (ash tnum -1)
        for config = (if (= (logand 1 tnum) 0) :s :r)
        do (chem:set-configuration atom config)
        do (format t "~a -> ~a   num: ~a~%" atom config tnum)))

;;; Set all of the stereocenters to config (:S or :R)
(defun set-all-stereocenters-to (list-of-centers config &key show)
  (set-stereoisomer-func list-of-centers (constantly config) :show show)
  (format t "~a stereocenters set~%" (length list-of-centers)))

;;; Return the number of stereoisomers 
(defun number-of-stereoisomers (vector-of-centers)
  (expt 2 (length vector-of-centers)))

(defun stereocenters-sorted-by-name (agg)
  (let ((stereocenters (cando:gather-stereocenters agg)))
    (sort stereocenters #'string< :key #'chem:get-name)))

(defmacro with-handle-linear-angles-dihedrals ((&key (max-times 3) verbose) &body body)
  `(let ((remaining-linear-handles ,max-times))
     (handler-bind
         ((chem::linear-atoms-error
            (lambda (condition)
              (let* ((width 2.0)
                     (half-width (/ width 2.0)))
                (loop for atom in (chem:atoms condition)
                      do (let* ((cp (chem:get-position atom))
                                (pos (geom:vec
                                      (+ (- (random width) half-width) (geom:get-x cp))
                                      (+ (- (random width) half-width) (geom:get-y cp))
                                      (+ (- (random width) half-width) (geom:get-z cp)))))
                           (chem:set-position atom pos))))
              (when ,verbose (format t "Encountered linear-atoms-error restarting~%"))
              (decf remaining-linear-handles)
              (when (> remaining-linear-handles 0)
                (chem:restart-minimizer)))))
       (progn
         ,@body))))

;; Recover from minimization problems using Common Lisp restarts
(defun minimize-no-fail (minimizer &key resignal-error verbose active-atoms-mask coords)
  (if verbose
      (progn
        (chem:enable-print-intermediate-results minimizer)
        (format t "enable-print-intermediate-results ~%"))
      (chem:disable-print-intermediate-results minimizer))
  (restart-case
      (handler-bind
          ((chem:minimizer-error (lambda (err)
                                   (warn "In minimize-no-fail - the minimizer reported: ~a" err)
                                   (invoke-restart 'cando:skip-rest-of-minimization err))))
        (with-handle-linear-angles-dihedrals (:max-times 3 :verbose verbose)
          (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
            (chem:minimize minimizer :active-atom-mask active-atoms-mask :coords coords :verbose verbose))))
    ;; skip-rest-of-minimization can also be triggered by the user from the debugger
    (skip-rest-of-minimization (err)
      :report "Skip the rest of the current minimization - continue processing"
      (chem:write-intermediate-results-to-energy-function minimizer)
      (when resignal-error (error err)))))

(defmacro with-ignore-minimizer-errors (&body body)
  `(handler-bind
       ((chem:minimizer-error (lambda (err)
                                (warn "In with-ignore-minimizer-errors - the minimizer reported: ~a" err)
                                (invoke-restart 'cando:skip-rest-of-minimization err))))
     (progn
       ,@body)))

;; Recover from minimization problems using Common Lisp restarts
(defun minimize-with-restarts (minimizer &key resignal-error verbose)
  (if verbose
      (progn
        (chem:enable-print-intermediate-results minimizer)
        (format t "enable-print-intermediate-results ~%"))
      (chem:disable-print-intermediate-results minimizer))
  (restart-case
      (with-handle-linear-angles-dihedrals (:verbose verbose)
        (multiple-value-bind (pos total-energy)
            (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
              (chem:minimize minimizer))
          (unless total-energy
            (break "total-energy was NIL"))
          (values pos total-energy)
          ))
    ;; skip-rest-of-minimization can also be triggered by the user from the debugger
    (skip-rest-of-minimization (&optional err)
      :report "Skip the rest of the current minimization - continue processing"
      (chem:write-intermediate-results-to-energy-function minimizer)
      (when resignal-error (error err)))
    (save-and-skip-rest-of-minimization (pathname)
      :report "Save the current structure to a file and skip the rest of the current minimization - continue processing"
      (chem:write-intermediate-results-to-energy-function minimizer)
      (cando.serialize:save-cando (chem:get-matter (chem:get-energy-function minimizer)) pathname))))


(defun minimizer-obey-interrupt (minimizer)
  (restart-case
      (handler-bind
          ((ext:unix-signal-received (lambda (c)
                                       (declare (ignore c))
                                       (print "Done") (invoke-restart 'skip-rest))))
        (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
          (chem:minimize minimizer)))
    (skip-rest ()
      :report "Skip rest of minimization"
      (chem:write-intermediate-results-to-energy-function minimizer)
      (print "Skipping rest of minimization") )))


(defun configure-minimizer (minimizer
                            &key (max-sd-steps 50)
                              (max-cg-steps 100)
                              (max-tn-steps 100)
                              (sd-tolerance 2000.0)
                              (cg-tolerance 0.5)
                              (tn-tolerance 0.0001)
                              )
  (chem:set-maximum-number-of-steepest-descent-steps minimizer max-sd-steps)
  (chem:set-maximum-number-of-conjugate-gradient-steps minimizer max-cg-steps)
  (chem:set-maximum-number-of-truncated-newton-steps minimizer max-tn-steps)
  (chem:set-steepest-descent-tolerance minimizer sd-tolerance)
  (chem:set-conjugate-gradient-tolerance minimizer cg-tolerance)
  (chem:set-truncated-newton-tolerance minimizer tn-tolerance))

(defun validate-energy-components (components valid-components)
  (when (consp components)
    (loop for comp in components
          unless (member comp valid-components)
            do (error "~s is not a valid component name - must be one of ~s" comp valid-components))))

(defun optimize-structure-debug (matter &key (keep-interaction-factory t)
                                          (disable-components nil)
                                          (enable-components nil)
                                          (use-excluded-atoms nil)
                                          (electrostatic-scale 1.0 electrostatic-scale-p)
                                          (sd-tolerance 500.0)
                                          (cg-tolerance 0.5)
                                          (tn-tolerance 0.0001)
                                          (max-sd-steps 500)
                                          (max-cg-steps 500)
                                          (max-tn-steps 500)
                                          verbose)
  "Minimize energy of a structure with lots of control.
: matter - The matter to optimize.
: disable-components - Either T to disable all or a list of energy component class names to disable
Disabling happens before enabling - so you can disable all with T and then selectively enable.
: enable-componennts - A list of energy component class names to enable.
"
  (let* ((energy-function (chem:make-energy-function :matter matter
                                                     :use-excluded-atoms use-excluded-atoms
                                                     :assign-types t
                                                     :keep-interaction-factory keep-interaction-factory))
         (min (chem:make-minimizer energy-function))
         (valid-components nil))
    (when disable-components
      (loop for comp in (chem:all-components energy-function)
            for class = (class-of comp)
            for name = (class-name class)
            do (push name valid-components)
            when (or (eq disable-components t)
                     (and (listp disable-components) (member name disable-components)))
              do (progn
                   (when verbose (format t "Disabling ~a~%" name))
                   (chem:disable comp))))
    (when enable-components
      (loop for comp in (chem:all-components energy-function)
            for class = (class-of comp)
            for name = (class-name class)
            when (and (listp enable-components) (member name enable-components))
              do (progn
                   (when verbose (format t "Enabling ~a~%" name))
                   (chem:enable comp))))
    (validate-energy-components disable-components valid-components)
    (validate-energy-components enable-components valid-components)
    (configure-minimizer min
                         :sd-tolerance sd-tolerance
                         :cg-tolerance cg-tolerance
                         :tn-tolerance tn-tolerance
                         :max-sd-steps max-sd-steps
                         :max-cg-steps max-cg-steps
                         :max-tn-steps max-tn-steps)
    (when verbose (chem:enable-print-intermediate-results min))
    (minimize-no-fail min :verbose verbose)
    (finish-output t)
    (values matter energy-function)))

(defun optimize-structure (matter &key (keep-interaction-factory t)
                                    (disable-components nil)
                                    (enable-components nil)
                                    (use-excluded-atoms nil)
                                    (turn-off-nonbond t)
                                    (electrostatic-scale 1.0 electrostatic-scale-p)
                                    (sd-tolerance 500.0)
                                    (cg-tolerance 0.5)
                                    (tn-tolerance 0.0001)
                                    (max-sd-steps 1000)
                                    (max-cg-steps 2000)
                                    (max-tn-steps 5000)
                                    verbose)
  "Minimize energy of a structure with lots of control.
: matter - The matter to optimize.
: disable-components - Either T to disable all or a list of energy component class names to disable
Disabling happens before enabling - so you can disable all with T and then selectively enable.
: enable-componennts - A list of energy component class names to enable.
"
  (let* ((energy-function (chem:make-energy-function :matter matter
                                                     :use-excluded-atoms use-excluded-atoms
                                                     :assign-types t
                                                     :keep-interaction-factory keep-interaction-factory))
         (min (chem:make-minimizer energy-function))
         (valid-components nil))
    (when disable-components
      (loop for comp in (chem:all-components energy-function)
            for class = (class-of comp)
            for name = (class-name class)
            do (push name valid-components)
            when (or (eq disable-components t)
                     (and (listp disable-components) (member name disable-components)))
              do (progn
                   (when verbose (format t "Disabling ~a~%" name))
                   (chem:disable comp))))
    (when enable-components
      (loop for comp in (chem:all-components energy-function)
            for class = (class-of comp)
            for name = (class-name class)
            when (and (listp enable-components) (member name enable-components))
              do (progn
                   (when verbose (format t "Enabling ~a~%" name))
                   (chem:enable comp))))
    (validate-energy-components disable-components valid-components)
    (validate-energy-components enable-components valid-components)
    (configure-minimizer min
                         :sd-tolerance sd-tolerance
                         :cg-tolerance cg-tolerance
                         :tn-tolerance tn-tolerance
                         :max-sd-steps max-sd-steps
                         :max-cg-steps max-cg-steps
                         :max-tn-steps max-tn-steps)
    (when verbose (chem:enable-print-intermediate-results min))
    (when turn-off-nonbond
      (chem:set-option energy-function 'chem::nonbond-term nil)
      (finish-output t)
      (minimize-no-fail min :verbose verbose))
    (chem:set-option energy-function 'chem:nonbond-term t)
    (finish-output t)
    (minimize-no-fail min :verbose verbose)
    (finish-output t)
    (values matter energy-function)))

(defun optimize-structure-with-restarts (matter &key (keep-interaction-factory t)
                                                  (turn-off-nonbond t)
                                                  verbose
                                                  energy-function
                                                  (max-sd-steps 5000)
                                                  (max-cg-steps 5000)
                                                  (max-tn-steps 500)
                                                  (sd-tolerance 100.0)
                                                  (cg-tolerance 0.5)
                                                  (tn-tolerance 0.0001)
                                                  )
  (let* ((energy-function (if energy-function
                              energy-function
                              (chem:make-energy-function :matter matter
                                                         :assign-types t
                                                         :keep-interaction-factory keep-interaction-factory)))
         (min (chem:make-minimizer energy-function)))
    (configure-minimizer min
                         :sd-tolerance sd-tolerance
                         :cg-tolerance cg-tolerance
                         :tn-tolerance tn-tolerance
                         :max-sd-steps max-sd-steps
                         :max-cg-steps max-cg-steps
                         :max-tn-steps max-tn-steps
                         )
    (when verbose (chem:enable-print-intermediate-results min))
    (when turn-off-nonbond
      (chem:set-option energy-function 'chem::nonbond-term nil)
      (when verbose
        (format t "Optimization with nonbond term off~%")
        (finish-output t))
      (minimize-with-restarts min :verbose verbose))
    (chem:set-option energy-function 'chem:nonbond-term t)
    (when verbose
      (format t "Optimization with nonbond term on~%")
      (finish-output t))
    (multiple-value-bind (pos total-energy)
        (minimize-with-restarts min :verbose verbose)
      (finish-output t)
      (values matter energy-function total-energy))))

(defun indexed-pathname (template index)
  (let* ((filename (format nil "~a~3,'0d" (pathname-name template) index)))
    (make-pathname :name filename
                            :defaults template)))



;; Load aggregates from a collection of files with the
;; name PATHNAME-TEMPLATE using MAX-FILES
;; Ex: (array-stereoisomers #P"/tmp/input.mol2" 9)
;;    will load files /tmp/input000.mol2 to /tmp/input008.mol2, translate
;;    the first molecule in each file onto an X/Y grid and combine them
;;    all into a new Aggregate and write it out to /tmp/input_grid.mol2

(defun array-stereoisomers (pathname-template max-files &key (scale 20.0))
  (let* ((width (ceiling (sqrt max-files)))
         (accumulate-aggregate (chem:make-aggregate))
         (target-pn (make-pathname :name (format nil "~a_grid"
                                                 (pathname-name pathname-template))
                                   :defaults pathname-template)))
    (dotimes (i max-files)
      (multiple-value-bind (xi yi)
          (floor i width)
        (let ((xo (* scale xi))
              (yo (* scale yi)))
          (let ((pn (indexed-pathname pathname-template i)))
            (when (probe-file pn)
              (let* ((agg (chem:load-mol2 pn))
                     (mol (chem:first-molecule agg)))
                (format t "Loaded ~a~%" pn)
                (chem:translate-all-atoms mol (make-array 3 :initial-contents (list xo yo 0.0)))
                (chem:add-matter accumulate-aggregate mol)))))))
    (format t "Writing result to: ~a~%" target-pn)
    (save-mol2 accumulate-aggregate target-pn)))

(defun chimera (&rest objs)
  (let (files)
    (dolist (obj objs)
      (cond
        ((pathnamep obj)
         (push (namestring obj) files))
        ((typep obj 'chem:aggregate)
         (let ((pn (pathname (format nil "/tmp/temp~a.mol2" (length files)))))
           (cando:save-mol2 obj pn :use-sybyl-types t)
           (push pn files)))
        (t (error "You cannot run chimera on ~a" obj))))
    (let ((cmd (format nil "open -a Chimera.app ~{~a ~}" (mapcar (lambda (x) (namestring x)) (nreverse files)))))
      (format t "Executing: ~a~%" cmd)
      (ext:system cmd))))

(defun bad-geometry-p (agg)
  (error "I removed check-for-beyond-threshold-interactions from the C++ code and so I'll need to add something back to take its place - I only did this for stretch and angle terms and the code has been broken for a long time.  I may not even use this calling function bad-geometry-p - I do still use topology:bad-geometry-p - take a look at that.")
  (let ((coords (chem:matter/extract-coordinates agg)))
    (when (not (ext:array-no-nans-p coords))
      (return-from bad-geometry-p t))
    (let* ((energy-function (chem:make-energy-function :matter agg))
           (fails (progn
                    #+(or)(chem:check-for-beyond-threshold-interactions energy-function 0.2))))
      (if (> (length fails) 0)
          fails
          nil))))

(defclass sketch-nonbond-force-field () ())

(defmethod chem:find-atom-type-position ((force-field sketch-nonbond-force-field) type)
  (declare (ignore force-field type))
  0)

(defun randomize-coordinates (coordinates &key unfrozen from-zero (width 40.0) &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (index)
           (when from-zero
             (setf (elt coordinates (+ index 0)) 0.0
                   (elt coordinates (+ index 1)) 0.0
                   (elt coordinates (+ index 2)) 0.0))
           (let* ((xp (elt coordinates (+ index 0)))
                  (yp (elt coordinates (+ index 1)))
                  (zp (elt coordinates (+ index 2)))
                  (xn (+ (- (random width) half-width) xp))
                  (yn (+ (- (random width) half-width) yp))
                  (zn (+ (- (random width) half-width) zp)))
             (setf (elt coordinates (+ index 0)) xn
                   (elt coordinates (+ index 1)) yn
                   (elt coordinates (+ index 2)) zn))))
    (loop for coord-index below (length coordinates) by 3
          for index from 0
          do (if unfrozen
                 (when (= (elt unfrozen index) 1)
                   (jostle-atom coord-index))
                 (jostle-atom coord-index)))))

(defparameter *stage1-scale-sketch-nonbond* 0.04)
(defparameter *stage1-flatten-force-components* (list 0.85 0.85 0.85))
(defparameter *stage1-sketch-nonbond-force* 0.2)
(defparameter *stage1-bond-length* 1.5)
(defparameter *first-bond-force* 0.1)
(defparameter *stage1-nonbond-constant* 0.04)


(defun prepare-stage1-sketch-function (sketch-function &optional sketch)
  "Generate bond energy terms for a 2D sketch."
  (declare (ignore sketch))
  (let* ((molecule (chem:get-graph sketch-function))
         (atom-table (chem:node-table sketch-function))
         (nonbond-energy (chem:get-sketch-nonbond-component sketch-function))
         (bond-energy (chem:get-stretch-component sketch-function)))
    (chem:set-scale-sketch-nonbond nonbond-energy *stage1-scale-sketch-nonbond*)
    (loop for ia1 from 0 below (1- (chem:get-number-of-atoms atom-table))
          for atom1 = (chem:elt-atom atom-table ia1)
          for atom1-coordinate-index = (chem:get-coordinate-index-times3 atom-table atom1)
          for freeze1 = (or (= (chem:get-atomic-number atom1) 1)
                            (eq (chem:get-element atom1) :lp))
          do (loop for ia2 from (1+ ia1) below (chem:get-number-of-atoms atom-table)
                   for atom2 = (chem:elt-atom atom-table ia2)
                   for atom2-coordinate-index = (chem:get-coordinate-index-times3 atom-table atom2)
                   for freeze2 = (or (= (chem:get-atomic-number atom2) 1)
                                     (eq (chem:get-element atom2) :lp))
                   for freeze-flag = (if (or freeze1 freeze2) 1 0)
                   do (chem:add-sketch-nonbond-term nonbond-energy
                                                    atom1-coordinate-index
                                                    atom2-coordinate-index
                                                    freeze-flag
                                                    *stage1-sketch-nonbond-force*)))
    (flet
        ((add-stretch (bond-energy a1 a2 a1ci a2ci bond-force bond-length)
           (declare (ignore a1 a2))
           #+(or)(format t "Adding sketch-stretch-term ~a ~a ~a ~a~%" (chem:get-name a1) (chem:get-name a2) bond-force bond-length)
           (chem:add-sketch-stretch-term bond-energy a1ci a2ci bond-force bond-length)))
      (chem:map-bonds
       'nil
       (lambda (a1 a2 bond-order bond)
         (declare (ignore bond-order bond))
         (let ((atom1-coord-index (chem:get-coordinate-index-times3 atom-table a1))
               (atom2-coord-index (chem:get-coordinate-index-times3 atom-table a2))
               (bond-length *stage1-bond-length*))
           #+(or)
           (let ((valence1 (chem:matter-get-property-or-default a1 :valence 1))
                 (valence2 (chem:matter-get-property-or-default a2 :valence 1)))
             (format t "a1 ~a a2 ~a   valence1 ~a valence2 ~a~%" (chem:get-name a1) (chem:get-name a2) valence1 valence2))
           (add-stretch bond-energy
                        a1 a2
                        atom1-coord-index atom2-coord-index
                        *first-bond-force* bond-length)))
       molecule))
    #+(or)(when sketch
      (loop for double-bond-restraint in (double-bond-restraints sketch)
            do (with-slots (top-left-atom bottom-left-atom left-atom right-atom top-right-atom bottom-right-atom) double-bond-restraint
                 (let ((top-left-index (chem:get-coordinate-index-times3 atom-table top-left-atom))
                       (bottom-left-index (chem:get-coordinate-index-times3 atom-table bottom-left-atom))
                       (top-right-index (chem:get-coordinate-index-times3 atom-table top-right-atom))
                       (bottom-right-index (chem:get-coordinate-index-times3 atom-table bottom-right-atom)))
                   (add-stretch bond-energy
                                top-left-atom top-right-atom
                                top-left-index top-right-index
                                +first-bond-force+
                                (* 2.0 +stage1-bond-length+))
                   (add-stretch bond-energy
                                bottom-left-atom bottom-right-atom
                                bottom-left-index bottom-right-index
                                +first-bond-force+
                                (* 2.0 +stage1-bond-length+))))))
    sketch-function))

(defun advance-simulation (dynamics &key unfrozen)
  (dynamics:velocity-verlet-step
   dynamics
   :velocity-verlet-function #'chem:sketch-function-velocity-verlet-step
   :unfrozen unfrozen))

(defun starting-geometry (agg &key accumulate-coordinates verbose)
  "Rapidly calculate a starting geometry for the single molecule in the aggregate"
  (unless (= (chem:content-size agg) 1)
    (format t "The aggregate must have a single molecule.")
    (return-from starting-geometry nil))
  (let ((atom-types (make-hash-table)))
    (chem:map-atoms
     nil
     (lambda (atm)
       (setf (gethash atm atom-types) :sketch))
     agg)
    (let* ((mol (cando:mol agg 0))
           (dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function mol dummy-sketch-nonbond-ff atom-types))
           (dynamics (dynamics:make-atomic-simulation sketch-function 
                                                      :accumulate-coordinates accumulate-coordinates))
           (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
      (chem:disable (chem:get-point-to-line-restraint-component sketch-function))
      (apply #'chem:setf-velocity-scale sketch-function *stage1-flatten-force-components*)
      ;; stage 1  step 0 - 1999
      (progn
        (prepare-stage1-sketch-function sketch-function)
        (chem:set-scale-sketch-nonbond energy-sketch-nonbond *stage1-nonbond-constant*)
        ;; Everything interesting happens in 1000 steps
        (randomize-coordinates (dynamics:coordinates dynamics) :from-zero t :unfrozen nil)
        (dotimes (i 1000) (advance-simulation dynamics :unfrozen nil)))
      (dynamics:write-coordinates-back-to-matter dynamics)
      (optimize-structure agg :turn-off-nonbond nil :verbose verbose)
      dynamics
      )))


(defclass rapid-starting-geometry ()
  ((dynamics :initarg :dynamics :accessor dynamics)
   (energy-function :initarg :energy-function :accessor energy-function)))

(defun make-rapid-starting-geometry (agg &key energy-function accumulate-coordinates verbose (keep-interaction-factory t))
  "Rapidly calculate a starting geometry for the single molecule in the aggregate"
  (let ((atom-types (make-hash-table))
        (energy-function (if energy-function
                             energy-function
                             (chem:make-energy-function :matter agg
                                                        :assign-types t
                                                        :keep-interaction-factory keep-interaction-factory))))
    (unless (= (chem:content-size agg) 1)
      (format t "The aggregate must have a single molecule.")
      (return-from make-rapid-starting-geometry nil))
    (chem:map-atoms
     nil
     (lambda (atm)
       (setf (gethash atm atom-types) :sketch))
     agg)
    (let* ((mol (cando:mol agg 0))
           (dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function mol dummy-sketch-nonbond-ff atom-types))
           (dynamics (dynamics:make-atomic-simulation sketch-function 
                                                      :accumulate-coordinates accumulate-coordinates))
           (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
      (chem:disable (chem:get-point-to-line-restraint-component sketch-function))
      (apply #'chem:setf-velocity-scale sketch-function *stage1-flatten-force-components*)
      ;; stage 1  step 0 - 1999
      (when verbose
        (format t "Generating starting structure~%"))
      (progn
        (prepare-stage1-sketch-function sketch-function)
        (chem:set-scale-sketch-nonbond energy-sketch-nonbond *stage1-nonbond-constant*)
        ;; Everything interesting happens in 1000 steps
        (make-instance 'rapid-starting-geometry
                       :dynamics dynamics
                       :energy-function energy-function)))))


#+(or)
(defun next-rapid-starting-geometry (rapid-starting-geometry &key verbose)
  (with-slots (dynamics energy-function) rapid-starting-geometry
    (randomize-coordinates (dynamics:coordinates dynamics) :from-zero t :unfrozen nil)
    (dotimes (i 1000) (advance-simulation dynamics :unfrozen nil))
    (dynamics:write-coordinates-back-to-matter dynamics)
    (when verbose
      (format t "Optimizing structure~%"))
    (if energy-function
        (optimize-structure-with-restarts agg :turn-off-nonbond nil :verbose verbose :energy-function energy-function)
        (optimize-structure-with-restarts agg :turn-off-nonbond nil :verbose verbose))
    rapid-starting-geometry))

(defun starting-geometry-with-restarts (agg &key accumulate-coordinates verbose energy-function
                                              (max-sd-steps 5000)
                                              (max-cg-steps 5000)
                                              (max-tn-steps 500)
                                              (sd-tolerance 10.0)
                                              (cg-tolerance 1.0)
                                              (tn-tolerance 0.001)
                                              )
  "Rapidly calculate a starting geometry for the single molecule in the aggregate"
  (let ((atom-types (make-hash-table)))
    (unless (= (chem:content-size agg) 1)
      (format t "The aggregate must have a single molecule.")
      (return-from starting-geometry-with-restarts nil))
    (chem:map-atoms
     nil
     (lambda (atm)
       (setf (gethash atm atom-types) :sketch))
     agg)
    (let* ((mol (cando:mol agg 0))
           (dummy-sketch-nonbond-ff (make-instance 'sketch-nonbond-force-field))
           (sketch-function (chem:make-sketch-function mol dummy-sketch-nonbond-ff atom-types))
           (dynamics (dynamics:make-atomic-simulation sketch-function 
                                                      :accumulate-coordinates accumulate-coordinates))
           (energy-sketch-nonbond (chem:get-sketch-nonbond-component sketch-function)))
      (chem:disable (chem:get-point-to-line-restraint-component sketch-function))
      (apply #'chem:setf-velocity-scale sketch-function *stage1-flatten-force-components*)
      ;; stage 1  step 0 - 1999
      (when verbose
        (format t "Generating starting structure~%"))
      (prepare-stage1-sketch-function sketch-function)
      (chem:set-scale-sketch-nonbond energy-sketch-nonbond *stage1-nonbond-constant*)
      (randomize-coordinates (dynamics:coordinates dynamics) :from-zero t :unfrozen nil)
      ;; The time consuming part of this happens next 1000 steps of dynamics
      ;;  followed by structure optimization
      (dotimes (i 1000) (advance-simulation dynamics :unfrozen nil))
      (dynamics:write-coordinates-back-to-matter dynamics)
      (when verbose
        (format t "Optimizing structure~%"))
      (multiple-value-bind (matter energy-function total-energy)
          (if energy-function
              (optimize-structure-with-restarts agg :turn-off-nonbond nil :verbose verbose :energy-function energy-function
                                                    :max-sd-steps max-sd-steps
                                                    :max-cg-steps max-cg-steps
                                                    :max-tn-steps max-tn-steps
                                                    :sd-tolerance sd-tolerance
                                                    :cg-tolerance cg-tolerance
                                                    :tn-tolerance tn-tolerance
                                                    )
              (optimize-structure-with-restarts agg :turn-off-nonbond nil :verbose verbose
                                                    :max-sd-steps max-sd-steps
                                                    :max-cg-steps max-cg-steps
                                                    :max-tn-steps max-tn-steps
                                                    :sd-tolerance sd-tolerance
                                                    :cg-tolerance cg-tolerance
                                                    :tn-tolerance tn-tolerance
                                                    ))
        (unless energy-function
          (error "starting-geometry-with-restarts energy-function is nil"))
        (values dynamics energy-function total-energy)
        ))))

(defparameter *build-agg* nil)
(defun build-good-geometry-from-random (agg &key verbose (no-errors t))
  (let (bad-geom)
    (dotimes (i 10)
      (format t "Attempt ~a to build good geometry from a random starting point~%" i)
      (scramble-positions agg)
      (setf *build-agg* agg)
      (if no-errors
          (handler-case
              (optimize-structure agg :verbose verbose)
            (error (e) (warn " ~a" e))
            )
          (optimize-structure agg :verbose verbose))
      (setf bad-geom (topology:bad-geometry-p agg))
      (when (not bad-geom)
        (return-from build-good-geometry-from-random)))
    (error "Exceeded max number of tries to build good geometry~%Bad geometry:~%~a" bad-geom)))


(defun save-mol2 (matter pathname &key (use-sybyl-types t))
  "Save MATTER as a mol2 to the file PATHNAME.  USE-SYBYL-TYPES of nil requires atom types be set and then they will be used."
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*)))
    (format t "Saving matter to ~a~%" npn)
    (chem:save-mol2 matter npn use-sybyl-types)))

(defun load-mol2 (pathname)
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*)))
    (unless (probe-file npn)
      (error "Could not open ~a" npn))
    (chem:load-mol2 npn)))


(defun simple-vector-coordinate-for-atomspec (agg &optional atomspec)
  "Return coordinates of AGG.  If ATOMSPEC is not nil then use it to downselect atoms."
  (let ((atoms (if atomspec
                   (chem:atoms-with-chimera-specifications agg atomspec)
                   (chem:matter/all-atoms-as-list agg))))
    (chem:make-simple-vector-coordinate-from-atom-sequence atoms)))

(defun superpose-one (&key fixed-atoms-list moveable-matter moveable-atoms-list)
  "Calculate the transform to align the atoms in MOVEABLE-ATOMS-LIST to match the positions
in FIXED-ATOMS-LIST and apply it to MOVEABLE-MATTER"
  (let ((coords-fixed (chem:make-simple-vector-coordinate-from-atom-sequence fixed-atoms))
        (coords-moveable (chem:make-simple-vector-coordinate-from-atom-sequence moveable-atoms))
        (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer coords-fixed)
    (chem:set-moveable-all-points superposer coords-moveable)
    (let ((transform (chem:superpose superposer)))
      (chem:apply-transform-to-atoms moveable-matter transform))))

(defun superpose-all (aggs &optional atomspec)
  (let* ((fixed-agg (car aggs))
         (rest-aggs (cdr aggs))
         (fixed-coords (simple-vector-coordinate-for-atomspec fixed-agg atomspec))
         (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer fixed-coords)
    (dolist (moveable-agg rest-aggs)
      (let ((moveable-coords (simple-vector-coordinate-for-atomspec moveable-agg atomspec)))
        (chem:set-moveable-all-points superposer moveable-coords)
        (let ((transform (chem:superpose superposer)))
          (chem:apply-transform-to-atoms moveable-agg transform))))))


(defun superpose-against-template (&key fixed-aggregate fixed-atomspec moveable-aggregates moveable-atomspec key)
  (let* ((rms-results (make-hash-table :test #'eq))
         (fixed-coords (simple-vector-coordinate-for-atomspec fixed-aggregate fixed-atomspec))
         (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer fixed-coords)
    (dolist (moveable-agg moveable-aggregates)
      (let ((moveable-coords (simple-vector-coordinate-for-atomspec moveable-agg moveable-atomspec)))
        (chem:set-moveable-all-points superposer moveable-coords)
        (let ((transform (chem:superpose superposer)))
          (chem:apply-transform-to-atoms moveable-agg transform)
          (let ((k (if key
                       (funcall key moveable-agg)
                       moveable-agg)))
            (setf (gethash k rms-results) (chem:root-mean-square-difference superposer)))
          )))
    rms-results))
#|
(defun superpose-one-test (agg-fixed agg-moveable)
  (let ((selected-atoms (core:make-cxx-object 'chem:superpose-selected-atoms))
        (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-matter-with-selected-atoms selected-atoms agg-fixed)
    (let ((coords-fixed (chem:extract-coordinates selected-atoms agg-fixed)))
      (chem:set-matter-with-selected-atoms selected-atoms agg-moveable)
       (let ((coords-moveable (chem:extract-coordinates selected-atoms agg-moveable)))
        (chem:set-fixed-all-points superposer coords-fixed)
        (chem:set-moveable-all-points superposer coords-moveable)
        (format t "fixed-points ~a~%" (chem:get-number-of-fixed-points superposer))
        (format t "moveable-points ~a~%" (chem:get-number-of-moveable-points superposer))
         (let ((transform (chem:superpose superposer)))
          (chem:apply-transform-to-atoms agg-moveable transform))))))
|#
    
(defun superpose-one-test2 (agg-fixed agg-moveable)
  (let ((selected-atoms (core:make-cxx-object 'chem:superpose-selected-atoms))
        (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-matter-with-selected-atoms selected-atoms agg-fixed)
    (chem:copy-matter-coordinates-into-moveable-coordinates selected-atoms agg-fixed)
    (chem:set-matter-with-selected-atoms selected-atoms agg-moveable)
    (chem:copy-matter-coordinates-into-fixed-coordinates selected-atoms agg-moveable)
    (format t "fixed-points ~a~%" (chem:get-number-of-fixed-points superposer))
    (format t "moveable-points ~a~%" (chem:get-number-of-moveable-points superposer))
    (let* ((transform (chem:superpose superposer))
           (rms (chem:root-mean-square-difference superposer)))
      (chem:apply-transform-to-atoms agg-moveable transform)
      rms)))

(defun gnuplot-data (data)
  (let ((data-file "/tmp/gnuplot_temp.xy")
        (output-file "/tmp/out.png"))
    ;; Write data to temporary file
    (with-open-file (fout data-file :direction :output :if-exists :supersede)
      (let ((xi 0) x y)
        (dolist (one data)
          (if (consp one)
              (setf x (first one)
                    y (second one))
              (setf x (prog1 xi (incf xi))
                    y one))
          (format fout "~s ~s~%" x y))))
    
    ;; Run gnuplot command
    (uiop:run-program 
     (list "gnuplot" "-e" 
           (format nil "set terminal png; set output ~s; plot ~s w lp; quit"
                   output-file data-file))
     :output t :error-output t)
    
    ;; Open the output file (platform-specific)
    #+darwin (uiop:run-program (list "open" output-file))
    #+linux (uiop:run-program (list "xdg-open" output-file))
    #+windows (uiop:run-program (list "start" output-file) :shell t)
    
    #+(or)(swank::send-to-emacs (list :write-image output-file "x"))))

(defun find-aggregate (agg-list matcher)
  (dolist (a agg-list)
    (when (funcall matcher a)
      (return-from find-aggregate a))))



(defun jostle (matter &optional (width 20.0) from-zero &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (atom)
           (when from-zero
             (chem:set-position atom (geom:vec 0.0 0.0 0.0)))
           (let* ((cp (if (geom:is-defined (chem:get-position atom))
                          (chem:get-position atom)
                          (geom:vec 0.0 0.0 0.0)))
                  (pos (geom:vec
                       (+ (- (random width) half-width) (geom:get-x cp))
                       (+ (- (random width) half-width) (geom:get-y cp))
                       (+ (- (random width) half-width) (geom:get-z cp)))))
             (chem:set-position atom pos))))
    (cond
      ((typep matter 'chem:matter)
       (chem:map-atoms
        nil
        (lambda (a) (jostle-atom a)) matter))
      ((hash-table-p matter)
       (maphash (lambda (a v)
                  (declare (ignore v))
                  (jostle-atom a)) matter))
      ((listp matter)
       (mapc (lambda (a) (jostle-atom a)) matter))
      (t (error "Add support to jostle positions for ~s" matter)))))


    
(defun assign-configuration (agg configuration)
  (let* ((stereocenters (sort (cando:gather-stereocenters agg) #'string< :key #'chem:get-name)))
    (cando:set-stereoisomer stereocenters configuration :show t)))

(defun assign-configuration-simple (agg configuration)
  (let* ((stereocenters (sort (cando:gather-stereocenters agg) #'string< :key #'chem:get-name)))
    (cando:set-stereoisomer stereocenters configuration :show t)))

(defun atoms-to-fix (agg)
  "Return all atoms with the :fix property sorted by the value of that property"
  (let* ((fix-atoms (sort (select:atoms-with-property agg :fix) #'string<
                          :key (lambda (a) (string (getf (chem:properties a) :fix))))))
    fix-atoms))

(defun clear-all-atom-restraints (matter)
  (chem:map-atoms nil (lambda (a) (chem:clear-restraints a)) matter))

(defun anchor-atom (atom pos &optional (weight 1.0))
  "Clear any restraints on an atom and create an anchor restraint for it"
  (let ((restraint-anchor (core:make-cxx-object 'chem:restraint-anchor :atom atom :position pos :weight weight)))
    (chem:clear-restraints atom)
    (chem:add-restraint atom restraint-anchor)
    restraint-anchor))

(defun anchor-atoms (fix-atoms anchors &optional (weight 0.1))
  "Anchor a collection of atoms"
  (loop for atom in fix-atoms
     for pos in anchors
     do (anchor-atom atom pos weight)))

(defun circle-points (radius points &key (z 0.0))
  "Generate a list of (points) points around a circle centered on origin"
  (loop for aindex from 0 below points
     for angrad = (* 360.0 (/ (float aindex) points) 0.0174533)
     collect (geom:vec (* (cos angrad) radius)
                       (* (sin angrad) radius) z)))

(defun anchor-stereocenters-to-circle (agg)
  "Create an anchor restraint for every stereocenter to points on a circle"
  (let* ((centers (stereocenters-sorted-by-name agg))
         (points (circle-points 15.0 (length centers))))
    (format t "centers -> ~a~%" centers)
    (format t "points -> ~a~%" points)
    (chem:clear-restraints agg)
    (mapc (lambda (a p)
            (format t "anchoring ~a to ~a~%" a p)
            (anchor-atom a p))
          (coerce centers 'list) points)))
