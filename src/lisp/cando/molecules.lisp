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
                       (+ (- (random width) half-width) (geom:vx center))
                       (+ (- (random width) half-width) (geom:vy center))
                       (+ (- (random width) half-width) (geom:vz center)))))
             (chem:set-position atom pos))))
    (cond
      ((typep matter 'chem:matter)
       (chem:map-atoms
        nil
        (lambda (a) (scramble-atom a)) matter))
      ((hash-table-p matter)
       (maphash (lambda (a v) (scramble-atom a)) matter))
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
      (chem:set-configuration (elt atom-vec i) config)))
  (when show
    (let ((bits (ceiling (length atom-vec))))
      (format t "======== Stereoisomer ~3d/~v,'0b~%" index bits index)
      (dump-stereocenters atom-vec))))

(defun gather-stereocenters (matter)
  (let (chiral-atoms)
    (chem:map-atoms nil (lambda (a &aux p)
                          (when (eq (chem:get-stereochemistry-type a) :chiral)
                            (push a chiral-atoms)))
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

(defmacro with-handle-linear-angles-dihedrals (&body body)
  `(handler-bind
       ((chem::linear-atoms-error
          (lambda (condition)
            (let* ((width 2.0)
                   (half-width (/ width 2.0)))
              (loop for atom in (chem:atoms condition)
                    do (let* ((cp (chem:get-position atom))
                              (pos (geom:vec
                                    (+ (- (random width) half-width) (geom:vx cp))
                                    (+ (- (random width) half-width) (geom:vy cp))
                                    (+ (- (random width) half-width) (geom:vz cp)))))
                         (chem:set-position atom pos))))
            (chem:restart-minimizer))))
     (progn
       ,@body)))

;; Recover from minimization problems using Common Lisp restarts
(defun minimize-no-fail (minimizer)
  (format t "minimize-no-fail~%")
  (restart-case
      (handler-bind
          ((chem:minimizer-error (lambda (c)
                                  (warn "The minimizer reported: ~a" c)
                                  (invoke-restart 'skip-rest-of-minimization))))
        (with-handle-linear-angles-dihedrals
            (chem:minimize minimizer)))
    ;; skip-rest-of-minimization can also be triggered by the user from the debugger
    (skip-rest-of-minimization ()
      :report "Skip the rest of the current minimization - continue processing"
      (format t "Saving intermediate results~%")
      (chem:write-intermediate-results-to-energy-function minimizer))))


(defun minimizer-obey-interrupt (minimizer)
  (restart-case
      (handler-bind
          ((ext:unix-signal-received (lambda (c) (print "Done") (invoke-restart 'skip-rest))))
        (chem:minimize minimizer))
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
                              (tn-tolerance 0.000001))
  (chem:set-maximum-number-of-steepest-descent-steps minimizer max-sd-steps)
  (chem:set-maximum-number-of-conjugate-gradient-steps minimizer max-cg-steps)
  (chem:set-maximum-number-of-truncated-newton-steps minimizer max-tn-steps)
  (chem:set-steepest-descent-tolerance minimizer sd-tolerance)
  (chem:set-conjugate-gradient-tolerance minimizer cg-tolerance)
  (chem:set-truncated-newton-tolerance minimizer tn-tolerance))
  
                              
(defun optimize-structure (matter &optional active-atoms system)
  (error "This function is broken because make-energy-function uses keywords")
  (let* ((energy-function (chem:make-energy-function matter system active-atoms))
         (min (chem:make-minimizer energy-function)))
    (configure-minimizer min
                         :max-steepest-descent-steps 1000
                         :max-conjugate-gradient-steps 50000
                         :max-truncated-newton-steps 100)
    (chem:enable-print-intermediate-results min)
    (chem:set-option energy-function 'chem::nonbond-term nil)
    (format t "Starting minimization stage 1 nonbond=NIL~%")
    (minimize-no-fail min)
    (chem:set-option energy-function 'chem:nonbond-term t)
    (format t "Starting minimization stage 2 nonbond=T~%")
    (minimize-no-fail min))
  matter)

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
         (accumulate-aggregate (core:make-cxx-object 'chem:aggregate))
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

    
(defun bad-geometry-p (agg &optional system)
  "Return true if there are any beyond-threshold force field interactions"
  (error "This function is broken because make-energy-function uses keywords")
  (let ((energy-function (chem:make-energy-function agg system)))
    (let ((fails (chem:check-for-beyond-threshold-interactions system)))
      (if (> fails 0)
          fails
          nil))))

(defun build-good-geometry-from-random (agg force-field)
  (dotimes (i 20)
    (format t "Attempt ~a to build good geometry from a random starting point~%" i)
    (scramble-positions agg)
    (optimize-structure agg force-field)
    (when (not (bad-geometry-p agg force-field))
        (return-from build-good-geometry-from-random)))
  (error "Exceeded max number of tries to build good geometry"))


(defun save-mol2 (matter pathname &key use-sybyl-types)
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*)))
    (format t "Saving matter to ~a~%" npn)
    (chem:save-mol2 matter npn use-sybyl-types)))

(defun load-mol2 (pathname &key verbose)
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*))
        (chem:*verbose* verbose))
    (chem:load-mol2 npn)))


(defun simple-vector-coordinate-for-atomspec (agg atomspec)
  (let ((atoms (chem:atoms-with-chimera-specifications agg atomspec)))
    (chem:make-simple-vector-coordinate-from-atom-list atoms)))

(defun superpose-one (&key fixed-atoms moveable-matter moveable-atoms)
  (let ((coords-fixed (chem:make-simple-vector-coordinate-from-atom-list fixed-atoms))
        (coords-moveable (chem:make-simple-vector-coordinate-from-atom-list moveable-atoms))
        (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer coords-fixed)
    (chem:set-moveable-all-points superposer coords-moveable)
    (let ((transform (chem:superpose superposer)))
      (chem:apply-transform-to-atoms moveable-matter transform))))

(defun superpose-all (aggs atomspec)
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
        (let ((transform (chem:superpose superposer)))
           (chem:apply-transform-to-atoms agg-moveable transform))))

(defun gnuplot-data (data)
  (with-open-file (fout "/tmp/gnuplot_temp.xy" :direction :output)
    (let ((xi 0) x y)
      (dolist (one data)
        (if (consp one)
            (setf x (first one)
                  y (second one))
            (setf x (prog1 xi (incf xi))
                  y one))
        (format fout "~s ~s~%" x y))))
  (ext:system "gnuplot -e 'set terminal png; set output \"/tmp/out.png\"; plot \"/tmp/gnuplot_temp.xy\" w lp; quit'")
  (ext:system "open /tmp/out.png")
  #+(or)(swank::send-to-emacs (list :write-image "/tmp/out.png" "x")))

(defun find-aggregate (agg-list matcher)
  (dolist (a agg-list)
    (when (funcall matcher a)
      (return-from find-aggregate a))))



(defun jostle (matter &optional (width 20.0) from-zero &aux (half-width (/ width 2.0)))
  "Randomly jostle atoms from their current positions"
  (flet ((jostle-atom (atom)
           (when from-zero
             (chem:set-position atom (geom:vec 0.0 0.0 0.0)))
           (let* ((cp (chem:get-position atom))
                  (pos (geom:vec
                       (+ (- (random width) half-width) (geom:vx cp))
                       (+ (- (random width) half-width) (geom:vy cp))
                       (+ (- (random width) half-width) (geom:vz cp)))))
             (chem:set-position atom pos))))
    (cond
      ((typep matter 'chem:matter)
       (chem:map-atoms
        nil
        (lambda (a) (jostle-atom a)) matter))
      ((hash-table-p matter)
       (maphash (lambda (a v) (jostle-atom a)) matter))
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
