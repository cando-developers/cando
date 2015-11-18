(in-package :cando)
(defun scramble-positions (matter &key (center (geom:make-ovector3 0.0 0.0 0.0)) (width 10.0)
                           &aux (half-width (/ width 2.0)))
  (flet ((scramble-atom (atom)
           (let ((pos (geom:make-ovector3
                       (+ (- (random width) half-width) (geom:x center))
                       (+ (- (random width) half-width) (geom:y center))
                       (+ (- (random width) half-width) (geom:z center)))))
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
                            (chem:get-properties x)))
                  matter))

(defun atoms-with-property (matter prop-name)
  (let (atoms)
    (chem:map-atoms nil
                    (lambda (x)
                      (when (not (eq :eof (getf (chem:get-properties x) prop-name :eof)))
                        (push x atoms)))
                    matter)
    atoms))


(defun dump-stereocenters (centers)
  (dotimes (i (length centers))
    (let ((a (elt centers i)))
      (format t "~6a ~a~%" (chem:get-name a) (chem:get-configuration a)))))

(defun set-stereoisomer (atom-vec index &key show)
  (dotimes (i (length atom-vec))
    (let ((config (if (logbitp i index) 'chem:S 'chem:R))
          (atom (elt atom-vec i)))
      (chem:set-configuration (elt atom-vec i) config)))
  (when show
    (let ((bits (ceiling (length atom-vec))))
      (format t "======== Stereoisomer ~3d/~v,'0b~%" index bits index)
      (dump-stereocenters atom-vec))))

(defun gather-stereocenters (matter)
  (let (chiral-atoms)
    (chem:map-atoms nil (lambda (a &aux p)
                          (when (eq (chem:get-stereochemistry-type a) 'chem::chiral)
                            (push a chiral-atoms)))
                    matter)
    (make-array (length chiral-atoms) :initial-contents chiral-atoms)))


(defun number-of-stereoisomers (chiral-atoms)
  (expt 2 (length chiral-atoms)))

;; Recover from minimization problems using Common Lisp restarts
(defun minimize-no-fail (minimizer)
  (restart-case
      (handler-case
          (chem:minimize minimizer)
        (chem:minimizer-error (c)
          (warn "The minimizer reported: ~a" c)
          (invoke-restart 'skip-rest-of-minimization)))
    ;; skip-rest-of-minimization can also be triggered by the user from the debugger
    (skip-rest-of-minimization ()
      :report "Skip the rest of the current minimization - continue processing"
      (chem:write-intermediate-results-to-energy-function minimizer))))

(defun configure-minimizer (minimizer
                            &key max-steepest-descent-steps
                              max-conjugate-gradient-steps
                              max-truncated-newton-steps
                              (steepest-descent-tolerance 2000.0)
                              (conjugate-gradients-tolerance 0.5)
                              (truncated-newton-tolerance 0.000001)
                              )
  (when max-steepest-descent-steps
    (chem:set-maximum-number-of-steepest-descent-steps minimizer max-steepest-descent-steps))
  (when max-conjugate-gradient-steps
    (chem:set-maximum-number-of-conjugate-gradient-steps minimizer max-conjugate-gradient-steps))
  (when max-truncated-newton-steps
    (chem:set-maximum-number-of-truncated-newton-steps minimizer max-truncated-newton-steps))
  (chem:set-steepest-descent-tolerance minimizer steepest-descent-tolerance)
  (chem:set-conjugate-gradient-tolerance minimizer conjugate-gradients-tolerance)
  (chem:set-truncated-newton-tolerance minimizer truncated-newton-tolerance)

  )
  
                              
(defun optimize-structure (matter force-field &optional active-atoms)
  (let* ((energy-function (chem:make-energy-function matter force-field active-atoms))
         (min (chem:make-minimizer :energy-function energy-function)))
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

(defun chimera (obj)
  (cond
    ((pathnamep obj)
     (ext:system (format nil "open ~a" (namestring obj))))
    ((typep obj 'chem:aggregate)
     (chem:save-mol2 obj #P"/tmp/temp.mol2")
     (ext:system (format nil "open /tmp/temp.mol2")))
    (t (error "You cannot run chimera on ~a" obj)))
  (sleep 0.3))

(defun bad-geometry-p (agg force-field)
  (let ((energy-function (chem:make-energy-function agg force-field)))
    (let ((fails (chem:check-for-beyond-threshold-interactions energy-function)))
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


(defun save-mol2 (matter pathname)
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*)))
    (format t "Saving matter to ~a~%" npn)
    (chem:save-mol2 matter npn)))

(defun load-mol2 (pathname)
  (let ((npn (merge-pathnames pathname *default-pathname-defaults*)))
    (chem:load-mol2 npn)))


(defun coordinate-array-for-atomspec (agg atomspec)
  (let ((atoms (chem:atoms-with-chimera-specifications agg atomspec)))
    (chem:make-coordinate-array-from-atom-list atoms)))

(defun superpose-one (agg-fixed agg-moveable atomspec)
  (let ((coords-fixed (coordinate-array-for-atomspec agg-fixed atomspec))
        (coords-moveable (coordinate-array-for-atomspec agg-moveable atomspec))
        (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer coords-fixed)
    (chem:set-moveable-all-points superposer coords-moveable)
    (let ((transform (chem:superpose superposer)))
      (chem:apply-transform-to-atoms agg-moveable transform))))

(defun superpose-all (aggs atomspec)
  (let* ((fixed-agg (car aggs))
         (rest-aggs (cdr aggs))
         (fixed-coords (coordinate-array-for-atomspec fixed-agg atomspec))
         (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer fixed-coords)
    (dolist (moveable-agg rest-aggs)
      (let ((moveable-coords (coordinate-array-for-atomspec moveable-agg atomspec)))
        (chem:set-moveable-all-points superposer moveable-coords)
        (let ((transform (chem:superpose superposer)))
          (chem:apply-transform-to-atoms moveable-agg transform))))))


(defun superpose-against-template (&key fixed-aggregate fixed-atomspec moveable-aggregates moveable-atomspec key)
  (let* ((rms-results (make-hash-table :test #'eq))
         (fixed-coords (coordinate-array-for-atomspec fixed-aggregate fixed-atomspec))
         (superposer (core:make-cxx-object 'chem:superpose-engine)))
    (chem:set-fixed-all-points superposer fixed-coords)
    (dolist (moveable-agg moveable-aggregates)
      (let ((moveable-coords (coordinate-array-for-atomspec moveable-agg moveable-atomspec)))
        (chem:set-moveable-all-points superposer moveable-coords)
        (let ((transform (chem:superpose superposer)))
          (chem:apply-transform-to-atoms moveable-agg transform)
          (let ((k (if key
                       (funcall key moveable-agg)
                       moveable-agg)))
            (setf (gethash k rms-results) (chem:root-mean-square-difference superposer)))
          )))
    rms-results))



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



(defun jostle-atoms (matter &optional (width 0.1) &aux (half-width (/ width 2.0)))
  (flet ((jostle-atom (atom)
           (let* ((cp (chem:get-position atom))
                  (pos (geom:make-ovector3
                       (+ (- (random width) half-width) (geom:x cp))
                       (+ (- (random width) half-width) (geom:y cp))
                       (+ (- (random width) half-width) (geom:z cp)))))
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
