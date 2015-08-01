(in-package :cando)


(defun scramble-positions (matter)
  (chem:map-atoms
   nil
   (lambda (a) (chem:randomize-atom-position a) a)
   matter))
  
(defun dump-atoms (matter)
  (chem:map-atoms nil
                  (lambda (x)
                    (format t "~6a ~a ~a ~a~%"
                            (chem:get-name x)
                            (chem:get-type x)
                            (chem:get-configuration x)
                            (chem:get-position x)))
                  matter))

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
                              max-truncated-newton-steps)
  (when max-steepest-descent-steps
    (chem:set-maximum-number-of-steepest-descent-steps minimizer max-steepest-descent-steps))
  (when max-conjugate-gradient-steps
    (chem:set-maximum-number-of-conjugate-gradient-steps minimizer max-conjugate-gradient-steps))
  (when max-truncated-newton-steps
    (chem:set-maximum-number-of-truncated-newton-steps minimizer max-truncated-newton-steps)))
  
                              
(defun optimize-structure (matter force-field)
  (let* ((energy-function (chem:make-energy-function matter force-field))
         (min (chem:make-minimizer :energy-function energy-function)))
    (configure-minimizer min
                         :max-steepest-descent-steps 100
                         :max-conjugate-gradient-steps 500
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
    (t (error "You cannot run chimera on ~a" obj))))

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


(defun save-mol2 (matter pn)
  (let ((npn (merge-pathnames pn *default-pathname-defaults*)))
    (format t "Saving matter to ~a~%" npn)
    (chem:save-mol2 matter npn)))
