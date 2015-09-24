(in-package :cando)

(progn
  (require :asdf)
  (let ((central-registry (find-symbol "*CENTRAL-REGISTRY*" :asdf))
        (load-system (find-symbol "LOAD-SYSTEM" :asdf)))
    (setf (symbol-value central-registry)
          (cons (translate-logical-pathname "cando:cando;")
                (symbol-value central-registry)))
    (funcall load-system "cando"))
  (format t "Done initialization pid = ~a~%" (getpid)))

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

(defun build-starting-structure-from-random (matter)
  (scramble-positions matter)
  (let* ((energy-function (chem:make-energy-function matter *ff*))
         (min (chem:make-minimizer :energy-function energy-function)))
    (chem:enable-print-intermediate-results min)
    (chem:set-option energy-function 'chem::nonbond-term nil)
    (format t "Starting minimization stage2~%")
    (minimize-no-fail min)
    (format t "Saving stage2~%")
    (chem:set-option energy-function 'chem:nonbond-term t)
    (format t "Starting minimization stage3~%")
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
         (accumulate-aggregate (core:make-cxx-object 'chem:aggregate)))
    (dotimes (i max-files)
      (multiple-value-bind (xi yi)
          (floor i width)
        (let ((xo (* scale xi))
              (yo (* scale yi)))
          (let ((pn (indexed-pathname pathname-template i))
                (target-pn (make-pathname :name (format nil "~a_grid" (pathname-name pathname-template))
                                          :defaults pathname-template)))
            (when (probe-file pn)
              (let* ((agg (chem:load-mol2 pn))
                     (mol (chem:first-molecule agg)))
                (chem:translate-all-atoms
                 mol
                 (make-array 3 :initial-contents (list xo yo 0.0)))
                (chem:add-matter accumulate-aggregate mol)))
            (chem:save-mol2 accumulate-aggregate target-pn)))))))

(defun chimera (obj)
  (cond
    ((pathnamep obj)
     (ext:system (format nil "open ~a" (namestring obj))))
    ((typep obj 'chem:aggregate)
     (chem:save-mol2 obj #P"/tmp/temp.mol2")
     (ext:system (format nil "open /tmp/temp.mol2")))
    (t (error "You cannot run chimera on ~a" obj))))

(defun bad-geometry-p (agg)
  (let ((energy-function (chem:make-energy-function agg *ff*)))
    (let ((fails (chem:check-for-beyond-threshold-interactions energy-function)))
      (if (> fails 0)
          fails
          nil))))

(defun build-good-geometry-from-random (&optional (agg *agg*))
  (dotimes (i 20)
    (format t "Attempt ~a to build good geometry from a random starting point~%" i)
    (scramble-positions agg)
    (build-starting-structure agg)
    (when (not (bad-geometry-p agg))
        (return-from build-good-geometry-from-random)))
  (error "Exceeded max number of tries to build good geometry"))

(defun build-good-geometry-from-random (&optional (agg *agg*))
  (dotimes (i 20)
    (format t "Attempt ~a to build good geometry from a random starting point~%" i)
    (let ((new-agg (build-starting-structure-from-random agg)))
      (when (not (bad-geometry-p new-agg))
        (return-from build-good-geometry-from-random new-agg))))
  (error "Exceeded max number of tries to build good geometry"))
