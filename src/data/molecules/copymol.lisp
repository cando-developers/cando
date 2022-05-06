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
                    (format t "~6a ~a ~a ~a~%" (chem:get-name x) (chem:get-type x) (chem:get-configuration x) (chem:get-position x)))
                  matter))

(defun dump-stereocenters (centers)
  (format t "-------- Stereocenters~%")
  (dotimes (i (length centers))
    (let ((a (elt centers i)))
      (format t "~6a ~a~%" (chem:get-name a) (chem:get-configuration a)))))

(defun set-stereoisomer (atom-vec index)
  (dotimes (i (length atom-vec))
    (let ((config (if (logbitp i index) 'chem:S 'chem:R))
          (atom (elt atom-vec i)))
      (chem:set-configuration (elt atom-vec i) config))))

(defun gather-stereocenters (matter)
  (let (chiral-atoms)
    (chem:map-atoms nil (lambda (a &aux p)
                          (when (eq (chem:get-stereochemistry-type a) 'chem::chiral)
                            (push a chiral-atoms)))
                    matter)
    (make-array (length chiral-atoms) :initial-contents chiral-atoms)))


(defun number-of-stereoisomers (chiral-atoms)
  (expt 2 (length chiral-atoms)))


(defun safe-minimize (minimizer)
  (restart-case (chem:minimize minimizer)
    (skip-rest-of-minimization (something)
      :report "Skip the rest of the current minimization - continue processing"
      (chem:write-intermediate-results-to-energy-function minimizer))))

(defun build-starting-structure-from-random (matter)
  (scramble-positions matter)
  (let* ((energy-function (chem:make-energy-function matter *ff*))
         (min (chem:make-minimizer :energy-function energy-function)))
    (chem:enable-print-intermediate-results min)
    (chem:save-mol2 matter #P"/tmp/before.mol2")
    ;;(defparameter min (chem:make-minimizer :energy-function *energy-function*))
    (format t "Starting minimization stage1~%")
    (safe-minimize min)
    (format t "Saving stage1~%")
    (chem:save-mol2 matter #P"/tmp/after1.mol2")
    (chem:set-option energy-function 'chem::nonbond-term nil)
    (format t "Starting minimization stage2~%")
    (safe-minimize min)
    (format t "Saving stage2~%")
    (chem:save-mol2 matter #P"/tmp/after2.mol2")
    (chem:set-option energy-function 'chem:nonbond-term t)
    (format t "Starting minimization stage3~%")
    (safe-minimize min)
    (format t "Saving stage3~%")
    (chem:save-mol2 matter #P"/tmp/after3.mol2")))

(format t "Function definitions done~%")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Work with a molecule

(getpid)35582
(setq *default-pathname-defaults*
      #P"~/Development/Clasp/projects/cando/work/molecules/")

;; Load AMBER parameters and Antechamber type rules
(let ((*default-pathname-defaults*
       #P"~/Development/clasp/projects/cando/work/forceField/"))
  (defparameter *parms*
    (let ((parms (chem:make-read-amber-parameters)))
      (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
        (chem:read-types parms fin))
      (with-open-file (fin "gaff.dat" :direction :input)
        (chem:read-parameters parms fin)
        parms)))
  (defparameter *ff* (chem:get-force-field *parms*)))

;; Load the chemdraw file containing a molecule
(let ((*default-pathname-defaults*
       #P"~/Development/Clasp/projects/cando/work/molecules/"))
  (defparameter *cd* (with-open-file (fin "onemol2.cdxml" :direction :input)
                       (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*))
  (scramble-positions *agg*)
  (chem:save-mol2 *agg* #P"/tmp/initial.mol2"))

;; Gather the atoms that are stereocenters into a vector
(defparameter *stereocenters* (gather-stereocenters *agg*))
(format t "~a~%" *stereocenters*)

;; Pick a particular stereoisomer
(set-stereoisomer *stereocenters* 255)
(dump-stereocenters *stereocenters*)

(chem:matter-copy *agg*)
