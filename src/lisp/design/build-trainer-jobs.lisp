(in-package :design)


(defmacro save-values (&body body)
  `(progn
     ,@body))

(defun trainer-name (trainer-index sequence)
  (format nil "trainer-~a-~a" trainer-index sequence))

(defvar *ef*)
(defun build-trainer-job (&key design
                            (trainer-index 0)
                            (sequence-index 0)
                            pathname-template)
  (format t "*default-pathname-defaults* -> ~s~%" *default-pathname-defaults*)
  (let* ((trainer (elt (design:trainers design) trainer-index))
         (oligomer (design:oligomer trainer)))
    (leap:easy-gaff)
    (chem:goto-sequence oligomer sequence-index)
    (let* ((molecule (chem:get-molecule oligomer))
           (tname (pathname-name pathname-template))
           (aggregate (chem:make-aggregate (intern (pathname-name pathname-template) :keyword))))
      (chem:add-molecule aggregate molecule)
      (cando:jostle aggregate)
      (let ((sequence-filename (make-pathname :type "mol2" :name tname :defaults pathname-template)))
        (ensure-directories-exist sequence-filename)
        (format t "Minimizing ~a~%" tname)
        (let ((energy-fn (chem:make-energy-function aggregate :assign-types t)))
          (loop named min
                for iteration below 10
                do (format t "Trying build ~s~%" iteration)
                   (energy:minimize-energy-function-from-bad-geometry energy-fn
                                                                      :cg-tolerance 0.1
                                                                      :max-cg-steps 10000
                                                                      :tn-tolerance 0.001
                                                                      :max-tn-steps 100)
                if (null (chem:check-for-beyond-threshold-interactions energy-fn 0.2))
                  do (return-from min nil)
                else
                  do (format t "Bad interactions -> ~d~%" (chem:check-for-beyond-threshold-interactions energy-fn 0.2))
                     (save-values (setf *ef* energy-fn))))
        (cando:save-mol2 aggregate sequence-filename :use-sybyl-types t)
        )
      molecule)))


(defun run-am1-calculation (aggregate pathname-template &key (maxcyc 9999))
  (let ((sqm-input (make-pathname :type "sqm-in" :defaults pathname-template))
        (sqm-output (make-pathname :type "sqm-out" :defaults pathname-template)))
    (charges:calculate-am1-bcc-charges aggregate :maxcyc maxcyc
                                                 :input-pathname sqm-input
                                                 :output-pathname sqm-output)))

(defvar *mol*)
(defvar *pathname-template*)
(defconstant +sqm-maxcyc+ 10)

(defun build-charge-solvate (&rest args &key (design-pathname #P"./")
                                          (trainer-index 0)
                                          (sequence-index 0))
  (leap:setup-default-paths)
  (let* ((tname (trainer-name trainer-index sequence-index))
         (pathname-template (make-pathname :directory '(:relative)
                                           :name (trainer-name trainer-index sequence-index)))
         (design (design:load-design design-pathname))
         (mol (build-trainer-job :design design
                                 :pathname-template pathname-template
                                 :trainer-index trainer-index
                                 :sequence-index sequence-index))
         (aggregate (chem:make-aggregate (chem:get-name mol))))
    (chem:add-molecule aggregate mol)
    (run-am1-calculation aggregate design-pathname :maxcyc +sqm-maxcyc+)
    (leap:load-off "solvents.lib")
    (leap:solvate-box aggregate (leap.core:lookup-variable :tip3pbox) 12.0 0.75)
    (cando:save-mol2 aggregate (make-pathname :name (format nil "~a-solvated" tname) :type "mol2"))
    (save-values
     (setq *mol* mol)
     (setq *pathname-template* pathname-template))))

(defparameter *build-charge-solvate*
  (cl:format nil "~{~s~%~}"
             '((in-package :cando-user)
               (defparameter *design-pathname* #P":%INPUT%")
               (defparameter *trainer-index* :%TRAINER-INDEX%)
               (defparameter *sequence-index* :%SEQUENCE-INDEX%)
               (design:build-charge-solvate :design-pathname *design-pathname*
                                            :trainer-index *trainer-index*
                                            :sequence-index *sequence-index*)
               )))


(defclass trainer-node (amber:node-file)
  ((trainer-index :initarg :trainer-index :accessor trainer-index)
   (sequence-index :initarg :sequence-index :accessor sequence-index)))

(defclass trainer-script (trainer-node amber:script-mixin)
  ())

(defclass trainer-job (amber:job)
  ((trainer-index :initarg :trainer-index :accessor trainer-index)
   (sequence-index :initarg :sequence-index :accessor sequence-index)))


(defmethod amber:substitutions ((job trainer-job) node-file)
  (list* (cons :%TRAINER-INDEX% (format nil "~a" (trainer-index job)))
         (cons :%SEQUENCE-INDEX% (format nil "~a" (sequence-index job)))
         (call-next-method)))

(defmethod amber:substitutions (job (node-file trainer-script))
  (list* (cons :%TRAINER-INDEX% (format nil "~a" (trainer-index job)))
         (cons :%SEQUENCE-INDEX% (format nil "~a" (sequence-index job)))
         (call-next-method)))

  
(defun build-all-trainers (design top-directory)
  (let ((*default-pathname-defaults* (merge-pathnames top-directory))
        (design-pathname (merge-pathnames #P"design.dat")))
    (ensure-directories-exist *default-pathname-defaults*)
    (format t "design-pathname -> ~s~%" design-pathname)
    (format t "Number of trainers = ~s~%" (length (design:trainers design)))
    (design:save-design design design-pathname)
    (let* ((design-node (amber:make-node-file design-pathname))
           (jupyter-job (amber:connect-graph (make-instance 'amber:jupyter-job :outputs (amber:arguments :output #P"design.dat")))))
      (loop for trainer-index from 0 below (length (design:trainers design))
            for trainer in (design:trainers design)
            for oligomer = (design:oligomer trainer)
            do (loop for sequence-index below (chem:number-of-sequences oligomer)
                     for tname = (trainer-name trainer-index sequence-index)
                     for trainer-script = (make-instance 'trainer-script
                                                         :trainer-index trainer-index
                                                         :sequence-index sequence-index
                                                         :script *build-charge-solvate*
                                                         :filename (make-pathname :name tname :type "lisp"))
                     for parm-file = (make-instance 'trainer-node
                                                    :trainer-index trainer-index
                                                    :sequence-index sequence-index
                                                    :filename (make-pathname :type "parm"))
                     for crd-file = (make-instance 'trainer-node
                                                   :trainer-index trainer-index
                                                   :sequence-index sequence-index
                                                   :filename (make-pathname :type "crd"))
                     for trainer-job = (amber:connect-graph
                                        (make-instance 'trainer-job
                                                       :trainer-index trainer-index
                                                       :sequence-index sequence-index
                                                       :script trainer-script
                                                       :inputs (amber:arguments :%INPUT% design-node)
                                                       :outputs (amber:arguments :parm parm-file
                                                                                 :crd crd-file)))))
      (amber:generate-all-code jupyter-job nil :pathname-defaults *default-pathname-defaults*))))
