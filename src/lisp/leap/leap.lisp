(in-package :leap)

(defun setup-amber-paths ()
  (leap.core::clear-path)
  (leap.core:add-path "amber:dat;leap;prep;")
  (leap.core:add-path "amber:dat;leap;lib;")
  (leap.core:add-path "amber:dat;leap;parm;")
  (leap.core:add-path "amber:dat;leap;cmd;")
  (leap.core:add-path "amber:dat;antechamber;"))

(defun setup-default-paths ()
  (setup-amber-paths)
  (leap.core:add-path *default-pathname-defaults*)
  (show-paths))

(defun show-paths ()
  (format t "Paths: ~a~%" leap.core::*path*))

(defun assign-atom-types (matter)
  (chem:map-molecules
   nil
   (lambda (molecule)
     (let* ((force-field-name (chem:force-field-name molecule))
            (combined-force-field (chem:find-force-field force-field-name))
            (chem:*current-rings* (chem:identify-rings molecule)))
       (chem:assign-force-field-types combined-force-field molecule)))
   matter))


(defun list-force-fields ()
  (maphash (lambda (k v)
             (format t "~s ~s~%" k (mapcar #'identity (chem:force-fields-as-list v))))
           leap.core:*force-fields*))

(defun source (filename)
"    source filename
      STRING                       _filename_

This command executes LEaP commands within a text file.  To display the
commands as they are read, see the verbosity command.  The text within
the source file must be formatted exactly like the text the user types
into LEaP. This command temporarily adds the path of _filename_ to the 
leap path list so that files local to the script can be loaded.
"
  (let* ((path (let ((p (leap.core:ensure-path filename)))
                 (unless p (error "Could not find file ~a" filename))
                 p))
         (entire-file (alexandria:read-file-into-string path))
         (ast (architecture.builder-protocol:with-builder ('list)
                (esrap:parse 'leap.parser:leap entire-file))))
    (leap.core:with-path (namestring (make-pathname :directory (pathname-directory path)))
      (leap.core:evaluate 'list ast leap.core:*leap-env*)))
  t)

(defun save-amber-parm (aggregate topology-file-name &optional crd-pathname show-progress)
  (let* ((top-pathname (merge-pathnames (pathname topology-file-name)))
         (crd-pathname (if crd-pathname
                           (merge-pathnames crd-pathname)
                           (make-pathname :type "crd" :defaults top-pathname))))
    (let ((chem:*verbose* show-progress))
      (leap.topology:save-amber-parm-format aggregate top-pathname crd-pathname :assign-types t))))


(defun load-amber-params (filename &optional (force-field :default))
  "    variable = loadAmberParams filename
      PARMSET                      _variable_
      STRING                       _filename_

Load an AMBER format parameter set file and place it in _variable_.
All interactions defined in the parameter set will be contained within
_variable_.   This command causes the loaded parameter set to be included
in LEaP's list of parameter sets that are searched when parameters are
required.  General proper and improper torsion parameters are modified,
the AMBER general type \"X\" is replaced with the LEaP general type \"?\".
"
  (let* ((fn (leap.core:ensure-path filename))
         (ff (let ((parmreader (chem:make-read-amber-parameters)))
               (with-open-file (fin fn :direction :input)
                 (chem:read-parameters parmreader fin leap:*amber-system*)
                 (chem:get-force-field parmreader)))))
    (chem:set-title ff filename)
    (format *standard-output* "Adding force field ~a to ~a~%" filename force-field)
    (leap.core:add-force-field-or-modification ff
                                               :force-field-name force-field
                                               :force-field-info filename
                                               :combined-force-field-class-name 'chem:combined-force-field)))

(defun load-smirnoff-params (filename)
  (let* ((path (leap.core:ensure-path filename))
         (force-field (smirnoff:load-smirnoff path)))
    (leap.core:add-force-field-or-modification force-field
                                               :force-field-name :smirnoff
                                               :force-field-info path
                                               :combined-force-field-class-name 'smirnoff:combined-force-field)))

(defun set-force-field (matter force-field-name)
  (chem:find-force-field force-field-name)
  (unless (or (typep matter 'chem:aggregate) (typep matter 'chem:molecule))
    (error "You can only set the force-field-name of an aggregate or a molecule"))
  (chem:setf-force-field-name matter force-field-name))
  
(defun load-atom-type-rules (filename &optional (force-field :default))
  (let* ((path (leap.core:ensure-path filename))
         (fftypedb (with-open-file (fin path)
                       (leap.antechamber-type-definition-parser:read-antechamber-type-rules fin)))
         (ff (core:make-cxx-object 'chem:force-field)))
    (chem:set-title ff filename)
    (chem:set-type-db ff fftypedb)
    (leap.core:add-force-field-or-modification ff :force-field-name force-field
                                                  :force-field-info filename
                                                  :combined-force-field-class-name 'chem:combined-force-field)))

(defun solvate-box (solute solvent buffer &key isotropic closeness)
  (when (numberp buffer)
    (setf buffer (list buffer buffer buffer)))
  (leap.solvate:tool-solvate-and-shell solute solvent buffer :closeness closeness :isotropic isotropic))

(defun solvate-oct (solute solvent buffer &key isotropic closeness)
  (when (numberp buffer)
    (setf buffer (list buffer buffer buffer)))
  (leap.solvate:tool-solvate-and-shell solute solvent buffer :oct t :isotropic isotropic :closeness closeness))


(defun solvate-shell (solute solvent farness &key isotropic closeness)
  (let ((buffer (list farness farness farness)))
    (leap.solvate:tool-solvate-and-shell solute solvent buffer :closeness closeness :farness farness :shell t :isotropic isotropic)))

(defun solvate-cap (solute solvent position radius &key closeness)
  (let ((count 0)
        (x-center 0.0)
        (y-center 0.0)
        (z-center 0.0)
        center)
    (if (listp position)
        (cond ((or (chem:is-aggregate (first position))
                   (chem:is-molecule (first position))
                   (chem:is-residue (first position))
                   (chem:is-atom (first position)))
               (loop for i from 0 below (list-length position) 
                     for x = (geom:vx (chem:geometric-center (nth i position)))
                     for y = (geom:vy (chem:geometric-center (nth i position)))
                     for z = (geom:vz (chem:geometric-center (nth i position)))
                     do (progn
                          (setf x-center (+ x-center x)
                                y-center (+ y-center y)
                                z-center (+ z-center z))
                          (incf count)))
               (setf x-center (/ x-center count)
                     y-center (/ y-center count)
                     z-center (/ z-center count))
               (setf (first center) x-center
                     (second center) y-center
                     (third center) z-center))
              ((numberp (first position))
               (geom:set-all3 center (first position) (second position) (third position)))
              (t (error "Position must be AGGREGATE, RESIDUE, MOLECULE, ATOM, LIST of atoms, or LIST of 3 double.")))
        (if (or (chem:is-aggregate position)
                (chem:is-molecule position)
                (chem:is-residue position)
                (chem:is-atom position))
            (setf center (chem:geometric-center position))
            (error "Position must be AGGREGATE, RESIDUE, MOLECULE, ATOM, LIST of atoms, or LIST of 3 double.")))
    (leap.solvate:tool-solvate-in-sphare solute solvent center radius :closeness closeness)))

(defstruct molecule-info name type force-field-name count)

(defun desc (matter)
  (flet ((describe-molecule (mi)
           (format t "~20s ~10s ~20s" (molecule-info-name mi)
                   (molecule-info-type mi)
                   (molecule-info-force-field-name mi))
           (when (> (molecule-info-count mi) 1)
             (format t "~5d" (molecule-info-count mi)))
           (format t "~%")))
    (cond
      ((typep matter 'chem:aggregate)
       (let ((mol-ht (make-hash-table :test #'equal)))
         (cando:do-molecules (mol matter)
           (let* ((key (list (chem:get-name mol) (chem:force-field-name mol) (chem:molecule-type mol)))
                  (entry (gethash key mol-ht)))
             (if entry
                 (progn
                   (incf (molecule-info-count entry)))
                 (setf (gethash key mol-ht) (make-molecule-info :name (chem:get-name mol)
                                                                :type (chem:molecule-type  mol)
                                                                :force-field-name (chem:force-field-name mol)
                                                                :count 1)))))
         (let ((molecule-list nil))
           (maphash (lambda (key value)
                      (push value molecule-list))
                    mol-ht)
           (loop for molecule-info in molecule-list
                 when (not (eq (molecule-info-type molecule-info) :solvent))
                   do (describe-molecule molecule-info))
           (loop for molecule-info in molecule-list
                 when (eq (molecule-info-type molecule-info) :solvent)
                   do (describe-molecule molecule-info)))))
      ((typep matter 'chem:molecule)
       (describe-molecule (make-molecule-info :name (chem:get-name matter)
                                              :type (chem:molecule-type matter)
                                              :force-field-name (chem:force-field-name matter)
                                              :count 1)))
      (t (format t "~s~%" matter)))))

