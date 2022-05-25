(in-package :leap.commands)


(defun leap-setup-amber-paths ()
  (leap:setup-amber-paths))

(defun leap.setup-default-paths ()
  (leap:setup-default-paths))

(defun leap.show-paths ()
  (leap:show-paths))

(defun leap.show (arg)
  (let ((val (leap.core:lookup-variable arg)))
    (funcall (find-symbol "SHOW" :cando-user) val)))


(defun leap.z-matrix (val list)
  "   Construct the external coordinates for the atoms.

      Arguments:
              [0]     - Container that contains the atoms.
              [1]     - A list of atoms and internal coordinates.

      The entries in the list of atoms and internal coordinates must
      look like:

      a1 a2 b12
      a1 a2 a3 b12 t123
      a1 a2 a3 a4 b12 t123 p1234
      a1 a2 a3 a4 b12 t123 t124 orientation

      Where a1,a2,a3,a4 can be an atom or an atom name which exists
      in the container.
"
  (let ((matter (leap.core:lookup-variable val)))
    (flet ((get-atom (name)
             (let* ((vname (leap.core:lookup-variable name nil :not-found))
                    (atm (if (eq vname :not-found)
                             (chem:first-atom-with-name matter name)
                             vname)))
               (unless (typep atm 'chem:atom)
                 (error "The object ~a -> ~a must represent an atom" name atm))
               atm)))
      (loop for entry in list
            for len = (length entry)
            do (case len
                 (3
                  (let* ((a1 (get-atom (first entry)))
                         (a2 (get-atom (second entry)))
                         (b12 (third entry))
                         (pa2 (chem:get-position a2))
                         (pa1 (geom:build-using-bond b12 pa2)))
                    (chem:set-position a1 pa1)))
                 (5
                  (let* ((a1 (get-atom (first entry)))
                         (a2 (get-atom (second entry)))
                         (a3 (get-atom (third entry)))
                         (b12 (fourth entry))
                         (t123 (* (fifth entry)#|degrees|# (/ PI 180.0)))#|rad|#
                         (pa3 (chem:get-position a3))
                         (pa2 (chem:get-position a2))
                         (pa1 (geom:build-using-bond-angle b12 pa2 t123 pa3)))
                    (chem:set-position a1 pa1)))
                 (7
                  (let* ((a1 (get-atom (first entry)))
                         (a2 (get-atom (second entry)))
                         (a3 (get-atom (third entry)))
                         (a4 (get-atom (fourth entry)))
                         (b12 (fifth entry))
                         (t123 (* (sixth entry)#|degrees|# (/ PI 180.0)))#|rad|#
                         (p1234 (* (seventh entry)#|degrees|# (/ PI 180.0)))#|rad|#
                         (pa4 (chem:get-position a4))
                         (pa3 (chem:get-position a3))
                         (pa2 (chem:get-position a2))
                         (pa1 (geom:build-using-bond-angle-dihedral b12 pa2 t123 pa3 p1234 pa4)))
                    (chem:set-position a1 pa1)))
                 (8
                  (let* ((a1 (get-atom (first entry)))
                         (a2 (get-atom (second entry)))
                         (a3 (get-atom (third entry)))
                         (a4 (get-atom (fourth entry)))
                         (b12 (fifth entry))
                         (t123 (* (sixth entry)#|degrees|# (/ PI 180.0)))#|rad|#
                         (t124 (* (seventh entry)#|degrees|# (/ PI 180.0)))#|rad|#
                         (orientation (eighth entry))
                         (pa4 (chem:get-position a4))
                         (pa3 (chem:get-position a3))
                         (pa2 (chem:get-position a2))
                         (pa1 (geom:build-using-bond-two-angles-orientation pa2 pa3 pa4 b12 t123 t124 orientation)))
                    (chem:set-position a1 pa1)))
                 (otherwise
                  (error "      The entries in the list of atoms and internal coordinates must
      look like:

      a1 a2 b12
      a1 a2 a3 b12 t123
      a1 a2 a3 a4 b12 t123 p1234
      a1 a2 a3 a4 b12 t123 t124 orientation

      Where a1,a2,a3,a4 can be an atom or an atom name which exists
      in the container.
")
                  )))))
  :no-output)

;;; ----------------------------------------------------------------------
;;;
;;; LEaP commands
;;;
(defvar *out* t)


(defun leap-easy-gaff ()
  " easyGaff

Setup the gaff-1 force-field."
  (leap:easy-gaff))

(defun set-variable (entry)
  (let ((variable (first entry))
        (value (third entry)))
    (leap.core:register-variable variable value)))


(defun leap-quickload (system)
  "quickload   system

Download a Cando package from the internet, compile it and make it available for use."
  (ql:quickload system))

(defun leap-log-file (entry)
  "    logFile filename

      STRING                       _filename_

This command opens the file with the file name _filename_ as a log
file.  User input and ALL output is written to the log file.  Output
is written to the log file as if the verbosity level were set to 2."
  (let ((filename (pathname (ensure-string entry))))
    (let* ((log-stream (open filename :direction :output :if-exists :supersede))
           (broadcast (make-broadcast-stream cl:*standard-output* #+(or)ext:+process-standard-output+ log-stream)))
      (setf cl:*standard-output* broadcast))))

(defmethod describe-object ((object chem:residue) stream)
  (cando:do-atoms (atom object)
    (format stream "~20s type: ~a~%" atom (chem:get-type atom))))

(defmethod describe-object ((object chem:topology) stream)
  (let* ((constitution (chem:get-constitution object))
         (constitution-atoms (chem:|Constitution_O::getConstitutionAtoms| constitution)))
    (loop for index from 0 below (chem:number-of-atoms constitution-atoms)
          for atom = (chem:atom-with-id constitution-atoms index)
          do (format t "~a~%" (core:encode atom)))))

(defmethod describe-object ((object chem:molecule) stream)
  (format stream "Force-field-name: ~a~%" (chem:force-field-name object))
  (format stream "   Molecule-type: ~a~%" (chem:molecule-type object))
  (cando:do-residues (res object)
    (format stream "~s~%" res)))

(defmethod describe-object ((object chem:aggregate) stream)
  (cando:do-molecules (mol object)
    (format stream "~s ~a atoms~%" mol (chem:number-of-atoms mol))))


(defgeneric set-property (object property value))

(defmethod set-property (object property value)
  (error "Provide a specialized version of this function for ~s ~s ~s" object property value))

(defun keyword-upcase (name)
  (intern (string-upcase (string name)) :keyword))

(defmethod set-property ((object chem:matter) (property (eql :name)) value)
  (chem:set-name object value))

(defmethod set-property ((object chem:matter) (property (eql :type)) value)
  (chem:set-type object value))

(defmethod set-property ((object chem:matter) (property (eql :element)) value)
  (chem:set-element object (intern (string value) :keyword)))

(defmethod set-property ((object chem:matter) (property (eql :charge)) value)
  (chem:set-charge object value))

(defmethod set-property ((object chem:matter) (property (eql :position)) value)
  (check-type value list)
  (let ((position (chem:get-position object)))
    (geom:set-all3 position (first value)
                   (second value)
                   (third value))))

(defmethod set-property ((object chem:residue) (property (eql :restype)) value)
  (chem:setf-residue-type object value))


(defmethod set-property ((object chem:molecule) (property (eql :force-field-name)) value)
  (let ((force-field (chem:find-force-field value)))
    (chem:setf-force-field-name object value)))

(defmethod set-property ((object chem:molecule) (property (eql :molecule-type)) value)
  (chem:setf-molecule-type object value))

(defun leap.set (object property value)
  (let ((object-fixed (leap.core:lookup-variable object))
        (property-fixed (keyword-upcase property))
        (value-fixed value))
    (set-property object-fixed property-fixed value-fixed))
  :no-output)

(defun leap-desc (name)
  "    desc variable
      object                       _variable_

Print a description of the object.
"
  (let ((val (leap.core:lookup-variable name)))
    (format t "Describing ~a~%" val)
    (cond
      ((null val)
       (format t "~a~%" val))
      (t (describe val)))
    :no-output))



(defgeneric do-remove (a b))

(defmethod do-remove ((a chem:matter) (b chem:matter))
  (error "~a does not directly contain ~a" a b))

(defmethod do-remove ((a chem:aggregate) (b chem:molecule))
  (cando:do-molecules (mol a)
    (when (eq mol b)
      (chem:remove-molecule a b)
      (return-from do-remove b)))
  (call-next-method))

(defmethod do-remove ((a chem:molecule) (b chem:residue))
  (cando:do-residues (res a)
    (when (eq res b)
      (chem:remove-residue a b)
      (return-from do-remove b)))
  (call-next-method))

(defmethod do-remove ((a chem:residue) (b chem:atom))
  (cando:do-atoms (atm a)
    (when (eq atm b)
      (chem:remove-atom-delete-bonds a b)
      (return-from do-remove b)))
  (call-next-method))

(defun leap-remove-matter (matter-a matter-b)
  "       remove a b

      UNIT/RESIDUE/ATOM          _a_
      UNIT/RESIDUE/ATOM          _b_

Remove the object _b_ from the object _a_.   If _b_ is not contained
by _a_ then an error message will be displayed.  This command is
used to remove ATOMs from RESIDUEs, and RESIDUEs from UNITs.   If
the object represented by _b_ is not referenced by some variable
name then it will be destroyed."
  (let ((matter-a (leap.core:lookup-variable matter-a))
        (matter-b (leap.core:lookup-variable matter-b)))
    (do-remove matter-a matter-b)))


(defun leap-combine (list)
  "
    variable = combine list
      object                       _variable_
      LIST                         _list_

Combine the contents of the UNITs within _list_ into a single UNIT.
The new UNIT is placed in _variable_.   This command is similar to
the sequence command except it does not link the ATOMs of the UNITs
together.

"
  (apply 'cando:combine list))

(defun leap-list-variables ()
  "    list

List the variables currently defined."
  (let* ((names (leap.core:all-variable-names))
         (sorted-names (sort names #'string<))
         (longest-name (loop for name in sorted-names
                             maximize (length name))))
    (loop for name in sorted-names
          for count from 1 by (+ 2 longest-name)
          do (format t "~va" (+ 2 longest-name) name)
          do (when (> count 60)
               (terpri)
               (setf count 0)))
    nil))



(defun leap.assign-atom-types (var-name)
  "    assignAtomTypes variable
     aggregate/molecule        _variable_

Assign atom types to all atoms in the object. The way that atom types
are assigned depends on the force-field-name of each molecule in the object.
The Amber protein force field assigns types based on residue/atom names.
The GAFFx force fields use GAFF type rules. Force fields like SMIRNOFF assign
Cando internally defined types that correspond to unique SMIRKS patterns."
  (let ((object (leap.core:lookup-variable var-name)))
    (unless (or (typep object 'chem:molecule)
                (typep object 'chem:aggregate))
      (error "The argument must be an aggregate or molecule - you passed ~a of class ~a" object (class-name (class-of object))))
    (leap:assign-atom-types object)
    (let ((total-atoms 0)
          (total-atoms-with-types 0))
      (cando:do-atoms (atm object)
        (incf total-atoms)
        (when (chem:get-type atm) (incf total-atoms-with-types)))
      (unless (= total-atoms total-atoms-with-types)
        (warn "There were ~a atoms that could not be assigned types" (- total-atoms total-atoms-with-types))
        (format t "Missing types for atoms in ~a...~%" object)
        (cond
          ((typep object 'chem:aggregate)
           (format t "Aggregate~%")
           (cando:do-molecules (mol object)
             (cando:do-residues (res mol)
               (cando:do-atoms (atm res)
                 (unless (chem:get-type atm)
                   (format t "~a ~a ~a~%" mol res atm))))))
          ((typep object 'chem:molecule)
           (cando:do-residues (res object)
             (cando:do-atoms (atm res)
               (unless (chem:get-type atm)
                 (format t "~a ~a~%" res atm)))))))
      (format t "Assigned ~a types of a total of ~a atoms.~%" total-atoms-with-types total-atoms))))

(defun leap.load-smirnoff-params (filename)
  (leap:load-smirnoff-params filename))

(defun leap-set-force-field (matter force-field-name)
  (leap:set-force-field matter force-field-name))

(defun object (name)
  "Return the object with name."
  (leap.core:lookup-variable name))

(defun leap.add-pdb-res-map (entry)
  "    addPdbAtomMap list
      LIST                         _list_

The atom Name Map is used to try to map atom names read from PDB files
to atoms within residue UNITs when the atom name in the PDB file does
not match an atom in the residue.  This enables PDB files to be read
in without extensive editing of atom names.  The LIST is a LIST of LISTs:
      { {sublist} {sublist} ... }
where each sublist is of the form
      { \"OddAtomName\" \"LibAtomName\" }
Many `odd' atom names can map to one `standard' atom name, but any single
odd atom name maps only to the last standard atom name it was mapped to.
"
  (valid-arguments entry 1)
  (let ((map (second entry)))
    (leap.pdb:add-pdb-res-map map)))

(defun create-atom (name &optional type (charge 0.0))
  "    variable = createAtom name type charge

      ATOM                         _variable_
      STRING                       _name_
      STRING                       _type_
      NUMBER                       _charge_

Return a new ATOM with _name_, _type_, and _charge_.
"
  (let* ((name-sym (intern (string name) :keyword))
         (element (chem:element-from-atom-name-string (string name-sym)))
         (atom (chem:make-atom name-sym element)))
    (chem:set-type atom type)
    (chem:set-charge atom charge)
    atom))


(defun leap-add-atom-types (info)
  (check-type info list)
  (loop for one in info
        for atom-type-string = (first one)
        for element-string = (second one)
        for hybridization-string = (third one)
        do (let ((atom-type (intern (string atom-type-string) :keyword))
                 (element (intern (string element-string) :keyword))
                 (hybridization (intern (string-upcase (string hybridization-string)) :keyword)))
             (leap.atom-types:add-one-atom-type atom-type element hybridization))))


(defun load-off (filename)
"    loadOff filename
      STRING                       _filename_

This command loads the Object File Format library within the file named
_filename_.  All UNITs and PARMSETs within the library will be loaded.
The objects are loaded into LEaP under the variable names the objects
had when they were saved.  Variables already in existence that have the
same names as the objects being loaded will be overwritten.  PARMSETs
loaded using this command are included in LEaP's library of PARMSETs
that is searched whenever parameters are required.
"
  (leap.off:load-off filename))

(defun leap-list-force-fields ()
  (let ((ffs (leap:list-force-fields)))
    (loop for ff in ffs
          do (format t "~s ~a~%" (car ff) (cdr ff)))))

(defun leap-source (entry)
"    source filename
      STRING                       _filename_

This command executes LEaP commands within a text file.  To display the
commands as they are read, see the verbosity command.  The text within
the source file must be formatted exactly like the text the user types
into LEaP. This command temporarily adds the path of _filename_ to the
leap path list so that files local to the script can be loaded.
"
  (let* ((filename (leap.core:ensure-path entry)))
    (leap:source filename)))


(defun leap.save-amber-parm (aggregate-name topology-file-name &optional (crd-pathname nil crd-pathname-p))
  "    saveAmberParm     unit topologyfilename coordinatefilename

      UNIT                         _unit_
      STRING                       _topologyfilename_
      STRING                       _coordinatefilename_

Save the AMBER topology and coordinate files for the UNIT into the
files named _topologyfilename_ and _coordinatefilename_ respectively.
This command will cause LEaP to search its list of PARMSETs for
parameters defining all of the interactions between the ATOMs within
the UNIT.  This command produces a topology file and a coordinate file
which are identical in format to those produced by the AMBER program
PARM, and which can be read into AMBER and SPASMS for energy minimization,
dynamics, or nmode calculations.
See also: saveAmberParmPol, saveAmberParmPert, and saveAmberParmPolPert
for including atomic polarizabilities and preparing free energy
perturbation calculations and saveAmberParmNetcdf for saving in a
binary format.

"
  (let ((aggregate (leap.core:lookup-variable aggregate-name))
        (crd-pathname (if crd-pathname-p
                          crd-pathname
                          (namestring (make-pathname :type "crd" :defaults topology-file-name)))))
    (leap.topology:save-amber-parm-format aggregate topology-file-name crd-pathname)))

(defun ensure-string (obj)
  (cond
    ((stringp obj) obj)
    ((symbolp obj)
     (let ((val (leap.core:lookup-variable obj nil :not-found)))
       (if (eq val :not-found)
           (string val)
           val)))))

(defun valid-arguments (entry num)
  (unless (= (length entry) (1+ num))
    (error "Bad arguments for ~a" (car entry))))



(defun process-iso-closeness (iso-closeness)
  (let ((iso nil)
        (closeness 1.0))
    (flet ((error-iso-closeness ()
             (error "You must provide either :iso or a number for closeness - you provided ~s" iso-closeness)))
      (cond
        ((= (length iso-closeness) 0))
        ((= (length iso-closeness) 1)
         (cond
           ((numberp (first iso-closeness))
            (setf closeness (first iso-closeness)))
           ((eq (first iso-closeness) :iso)
            (setf iso t))
           ((string-equal (first iso-closeness) "iso")
            (setf iso t))
           (t (error-iso-closeness))))
        ((and (= (length iso-closeness) 2)
              (or (find :iso iso-closeness)
                  (find "iso" iso-closeness :test #'string-equal))
              (find-if #'numberp iso-closeness))
         (setf iso t
               closeness (find-if #'numberp iso-closeness)))
        (t (error-iso-closeness))))
    (list :isotropic iso :closeness closeness)))

(defun leap.solvate-box (solute-name solvent-name buffer &rest iso-closeness)
  "    solvateBox solute solvent buffer [ \"iso\" ] [ closeness ]

      UNIT                         _solute_
      UNIT                         _solvent_
      object                       _buffer_
      NUMBER                       _closeness_

The solvateBox command creates a solvent box around the _solute_ UNIT.
The _solute_ UNIT is modified by the addition of _solvent_ RESIDUEs.

The user may want to first align long solutes that are not expected
to tumble using alignAxes, in order to minimize box volume.

The normal choice for a TIP3 _solvent_ UNIT is WATBOX216. Note that
constant pressure equilibration is required to bring the artificial box
to reasonable density, since Van der Waals voids remain due to the
impossibility of natural packing of solvent around the solute and at
the edges of the box.

The solvent box UNIT is copied and repeated in all three spatial directions
to create a box containing the entire solute and a buffer zone defined
by the _buffer_ argument. The _buffer_ argument defines the distance,
in angstroms, between the wall of the box and the closest ATOM in the
solute.

If the buffer argument is a single NUMBER, then the buffer distance is
the same for the x, y, and z directions, unless the \"iso\" option is used
to make the box isometric, with the shortest box clearance = buffer. If
\"iso\" is used, the solute is rotated to orient the principal axes,
otherwise it is just centered on the origin.

If the buffer argument is a LIST of three NUMBERS, then the NUMBERs are
applied to the x, y, and z axes respectively. As the larger box is created
and superimposed on the solute, solvent molecules overlapping the solute
are removed.

The optional _closeness_ parameter can be used to control the extent to
which _solvent_ ATOMs overlap _solute_ ATOMs.  The default value of
the _closeness_ argument is 1.0, which allows no overlap.  Smaller
values allow solvent ATOMs to overlap _solute_ ATOMs by (1 - closeness) *
R*ij, where R*ij is the sum of the Van der Waals radii of solute and
solvent atoms.  Values greater than 1 enforce a minimum gap between
solvent and solute of (closeness - 1) * R*ij.

This command modifies the _solute_ UNIT in several ways.  First, the
coordinates of the ATOMs are modified to move the center of a box
enclosing the Van der Waals radii of the atoms to the origin.  Secondly,
the UNIT is modified by the addition of _solvent_ RESIDUEs copied from
the _solvent_ UNIT. Finally, the box parameter of the new system (still
named for the _solute_) is modified to reflect the fact that a periodic,
rectilinear solvent box has been created around it.
"
  (let ((solute (leap.core:lookup-variable solute-name))
        (solvent (leap.core:lookup-variable solvent-name))
        (iso-closeness-args (process-iso-closeness iso-closeness)))
     (apply #'leap:solvate-box solute solvent buffer iso-closeness-args)))

(defun leap.solvate-oct (solute solvent farness &rest iso-closeness)
  "    solvateOct solute solvent buffer [ \"iso\" ] [ closeness ]

      UNIT                         _solute_
      UNIT                         _solvent_
      object                       _buffer_
      NUMBER                       _closeness_

The solvateOct command is the same as solvateBox, except the corners
of the box are sliced off, resulting in a truncated octahedron, which
typically gives a more uniform distribution of solvent around the
solute.

In solvateOct, when a LIST is given for the buffer argument, four
numbers are given instead of three, where the fourth is the diagonal
clearance. If 0.0 is given as the fourth number, the diagonal clearance
resulting from the application of the x,y,z clearances is reported. If
a non-0 value is given, this may require scaling up the other clearances,
which is also reported. Similarly, if a single number is given, any
scaleup of the x,y,z buffer to accommodate the diagonal clip is reported.

If the \"iso\" option is used, the isometric truncated octahedron is
rotated to an orientation used by the PME code, and the box and angle
dimensions output by the saveAmberParm* commands are adjusted for PME
code imaging.
"
  (let ((buffer (list farness farness farness))
        (iso-closeness-args (process-iso-closeness iso-closeness)))
    (apply #'leap:solvate-shell solute solvent buffer :farness farness :oct t iso-closeness-args)))

(defun leap.solvate-shell (solute solvent farness &rest iso-closeness)
  "    solvateShell solute solvent thickness [ closeness ]

      UNIT                         _solute_
      UNIT                         _solvent_
      NUMBER                       _thickness_
      NUMBER                       _closeness_

The solvateShell command creates a solvent shell around the _solute_ UNIT.
The _solute_ UNIT is modified by the addition of _solvent_ RESIDUEs.

The normal choice for a TIP3 _solvent_ UNIT is WATBOX216. The _solvent_ box
is repeated in all three spatial directions and _solvent_ RESIDUEs selected
to create a solvent shell with a radius of _thickness_ Angstroms around the
_solute_.

The _thickness_ argument defines the maximum distance a _solvent_ ATOM may
be from the closest _solute_ ATOM.

The optional _closeness_ parameter can be used to control overlap of _solvent_
with _solute_ ATOMs.   The default value of the _closeness_ argument is
1.0, which allows contact but no overlap.  Please see the solvateBox
command for more details on the _closeness_ parameter.
"
  (let ((buffer (list farness farness farness))
        (iso-closeness-args (process-iso-closeness iso-closeness)))
    (apply #'leap:solvate-shell solute solvent buffer iso-closeness-args)))

(defun leap.solvate-cap (solute-name solvent-name position radius &optional closeness)
  "    solvateCap solute solvent position radius [ closeness ]

      UNIT                         _solute_
      UNIT                         _solvent_
      object                       _position_
      NUMBER                       _radius_
      NUMBER                       _closeness_

The solvateCap command creates a solvent cap around the _solute_ UNIT
or a part thereof.  The _solute_ UNIT is modified by the addition of
_solvent_ RESIDUEs.
The normal choice for a TIP3 _solvent_ UNIT is WATBOX216. The _solvent_
box is repeated in all three spatial directions and _solvent_ RESIDUEs
selected to create a solvent sphere with a radius of _radius_ Angstroms.

The _position_ argument defines where the center of the solvent cap
is to be placed.   If _position_ is a UNIT, RESIDUE, ATOM, or a LIST
of UNITs, RESIDUEs, or ATOMs, then the geometric center of the ATOMs
within the object will be used as the center of the solvent cap sphere.
If _position_ is a LIST containing three NUMBERs then the _position_
argument will be treated as a vector that defines the position of
the solvent cap sphere center.

The optional _closeness_ parameter can be used to control the extent to
which _solvent_ ATOMs overlap _solute_ ATOMs.  The default value of
the _closeness_ argument is 1.0, which allows no overlap.  Smaller
values allow solvent ATOMs to overlap _solute_ ATOMs by (1 - closeness) *
R*ij, where R*ij is the sum of the Van der Waals radii of solute and
solvent atoms.  Values greater than 1 enforce a minimum gap between
solvent and solute of (closeness - 1) * R*ij.

This command modifies the _solute_ UNIT in several ways.  First,
the UNIT is modified by the addition of _solvent_ RESIDUEs copied
from the _solvent_ UNIT. Secondly, the cap parameter of the UNIT
_solute_ is modified to reflect the fact that a solvent cap has been
created around the solute.
"
  (let ((solute (leap.core:lookup-variable solute-name))
        (solvent (leap.core:lookup-variable solvent-name))
        (list-closeness (list :closeness closeness)))
    (apply #'leap:solvate-cap solute solvent position radius list-closeness)))


(defun leap.dir (&optional (path *default-pathname-defaults*))
  (let* ((name (or (pathname-name path) "*"))
         (dirpath (make-pathname :name name :defaults path)))
    (format t "~a~%" (directory dirpath))))

(defun load-chem-draw (filename)
  (handler-bind ((warning #'muffle-warning))
    (let ((pathname (leap.core:ensure-path filename)))
      (with-open-file (fin (open pathname :direction :input))
        (chem:make-chem-draw fin :add-hydrogens nil)))))

(defun leap-jostle-atoms (matter)
  "     jostleAtoms matter

Add a random displacement to every atom in the matter."
  (let ((matter (leap.core:lookup-variable matter)))
    (cando:jostle matter)))

(defun leap-minimize (matter
                      &rest args
                      &key (restraints-on t)
                        (max-sd-steps 1000)
                        (max-cg-steps 50000)
                        (max-tn-steps 0)
                        (sd-tolerance 5000.0)
                        (cg-tolerance 0.5)
                        (tn-tolerance 0.00001)
                        (use-excluded-atoms t)
                        (assign-types t)
                        (save-trajectory nil)
                        )
  "     minimize matter &key (restraints-on t)
                          (max-sd-steps 1000)
                          (max-cg-steps 50000)
                          (max-tn-steps 0)
                          (sd-tolerance 5000.0)
                          (cg-tolerance 0.5)
                          (tn-tolerance 0.00001)
                          (use-excluded-atoms t)
                          (assign-types t)
                          (save-trajectory nil))

      AGGREGATE/MOLECULE     _matter_
Minimize the energy of the given matter. If _assign-types_ is passed with T
then the atom types are assigned using the rules defined by the force-field
of each molecule. "
  (let ((matter (leap.core:lookup-variable matter)))
    (format t "args -> ~s~%" args)
    (apply 'energy:minimize matter args)))

(defun leap.load-chem-draw (filename)
  (load-chem-draw filename))

(defun leap.load-chem-draw-aggregate (filename)
  "   loadChemDrawAggregate filename

Load an aggregate from a ChemDraw cdxml file and return it."
  (handler-bind ((warning #'muffle-warning))
    (let ((pathname (leap.core:ensure-path filename)))
      (with-open-file (fin (open pathname :direction :input))
        (cando:load-chem-draw-aggregate pathname)))))

(defun leap.compile-smarts (smarts-string &optional tests)
  (unless (stringp smarts-string)
    (error "You must provide a smarts string"))
  (chem:compile-smarts smarts-string :tests tests))

(defun leap-quit ()
  (ext:quit))

(defun leap-cando ()
  "
This command brings up a cando command line.
Leap can be reentered by evaluating (leap).
"
  (core:top-level nil nil))


(defun leap-start-swank (&optional (port 4005))
  " startSwank &optional (port 4005)

Start a swank server for debugging cando."
  (funcall (find-symbol "START-SWANK" :cando-user) port))

(defun leap-help (&optional arg)
  "    help [string]

      STRING                       _string_

This command prints a description of the command in _string_.  If
the STRING is not given then a list of legal STRINGs is provided.
"
  (if arg
      (progn
        (let ((found (assoc (string arg) leap.parser:*function-names/alist* :test #'string-equal)))
          (if found
              (let ((cmd (cdr found)))
                (format t "~%~a~%" (documentation cmd 'function)))
              (format t "No help available for ~s package: ~s~%" arg (symbol-package arg)))))
      (progn
        (format t "Help is available on the following subjects: ~%")
        (let* ((commands (copy-seq leap.parser:*function-names/alist*))
               (sorted-commands (sort commands #'string< :key #'car)))
          (let ((longest-cmd (loop for cmd in sorted-commands maximize (length (car cmd)))))
            (loop for cmd in sorted-commands
                  for col from 1
                  do (format t "~v,a" (+ 3 longest-cmd) (car cmd))
                  when (and (> col 0) (= (rem col 4) 0))
                    do (terpri)
                  finally (when (/= (rem col 4) 0)
                            (terpri))))))))

(defun leap-add (object-a-name  object-b-name)
 "    add a b

      UNIT/RESIDUE/ATOM          _a_
      UNIT/RESIDUE/ATOM          _b_

Add the object _b_ to the object _a_.   This command is used to place
ATOMs within RESIDUEs, and RESIDUEs within UNITs.
"
  (let ((object-a (leap.core:lookup-variable object-a-name))
        (object-b (leap.core:lookup-variable object-b-name)))
    (chem:add-matter object-a object-b))
  :no-output)

(defun leap-bond (atom1-name atom2-name &rest order)
"    bond atom1 atom2 [ order ]

      ATOM                         _atom1_
      ATOM                         _atom2_
      STRING                       _order_

Create a bond between _atom1_ and _atom2_.  Both of these ATOMs must
be contained by the same UNIT. By default, the bond will be a single
bond.  By specifying \"S\", \"D\", \"T\", or \"A\" as the optional argument
_order_ the user can specify a single, double, triple, or aromatic
bond.
"
  (let ((atom1 (leap.core:lookup-variable atom1-name))
        (atom2 (leap.core:lookup-variable atom2-name)))
    (case order
      (:s (chem:bond-to atom1 atom2 :single-bond nil))
      (:d (chem:bond-to atom1 atom2 :double-bond nil))
      (:t (chem:bond-to atom1 atom2 :triple-bond nil))
      (:a (chem:bond-to atom1 atom2 :aromatic-bond nil))
      (otherwise (chem:bond-to atom1 atom2 :single-bond nil))))
  :no-output)

(defun leap-add-h (object-name)
"    addH obj

      UNIT                       _obj_

Add missing hydrogens and build external coordinates for _obj_.
"
  (let ((object (leap.core:lookup-variable object-name)))
    (cando:build-unbuilt-hydrogens object)))

(defun leap-create-residue (name)
"    variable = createResidue name
      RESIDUE                      _variable_
      STRING                       _name_

Return a new and empty RESIDUE with the name _name_.
"
  (let ((residue (chem:make-residue name)))
    residue))

(defun leap-create-unit (name)
"    variable = createUnit name
      AGGREGATE                    _variable_
      STRING                       _name_

Return a new and empty AGGREGATE with the name _name_.
"
  (let ((aggregate (chem:make-aggregate name)))
    aggregate))

(defun leap-create-aggregate (name)
"    variable = createAggregate name
      AGGREGATE                    _variable_
      STRING                       _name_

Return a new and empty AGGREGATE with the name _name_.
"
  (let ((aggregate (chem:make-aggregate name)))
    aggregate))

(defun leap-create-molecule (name)
"    variable = createMolecule name
      MOLECULE                    _variable_
      STRING                       _name_

Return a new and empty MOLECULE with the name _name_.
"
  (let ((molecule (chem:make-molecule name)))
    molecule))

(defun leap-translate (aggregate-name direction)
"      translate atoms direction

      UNIT/RESIDUE/ATOM            _atoms_
      LIST                         _direction_

Translate all of the ATOMs within _atoms_ by the vector defined by
the three NUMBERs in the LIST _ direction_.
"
  (check-type direction list)
  (let ((aggregate (leap.core:lookup-variable aggregate-name)))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (let* ((pos (chem:get-position atom))
              (x-position (geom:vx pos))
              (y-position (geom:vy pos))
              (z-position (geom:vz pos)))
         (geom:set-all3 pos
                        (+ x-position (first direction))
                        (+ y-position (second direction))
                        (+ z-position (third direction)))
         (chem:set-position atom pos)))
   aggregate)))

(defun leap-copy (object-name)
"    newvariable = copy variable

      object                       _newvariable_
      object                       _variable_

Create an exact duplicate of the object _variable_. Changing the object
_variable_ will not affect the object _newvariable_.
This is in contrast to the situation created by \"newvariable = variable\"
in which both names reference the same object.
"
  (let* ((object (leap.core:lookup-variable object-name))
         (new-object (chem:matter-copy object)))
    new-object))

(defun leap-center (container-name)
"    center container
      UNIT/RESIDUE/ATOM          _container_

Display the coordinates of the geometric center of the ATOMs within
_container_.
"
  (let* ((container (leap.core:lookup-variable container-name))
         (position (chem:geometric-center container)))
    (format t "The center is at: ~,2f, ~,2f, ~,2f~%" (geom:vx position) (geom:vy position) (geom:vz position))))

(defun leap-measure-geom (atom1 atom2 &optional atom3 atom4)
"    measureGeom atom1 atom2 [ atom3 [ atom4 ] ]
      ATOM                _atom1_ _atom2_ _atom3_ _atom4_

Measure the distance, angle, or torsion between two, three, or four ATOMs,
respectively.
"
  (let ((atom1-position (chem:get-position atom1))
        (atom2-position (chem:get-position atom2)))
    (if atom3
        (let ((atom3-position (chem:get-position atom3)))
          (if atom4
              (let ((atom4-position (chem:get-position atom4)))
                (format t "Torsion angle: ~,2f degrees~%"
                        (* 57.2958
                           (geom:calculate-dihedral atom1-position atom2-position atom3-position atom4-position))))
              (format t  "Angle: ~,2f degrees~%"
                      (* 57.2958 (geom:calculate-angle atom1-position atom2-position atom3-position)))))
        (format t "Distance: ~,2f angstroms~%"
                (geom:calculate-distance
                 atom1-position atom2-position)))))

(defun leap-add-path (path)
"    addPath path

      STRING                       _path_

Add the directory in _path_ to the list of directories that are searched
for files specified by other commands.
"
  (setf *default-pathname-defaults* path))

(defun leap-align-axes (unit-name)
"     alignAxes unit

      UNIT                         _unit_

Translate the geometric center of _unit_ to the origin and align
the principle axes of the ATOMs within _unit_ along the coordinate
axes.  This is done by calculating the moment of inertia of the UNIT
using unit mass for each ATOM, and then diagonalizing the resulting
matrix and aligning the eigenvectors along the coordinate axes.
This command modifies the coordinates of the UNIT. It may be
especially useful for preparing long solutes such as nucleic acids
for solvation.
"
  (let ((unit (leap.core:lookup-variable unit-name)))
  (leap.align-axes:tool-orient-principle-axis-along-coordinate-axis unit)))

(defun leap-transform (atoms-name matrix)
"    transform atoms matrix

      CONTAINER/LIST               _atoms_
      LIST                         _matrix_

Translate all of the ATOMs within _atoms_ by the (3X3) or (4X4) matrix
defined by the 9 or 16 NUMBERs in the LIST of LISTs _matrix_.
    { { r11  r12  r13  -tx }
      { r21  r22  r23  -ty }
      { r31  r32  r33  -tz }
      { 0    0    0     1  } }
The diagonal upper left elements, rII can be used for symmetry operations,
e.g. a reflection in the XY plane can be produced with r11=1, r22=1, r33=-1
where the other rIJ elements are 0. The -t column is used to specify
translations along the appropriate axes (0 for no translation).
"
  (let ((atoms (leap.core:lookup-variable atoms-name))
        (transform (geom:make-matrix-identity))
        (ix 0)
        (iy 0))
    (check-type matrix list)
    (dolist (elements-y matrix)
      (format t "elements-y ~s~%" elements-y)
      (check-type elements-y list)
      (if (< iy 4)
          (dolist (elements-x elements-y)
            (if (< ix 4)
                (progn
                  (geom:at-row-col-put transform ix iy elements-x)
                  (incf ix)))))
      (setf ix 0)
      (incf iy))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (let* ((atom-position (chem:get-position atom))
              (x-position (geom:vx atom-position))
              (y-position (geom:vy atom-position))
              (z-position (geom:vz atom-position)))
         (geom:set-all3 atom-position
                        (+ (* (geom:at-row-col-get transform 0 0) (geom:vx atom-position))
                           (* (geom:at-row-col-get transform 1 0) (geom:vy atom-position))
                           (* (geom:at-row-col-get transform 2 0) (geom:vz atom-position)))
                        (+ (* (geom:at-row-col-get transform 0 1) (geom:vx atom-position))
                           (* (geom:at-row-col-get transform 1 1) (geom:vy atom-position))
                           (* (geom:at-row-col-get transform 2 1) (geom:vz atom-position)))
                        (+ (* (geom:at-row-col-get transform 0 2) (geom:vx atom-position))
                           (* (geom:at-row-col-get transform 1 2) (geom:vy atom-position))
                           (* (geom:at-row-col-get transform 2 2) (geom:vz atom-position))))
         #+(or)(geom:set-all3 atom-position
                              (geom:vx (geom:m*v transform atom-position))
                              (geom:vy (geom:m*v transform atom-position))
                              (geom:vz (geom:m*v transform atom-position)))
         (chem:set-position atom atom-position)))
     atoms)))

(defun leap-charge (container-name)
"    charge container

      UNIT/RESIDUE/ATOM          _container_

This command calculates the total charge of the ATOMs within _container_.
The unperturbed and perturbed total charge are displayed.
"
  (let ((container (leap.core:lookup-variable container-name))
        (total-charge 0.0)
        (pert-charge 0.0))
    ;; We don't handle pert-charge - and it would be a delta
    (cando:do-atoms (atm container)
      (incf total-charge (chem:get-charge atm)))
    (format t "Total unperturbed charge:    ~,6f~%" total-charge)
    (format t "Total perturbed charge:      ~,6f~%" (+ total-charge pert-charge))))

(defun leap-add-ions (mol-name ion1-name ion1-number &optional ion2-name ion2-number)
"
    addIons unit ion1 #ion1 [ion2 #ion2]
      UNIT   U                   _unit_
      UNIT                      _ion1_
      NUMBER                    _#ion1_
      UNIT                      _ion2_
      NUMBER                    _#ion2_
Adds counterions in a shell around _unit_ using a Coulombic potential
on a grid. If _#ion1_ is 0, the _unit_ is neutralized (_ion1_ must be
opposite in charge to _unit_, and _ion2_ cannot be specified). Otherwise,
the specified numbers of _ion1_ [_ion2_] are added [in alternating order].
If solvent is present, it is ignored in the charge and steric calculations,
and if an ion has a steric conflict with a solvent molecule, the ion is
moved to the center of said molecule, and the latter is deleted. (To
avoid this behavior, either solvate _after_ addIons, or use addIons2.)
Ions must be monoatomic. Note that the one-at-a-time procedure is not
guaranteed to globally minimize the electrostatic energy. When neutralizing
regular-backbone nucleic acids, the first cations will generally be added
between phosphates, leaving the final two ions to be placed somewhere around
the middle of the molecule.
The default grid resolution is 1 Angstrom, extending from an inner radius
of (max ion size + max solute atom size) to an outer radius 4 Angstroms
beyond. A distance-dependent dielectric is used for speed.
"
  (let ((mol (leap.core:lookup-variable mol-name))
        (ion1 (intern (string ion1-name) :keyword)))
        (if (and ion2-name ion2-number)
            (let ((ion2 (intern (string ion2-name) :keyword)))
              (leap.add-ions:add-ions mol ion1 ion1-number ion2 ion2-number))
            (leap.add-ions:add-ions mol ion1 ion1-number))))

(defun leap-add-ions-2 (mol-name ion1-name ion1-number &optional ion2-name ion2-number)
"    addIons2 unit ion1 #ion1 [ion2 #ion2]
        UNIT                    _unit_
        UNIT                    _ion1_
        NUMBER                  _#ion1_
        UNIT                    _ion2_
        NUMBER                  _#ion2_
Same as addIons, except solvent and solute are treated the same.
"
  (let ((mol (leap.core:lookup-variable mol-name))
        (ion1 (intern (string ion1-name) :keyword)))
        (if (and ion2-name ion2-number)
            (let ((ion2 (intern (string ion2-name) :keyword)))
              (leap.add-ions:add-ions mol ion1 ion1-number ion2 ion2-number))
            (leap.add-ions:add-ions mol ion1 ion1-number))))

(defun leap-add-ions-rand (mol-name ion1-name ion1-number &rest ion2-separation)
"    addIonsRand unit ion1 #ion1 [ion2 #ion2] [separation]

        UNIT      _unit_
        UNIT      _ion1_
        NUMBER    _#ion1_
        UNIT      _ion2_
        NUMBER    _#ion2_
        NUMBER    _separation_

Adds counterions in a shell around _unit_ by replacing random solvent
molecules. If _#ion1_ is 0, the _unit_ is neutralized (_ion1_ must be
opposite in charge to _unit_, and _ion2_ cannot be specified). Otherwise,
the specified numbers of _ion1_ [_ion2_] are added [in alternating order].
If _separation_ is specified, ions will be guaranteed to be more than that
distance apart in Angstroms.

Ions must be monoatomic. This procedure is much faster than addIons, as
it does not calculate charges. Solvent must be present. It must be possible
to position the requested number of ions with the given separation in the
solvent.
"
    (let ((mol (leap.core:lookup-variable mol-name))
          (ion1 (intern (string ion1-name) :keyword)))
      (cond
        ((= (length ion2-separation) 0)
         (leap.add-ions:add-ions-rand mol ion1 ion1-number))
        ((= (length ion2-separation) 1)
         (leap.add-ions:add-ions-rand mol ion1 ion1-number :separation (first ion2-separation)))
        ((= (length ion2-separation) 2)
         (let ((ion2 (intern (string (first ion2-separation)) :keyword))
               (ion2-number (second ion2-separation)))
           (leap.add-ions:add-ions-rand mol ion1 ion1-number :ion2 ion2 :ion2-number ion2-number)))
        ((= (length ion2-separation) 3)
         (let ((ion2 (intern (string (first ion2-separation)) :keyword))
               (ion2-number (second ion2-separation))
               (separation (third ion2-separation)))
           (leap.add-ions:add-ions-rand mol ion1 ion1-number :ion2 ion2 :ion2-number ion2-number :separation separation))))))

(defun leap-save-mol2 (aggregate-name path-name option)
"    saveMol2 unit filename option

      UNIT                         _unit_
      STRING                       _filename_
      NUMBER                       _option_

Write UNIT to the file _filename_ as a Mol2 format file.
option = 0 for Default atom types
option = 1 for AMBER atom types

More information:
http://q4md-forcefieldtools.org/Tutorial/leap.php
"
  (let ((aggregate (leap.core:lookup-variable aggregate-name)))
    (cond
      ((= option 0)
       (cando:save-mol2 aggregate path-name :use-sybyl-types t))
      ((= option 1)
       (cando:save-mol2 aggregate path-name :use-sybyl-types nil))
      (t (error "Option must be 0 or 1")))))

(defun leap-save-pdb (aggregate-name path-name)
"    savePdb unit filename

      UNIT                         _unit_
      STRING                       _filename_

Write UNIT to the file _filename_ as a PDB format file.
"
  (let ((aggregate (leap.core:lookup-variable aggregate-name)))
    (format t "Writing pdb file: ~a~%" path-name)
    (chem:save-pdb aggregate path-name)))

(defun leap-delete-bond (atom1-name atom2-name)
"    deleteBond atom1 atom2

      ATOM                         _atom1_
      ATOM                         _atom2_

Remove the bond between the ATOMs _atom1_ and _atom2_.   If no bond
exists, an error will be displayed.
"
  (let ((atom1 (leap.core:lookup-variable atom1-name))
        (atom2 (leap.core:lookup-variable atom2-name)))
    (if (chem:is-bonded-to atom1 atom2)
        (chem:remove-bond-to atom1 atom2)
        (error "Atoms are not bonded."))))

(defun leap-select (matter-name)
"    select obj

      UNIT/RESIDUE/ATOM          _obj_

Sets the SELECT flag on all ATOMs within _obj_.  See the deSelect command.
"
  (let ((matter (leap.core:lookup-variable matter-name)))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (chem:set-property atom :select t))
     matter)))

(defun leap-de-select (matter-name)
"    deSelect obj

      UNIT/RESIDUE/ATOM          _obj_

Clears the SELECT flag on all ATOMs within _obj_.  See the select command.
"
  (let ((matter (leap.core:lookup-variable matter-name)))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (if (chem:has-property atom :select)
           (chem:clear-property atom :select)))
     matter)))

(defun leap-bond-by-distance (matter-name &optional (max-bond 2.0))
"    bondByDistance container [ maxBond ]

      UNIT/RESIDUE/ATOM         _container_
      NUMBER                    _maxBond_

Create single bonds between all ATOMs in _container_ that are within
_maxBond_ angstroms of each other.   If _maxBond_ is not specified,
a default distance of 2 angstroms used.
"
  (let ((matter (leap.core:lookup-variable matter-name))
        (atom-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (count 0))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (vector-push-extend atom atom-vec)
       (incf count))
     matter)
    (loop for i from 0 below count
          do (loop for j from (+ i 1) below count
                   for atom1 = (aref atom-vec i)
                   for atom2 = (aref atom-vec j)
                   do (if (and (not (chem:is-bonded-to atom1 atom2))
                               (> (* max-bond max-bond) (chem:distance-squared-between-two-atoms atom1 atom2)))
                          (chem:bond-to atom1 atom2 :single-bond))))))

(defun distance-search (matter closeness absolute-distance)
  (let ((atom-vec (make-array 256 :fill-pointer 0 :adjustable t))
        (distance 0.0)
        (dx 0.0)
        (dy 0.0)
        (dz 0.0)
        (radius 0.0))
    (chem:map-atoms
     'nil
     (lambda (atom)
       (vector-push-extend atom atom-vec))
     matter)
    (loop for i from 0 below (length atom-vec)
          for atom1 = (aref atom-vec i)
          for atom1-pos = (chem:get-position atom1)
          for atom1-rad = (chem:get-vdw-radius atom1)
          for atom1-element = (chem:get-element atom1)
          do  (if (< atom1-rad 0.1)
                  (if (eq atom1-element :H)
                      (setf atom1-rad 1.0)
                      (setf atom1-rad 1.5)))
          do (setf atom1-rad (* atom1-rad closeness))
          do (loop for j from (+ i 1) below (length atom-vec)
                   for atom2 = (aref atom-vec j)
                   for atom2-pos = (chem:get-position atom2)
                   for atom2-rad = (chem:get-vdw-radius atom2)
                   do (if (chem:is-bonded-to atom1 atom2)
                          (progn
                            (if absolute-distance
                                (setf radius closeness)
                                (setf radius (* 0.5 (+ atom1-rad atom2-rad))))
                            (if (or (> (geom:vx atom1-pos) (+ (geom:vx atom2-pos) radius))
                                    (> (geom:vy atom1-pos) (+ (geom:vy atom2-pos) radius))
                                    (> (geom:vz atom1-pos) (+ (geom:vz atom2-pos) radius))
                                    (< (geom:vx atom1-pos) (- (geom:vx atom2-pos) radius))
                                    (< (geom:vy atom1-pos) (- (geom:vy atom2-pos) radius))
                                    (< (geom:vz atom1-pos) (- (geom:vz atom2-pos) radius)))
                                (progn
                                  (setf dx (- (geom:vx atom1-pos) (geom:vx atom2-pos))
                                        dy (- (geom:vy atom1-pos) (geom:vy atom2-pos))
                                        dz (- (geom:vz atom1-pos) (geom:vz atom2-pos)))
                                  (setf radius (* radius radius))
                                  (setf distance (+ (* dx dx) (* dy dy) (* dz dz)))
                                  (if (< distance radius)
                                      (warn "Close contact of ~a angstroms between ~a and ~a~%"
                                             (sqrt distance) (chem:get-name (aref atom-vec i))
                                             (chem:get-name (aref atom-vec j)))))))))
                   )))



(defun leap-check (matter-name)
"    check unit [ parmset ]
      UNIT                      _unit_
      PARMSET/STRING            _parmset_

This command can be used to check the UNIT for internal inconsistencies
that could cause problems when performing calculations.  This is
a very useful command that should be used before a UNIT is saved
with saveAmberParm or its variations.
With the optional parmset, all missing parameters are placed in the
PARMSET to allow for easy editing of those parameters.  If a string is
passed, a PARMSET will be created with that name.
Currently it checks for the following possible problems:

   - Long bonds.

   - Short bonds.

   - Non-integral total charge of the UNIT.

   - Missing types.

   - Close contacts between non-bonded ATOMs.  A close contact is
     less than 1.5 angstroms.
"
  (let ((matter (leap.core:lookup-variable matter-name)))
    (distance-search matter 1.5 t)
    (format t "Checking parameters for unit ~a~%" (chem:get-name matter))
    (let* ((energy-function (chem:make-energy-function :matter matter
                                                       :use-excluded-atoms t
                                                       :assign-types t))
           (atom-table (chem:atom-table energy-function))
           (energy-stretch (chem:get-stretch-component energy-function))
           (energy-angle (chem:get-angle-component energy-function))
           (energy-dihedral (chem:get-dihedral-component energy-function)))
      (format t "Checking for stretch parameters. ~%")
      (let* ((stretch-vectors (chem:extract-vectors-as-alist energy-stretch))
             (kb-vector (cdr (assoc :kb stretch-vectors)))
             (r0-vector (cdr (assoc :r0 stretch-vectors)))
             (i1-vector (cdr (assoc :i1 stretch-vectors)))
             (i2-vector (cdr (assoc :i2 stretch-vectors)))
             (atom1-vector (cdr (assoc :atom1 stretch-vectors)))
             (atom2-vector (cdr (assoc :atom2 stretch-vectors))))
        (loop for i from 0 below (length kb-vector)
              for kb = (aref kb-vector i)
              for r0 = (aref r0-vector i)
              for i1 = (aref i1-vector i)
              for i2 = (aref i2-vector i)
              for atom1 = (aref atom1-vector i)
              for atom2 = (aref atom2-vector i)
              do (if (or (not kb)
                         (not r0)
                         (not i1)
                         (not i2))
                     (warn "Could not find stretch parameter: ~a - ~a~%" atom1 atom2))))
      (format t "Checking for angle parameters. ~%")
      (let* ((angle-vectors (chem:extract-vectors-as-alist energy-angle))
             (kt-vector (cdr (assoc :kt angle-vectors)))
             (t0-vector (cdr (assoc :t0 angle-vectors)))
             (i1-vector (cdr (assoc :i1 angle-vectors)))
             (i2-vector (cdr (assoc :i2 angle-vectors)))
             (i3-vector ( cdr (assoc :i3 angle-vectors)))
             (atom1-vector (cdr (assoc :atom1 angle-vectors)))
             (atom2-vector (cdr (assoc :atom2 angle-vectors)))
             (atom3-vector (cdr (assoc :atom3 angle-vectors))))
        (loop for i from 0 below (length kt-vector)
              for kt = (aref kt-vector i)
              for t0 = (aref t0-vector i)
              for i1 = (aref i1-vector i)
              for i2 = (aref i2-vector i)
              for i3 = (aref i3-vector i)
              for atom1 = (aref atom1-vector i)
              for atom2 = (aref atom2-vector i)
              for atom3 = (aref atom3-vector i)
              do (if (or (not kt)
                         (not t0)
                         (not i1)
                         (not i2)
                         (not i3))
                     (warn "Could not find angle parameter: ~a - ~a - ~a~%" atom1 atom2 atom3))))
      (format t "Checking for dihedral parameters. ~%")
      (let* ((energy-dihedral (chem:get-dihedral-component energy-function))
             (dihedral-vectors (chem:extract-vectors-as-alist energy-dihedral))
             (v-vector (cdr (assoc :v dihedral-vectors)))
             (in-vector (cdr (assoc :in dihedral-vectors)))
             (phase-vector (cdr (assoc :phase dihedral-vectors)))
             (i1-vector (cdr (assoc :i1 dihedral-vectors)))
             (i2-vector (cdr (assoc :i2 dihedral-vectors)))
             (i3-vector (cdr (assoc :i3 dihedral-vectors)))
             (i4-vector (cdr (assoc :i4 dihedral-vectors)))
             (proper-vector (cdr (assoc :proper dihedral-vectors)))
             (atom1-vector (cdr (assoc :atom1 dihedral-vectors)))
             (atom2-vector (cdr (assoc :atom2 dihedral-vectors)))
             (atom3-vector (cdr (assoc :atom3 dihedral-vectors)))
             (atom4-vector (cdr (assoc :atom4 dihedral-vectors))))
        (loop for i from 0 below (length v-vector)
              for v = (aref v-vector i)
              for in = (aref in-vector i)
              for phase = (aref phase-vector i)
              for i1 = (aref i1-vector i)
              for i2 = (aref i2-vector i)
              for i3 = (aref i3-vector i)
              for i4 = (aref i4-vector i)
              for atom1 = (aref atom1-vector i)
              for atom2 = (aref atom2-vector i)
              for atom3 = (aref atom3-vector i)
              for atom4 = (aref atom4-vector i)
              do (if (or (not v)
                         (not in)
                         (not phase)
                         (not i1)
                         (not i2)
                         (not i3)
                         (not i4))
                         (warn "Could not find dihedral parameter: ~a - ~a - ~a - ~a~%" atom1 atom2 atom3 atom4)))))))



(defparameter *leap-commands* (list       "add" "addAtomTypes"
    "addH" "addIons" "addIons2" "addIonsRand" "addPath" "addPdbAtomMap" "addPdbResMap" "alias" "alignAxes"
    "bond" "bondByDistance" "center" "charge" "check" "clearPdbAtomMap" "clearPdbResMap" "clearVariables" "combine"
    "copy" "createAtom" "createParmset" "createResidue" "createUnit" "createAggregate" "createMolecule"
    "crossLink" "debugOff" "debugOn" "debugStatus"
    "deleteBond" "deleteOffLibEntry" "deleteRestraint" "desc" "deSelect" "displayPdbAtomMap" "displayPdbResMap" "edit"
    "flip" "groupSelectedAtoms" "help" "impose" "list" "listOff" "loadAmberParams" "loadAmberPrep" "loadMol2"
    "loadMol3" "loadOff" "loadPdb" "loadPdbUsingSeq" "logFile" "matchVariables" "measureGeom" "quit" "relax"
    "remove" "restrainAngle" "restrainBond" "restrainTorsion" "saveAmberParm" "saveAmberParmNetcdf"
    "saveAmberParmPert" "saveAmberParmPol" "saveAmberParmPolPertsaveAmberPrep" "saveMol2" "saveMol3" "saveOff" "saveOffParm"
    "savePdb" "scaleCharges" "select" "sequence" "set" "set_default" "setBox" "showdefault" "solvateBox"
    "solvateCap" "solvateDontClip" "solvateOct" "solvateShell" "source" "transform" "translate" "verbosity" "zMatrix"
    ))

(defun leap-commands-missing ()
  "   leapCommandsMissing

Provide a list of commands that cleap is missing when compared to tleap."
  (loop for cmd-name in *leap-commands*
        unless (member cmd-name leap.parser:*function-names/alist* :key #'car :test #'string=)
          do (format t "~a~%" cmd-name)))


(defun leap-commands-available ()
  "   leapCommandsAvailable

Provide a list of commands that cleap has available to mimic tleap."
  (loop for cmd-name in *leap-commands*
        when (member cmd-name leap.parser:*function-names/alist* :key #'car :test #'string=)
          do (format t "~a~%" cmd-name)))







(defun leap-load-moe (filename)
  (let ((pathname (leap.core:ensure-path filename)))
    (unless pathname
      (error "Could not file file ~a" filename))
    (chem:moe-read-aggregate (namestring pathname))))

(defvar *error-on-bad-alias* t)

(eval-when (:load-toplevel :execute)
  ;; Define a list of names that do not get exported as lisp aliases
  (defvar *ignore-lisp-aliases* '( "list" "set" "remove"
                                  "source" "quit" "desc" "minimize" "show" "combine"))
  (setf leap.parser:*function-names/alist*
        '(("logFile" . leap-log-file)
          ("quickload" . leap-quickload)
          ("combine" . leap-combine)
          ("leapCommandsMissing" . leap-commands-missing)
          ("leapCommandsAvailable" . leap-commands-available)
          ("setupAmberPaths" . leap-setup-amber-paths)
          ("easyGaff" . leap-easy-gaff)
          ("jostleAtoms" . leap-jostle-atoms)
          ("desc" . leap-desc)
          ("listForceFields" . leap-list-force-fields)
          ("loadOff" . load-off)
          ("loadMol2" . cando:load-mol2)
          ("loadMoe" . leap-load-moe )
          ("minimize" . leap-minimize)
          ("dir" . leap.dir )
          ("ls" . leap.dir)
          ("list" . leap-list-variables)
          ("listVariables" . leap-list-variables) ; alternative to "list"
          ("loadPdb" . leap.pdb:load-pdb)
          ("source" . leap-source)
          ("set" . leap.set )
          ("loadChemDraw" . leap.load-chem-draw)
          ("loadChemDrawAggregate" . leap.load-chem-draw-aggregate)
          ("compileSmarts" . leap.compile-smarts)
          ("setForceField" . leap-set-force-field)
          ("loadSmirnoffParams" . leap.load-smirnoff-params)
          ("loadAmberParams" . leap:load-amber-params)
          ("assignAtomTypes" . leap.assign-atom-types)
          ("addPdbResMap" . leap.pdb:add-pdb-res-map)
          ("addPdbAtomMap" . leap.pdb:add-pdb-atom-map)
          ("addAtomTypes" . leap-add-atom-types)
          ("remove" . leap-remove-matter)
          ("removeMatter" . leap-remove-matter)
          ("saveAmberParm" . leap.save-amber-parm)
          ("solvateBox" . leap.solvate-box)
          ("solvateOct" . leap.solvate-oct)
          ("solvateShell" . leap.solvate-shell)
          ("solvateCap" . leap.solvate-cap)
          ("setupDefaultPaths" . leap.setup-default-paths )
          ("addIons" . leap-add-ions)
          ("addIons2" . leap-add-ions-2)
          ("addIonsRand" . leap-add-ions-rand)
          ("setBox" . leap.set-box:set-box)
          ("showPaths" . leap.show-paths)
          ("show" . leap.show )
          ("zMatrix" . leap.z-matrix )
          ("createAtom" . create-atom )
          ("help" . leap-help)
          ("cando" . leap-cando )
          ("startSwank" . leap-start-swank)
          ("quit" . leap-quit)
          ("add" . leap-add)
          ("bond" . leap-bond)
          ("addH" . leap-add-h)
          ("createResidue" . leap-create-residue)
          ("createUnit" . leap-create-unit)
          ("createAggregate" . leap-create-aggregate)
          ("createMolecule" . leap-create-molecule)
          ("translate" . leap-translate)
          ("transform" . leap-transform)
          ("copy" . leap-copy)
          ("center" . leap-center)
          ("measureGeom" . leap-measure-geom)
          ("addPath" . leap-add-path)
          ("alignAxes" . leap-align-axes)
          ("charge" . leap-charge)
          ("saveMol2" . leap-save-mol2)
          ("savePdb" . leap-save-pdb)
          ("deleteBond" . leap-delete-bond)
          ("select" . leap-select)
          ("deSelect" . leap-select)
          ("bondByDistance" . leap-bond-by-distance)
          ("check" . leap-check)
          ))
  ;; register all of the leap command and export them from the leap package
  (loop for command-pair in leap.parser:*function-names/alist*
        for command-name = (car command-pair)
        for command-symbol = (cdr command-pair)
        unless (fboundp command-symbol)
          do (error "~a must be fbound to a function" command-symbol)
        do (progn
             (setf (leap.core:function-lookup command-name leap.core:*leap-env*) command-symbol)
             ;; Bind the leap name to an alias in the LEAP package.
             (when (not (member command-name *ignore-lisp-aliases* :test #'string=))
               (let ((leap-sym (intern (string-upcase command-name) :leap.commands)))
                 (if (not (fboundp leap-sym))
                   (progn
                     (setf (fdefinition leap-sym) (fdefinition command-symbol))
                     (export leap-sym))
                   (when *error-on-bad-alias*
                     (error "Cannot create a lisp alias for function name ~s because it already defines a function ~s~%" leap-sym (fdefinition leap-sym))))))))
  (setf *error-on-bad-alias* nil))

(defun parse-leap-code (code)
  "Parse the CODE and return the result of the parse. This is for debugging the parser."
  (let ((ast (architecture.builder-protocol:with-builder ('list)
               (esrap:parse 'leap.parser:leap code))))
    ast))

(defun process-command-line-options ()
  (let ((includes (core:leap-command-line-includes))
        (scripts (core:leap-command-line-scripts)))
    (loop for include-path in includes
          do (format t "Adding include path: ~s~%" include-path)
          do (leap.core:add-path include-path))
    (loop for script in scripts
          do (format t "Sourcing script: ~s~%" script)
          do (leap:source script))
    ))

(defun parse-evaluate-leap-command (code)
  (let* ((*debugger-hook*
           #'(lambda (condition &rest argss)
               (format t "*debugger-hook* Encountered error ~s while parsing ~s~%" condition code)))
         (ast (architecture.builder-protocol:with-builder
                  ('list)
                (handler-bind ((esrap:esrap-parse-error
                                 (lambda (c)
                                   (format t "hander-bind Encountered error ~s while parsing ~s~%" c code))))
                  (esrap:parse 'leap.parser:leap code)))))
    (leap.core:evaluate 'list ast leap.core:*leap-env*)))

(defun leap-repl ()
  (process-command-line-options)
  (format t "Welcome to Cando-LEaP!~%")
  (clear-input)
  (unwind-protect
       (loop for code = (progn
                          (format t "> ") (finish-output)
                          (let ((line (read-line *standard-input* nil :eof)))
                            line))
             do (handler-case
                    (if (eq code :eof)
                        (progn
                          (clear-input *standard-input*))
                        (parse-evaluate-leap-command code))
                  (error (var)
                    (format t "Error ~a - while evaluating: ~a~%" var code))))
    (progn
      (format t "Goodbye~%"))))

(defun leap ()
  (leap-repl)
  (format t "Entering cando repl - use (leap) to re-enter leap.~%"))

(defun leap-repl-then-exit ()
  (if (leap-repl)
      (progn
        (format t "Entering cando repl - use (leap) to re-enter leap.~%"))
      (progn
        (format t "Leaving leap~%")
        (ext:quit 0))))






;;(defun solvate-box (solute solvent width-list
;;
;;)
;;; ----------------------------------------------------------------------
;;;
;;; Implement a simple leap script parser
;;;
;;; This uses the Common Lisp reader and it would probably be better
;;; to develop a proper parser
;;;
#++
(progn
  (defun leap-list-reader (stream char)
    "Read a Leap list delimited by #\} and return it."
    (let ((*readtable* (copy-readtable)))
      (set-syntax-from-char #\newline #\space)
      (read-delimited-list #\} stream t)))

  (defun leap-parse-entry (stream &optional (eof-error-p t) eof)
    "Parse a leap command from the stream and return it.
A leap command can have the form:
[command] [arg1] [arg2] [arg3] ... #\NEWLINE  - execute command
[variable] = [expression] #\NEWLINE           - set variable
[variable] = #\NEWLINE                        - clear variable
# ... #\NEWLINE                               - comment
An [expression] is either a command or a number or variable."
    (if (eq (peek-char nil stream eof-error-p eof) eof)
        eof
        (let ((*readtable* (copy-readtable))
              (*package* (find-package :keyword)))
          (setf (readtable-case *readtable*) :preserve)
          (set-macro-character #\{ #'leap-list-reader nil)
          (set-syntax-from-char #\# #\;)
          (set-syntax-from-char #\newline #\))
          (read-delimited-list #\newline stream nil))))

  (defun interpret-leap-entry (entry)
    "Interpret a leap command or variable assignment.
Nothing is returned, it's all side effects."
    (cond
      ((eq (second entry) :=)
       (set-variable entry))
      (t (let* ((cmd (intern (string-upcase (string (first entry))) :keyword))
                (cmd-assoc (assoc cmd leap.parser:*function-names/alist*))
                (cmd-func (cdr cmd-assoc)))
           (format *out* "cmd cmd-assoc cmd-func: ~a ~a ~a~%" cmd cmd-assoc cmd-func)
           (if (null cmd-func)
               (error "Illegal function ~a cmd-assoc: ~a  cmd-func: ~a" cmd cmd-assoc cmd-func)
               (progn
                 (funcall cmd-func entry)))))))

  (defun leap-interpreter (fin &optional eof-error-p eof)
    "Load the file of leap commands and execute them one by one.
Nothing is returned."
    (loop for entry = (leap-parse-entry fin eof-error-p eof)
       do (format *out* "entry = ~a~%" entry)
       unless (or (null entry) (eq entry eof))
       do (interpret-leap-entry entry)
       until (eq entry eof)))

  )
