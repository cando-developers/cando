(in-package :leap.commands)


(defun leap.setup-amber-paths ()
  (leap:setup-amber-paths))

(defun leap.setup-default-paths ()
  (leap:setup-default-paths))

(defun leap.show-paths ()
  (leap:show-paths))

(defun leap.show (arg)
  (let ((val (leap.core:lookup-variable arg)))
    (funcall (find-symbol "SHOW" :cando-user) val)))


;;; ----------------------------------------------------------------------
;;;
;;; Quick and dirty setup-gaff
;;;

(defun setup-gaff ()
  (let ((*default-pathname-defaults*
          (translate-logical-pathname #P"source-dir:extensions;cando;src;data;force-field;"))
        (parms (chem:make-read-amber-parameters)))
    (with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
      (chem:read-types parms fin))
    (with-open-file (fin "gaff.dat" :direction :input)
      (chem:read-parameters parms fin (symbol-value (find-symbol "*AMBER-SYSTEM*" :leap))))
    (setf energy:*ff* (chem:get-force-field parms))
    energy:*ff*))

;;; ----------------------------------------------------------------------
;;;
;;; LEaP commands
;;;
(defvar *out* t)

(defun set-variable (entry)
  (let ((variable (first entry))
        (value (third entry)))
    (leap.core:register-variable variable value)))


(defun log-file (filename)
  #+(or)(format t "log-file doesn't do anything at the moment - to work with cl4py~%")
  (let* ((log-stream (open filename :direction :output :if-exists :supersede))
         (broadcast (make-broadcast-stream cl:*standard-output* #+(or)ext:+process-standard-output+ log-stream)))
    (setf cl:*standard-output* broadcast)))
  
(defun leap.log-file (entry)
  (valid-arguments entry 1)
  (let ((filename (pathname (ensure-string (second entry)))))
    (log-file filename)))

(defmethod describe-object ((object chem:aggregate) stream)
  (format stream "Aggregate name: ~a~%" (chem:get-name object))
  (let ((ht (make-hash-table :test #'equalp)))
    (cando:do-molecules (mol object)
      (let ((key (list (chem:get-name mol) (chem:number-of-atoms mol) (chem:force-field-name mol))))
        (push mol (gethash key ht))))
    (let* ((alist (alexandria:hash-table-alist ht))
           (sorted (sort alist #'> :key (lambda (entry) (second (first entry))))))
      (loop for (key . value) in sorted
            do (destructuring-bind (name number-of-atoms force-field-name)
                   key
                 (format stream "~a cop~:@P of molecule ~a with ~a atoms - force-field ~a~%" (length value) name number-of-atoms force-field-name))))))

(defun desc (&optional (name nil namep))
  "    desc variable
      object                       _variable_

Print a description of the object.
"
  (if namep
      (let ((val (leap.core:lookup-variable name)))
        (describe val))
      (format t "~a" (leap.core:all-variables))))

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
        (warn "There were ~a atoms that could not be assigned types" (- total-atoms total-atoms-with-types)))
      (format t "Assigned ~a types of a total of ~a atoms.~%" total-atoms-with-types total-atoms))))
    
(defun leap.load-smirnoff-params (filename)
  (leap:load-smirnoff-params filename))

(defun leap.set-force-field (matter force-field-name)
  (leap:set-force-field matter force-field-name))

(defun object (name)
  "Return the object with name."
  (leap.core:lookup-variable name))

(defun leap.desc (entry)
  (leap:desc entry))

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

(defun add-atom-types (info)
  (warn "add-atom-types doesn't do anything"))

(defun leap.add-atom-types (entry)
  (valid-arguments entry 1)
  (let ((info (second entry)))
    (add-atom-types info)))

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

(defun leap.list-force-fields ()
  (let ((ffs (leap:list-force-fields)))
    (loop for ff in ffs
          do (format t "~s ~a~%" (car ff) (cdr ff)))))

(defun leap.source (entry)
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
    (apply #'leap:solvate-shell solute solvent buffer :farness farness :shell t iso-closeness-args)))

(defun leap.solvate-cap (solute-name solvent-name position radius &rest iso-closeness)
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
        (iso-closeness-args (process-iso-closeness iso-closeness)))
    (apply #'leap:solvate-cap solute solvent position radius iso-closeness-args)))


(defun leap.dir (&optional (path *default-pathname-defaults*))
  (let* ((name (or (pathname-name path) "*"))
         (dirpath (make-pathname :name name :defaults path)))
    (format t "~a~%" (directory dirpath))))

(defun load-chem-draw (filename)
  (handler-bind ((warning #'muffle-warning))
    (let ((pathname (leap.core:ensure-path filename)))
      (with-open-file (fin (open pathname :direction :input))
        (chem:make-chem-draw fin :add-hydrogens nil)))))

(defun leap.load-chem-draw (filename)
  (load-chem-draw filename))

(defun leap.load-chem-draw-aggregate (filename)
  (handler-bind ((warning #'muffle-warning))
    (let ((pathname (leap.core:ensure-path filename)))
      (with-open-file (fin (open pathname :direction :input))
        (cando:load-chem-draw-aggregate pathname)))))

(defun leap.compile-smarts (smarts-string &optional tests)
  (unless (stringp smarts-string)
    (error "You must provide a smarts string"))
  (chem:compile-smarts smarts-string :tests tests))
    
(defun leap-quit ()
  (core:quit))

(defun leap-cando ()
  "
This command brings up a cando command line. 
Leap can be reentered by evaluating (leap).
"
  (core:top-level nil nil))


(defun leap-start-swank ()
  (funcall (find-symbol "START-SWANK" :cando-user)))

(defun leap-help (&optional arg)
  "    help [string]

      STRING                       _string_

This command prints a description of the command in _string_.  If
the STRING is not given then a list of legal STRINGs is provided.
"
  (if arg
      (progn
        (let ((found (assoc arg leap.parser:*function-names/alist* :test #'string=)))
          (if found
              (let ((cmd (cdr found)))
                (format t "~%~a~%" (documentation cmd 'function)))
              (format t "No help available for ~a~%" arg))))
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
  (chem:add-matter object-a object-b)))

(defun leap-bond (atom1 atom2 &rest order)
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
  (case order
    (:s (chem:bond-to atom1 atom2 :single-bond))
    (:d (chem:bond-to atom1 atom2 :double-bond))
    (:t (chem:bond-to atom1 atom2 :triple-bond))
    (:a (chem:bond-to atom1 atom2 :aromatic-bond))
    (otherwise (chem:bond-to atom1 atom2 :single-bond))))

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
      UNIT                         _variable_
      STRING                       _name_

Return a new and empty UNIT with the name _name_.
"
  (let ((aggregate (chem:make-aggregate name)))
    aggregate))

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
   (let ((container (leap.core:lookup-variable container-name)))
     (format t "Total unperturbed charge:    ~,6f~%" (chem:get-charge container))
     (format t "Total perturbed charge:      ~,6f~%" (chem:get-charge container))))

(defun leap-add-ions (mol-name ion1-name ion1-number &optional ion2-name ion2-number)
"
    addIons unit ion1 #ion1 [ion2 #ion2]
      UNIT                      _unit_
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

(defun leap-load-moe (filename)
  (let ((pathname (leap.core:ensure-path filename)))
    (unless pathname
      (error "Could not file file ~a" filename))
    (chem:moe-read-aggregate (namestring pathname))))
    

(eval-when (:load-toplevel :execute)
  (setf leap.parser:*function-names/alist*
        '(("logFile" . log-file)
          ("setupAmberPaths" . leap.setup-amber-paths)
          ("desc" . desc)
          ("listForceFields" . leap.list-force-fields)
          ("loadOff" . load-off)
          ("loadMol2" . cando:load-mol2)
          ("loadMoe" . leap-load-moe )
          ("dir" . leap.dir )
          ("ls" . leap.dir)
          ("loadPdb" . leap.pdb:load-pdb)
          ("source" . leap.source)
          ("loadChemDraw" . leap.load-chem-draw)
          ("loadChemDrawAggregate" . leap.load-chem-draw-aggregate)
          ("compileSmarts" . leap.compile-smarts)
          ("setForceField" . leap.set-force-field)
          ("loadSmirnoffParams" . leap.load-smirnoff-params)
          ("loadAmberParams" . leap:load-amber-params)
          ("assignAtomTypes" . leap.assign-atom-types)
          ("addPdbResMap" . leap.pdb:add-pdb-res-map)
          ("addPdbAtomMap" . leap.pdb:add-pdb-atom-map)
          ("addAtomTypes" . add-atom-types)
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
          ("translate" . leap-translate)
          ("transform" . leap-transform)
          ("copy" . leap-copy)
          ("center" . leap-center)
          ("measureGeom" . leap-measure-geom)
          ("addPath" . leap-add-path)
          ("alignAxes" . leap-align-axes)
          ("charge" . leap-charge)
          ))
  (dolist (command leap.parser:*function-names/alist*)
    (if (fboundp (cdr command))
        (setf (leap.core:function-lookup (car command) leap.core:*leap-env*) (cdr command))
        (error "~a is not a function" (cdr command)))))


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
    (format t "ast -> ~a~%" ast)
    (core:call-with-stack-top-hint
     (lambda ()
       (leap.core:evaluate 'list ast leap.core:*leap-env*)))))

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
        (core:exit 0))))






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



