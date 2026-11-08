(in-package :topology)

(defclass joint-template (cando.serialize:serializable)
  ((parent :initform nil :initarg :parent :accessor parent)
   (atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atoms-index :initarg :constitution-atoms-index :accessor constitution-atoms-index)
   (children :initform nil :initarg :children :accessor children)))


(defmethod print-object ((obj joint-template) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (format stream "~a" (atom-name obj)))))

(defclass bonded-joint-template (joint-template) ())

(defclass adjustable-bonded-joint-mixin ()
  ((adjustment :initarg :adjustment :accessor adjustment)))

(defclass adjustment ()
  ((joint :initarg :joint :accessor joint)))

(defclass internal-adjustment (adjustment) ()
  (:documentation "Adjust the internal coordinates of the joint when you can use internal coordinates. This is the best way to adjust coordinates.  If you need to adjust an internal coordinate to be planar relative to the other child of the parent of joing use topology:internal-planar-adjustment"))

(defclass external-adjustment (adjustment) ()
  (:documentation "Adjust the external coordinates of the joint when you cannot use internal coordinates but can only use external coordinates"))

(defclass adjustments ()
  ((internal-adjustments :initform (make-hash-table) :initarg :internal-adjustments :accessor internal-adjustments)
   (external-adjustments :initform (make-hash-table) :initarg :external-adjustments :accessor external-adjustments)))

(defgeneric add-to-adjustments (adjustments adjustment atresidue))

(defmethod add-to-adjustments (adjustments (adjustment internal-adjustment) atresidue)
  (push adjustment (gethash atresidue (internal-adjustments adjustments) nil)))

(defmethod add-to-adjustments (adjustments (adjustment external-adjustment) atresidue)
  (push adjustment (gethash atresidue (external-adjustments adjustments) nil)))

(defgeneric initialize-adjustment (adjustment assembler))

(defgeneric external-adjust (adjustment assembler assembler-internals coordinates))
(defmethod external-adjust ((adjustment internal-adjustment) assembler assembler-internals coordinates)
  "default - don't do anything"
  nil)

(defgeneric internal-adjust (adjustment assembler internals))

(defun validate-new-adjustment (spec)
  "Validate the declarative interface independently of assembled atom lookup."
  (unless (and (consp spec) (member (car spec) '(:tetrahedral-hydrogens :planar))
               (listp (cdr spec)) (evenp (length (cdr spec))))
    (error "Invalid :new-adjustment declaration: ~s" spec))
  (let ((options (cdr spec)) (seen nil)
        (allowed (if (eq (car spec) :planar)
                     '(:references)
                     '(:other-hydrogen :references :orientation))))
    (loop for (key value) on options by #'cddr
          do (unless (member key allowed)
               (error "Unknown :new-adjustment option ~s" key))
             (when (member key seen) (error "Duplicate adjustment option ~s" key))
             (push key seen))
    (when (eq (car spec) :tetrahedral-hydrogens)
      (unless (and (getf options :other-hydrogen)
                   (symbolp (getf options :other-hydrogen)))
        (error "Specify :other-hydrogen as an atom name: ~s" spec))
      (unless (member (getf options :orientation) '(:positive :negative))
        (error "Specify :orientation :positive or :negative: ~s" spec)))
    (let ((references (getf options :references)))
      (unless (and (listp references) (= (length references) 2)
                   (every (lambda (ref)
                            (or (and ref (symbolp ref))
                                (and (listp ref) (= (length ref) 2)
                                     (eq (first ref) :plug)
                                     (keywordp (second ref))))) references))
        (error "Specify two atom names or (:plug keyword) references: ~s" spec))
      (when (equal (first references) (second references))
        (error "Adjustment references must be distinct: ~s" spec))))
  spec)

(defclass named-adjustment-definition ()
  ((specification :initarg :specification :reader specification)))

(defvar *named-adjustment-definitions* (make-hash-table :test #'eq))

(defun register-adjustment-definition (name specification)
  (unless (and name (symbolp name) (not (keywordp name)) (not (eq name t)))
    (error "Adjustment name must be a non-keyword, non-NIL symbol: ~s" name))
  (setf (gethash name *named-adjustment-definitions*)
        (make-instance 'named-adjustment-definition
                       :specification (copy-tree (validate-new-adjustment specification))))
  name)

(defmacro define-adjustment (name specification)
  "Define a reusable, quoted declaration for :NEW-ADJUSTMENT NAME.
Names are package-qualified symbols, independent of legacy adjustment classes.
Define before DEFINE-TOPOLOGY. Redefinition affects subsequently built topology
templates; reload dependent topology definitions to update existing templates."
  `(register-adjustment-definition ',name ',specification))

(defun resolve-new-adjustment (specification)
  "Resolve a named or inline declaration to a fresh validated specification."
  (copy-tree
   (validate-new-adjustment
    (if (and specification (symbolp specification))
        (let ((definition (gethash specification *named-adjustment-definitions*)))
          (unless definition
            (error "Unknown :new-adjustment name ~s; define it with topology:define-adjustment first"
                   specification))
          (specification definition))
        specification))))

(defun atom-adjustment-specification (properties)
  "Keep legacy class names and new declarations on the same template path."
  (let ((keys (loop for (key value) on properties by #'cddr collect key)))
    (when (and (member :adjust keys) (member :new-adjustment keys))
      (error "An atom cannot specify both :adjust and :new-adjustment: ~s" properties))
    (if (member :new-adjustment keys)
        (resolve-new-adjustment (getf properties :new-adjustment))
        (getf properties :adjust))))

(defclass tetrahedral-hydrogens-adjustment (internal-adjustment)
  ((specification :initarg :specification :reader specification)
   (hydrogen-pair :accessor hydrogen-pair)
   (reference-joints :accessor reference-joints)))

(defclass planar-adjustment (internal-adjustment)
  ((specification :initarg :specification :reader specification)
   (reference-joints :accessor reference-joints))
  (:documentation "Place the adjusted atom opposite the bisector of two bonded
references, preserving its bond length. Reference order does not matter."))

(defun instantiate-joint-adjustment (spec joint)
  (if (consp spec)
      (progn
        (validate-new-adjustment spec)
        (make-instance (ecase (car spec)
                         (:tetrahedral-hydrogens 'tetrahedral-hydrogens-adjustment)
                         (:planar 'planar-adjustment))
                       :specification spec :joint joint))
      (make-instance spec :joint joint)))

(defun adjustment-plug-reference (assembler center plug-name)
  "Resolve the bond through PLUG-NAME at CENTER, including incoming and ring
couplings. Match plug-bond indices just as CONNECT-RESIDUES does; joint-tree
children and :OUT-PLUG-NAME alone do not describe ring-closing bonds."
  (destructuring-bind (mol-index res-index atom-index) (kin:id center)
    (declare (ignore atom-index))
    (let* ((monomer
             (or (loop for monomer being the hash-keys of (monomer-positions assembler)
                         using (hash-value position)
                       when (and (= mol-index (molecule-index position))
                                 (= res-index (residue-index position)))
                         return monomer)
                 (error "No monomer owns adjustment center ~s" center)))
           (coupling (or (monomer-plug-named monomer plug-name)
                         (error "Adjustment center ~s has no coupling at ~s" center plug-name)))
           (other (other-monomer coupling monomer))
           (other-plug-name
             (etypecase coupling
               (directional-coupling
                (if (eq monomer (source-monomer coupling))
                    (target-plug-name coupling) (source-plug-name coupling)))
               (ring-coupling
                (if (eq monomer (monomer1 coupling)) (plug2 coupling) (plug1 coupling)))))
           (atres (assembler-atresidue assembler monomer))
           (other-atres (assembler-atresidue assembler other))
           (plug (plug-named (topology atres) plug-name))
           (other-plug (plug-named (topology other-atres) other-plug-name)))
      (unless (and plug other-plug
                   (= (length (plug-bonds plug)) (length (plug-bonds other-plug))))
        (error "Invalid adjustment coupling from ~s through ~s" center plug-name))
      (let ((matches
              (loop for bond across (plug-bonds plug)
                    for other-bond across (plug-bonds other-plug)
                    when (eq (atom-name bond) (kin:joint/name center))
                      collect (joint-with-name other-atres (atom-name other-bond)))))
        (unless (= (length matches) 1)
          (error "Adjustment center ~s: plug ~s resolves to ~s" center plug-name matches))
        (first matches)))))

(defun resolve-adjustment-reference (assembler center neighbors reference)
  "Resolve exactly one chemical neighbor by atom name or coupling endpoint."
  (let ((matches
          (if (symbolp reference)
              (remove-if-not (lambda (neighbor)
                               (string= (symbol-name reference)
                                        (symbol-name (kin:joint/name neighbor))))
                             neighbors)
              (let ((neighbor (adjustment-plug-reference assembler center (second reference))))
                (when (member neighbor neighbors :test #'eq) (list neighbor))))))
    (unless (= (length matches) 1)
      (error "Adjustment center ~s: reference ~s must identify one bonded atom, found ~s"
             center reference matches))
    (first matches)))

(defmethod initialize-adjustment ((adjustment tetrahedral-hydrogens-adjustment) assembler)
  (let* ((hydrogen (joint adjustment))
         (center (kin:parent hydrogen))
         (neighbors (mapcar (lambda (atom) (find-joint-for-atom assembler atom))
                           (chem:atom/bonded-atoms-as-list
                            (find-atom-for-joint assembler center))))
         (options (cdr (specification adjustment))))
    (labels ((resolve (ref)
               (resolve-adjustment-reference assembler center neighbors ref)))
      (let* ((other (resolve (getf options :other-hydrogen)))
             (references (mapcar #'resolve (getf options :references)))
             (all (list* hydrogen other references)))
        (unless (and (= (length neighbors) 4)
                     (= (length (remove-duplicates all :test #'eq)) 4)
                     (every (lambda (h)
                              (and (typep h 'kin:bonded-joint)
                                   (eq (kin:parent h) center)
                                   (eq (chem:get-element (find-atom-for-joint assembler h)) :H)))
                            (list hydrogen other)))
          (error "Tetrahedral adjustment requires two bonded hydrogen children and two distinct references: ~s" all))
        (setf (hydrogen-pair adjustment) (list hydrogen other)
              (reference-joints adjustment) references)))))

(defmethod initialize-adjustment ((adjustment planar-adjustment) assembler)
  (let ((target (joint adjustment)))
    (unless (and (typep target 'kin:bonded-joint)
                 (kin:joint/parent-bound-p target)
                 (kin:joint/parent-bound-p (kin:parent target))
                 (kin:joint/parent-bound-p (kin:parent (kin:parent target))))
      (error "Planar adjustment needs a bonded joint with three ancestors: ~s" target))
    (let* ((center (kin:parent target))
           (neighbors (mapcar (lambda (atom) (find-joint-for-atom assembler atom))
                              (chem:atom/bonded-atoms-as-list
                               (find-atom-for-joint assembler center))))
           (references
             (mapcar (lambda (ref) (resolve-adjustment-reference assembler center neighbors ref))
                     (getf (cdr (specification adjustment)) :references))))
      (unless (and (= (length neighbors) 3)
                   (member target neighbors :test #'eq)
                   (= (length (remove-duplicates (cons target references) :test #'eq)) 3))
        (error "Planar adjustment on ~s requires a three-coordinate center and two other distinct bonded references"
               target))
      (setf (reference-joints adjustment) references))))

(defun planar-adjustment-direction (first-bond second-bond)
  "Unit vector opposite the two normalized reference bonds' bisector."
  (flet ((unit (v)
           (let ((norm (sqrt (reduce #'+ (map 'vector (lambda (x) (* x x)) v)))))
             (when (< norm 1.0d-10)
               (error "Degenerate reference geometry in planar adjustment"))
             (map 'vector (lambda (x) (/ x norm)) v))))
    (let* ((a (unit first-bond)) (b (unit second-bond))
           (dot (reduce #'+ (map 'vector #'* a b))))
      (when (> (abs dot) (- 1.0d0 1.0d-10))
        (error "Collinear references cannot define a planar adjustment"))
      (map 'vector #'- (unit (map 'vector #'+ a b))))))

(defmethod internal-adjust ((adjustment planar-adjustment) assembler internals)
  (let* ((target (joint adjustment))
         (center (kin:parent target))
         (grandparent (kin:parent center))
         (great-grandparent (kin:parent grandparent))
         (root (loop for ancestor = center then (kin:parent ancestor)
                     unless (kin:joint/parent-bound-p ancestor) return ancestor))
         (coords (make-coordinates-for-assembler assembler))
         (ic (kin:joint/position-index-x3 center))
         (itarget (kin:joint/position-index-x3 target))
         (distance (kin:bonded-joint/get-distance target internals)))
    ;; Reconstruct current coordinates so a ring-closure reference uses its
    ;; actual geometry, regardless of where it occurs in the atom tree.
    (kin:update-xyz-coords root internals coords)
    (flet ((bond (reference)
             (let ((index (kin:joint/position-index-x3 reference)))
               (map 'vector #'- (subseq coords index (+ index 3))
                               (subseq coords ic (+ ic 3))))))
      (let ((direction (planar-adjustment-direction
                        (bond (first (reference-joints adjustment)))
                        (bond (second (reference-joints adjustment))))))
        (unless (> distance 0.0d0) (error "Invalid planar-adjustment bond length ~s" distance))
        (dotimes (axis 3)
          (setf (aref coords (+ itarget axis))
                (+ (aref coords (+ ic axis)) (* distance (aref direction axis)))))
        (kin:bonded-joint/set-theta target internals
          (geom:calculate-angle-array (kin:joint/position-index-x3 grandparent) ic itarget coords))
        (kin:bonded-joint/set-phi target internals
          (geom:calculate-dihedral-array (kin:joint/position-index-x3 great-grandparent)
                                        (kin:joint/position-index-x3 grandparent) ic itarget coords))))))

(defun tetrahedral-hydrogen-directions (first-bond second-bond orientation)
  "Return unit vectors. Positive means (FIRST-BOND x SECOND-BOND) dot H1 > 0."
  (labels ((unit (v)
             (let ((norm (sqrt (reduce #'+ (map 'vector (lambda (x) (* x x)) v)))))
               (when (< norm 1.0d-10)
                 (error "Degenerate heavy-atom geometry in tetrahedral adjustment"))
               (map 'vector (lambda (x) (/ x norm)) v)))
           (cross (a b)
             (vector (- (* (aref a 1) (aref b 2)) (* (aref a 2) (aref b 1)))
                     (- (* (aref a 2) (aref b 0)) (* (aref a 0) (aref b 2)))
                     (- (* (aref a 0) (aref b 1)) (* (aref a 1) (aref b 0))))))
    (unless (member orientation '(:positive :negative))
      (error "Invalid tetrahedral orientation ~s" orientation))
    (let* ((a (unit first-bond)) (b (unit second-bond))
           (bisector (unit (map 'vector #'+ a b)))
           (normal (unit (cross a b)))
           (sign (if (eq orientation :positive) 1.0d0 -1.0d0)))
      (values (map 'vector (lambda (x n) (+ (- (/ x (sqrt 3.0d0)))
                                           (* sign (sqrt (/ 2.0d0 3.0d0)) n))) bisector normal)
              (map 'vector (lambda (x n) (- (- (/ x (sqrt 3.0d0)))
                                           (* sign (sqrt (/ 2.0d0 3.0d0)) n))) bisector normal)))))

(defmethod internal-adjust ((adjustment tetrahedral-hydrogens-adjustment) assembler internals)
  (let* ((hydrogen (joint adjustment))
         (center (kin:parent hydrogen))
         (grandparent (kin:parent center))
         (great-grandparent (kin:parent grandparent))
         (root (loop for ancestor = center then (kin:parent ancestor)
                     unless (kin:joint/parent-bound-p ancestor) return ancestor))
         (coords (make-coordinates-for-assembler assembler))
         (ic (kin:joint/position-index-x3 center)))
    ;; Use current internals, never stale atom positions. Do not recurse into adjustments.
    (kin:update-xyz-coords root internals coords)
    (flet ((bond (reference)
             (let ((index (kin:joint/position-index-x3 reference)))
               (map 'vector #'- (subseq coords index (+ index 3))
                               (subseq coords ic (+ ic 3))))))
      (multiple-value-bind (h1 h2)
          (tetrahedral-hydrogen-directions
           (bond (first (reference-joints adjustment)))
           (bond (second (reference-joints adjustment)))
           (getf (cdr (specification adjustment)) :orientation))
        (loop for h in (hydrogen-pair adjustment) for direction in (list h1 h2)
              for ih = (kin:joint/position-index-x3 h)
              for distance = (kin:bonded-joint/get-distance h internals)
              do (unless (> distance 0.0d0) (error "Invalid hydrogen bond length ~s" distance))
                 (dotimes (axis 3)
                   (setf (aref coords (+ ih axis))
                         (+ (aref coords (+ ic axis)) (* distance (aref direction axis)))))
        ;; Both positions are set before either internal coordinate is changed.
        (dolist (h (hydrogen-pair adjustment))
          (let ((ih (kin:joint/position-index-x3 h)))
            (kin:bonded-joint/set-theta h internals
              (geom:calculate-angle-array (kin:joint/position-index-x3 grandparent) ic ih coords))
            (kin:bonded-joint/set-phi h internals
              (geom:calculate-dihedral-array (kin:joint/position-index-x3 great-grandparent)
                                            (kin:joint/position-index-x3 grandparent) ic ih coords)))))))))

(defmethod internal-adjust ((adjustment external-adjustment) assembler internals)
  "default - don't do anything"
  nil)

(defmethod print-object ((obj adjustment) stream)
  "Name the joint this adjustment acts on.

Subclasses add whatever state decides whether they FIRE - see the INTERNAL-PLANAR-ADJUSTMENT
method.  Reporting an adjustment is the class's own job: an outside caller cannot know which slot
gates which subclass, and one that reads OTHER directly will describe BETA-CARBON - which has no
OTHER and is gated on nothing - as though it were disabled."
  (print-unreadable-object (obj stream :type t)
    (format stream "joint ~a"
            (if (slot-boundp obj 'joint) (kin:joint/name (joint obj)) :unbound))))

(defclass internal-planar-adjustment (internal-adjustment)
  ((other :initarg :other :accessor other))
  (:documentation "Rotate the joint 180 degrees from a sibling or ring-closure
reference, measured in the adjusted joint's internal-coordinate frame."))

(defmethod print-object ((obj internal-planar-adjustment) stream)
  "OTHER is included because INTERNAL-ADJUST is a no-op unless it is bound, so it is the difference
between an adjustment that fires and one that silently does nothing."
  (print-unreadable-object (obj stream :type t)
    (format stream "joint ~a other ~a"
            (if (slot-boundp obj 'joint) (kin:joint/name (joint obj)) :unbound)
            (if (slot-boundp obj 'other)
                (let ((o (other obj))) (if o (kin:joint/name o) nil))
                :unbound))))

(defun joint-plug-name (joint)
  "The :OUT-PLUG-NAME property WRITE-INTO-JOINT-TREE stamps on an in-plug joint
  (joint-templates.lisp:566-568), or NIL for a joint that does not start a coupling.
                                                    
  This is what says WHICH plug a child arrived through, and it is the only thing that separates a
  sidechain root from the backbone continuation without relying on atom names - which vary."
  (kin:get-property-or-default joint :out-plug-name nil))

(defun planar-adjustment-bonded-reference (joint assembler)
  "Find the chemical neighbor that supplies a missing sibling reference.
Exclude the adjusted atom, its grandparent, and alternate sidechains on the
same outgoing plug. The remaining bond can be an incoming ring closure."
  (let* ((parent (kin:parent joint))
         (grandparent (kin:parent parent))
         (plug (joint-plug-name joint))
         (candidates
           (loop for atom in (chem:atom/bonded-atoms-as-list
                              (find-atom-for-joint assembler parent))
                 for neighbor = (find-joint-for-atom assembler atom)
                 unless (or (eq neighbor joint)
                            (eq neighbor grandparent)
                            (and (kin:joint/parent-bound-p neighbor)
                                 (eq (kin:parent neighbor) parent)
                                 (eql (joint-plug-name neighbor) plug)))
                   collect neighbor)))
    (unless (= (length candidates) 1)
      (error "Planar adjustment on ~s: expected one bonded reference beyond its parent's tree parent and sidechain children, found ~s"
             joint candidates))
    (first candidates)))

(defmethod initialize-adjustment ((adjustment internal-planar-adjustment) assembler)
  ;; Cartesian/root joints do not run INTERNAL-ADJUST. In trainers these can
  ;; have several children; there is no internal dihedral or sibling to resolve.
  (unless (typep (joint adjustment) 'kin:bonded-joint)
    (return-from initialize-adjustment adjustment))
  (let* ((joint (joint adjustment))
         (jparent (kin:parent joint))
         (my-plug-name (joint-plug-name joint))
         (jother (if (null my-plug-name)
                     (if (= (kin:number-of-children jparent) 2)
                         (kin:joint/only-other-child jparent joint)
                         (error "~s on joint ~a: parent has ~d children and this joint has no :out-plug-name, ~
         so the sibling cannot be identified.  Expected exactly 2 children."
                                (class-name (class-of adjustment)) (kin:joint/name joint)
                                (kin:number-of-children jparent)))
                     (let ((candidates (remove-if (lambda (child)
                                                    (or (eq child joint)
                                                        (eql (joint-plug-name child) my-plug-name)))
                                                  (kin:joint/joint-children jparent))))
                       (cond
                         ((= (length candidates) 1) (first candidates))
                         ((null candidates)
                          (planar-adjustment-bonded-reference joint assembler))
                         (t (error "~s on joint ~a: expected exactly one sibling arriving through a plug other than ~s, ~
          found ~d among ~d children.  Under a fan-out every alternate shares this joint's plug ~
          and is excluded, so ~:[none remained~;the parent has more than one continuation~]."
                                  (class-name (class-of adjustment)) (kin:joint/name joint) my-plug-name
                                  (length candidates) (kin:number-of-children jparent) (plusp (length candidates)))))))))
    (setf (other adjustment) jother)))

(defmethod print-object ((adjustment internal-planar-adjustment) stream)
  (print-unreadable-object (adjustment stream :type t)
    (format stream "other ~s" (other adjustment))))

(defparameter *internal-planar-minimum-separation* 1.0d0
  "Minimum distance in Angstroms from an adjusted atom to either other neighbor of its parent.")

(defun planar-adjustment-neighbor-distances (distance theta phi
                                            other-distance other-theta other-phi
                                            parent-distance)
  "Return distances to the sibling and grandparent from current internal coordinates."
  (let* ((cos-separation (+ (* (cos theta) (cos other-theta))
                            (* (sin theta) (sin other-theta)
                               (cos (- phi other-phi)))))
         (sibling-squared (- (+ (* distance distance)
                                (* other-distance other-distance))
                             (* 2 distance other-distance cos-separation)))
         (grandparent-squared (- (+ (* distance distance)
                                    (* parent-distance parent-distance))
                                 (* 2 distance parent-distance (cos theta)))))
    (values (sqrt (max 0.0d0 sibling-squared))
            (sqrt (max 0.0d0 grandparent-squared)))))

(defun planar-adjustment-reference-internals (adjustment internals assembler)
  "Return reference distance, angle and dihedral in the adjusted joint's frame.
Siblings already use this frame. Ring-closure references must be measured from
current coordinates rather than using their own, unrelated internal frame."
  (let* ((joint (joint adjustment))
         (reference (other adjustment))
         (parent (kin:parent joint)))
    (if (and (kin:joint/parent-bound-p reference)
             (eq (kin:parent reference) parent)
             (typep reference 'kin:bonded-joint))
        (values (kin:bonded-joint/get-distance reference internals)
                (kin:bonded-joint/get-theta reference internals)
                (kin:bonded-joint/get-phi reference internals))
        (let* ((grandparent (kin:parent parent))
               (great-grandparent (kin:parent grandparent))
               (root (loop for ancestor = joint then (kin:parent ancestor)
                           unless (kin:joint/parent-bound-p ancestor)
                             return ancestor))
               (coords (make-coordinates-for-assembler assembler))
               (ia (kin:joint/position-index-x3 great-grandparent))
               (ib (kin:joint/position-index-x3 grandparent))
               (ic (kin:joint/position-index-x3 parent))
               (id (kin:joint/position-index-x3 reference)))
          ;; This primitive rebuilds coordinates without invoking adjustments.
          (kin:update-xyz-coords root internals coords)
          (values (geom:calculate-distance-array ic id coords)
                  (geom:calculate-angle-array ib ic id coords)
                  (geom:calculate-dihedral-array ia ib ic id coords))))))

(defun assert-planar-adjustment-separation (adjustment internals &optional assembler)
  (let* ((joint (joint adjustment))
         (other (other adjustment))
         (parent (kin:parent joint))
         (grandparent (kin:parent parent))
         (parent-distance
           (etypecase parent
             (kin:bonded-joint (kin:bonded-joint/get-distance parent internals))
             (kin:xyz-joint
              (geom:calculate-distance (kin:xyz-joint/transformed-pos parent)
                                       (kin:xyz-joint/transformed-pos grandparent))))))
    (multiple-value-bind (other-distance other-theta other-phi)
        (planar-adjustment-reference-internals adjustment internals assembler)
      (multiple-value-bind (sibling-distance grandparent-distance)
          (planar-adjustment-neighbor-distances
           (kin:bonded-joint/get-distance joint internals)
           (kin:bonded-joint/get-theta joint internals)
           (kin:bonded-joint/get-phi joint internals)
           other-distance other-theta other-phi
           parent-distance)
        (loop for neighbor in (list other grandparent)
              for distance in (list sibling-distance grandparent-distance)
              do (assert (>= distance *internal-planar-minimum-separation*) ()
                         "Internal planar adjustment of ~s bonded to ~s puts it only ~,4f A from ~s (minimum ~,4f A)."
                         joint parent distance neighbor *internal-planar-minimum-separation*))))))

(defmethod internal-adjust ((adjustment internal-planar-adjustment) assembler internals)
  ;; Only when the "other" slot is bound do we make adjustment
  (when (and (slot-boundp adjustment 'other)
             (typep (joint adjustment) 'kin:bonded-joint)
             (other adjustment))
    (let* ((joint (joint adjustment))
           (phi-original (kin:bonded-joint/get-phi joint internals))
           (phi-other (nth-value 2
                       (planar-adjustment-reference-internals adjustment internals assembler)))
           (phi-adjust (radians-add phi-other PI)))
      #+(or)(format t "internal-adjust ~s phi-original ~8,3f  phi-adjust ~8,3f~%"
                    joint
                    (topology:rad-to-deg phi-original) (topology:rad-to-deg phi-adjust))
      #+(or)(kin:bonded-joint/set-distance joint 1.47)
      #+(or)(kin:bonded-joint/set-theta joint (topology:deg-to-rad 120.0))
      (fill-joint-phi joint phi-adjust assembler internals)
      (assert-planar-adjustment-separation adjustment internals assembler)
      )))

(defun make-bonded-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defun make-adjustable-bonded-joint-template (constitution-atoms-index &key atom-name parent adjustment)
  (make-instance 'adjustable-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :adjustment adjustment))
(defclass in-plug-bonded-joint-template (bonded-joint-template)
  ((in-plug :initarg :in-plug :accessor in-plug)))

(defclass adjustable-in-plug-bonded-joint-template (in-plug-bonded-joint-template
                                                    adjustable-bonded-joint-mixin)
  ())

(defclass adjustable-bonded-joint-template (bonded-joint-template
                                            adjustable-bonded-joint-mixin)
  ())

(defun make-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug)
  (make-instance 'in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug))

(defun make-adjustable-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug adjustment)
  (make-instance 'adjustable-in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug
                 :adjustment adjustment))

(defclass complex-bonded-joint-template (bonded-joint-template)
  ((input-stub-joints :initform (make-array 2) :initarg :input-stub-joints :accessor input-stub-joints)))

(defmethod print-object ((obj complex-bonded-joint-template) stream)
  (if *print-readably*
      (call-next-method)
      (print-unreadable-object (obj stream :type t)
        (let ((stubs (loop for stub across (input-stub-joints obj)
                           when (not (null stub))
                             count 1)))
          (format stream "~a :stubs ~a" (atom-name obj) (1+ stubs))))))

(defun make-complex-bonded-joint-template (constitution-atoms-index &key atom-name stub-joints)
  (cond
    ((eq (length stub-joints) 1)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-element nil)))
    ((eq (length stub-joints) 2)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-contents (list (second stub-joints)
                                                                             nil))))
    ((eq (length stub-joints) 3)
     (make-instance 'complex-bonded-joint-template
                    :constitution-atoms-index constitution-atoms-index
                    :atom-name atom-name
                    :parent (first stub-joints)
                    :input-stub-joints (make-array 2 :initial-contents (list (second stub-joints)
                                                                             (third stub-joints)))))
    (t (error "Illegal number of stub-joints: ~s" stub-joints))))

(defclass jump-joint-template (joint-template) ())

(defclass xyz-joint-template (joint-template) ())

(defun make-jump-joint-template (constitution-atoms-index &key atom-name)
  (make-instance 'jump-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name))

(defun make-xyz-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'xyz-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defun add-child (joint-template child-template)
  (setf (parent child-template) joint-template)
  (setf (children joint-template) (append (children joint-template) (list child-template))))

(defun sibling (joint-template index)
  (nth index (children joint-template)))


(defun walk-joint-template (joint-template callback)
  "Walk a joint-template tree"
  (funcall callback joint-template)
  (loop for child in (children joint-template)
        do (walk-joint-template child callback)))

(defun new-joint-template-factory (parent-template node in-plug xyz-joints)
  (let* ((atom-name (topology:name node))
         (child-indexes (mapcar #'topology:constitution-atom-index (topology:children node)))
         (constitution-atoms-index (topology:constitution-atom-index node))
         (gparent-template (if parent-template
                               (parent parent-template)
                               nil))
         (ggparent-template (if gparent-template
                                (parent gparent-template)
                                nil)))
    (cond
      ((and xyz-joints
            (or (eq xyz-joints :all)
                (member atom-name xyz-joints)))
       (make-xyz-joint-template constitution-atoms-index
                                :atom-name atom-name
                                :parent parent-template))
      ((and (null parent-template) (typep in-plug 'topology:in-plug))
       (let ((adjust (atom-adjustment-specification (property-list node))))
         (if adjust
             (make-adjustable-in-plug-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent nil
              :in-plug in-plug
              :adjustment adjust)
             (make-in-plug-bonded-joint-template constitution-atoms-index
                                                 :atom-name atom-name
                                                 :parent nil
                                                 :in-plug in-plug))))
      ((typep in-plug 'topology:in-plug)
       (let ((adjust (atom-adjustment-specification (property-list node))))
         (if adjust
             (make-adjustable-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent nil
              :adjustment adjust)
             (make-bonded-joint-template constitution-atoms-index
                                         :atom-name atom-name
                                         :parent parent-template))))
      ((null parent-template)
       (make-jump-joint-template constitution-atoms-index
                                 :atom-name atom-name
                                 ))
      ((null gparent-template)
       (let ((stub-joints (cond
                            ((= 1 constitution-atoms-index)
                             (list parent-template))
                            ((= 2 constitution-atoms-index)
                             (list parent-template (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   (sibling parent-template 0)
                                   (sibling parent-template 1))))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      ((null ggparent-template)
       (let* ((sibling (when (> (length (children parent-template)) 0)
                         (elt (children parent-template) 0)))
              (aunt (when (> (length (children gparent-template)) 0)
                      (elt (children gparent-template) 0)))
              (stub-joints (cond
                            ((= constitution-atoms-index 2)
                             (list parent-template
                                   gparent-template))
                            (sibling
                             (list parent-template
                                   sibling
                                   gparent-template))
                            (aunt
                             (list parent-template
                                   gparent-template
                                   aunt))
                            (t
                             (error "One of the previous cases must be true")))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      (t
       (let ((adjust (atom-adjustment-specification (property-list node))))
         (if adjust
             (make-adjustable-bonded-joint-template
              constitution-atoms-index
              :atom-name atom-name
              :parent parent-template
              :adjustment adjust)
             (make-bonded-joint-template constitution-atoms-index
                                         :atom-name atom-name
                                         :parent parent-template)))))))


(defun build-joint-template-recursively (parent root in-plug xyz-joints)
  (let ((root-template (new-joint-template-factory parent root in-plug xyz-joints))
        (children (topology:children root)))
    (loop for child in children
          for sub-child-index from 0
          do (let ((child-template (build-joint-template-recursively root-template
                                                                     child
                                                                     in-plug
                                                                     xyz-joints)))
               (add-child root-template child-template)))
    root-template))

(defun build-joint-template (graph xyz-joints)
  (let ((root-node (topology:root-node graph)))
    (build-joint-template-recursively nil root-node (topology:in-plug graph) xyz-joints)))

(defun kin:joint-calculate-position-index (joint atom-table)
  (let ((atom-id (kin:joint/id joint)))
    (destructuring-bind (molecule-id residue-id atom-id)
        atom-id
      (let* ((mol (if (< molecule-id (length (chem:atom-table/atom-table-molecules atom-table)))
                      (aref (chem:atom-table/atom-table-molecules atom-table) molecule-id)
                      (error "The atom-table/atom-table-molecules is wrong")))
             (res (chem:content-at mol residue-id))
             (atm (chem:content-at res atom-id)))
        #+(or)(format t "joint-calculate-position-index number-of-atoms: ~a~%" (chem:atom-table/get-number-of-atoms atom-table))
        (unless (eq (chem:get-name atm) (kin:joint/name joint))
          (break "The joint ~a does not match the atom-table ~a" joint atom-table))
        (let* ((index*3 (let ((res (gethash atm (chem:atom-table/get-atom-table-indexes atom-table))))
                          (unless res
                            (error "Could not find atm ~a in atom-table ~a for joint ~a" atm atom-table joint))
                          res)))
          index*3)))))

(defun residue-shape-atom (nil-or-monomer-shape atom-name)
  (let* ((residue (residue nil-or-monomer-shape))
         (atm (chem:matter/first-atom-with-name residue atom-name)))
    (unless atm (error "Could not find atom with name ~s in ~s" atom-name residue))
    atm))

(defgeneric write-into-joint-tree (joint-template
                                   parent-joint
                                   atresidue
                                   atmolecule-index
                                   atresidue-index
                                   atom-table
                                   adjustments
                                   monomer
                                   nil-or-monomer-shape
                                   monomer-subset))

(defmethod write-into-joint-tree ((joint-template t)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset)
  (error "write-into-joint-tree - handle joint-template ~a" joint-template))

(defmethod write-into-joint-tree ((joint-template jump-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (etypecase nil-or-monomer-shape
                  (null (kin:make-jump-joint atomid atom-name atom-table))
                  (rotamer-shape (kin:make-jump-joint atomid atom-name atom-table))
                  (residue-shape (kin:make-xyz-joint atomid atom-name atom-table))))
         )
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defmethod write-into-joint-tree ((joint-template complex-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (input-stub-joints (input-stub-joints joint-template)))
    (flet ((make-complex-bonded-joint ()
             (let ((joint (kin:make-complex-bonded-joint atomid atom-name atom-table)))
               (put-joint atresidue joint constitution-atoms-index)
               (let ((input-stub0-template (aref input-stub-joints 0))
                     (input-stub1-template (aref input-stub-joints 1)))
                 (cond
                   ((and (null input-stub0-template) (null input-stub1-template))
                    ;; Do nothing
                    )
                   ((and input-stub0-template (null input-stub1-template))
                    (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                           (input-stub0 (aref (joints atresidue) input-stub0-index)))
                      (kin:set-input-stub-joint1 joint input-stub0)))
                   ((and (null input-stub0-template) input-stub1-template)
                    (let* ((input-stub1-index (constitution-atoms-index input-stub1-template))
                           (input-stub1 (aref (joints atresidue) input-stub1-index)))
                      (kin:set-input-stub-joint1 joint input-stub1)))
                   (t
                    (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                           (input-stub0 (aref (joints atresidue) input-stub0-index))
                           (input-stub1-index (constitution-atoms-index input-stub1-template))
                           (input-stub1 (aref (joints atresidue) input-stub1-index)))
                      (kin:set-input-stub-joint1 joint input-stub0)
                      (kin:set-input-stub-joint2 joint input-stub1))))
                 (when parent-joint (kin:joint/add-child parent-joint joint))
                 joint))))
      (etypecase nil-or-monomer-shape
        (null (make-complex-bonded-joint))
        (rotamer-shape (make-complex-bonded-joint))
        (residue-shape
         (let ((joint (kin:make-xyz-joint atomid atom-name atom-table)))
           (put-joint atresidue joint constitution-atoms-index)
           (when parent-joint (kin:joint/add-child parent-joint joint))
           joint))))))

(defgeneric make-appropriate-joint (joint-template nil-or-monomer-shape parent-joint atomid atom-name atom-table))

(defmethod write-into-joint-tree ((joint-template bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (when (and monomer-subset (null parent-joint))
    (error "Handle ~s joint if we have a monomer-subset" joint-template))
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (make-appropriate-joint joint-template nil-or-monomer-shape parent-joint atomid atom-name atom-table)))
    (put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

#+(or) ;; use create-anchored-joint
(defun create-stub-joint (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (let* ((base-assembler (backbone-assembler monomer-subset))
         (base-ataggregate (ataggregate base-assembler))
         (base-monomer-pos (gethash monomer (monomer-positions base-assembler)))
         (base-atresidue (at-position base-ataggregate base-monomer-pos))
         (base-joint0 (aref (joints base-atresidue) 0))
         (base-parent-joint (kin:parent base-joint0))
         (base-grand-parent-joint (kin:parent base-parent-joint))
         (base-great-grand-parent-joint (kin:parent base-grand-parent-joint))
         (constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-stub-joint atomid atom-name atom-table
                                     (kin:joint/position-index-x3 base-joint0)
                                     (kin:joint/position-index-x3 base-parent-joint)
                                     (kin:joint/position-index-x3 base-grand-parent-joint)
                                     (kin:joint/position-index-x3 base-great-grand-parent-joint)))
         (backbone-atom-tree-spliced-parent (gethash monomer (monomers monomer-subset))))
    (put-joint atresidue joint constitution-atoms-index)
    (kin:joint/set-parent joint backbone-atom-tree-spliced-parent)
    joint))

(defun create-anchored-joint (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (let* ((base-assembler (backbone-assembler monomer-subset))
         (base-ataggregate (ataggregate base-assembler))
         (base-monomer-pos (gethash monomer (monomer-positions base-assembler)))
         (base-atresidue (at-position base-ataggregate base-monomer-pos))
         (base-joint0 (aref (joints base-atresidue) 0))
         (base-parent-joint (kin:parent base-joint0))
         (base-grand-parent-joint (kin:parent base-parent-joint))
         (base-great-grand-parent-joint (kin:parent base-grand-parent-joint))
         (constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         ;; AnchoredBondedJoint computes CB from internals against the 3 anchors -
         ;; it does NOT need base-joint0's own index (StubJoint's fixed-cartesian source).
         (joint (kin:make-anchored-bonded-joint atomid atom-name atom-table
                                                (kin:joint/position-index-x3 base-parent-joint)
                                                (kin:joint/position-index-x3 base-grand-parent-joint)
                                                (kin:joint/position-index-x3 base-great-grand-parent-joint)))
         (backbone-atom-tree-spliced-parent (gethash monomer (monomers monomer-subset))))
    (put-joint atresidue joint constitution-atoms-index)
    (kin:joint/set-parent joint backbone-atom-tree-spliced-parent)
    joint))


(defun create-xyz-joint-connect-to-backbone (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-xyz-joint atomid atom-name atom-table))
         (backbone-atom-tree-spliced-parent (gethash monomer (monomers monomer-subset))))
    (put-joint atresidue joint constitution-atoms-index)
    (kin:joint/set-parent joint backbone-atom-tree-spliced-parent)
    joint))


(defun deal-with-pair-scan-sidechain (monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  #+(or)(create-xyz-joint-connect-to-backbone monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  #+(or)(create-stub-joint monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  (create-anchored-joint   monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
  )

(defmethod write-into-joint-tree ((joint-template xyz-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
             (atom-name (atom-name joint-template))
             (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
             (new-joint (kin:make-xyz-joint atomid atom-name atom-table)))
        (put-joint atresidue new-joint constitution-atoms-index)
        (when parent-joint (kin:joint/add-child parent-joint new-joint))
        new-joint)))

(defmethod write-into-joint-tree ((joint-template adjustable-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (when (and monomer-subset (null parent-joint))
    (error "Handle ~s joint if we have a monomer-subset" joint-template))
  (let ((joint (call-next-method)))
    (when adjustments
      (add-to-adjustments adjustments
                          (instantiate-joint-adjustment (adjustment joint-template) joint)
                          atresidue))
    joint))

(defmethod write-into-joint-tree ((joint-template in-plug-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let* ((joint (call-next-method))
             (in-plug (in-plug joint-template))
             (in-plug-name (name in-plug))
             (coupling-name (coupling-name in-plug-name))
             (out-plug-name (out-plug-name coupling-name)))
        (if (= (length (plug-bonds in-plug)) 1)
            (kin:joint/set-property joint :out-plug-name out-plug-name)
            (let ((indexed-out-plug-name (intern (format nil "~a.0" out-plug-name) :keyword)))
              (kin:joint/set-property joint :out-plug-name indexed-out-plug-name)))
        joint)))

(defmethod write-into-joint-tree ((joint-template adjustable-in-plug-bonded-joint-template)
                                  parent-joint
                                  atresidue
                                  atmolecule-index
                                  atresidue-index
                                  atom-table
                                  adjustments
                                  monomer
                                  nil-or-monomer-shape
                                  monomer-subset
                                  )
  (if (and monomer-subset (null parent-joint))
      (deal-with-pair-scan-sidechain monomer-subset monomer joint-template atmolecule-index atresidue-index atom-table atresidue)
      (let ((joint (call-next-method)))
        (when adjustments
          (add-to-adjustments adjustments
                              (instantiate-joint-adjustment (adjustment joint-template) joint)
                              atresidue))
        joint)))

(defclass topology-graph ()
  ((topology :initarg :topology :reader topology)
   (nodes :initarg :nodes :reader nodes)
   (nodes-to-ids :initarg :nodes-to-ids :reader nodes-to-ids)
   (names-to-nodes :initarg :names-to-nodes :reader names-to-nodes)
   (undirected-edges :initarg :undirected-edges :reader undirected-edges)))

(defclass topology-edge ()
  ((from :initarg :from :reader from)
   (to :initarg :to :reader to)
   (label :initarg :label :initform nil :reader label)))

(defmethod topology.graphviz:make-graph ((obj topology:topology))
  (let* ((nodes (labels ((accumulate-nodes (joint-template)
                           (list* joint-template
                                  (loop for child in (topology:children joint-template)
                                        append (accumulate-nodes child)))))
                  (accumulate-nodes (topology::joint-template obj))))
         (nodes-to-ids (make-hash-table))
         (names-to-nodes (make-hash-table)))
    (loop for node in nodes
          do (setf (gethash node nodes-to-ids) (gensym))
          do (setf (gethash (atom-name node) names-to-nodes) node))
    (let* ((constitution (topology:constitution obj))
           (constitution-atoms (topology:constitution-atoms constitution))
           (constitution-bond-edges nil))
      (loop for catom across constitution-atoms
            for cbond-from-index from 0
            for from-atom-name = (atom-name catom)
            do (loop for cbond in (topology:bonds catom)
                     for cbond-to-index = (topology:to-atom-index cbond)
                     for catom-to = (aref constitution-atoms cbond-to-index)
                     for to-atom-name = (atom-name catom-to)
                     when (< cbond-from-index cbond-to-index)
                       do (let ((edge (make-instance 'topology-edge
                                                     :from from-atom-name
                                                     :to to-atom-name
                                                     :label (case (order cbond)
                                                              (:single-bond nil)
                                                              (:double-bond "=")
                                                              (otherwise (string (order cbond)))))))
                            (push edge constitution-bond-edges))))
      (make-instance 'topology-graph
                     :topology obj
                     :nodes nodes
                     :nodes-to-ids nodes-to-ids
                     :names-to-nodes names-to-nodes
                     :undirected-edges constitution-bond-edges
                     ))))

(defmethod topology.graphviz:graph-label ((graph topology-graph))
  (string (name (topology graph))))

(defmethod topology.graphviz:graph-name ((graph topology-graph))
  (string (name (topology graph))))

(defmethod topology.graphviz:nodes ((obj topology-graph))
  (nodes obj))

(defmethod topology.graphviz:node-label ((node topology:joint-template) graph)
  (declare (ignore graph))
  (format nil "~a" (symbol-name (topology::atom-name node))))

(defmethod topology.graphviz:node-id ((node topology:joint-template) graph)
  (format nil "~a" (symbol-name (gethash node (nodes-to-ids graph)))))


(defmethod topology.graphviz:directed-edges ((obj topology-graph))
  (labels ((accumulate-directed-edges (node)
             (loop for child in (children node)
                   append (list*
                            (make-instance 'topology-edge
                                           :from node
                                           :to child)
                            (accumulate-directed-edges child)))))
    (accumulate-directed-edges (joint-template (topology obj)))))

(defmethod topology.graphviz:directed-edge-from (edge) (from edge))
(defmethod topology.graphviz:directed-edge-to (edge) (to edge))
(defmethod topology.graphviz:directed-edge-label (edge) (declare (ignore edge)) nil)

(defmethod topology.graphviz:undirected-edges ((obj topology-graph))
  (undirected-edges obj))

(defmethod topology.graphviz:undirected-edge-from (edge (graph topology-graph))
  (let ((node (gethash (from edge) (names-to-nodes graph))))
    node))
(defmethod topology.graphviz:undirected-edge-to (edge (graph topology-graph))
  (let ((node (gethash (to edge) (names-to-nodes graph))))
    node))

(defmethod topology.graphviz:undirected-edge-label (edge (graph topology-graph))
  (declare (ignore graph))
  (label edge))
