(in-package :design.joint-tree)

(defclass joint-template ()
  ((parent :initform nil :initarg :parent :accessor parent)
   (atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atoms-index :initarg :constitution-atoms-index :accessor constitution-atoms-index)
   ))

(defclass bonded-joint-template (joint-template)
  ((children :initform nil :initarg :children :accessor children)
   ))

(defun make-bonded-joint-template (constitution-atoms-index &key atom-name parent)
  (make-instance 'bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent))

(defclass in-plug-bonded-joint-template (bonded-joint-template)
  ((in-plug :initarg :in-plug :accessor in-plug)))

(defun make-in-plug-bonded-joint-template (constitution-atoms-index &key atom-name parent in-plug)
  (make-instance 'in-plug-bonded-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name
                 :parent parent
                 :in-plug in-plug))

(defclass complex-bonded-joint-template (bonded-joint-template)
  ((input-stub-joints :initform (make-array 2) :initarg :input-stub-joints :accessor input-stub-joints)))

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

(defclass jump-joint-template (joint-template)
  ((children :initform nil :initarg :children :accessor children)))

(defun make-jump-joint-template (constitution-atoms-index &key atom-name)
  (make-instance 'jump-joint-template
                 :constitution-atoms-index constitution-atoms-index
                 :atom-name atom-name))

(defun add-child (joint-template child-template)
  (setf (parent child-template) joint-template)
  (setf (children joint-template) (append (children joint-template) (list child-template))))

(defun sibling (joint-template index)
  (nth index (children joint-template)))

(defun new-joint-template-factory (parent-template atom child-indexes in-plug constitution-atoms constitution-name topology-name)
  (declare (ignore atom constitution-atoms topology-name constitution-name))
  (let* ((atom-name (chem:get-name atom))
         (constitution-atoms-index (chem:constitution-atoms/index constitution-atoms atom-name))
         (gparent-template (if parent-template
                               (parent parent-template)
                               nil))
         (ggparent-template (if gparent-template
                                (parent gparent-template)
                                nil)))
    (cond
      ((and (null parent-template) (typep in-plug 'chem:in-plug))
       (make-in-plug-bonded-joint-template constitution-atoms-index
                                           :atom-name atom-name
                                           :parent nil
                                           :in-plug in-plug))
      ((typep in-plug 'chem:in-plug)
       (make-bonded-joint-template constitution-atoms-index
                                   :atom-name atom-name
                                   :parent parent-template))
      ((null parent-template)
       (make-jump-joint-template constitution-atoms-index
                                 :atom-name atom-name
                                 ))
      ((null gparent-template)
       (let ((stub-joints (cond
                            ((eql 0 (first child-indexes))
                             (list parent-template))
                            ((eql 1 (first child-indexes))
                             (list parent-template (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   (sibling parent-template 1)
                                   (sibling parent-template 0))))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      ((null ggparent-template)
       (let ((stub-joints (cond
                            ((and (eql 0 (first child-indexes)) (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template))
                            ((and (> 0 (first child-indexes)) (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template
                                   (sibling parent-template 0)))
                            (t
                             (list parent-template
                                   gparent-template
                                   (sibling gparent-template 0))))))
         (make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      (t (make-bonded-joint-template constitution-atoms-index
                                     :atom-name atom-name
                                     :parent parent-template)))))

;;; ------------------------------------------------------------
;;; ------------------------------------------------------------
;;; ------------------------------------------------------------
;;; ------------------------------------------------------------
;;; ------------------------------------------------------------
;;; ------------------------------------------------------------
;;; ------------------------------------------------------------

#+(or)
(progn
  (defclass joint ()
    ((parent :initarg :parent :accessor parent)
     (children :initarg :children :initform (make-array 4 :fill-pointer 0 :adjustable t) :accessor children)))

  (defun add-child (joint-template child-template)
    (setf (parent child-template) joint-template)
    (vector-push-extend child-template (children joint-template)))

  (defclass checkpoint-out-plug-joint (joint)
    ((constitution-name :initarg :constitution-name :accessor constitution-name)
     (topology-name :initarg :topology-name :accessor topology-name)
     (out-plug :initarg :out-plug :accessor out-plug)
     ))

  (defun make-checkpoint-out-plug-joint (&rest args)
    (apply 'make-instance 'checkpoint-out-plug-joint args))

  (defclass checkpoint-joint (joint)
    ((constitution-name :initarg :constitution-name :accessor constitution-name)
     (topology-name :initarg :topology-name :accessor topology-name)
     (out-plug :initarg :out-plug :accessor out-plug)
     ))

  (defun make-checkpoint-joint (&rest args)
    (apply 'make-instance 'checkpoint-joint args))


  (defclass delayed-bonded-joint-template (joint)
    (
     (id :initarg :id :accessor id )
     (name :initarg :name :accessor name )
     (parent :initarg :parent :accessor parent )
     (children :initarg :children :accessor children )
     (checkpoint :initarg :checkpoint :accessor checkpoint )
     (comment :initarg :comment :accessor comment )
     (out-plug :initarg :out-plug :accessor out-plug )
     ))

  (defun make-delayed-bonded-joint-template (&rest args)
    (apply 'make-instance 'delayed-bonded-joint-template args))


  (defclass root-bonded-joint-template (joint)
    (
     (id :initarg :id :accessor id )
     (parent :initarg :parent :accessor parent )
     (name :initarg :name :accessor name )
     (children :initarg :children :accessor children )
     (constitution-name :initarg :constitution-name :accessor constitution-name )
     (topology-name :initarg :topology-name :accessor topology-name )
     (in-plug :initarg :in-plug :accessor in-plug )
     (comment :initarg :comment :accessor comment )
     (out-plug :initarg :out-plug :accessor out-plug )
     ))

  (defmethod print-object ((obj root-bonded-joint-template) stream)
    (print-unreadable-object (obj stream :type t)
      (format stream ":name ~a" (name obj))))

  (defun make-root-bonded-joint-template (&rest args)
    (apply 'make-instance 'root-bonded-joint-template args))

  (defclass bonded-joint-template (joint)
    (
     (id :initarg :id :accessor id )
     (parent :initarg :parent :accessor parent )
     (name :initarg :name :accessor name )
     (children :initarg :children :accessor children )
     (comment :initarg :comment :accessor comment )
     (out-plug :initarg :out-plug :accessor out-plug )
     ))

  (defun make-bonded-joint-template (&rest args)
    (apply 'make-instance 'bonded-joint-template args))
  )
  
  #+(or)
(defun interpret-builder-info-instruction (instr residue)
  (error "This is currently not used")
  (cond
    ((eq (car instr) 'chem:origin-plug)
     (destructuring-bind (cmd atom-name plug-name) instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-origin-plug plug-name atom-name)))))
    ((eq (car instr) 'chem:jum*p-plug)
     (destructuring-bind (cmd atom-name plug-name) instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-jump-plug plug-name atom-name)))))
    ((eq (car instr) 'chem:in-plug)
     (destructuring-bind (cmd
                          atom-name
                          plug-name
                          &key
                            frame
                            (bond-order :single-bond)
                            atom2-name
                            (bond-order2 :single-bond))
         instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-in-plug plug-name nil atom-name bond-order atom2-name bond-order2)))))
    ((eq (car instr) 'chem:out-plug)
     (destructuring-bind (cmd
                          atom-name
                          plug-name
                          &key
                            (bond-order :single-bond)
                            atom2-name
                            (bond-order2 :single-bond))
         instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (chem:make-out-plug plug-name nil nil atom-name bond-order atom2-name bond-order2)))))
    ((error "handle (eq (car instr) 'chem:ring-closing-plug)")
     (destructuring-bind (cmd
                          atom-name
                          plug-name
                          &key
                            (bond-order :single-bond)
                            atom2-name
                            (bond-order2 :single-bond))
         instr
       (declare (ignore cmd))
       (let* ((atom (when (chem:has-atom-with-name residue atom-name)
                      (chem:first-atom-with-name residue atom-name))))
         (when atom
           (error "(chem:make-ring-closing-plug plug-name nil nil atom-name bond-order atom2-name bond-order2)")
           ))))
    (t (error "Unknown instruction ~a" instr))))

(defvar *valid-atom-properties*
  '(:in :origin :in/origin/out :in/out :out))

(defclass plug-factory ()
  ((in-reference-atom :initarg :in-reference-atom :initform nil :reader in-reference-atom)
   (in-name :initarg :in-name :reader in-name)
   (in-plug :initform nil :accessor in-plug)
   (origin-atom :initarg :origin-atom :initform nil :reader origin-atom)
   (origin-name :initarg :origin-name :reader origin-name)
   (origin-plug :initform nil :accessor origin-plug)
   (out-reference-atom :initarg :out-reference-atom :initform nil :reader out-reference-atom)
   (out-name :initarg :out-name :reader out-name)
   (out-plug :initform nil :accessor out-plug))
  (:documentation "Create plugs of different kinds" ))

(defun maybe-make-origin-plug (origin-plug-factory)
  "Either return the cached plug or create and cache one"
  (if (origin-plug origin-plug-factory)
      (origin-plug origin-plug-factory)
      (setf (origin-plug origin-plug-factory)
            (chem:make-origin-plug (origin-name origin-plug-factory)
                                   (chem:get-name (origin-atom origin-plug-factory))))))

(defun maybe-make-in-plug (in-plug-factory complex-plugs)
  "Either return the cached plug or create and cache one"
  (if (in-plug in-plug-factory)
      (in-plug in-plug-factory)
      (let ((reference-atom (in-reference-atom in-plug-factory))
            (maybe-complex-plug (design:complex-plug-or-nil (in-name in-plug-factory) complex-plugs)))
        (if maybe-complex-plug
            (let* ((compiled-smarts (design:compiled-smarts maybe-complex-plug))
                   (matches (chem:matches compiled-smarts reference-atom)))
              (unless matches
                (error "No match could be found for in-plug ~a reference-atom ~a smarts ~a"
                       (name maybe-complex-plug)
                       reference-atom
                       (design:smarts maybe-complex-plug)))
              (let ((atom0 (chem:get-atom-with-tag matches 0))
                    (atom1 (chem:get-atom-with-tag matches 1)))
                (format t "matches: ~s ~s ~s~%" matches atom0 atom1)
                (setf (in-plug in-plug-factory)
                      (chem:make-in-plug (in-name in-plug-factory)
                                         nil
                                         (chem:get-name atom0)
                                         :single-bond
                                         (chem:get-name atom1)
                                         :single-bond))))
            (setf (in-plug in-plug-factory)
                  (chem:make-in-plug (in-name in-plug-factory)
                                     nil
                                     (chem:get-name reference-atom)
                                     :single-bond
                                     nil ; Handle two bond plugs
                                     :no-bond))))))

(defun maybe-make-out-plug (out-plug-factory complex-plugs)
  "Either return the cached plug or create and cache one"
  (if (out-plug out-plug-factory)
      (out-plug out-plug-factory)
      (let ((reference-atom (out-reference-atom out-plug-factory))
            (maybe-complex-plug (design:complex-plug-or-nil (out-name out-plug-factory) complex-plugs)))
        (if maybe-complex-plug
            (let* ((compiled-smarts (design:compiled-smarts maybe-complex-plug))
                   (matches (chem:matches compiled-smarts reference-atom)))
              (unless matches
                (error "No match could be found for out-plug ~a reference-atom ~a smarts ~a"
                       (name maybe-complex-plug)
                       reference-atom
                       (design:smarts maybe-complex-plug)))
              (let ((atom0 (chem:get-atom-with-tag matches 0))
                    (atom1 (chem:get-atom-with-tag matches 1)))
                (format t "matches: ~s ~s ~s~%" matches atom0 atom1)
                (setf (out-plug out-plug-factory)
                      (chem:make-out-plug (out-name out-plug-factory)
                                          nil
                                          nil
                                          (chem:get-name atom0)
                                          :single-bond
                                          (chem:get-name atom1)
                                          :single-bond))))
            (setf (out-plug out-plug-factory)
                  (chem:make-out-plug (out-name out-plug-factory)
                                      nil
                                      nil
                                      (chem:get-name reference-atom)
                                      :single-bond
                                      nil ; Handle two bond plugs
                                      :no-bond))))))

(defun gather-plugs (plug-factory)
  "Return a list of plugs that the plug-factory created"
  (let (plugs)
    (when (out-plug plug-factory) (push (out-plug plug-factory) plugs))
    (when (in-plug plug-factory) (push (in-plug plug-factory) plugs))
    (when (origin-plug plug-factory) (push (origin-plug plug-factory) plugs))))
    
(defun in-plug-factory-p (plug-factory)
  "Return generalized boolean for whether this plug-factory is capable of generating an in-plug"
  (in-reference-atom plug-factory))

(defun origin-plug-factory-p (plug-factory)
  "Return generalized boolean for whether this plug-factory is capable of generating an origin-plug"
  (origin-atom plug-factory))

(defun out-plug-factory-p (plug-factory)
  "Return generalized boolean for whether this plug-factory is capable of generating an out-plug"
  (out-reference-atom plug-factory))

(defun gather-plug-factories (residue)
  "Loop over the atoms in the residue and use the properties to create plug-factory(s).
The properties can look like:
:in name (default name is :-default)
:out name (default name is :+default)
:origin name (default name is :-origin)
:in/out (default name is (:-default :+out)
:in/origin/out (in-name origin-name out-name) default names are (:-default :-origin :+origin)"
  (flet ((plug-factory-argument-error (atom property value type msg)
           (declare (ignore atom))
           (error "The plug argument was incorrect for ~a - the ~a plug was ~a but it must be a ~a - ~a"
                  (chem:get-name residue) property value type msg))
         (in-plug-name (name &optional (default :-default)) (if name name default))
         (origin-plug-name (name &optional (default :+origin)) (if name name default))
         (out-plug-name (name &optional (default :+default)) (or name default))
         )
    (let (plug-factories)
      (cando:do-atoms (atom residue)
        (when (chem:properties atom)
          (loop for cur = (chem:properties atom) then (cddr cur)
                for save-cur = cur
                for property = (first cur)
                for value = (second cur)
                until (null cur)
                do (let ((maybe-plug-factory
                           (case property
                             (:in
                              (unless (symbolp value) (plug-factory-argument-error atom property value 'symbol "Expected plug-name"))
                              (make-instance 'plug-factory :in-reference-atom atom
                                                           :in-name (in-plug-name value)))
                             (:out
                              (unless (symbolp value) (plug-factory-argument-error atom property value 'symbol "Expected plug-name)"))
                              (make-instance 'plug-factory :out-reference-atom atom
                                                           :out-name (out-plug-name value)))
                             (:in/out
                              (unless (listp value) (plug-factory-argument-error atom property value 'symbol "Expected plug-name)"))
                              (make-instance 'plug-factory :in-reference-atom atom
                                                           :in-name (in-plug-name (first value) :-default)
                                                           :out-reference-atom atom
                                                           :out-name (out-plug-name value :+out)))
                             (:origin
                              (make-instance 'plug-factory :origin-atom atom
                                                           :origin-name (origin-plug-name value :-origin)))
                             (:in/origin/out
                              (unless (listp value) (plug-factory-argument-error atom property value 'symbol "Expected plug-name)"))
                              (make-instance 'plug-factory :in-reference-atom atom
                                                           :in-name (in-plug-name (first value))
                                                           :origin-atom atom
                                                           :origin-name (origin-plug-name (second value) :-origin)
                                                           :out-reference-atom atom
                                                           :out-name (out-plug-name (third value) :+origin))))))
                     (when maybe-plug-factory (push maybe-plug-factory plug-factories))))))
      plug-factories)))

(defclass prepare-topology ()
  ((name :initarg :name :reader name)
   (constitution :initarg :constitution :reader constitution)
   (constitution-atoms :initarg :constitution-atoms :reader constitution-atoms)
   (plugs :initarg :plugs :reader plugs)
   (stereo-informaton :initarg :stereo-information :reader stereo-information)
   (residue :initarg :residue :reader residue))
  (:documentation "Store information needed to create a new topology"))

(defmethod print-object ((obj prepare-topology) stream)
  (print-unreadable-object (obj stream)
    (format stream "~a ~a" (class-name (class-of obj)) (name obj))))


(defun extract-residue-plug-factory-alist (chem-draw-object)
  "Loop through all of the fragments in the **chem-draw-object** and then loop through all of the residues
in each fragment. Those residues that have the property :topology are topology's that we need to extract.
Return an alist of (residue . plug-factory-list)."
  (declare (debug 3))
  (let ((residue-plug-factory-alist))
    (loop for fragment in (chem:all-fragments-as-list chem-draw-object)
          do (cando:do-residues (residue (chem:get-molecule fragment))
               (when (chem:has-property residue :topology)
                 (chem:set-name residue (chem:matter-get-property residue :topology))
                 (let ((plug-factories (gather-plug-factories residue)))
                   (push (cons residue plug-factories) residue-plug-factory-alist)))))
    residue-plug-factory-alist))

(defun identify-root-atom-name (residue plug-factories)
  "Return the name of the atom that will serve to organize the stereocenters.
It's either the origin atom or the first in atom found."
  (loop for plug-factory in plug-factories
        when (origin-plug-factory-p plug-factory)
          do (return-from identify-root-atom-name (chem:get-name (origin-atom plug-factory)))
        when (in-plug-factory-p plug-factory)
          do (return-from identify-root-atom-name (chem:get-name (in-reference-atom plug-factory))))
  (error "There must either be an origin-plug or an in-plug in the residue ~a" residue))

(defun build-stereoisomer (name chiral-atoms stereoisomer-index)
  "Build a stereoisomer. Assemble a name for the stereiosomer using the name and
the stereoisomer-index.   The stereoisomer-index is like a bit-vector that indicates
the stereochemistry of this stereoisomer (0 is :S and 1 is :R).
The name is constructed by appending \"{chiral-atom-name/(S|R),...}\".
So if name is \"ALA\" and stereoisomer-index is 1 the name becomes ALA{CA/S}."
  (if chiral-atoms
      (let* ((configurations (loop for chiral-index below (length chiral-atoms)
                                   for config = (if (logbitp chiral-index stereoisomer-index) :S :R)
                                   collect config))
             (stereo-configurations (mapcar (lambda (atom config) (chem:make-stereo-configuration (chem:get-name atom) config))
                                            chiral-atoms configurations))
             (new-name-string (format nil "~A[~{~A~^,~}]"
                                      name
                                      (mapcar (lambda (atom config)
                                                (format nil "~A=~A"
                                                        (string (chem:get-name atom))
                                                        (string config))) chiral-atoms configurations)))
             (new-name (intern new-name-string :keyword)))
        #+(or)
        (progn
          (format t "name -> ~a~%" name)
          (format t "chiral-atoms -> ~a~%" chiral-atoms)
          (format t "stereoisomer-index -> ~a~%" stereoisomer-index)
          (format t "new-name -> ~a~%" new-name))
        (chem:make-stereoisomer new-name nil stereoisomer-index stereo-configurations))
      (progn
        ;;#+(or)(format t "Single stereoisomer name -> ~a~%" name)
        (chem:make-stereoisomer name name 0 nil))))

(defun build-stereoinformation (name constitution-atoms root-atom-name)
  "Build a residue from the constitution-atoms in **prep-topology**.
Then identify all stereocenters sorted by distance from the **root-atom-name**.
Build stereoinformation from all of this and return it."
  (let* ((residue (chem:build-residue constitution-atoms)))
;;;      (chem:define-stereochemical-configurations-for-all-atoms residue)
    (let* ((root-atom (chem:atom-with-name residue root-atom-name))
           (spanning-loop (chem:spanning-loop/make root-atom))
           (atoms (chem:all-atoms spanning-loop))
           (unsorted-chiral-atoms (remove-if (lambda (x) (null (eq (chem:get-stereochemistry-type x) :chiral))) atoms))
           (chiral-atoms (sort unsorted-chiral-atoms #'string<
                               :key (lambda (atm) (string (chem:get-name atm)))))
           (number-of-chiral-atoms (length chiral-atoms))
           (number-of-stereoisomers (expt 2 number-of-chiral-atoms)))
      #+(or)
      (progn
        (format t "residue: ~s~%" residue)
        (format t "root atom: ~s~%" root-atom)
        (format t "chiral atoms: ~s~%" chiral-atoms)
        (terpri)
        (terpri))
      (let ((stereoisomers (loop for stereoisomer-index below number-of-stereoisomers
                                 for stereoisomer = (build-stereoisomer name chiral-atoms stereoisomer-index)
                                 collect stereoisomer)))
        (chem:make-stereoinformation stereoisomers nil)))))



(defun prepare-topologys-from-plug-factories (name residue plug-factories constitution constitution-atoms
                                              stereo-information complex-plugs)
  "Create prepare-topology objects using the list of plug-factory objects.
If there is a plug-factory that can create an origin that is one topology.
For every plug-factory that can generate an in-plug use that to generate another topology.
All other plug-factory objects must be capable of generating out-plugs."
  (let (result
        (origin-plug-factory (find-if #'origin-plug-factory-p plug-factories)))
    (when origin-plug-factory
      (let ((counter 0))
        (loop for plug-factory in plug-factories
              if (origin-plug-factory-p plug-factory)
                count counter)
        (when (> counter 1) (error "There can only be one plug factory capable of generating an origin-plug"))))
      (when origin-plug-factory
        (let* ((origin-name (origin-name origin-plug-factory))
               (origin-atom-name (chem:get-name (origin-atom origin-plug-factory)))
               (origin-plug (maybe-make-origin-plug origin-plug-factory))
               (out-plugs (loop for plug-factory in plug-factories
                                if (out-plug-factory-p plug-factory)
                                  collect (maybe-make-out-plug plug-factory complex-plugs)
                                else
                                  do (when (and (in-plug-factory-p plug-factory)
                                                (null (out-plug-factory-p plug-factory)))
                                       (error "In residue ~a you have specified an origin-plug ~a on atom ~a and elsewhere there is atom ~a indicated as an in-plug ~a - it must also be convertable to an out-plug"
                                              origin-name
                                              origin-atom-name
                                              (if (slot-boundp plug-factory 'out-name) (out-name plug-factory) :unbound)
                                              (chem:get-name (out-reference-atom plug-factory))))))
                (plugs (cons origin-plug out-plugs)))
           ;; If the topology has an origin-plug then name it XXX.ORIGIN
           (push (make-instance 'prepare-topology
                                :name (intern (format nil "~a.ORIGIN" (string name)) :keyword)
                                :constitution constitution
                                :constitution-atoms constitution-atoms
                                :stereo-information stereo-information
                                :plugs plugs
                                :residue residue)
                 result)))
    ;; Get all of the plug-factories that can generate in-plugs
    (let ((in-plug-factories (remove-if (lambda (x) (null (in-plug-factory-p x))) plug-factories)))
      (loop for in-plug-factory in in-plug-factories
            for in-plug = (maybe-make-in-plug in-plug-factory complex-plugs)
            for out-plugs = (loop for out-plug-factory in (remove-if (lambda (x) (eq x in-plug-factory)) plug-factories)
                                  collect (maybe-make-out-plug out-plug-factory complex-plugs))
            for plugs = (cons in-plug out-plugs)
            do (push (make-instance 'prepare-topology
                                    :name (if (eq (chem:get-name in-plug) :-default)
                                              name
                                              (intern (format nil "~a~a" (string name) (string (chem:get-name in-plug))) :keyword))
                                    :constitution constitution
                                    :constitution-atoms constitution-atoms
                                    :stereo-information stereo-information
                                    :plugs plugs
                                    :residue residue)
                     result)))
    result))

(defun extract-prepare-topologys (design)
  "Loop through all of the fragments in the chemdraw object in **design** and then loop through all of the residues
in each fragment. Those residues that have the property :topology-name are topology's that we need to extract.
Return a list of prepare-topology objects - one for each residue that we need to extract."
  (let* ((chemdraw-object (design:chemdraw design))
         (extracted-prepare-topologys nil)
         (residue-plug-factory-alist (extract-residue-plug-factory-alist chemdraw-object)))
    ;; For each residue and set of plug-factories construct a constitution-atoms,
    ;; a stereo-information and a constitution
    (loop for (residue . plug-factories) in residue-plug-factory-alist
          for name = (chem:get-name residue)
          for root-atom-name = (identify-root-atom-name residue plug-factories)
          for constitution-atoms = (chem:make-constitution-atoms-from-residue residue)
          for stereo-information = (build-stereoinformation (chem:get-name residue) constitution-atoms root-atom-name)
          for constitution = (chem:make-constitution (chem:get-name residue)
                                                     "no-comment"
                                                     constitution-atoms
                                                     stereo-information
                                                     nil #| This should be a list of unique plugs - but I don't have them until prepare-topologys-from-plug-factories below |#
                                                     nil )
          ;; Collect a list of prepare-topology that we generate from the list of plug-factories
          do (loop for prepare-topology in (prepare-topologys-from-plug-factories name
                                                                                  residue
                                                                                  plug-factories
                                                                                  constitution
                                                                                  constitution-atoms
                                                                                  stereo-information
                                                                                  (design:complex-plugs design))
                   do (push prepare-topology extracted-prepare-topologys))
             ;; For every plug that was generated by the plug-factories - add it to the constitution
          do (loop for plug-factory in plug-factories
                   do (loop for plug in (gather-plugs plug-factory)
                            do (chem:add-plug constitution (chem:get-name plug) plug))))
    extracted-prepare-topologys))

(defun find-in-plug (plugs)
  (loop for p in plugs
        when (chem:get-is-in p)
          do (return-from find-in-plug p)))

(defun find-out-plugs (plugs)
  (loop for p in plugs
        when (not (chem:get-is-in p))
          collect p))



(defun get-all-out-plug-bond0atoms-as-set (plugs)
  (error "This function is messed up - fix it")
  (let ((res nil))
    (loop for p in plugs
          do (let ((atom-name (chem:plug/get-b0 p)))
               (unless (not (contains res atom-name))
                 (error "The atom name[~s] is already in the set[~s]" atom-name res))
               (push atom-name res)))
    res))

(defun get-all-out-plug-bond1atoms-as-set (plugs out-plug-bond0-atoms-set )
  (let ((res nil))
    (loop for p in plugs
          do (when (chem:plug/get-b1 p)
               (let ((atom-name (chem:plug/get-b1 p)))
                 (unless (not (member atom-name res))
                   (error "The atom name[~s] is already in the listset[~s]" atom-name res))
                 (unless (not (member atom-name out-plug-bond0-atoms-set))
                   (error "An atom can't be both a bond0 and bond1 atom!!! Something is wrong with the plugs"))
                 (push atom-name res))))
    res))

(defun get-all-stub-pivot-atoms-as-set (plugs out-plug-bond0-atoms-set )
  (error "There is a problem with ObjectSet")
  #+(or)
  (let ((res (ObjectSet))
        (plug-counter 0))
    (loop for p in plugs
          do (when (chem:has-stub-pivot-atom p)
               (let ((atom-name (chem:plug/get-b1 p)))
                 (unless (not (member atom-name res))
                   (error "The atom name[~s] is already in the set[~s]" atom-name res))
                 (unless (not (member atom-name out-plug-bond0-atoms-set))
                   (error "An atom can't be both a bond0 and bond1 atom!!! Something is wrong with the plugs"))
                 (push atom-name res)))
             (incf plug-counter))
    res))

;;
;; Define the stub atoms in terms of
;; 1) parent older-sibling even-older-sibling
;; or if no even-older-sibling 2) parent older-sibling grand-parent
;; or if no older-sibling      3) parent grand-parent great-grand-parent
;; unless your parent is a OddIcoorAtom then just use
;; parent grand-parent
;;
;; (defun setStubAtoms (buildAtom parent)
;;   (IMPLEMENT_ME))
;; 
;; 
;; (defun addAtomToICoorTableRecursively (icoorTable parent me)
;;   (let ((buildAtom (cond
;;                  ((getProperty me :bond1Atom) (OddICoorAtom (atom-name me)))
;;                  ( true (NormalICoorAtom (atom-name me))))))
;;     (setProperty me :buildAtom buildAtom)
;; ;    (setStubAtoms buildAtom parent)
;;     (appendEntry icoorTable buildAtom)
;;     (for child (getProperty me :children)
;;       (addAtomToIcoorTableRecursively icoorTable me child))))
;; 

(defun one-atom-as-string (me constitution-atoms)
  ;; kin:joint-template-atom-name is no longer available
  #+(or)(let ((target (make-string-output-stream))
        (result (make-string-output-stream)))
    (format target "~a" (kin:joint-template-atom-name me constitution-atoms))
    (format target "/~a" (class-name (class-of me)))
    (format result "~30a: " (get-output-stream-string target))
    (loop for child in (children me)
          do (let ((child-name (kin:joint-template-atom-name child constitution-atoms)))
               (format result "~a/~a" child-name (kin:joint-template-comment child))))
    (get-output-stream-string result)))

(defun dump-build-order-recursively (parent-template cur-template constitution-atoms sout)
  (format sout "~a~%" (one-atom-as-string cur-template constitution-atoms))
  (dolist (child (children cur-template))
    (dump-build-order-recursively cur-template child constitution-atoms sout)))


(defun joint-template-factory (parent-template atom constitution-atoms constitution-name topology-name)
  (let* ((out-plug-atom-prop (chem:matter-get-property-or-default atom :out-plug nil))
         (entity-to-delay-children-for (chem:matter-get-property-or-default atom :entity-to-delay-children-for nil))
         (root-atom-prop (chem:matter-get-property-or-default atom :root-atom nil))
         (atom-index (chem:constitution-atoms/index constitution-atoms (chem:atom-name atom)))
         (atom-name (chem:atom-name (chem:atom-with-id constitution-atoms atom-index)))
         (comment (format nil "~s [~a]" atom-name (chem:matter-get-property-or-default atom :weight -1))) )
    (cond
      (entity-to-delay-children-for
       (when root-atom-prop
         (error "joint-template-factory problem with atom ~a~% - atom has entity-to-delay-children-for ~a~% - and root-atom-prop ~a ~% -  constitutionName[~s] topologyName[~s]"
                atom
                entity-to-delay-children-for
                root-atom-prop
                atom
                constitution-name
                topology-name))
       (let ((checkpoint (if (typep entity-to-delay-children-for 'chem:plug)
                             (make-checkpoint-out-plug-joint
                              :constitution-name constitution-name
                              :topology-name topology-name
                              :out-plug entity-to-delay-children-for)
                             (make-checkpoint-joint
                              :constitution-name constitution-name
                              :topology-name topology-name
                              :atom-name (chem:atom-name entity-to-delay-children-for)))))
         (if out-plug-atom-prop
             (make-delayed-bonded-joint-template :id atom-index
                                                 :name atom-name
                                                 :parent parent-template
                                                 :checkpoint checkpoint
                                                 :comment comment
                                                 :out-plug out-plug-atom-prop )
             (make-delayed-bonded-joint-template :id atom-index
                                                 :name atom-name
                                                 :parent parent-template
                                                 :checkpoint checkpoint
                                                 :comment comment))))
      (root-atom-prop
       (if out-plug-atom-prop
           (make-root-bonded-joint-template :id atom-index
                                            :parent parent-template
                                            :name atom-name
                                            :constitution-name constitution-name
                                            :topology-name topology-name
                                            :in-plug root-atom-prop
                                            :comment comment
                                            :out-plug out-plug-atom-prop)
           (make-root-bonded-joint-template :id atom-index
                                            :parent parent-template
                                            :name atom-name
                                            :constitution-name constitution-name
                                            :topology-name topology-name
                                            :in-plug root-atom-prop
                                            :comment comment)))
      (t (if out-plug-atom-prop
             (make-bonded-joint-template :id atom-index
                                         :parent parent-template
                                         :name atom-name
                                         :comment comment
                                         :out-plug out-plug-atom-prop)
             (make-bonded-joint-template :id atom-index
                                         :parent parent-template
                                         :name atom-name
                                         :comment comment))))))



(defgeneric write-into-joint-tree (joint-template parent-joint atresidue atmolecule-index atresidue-index))

(defmethod write-into-joint-tree ((joint-template t) parent-joint atresidue atmolecule-index atresidue-index)
  (error "write-into-joint-tree - handle joint-template ~a" joint-template))


(defmethod write-into-joint-tree ((joint-template jump-joint-template) parent-joint atresidue atmolecule-index atresidue-index)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-jump-joint atomid atom-name)))
    (kin:put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

(defmethod write-into-joint-tree ((joint-template complex-bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-complex-bonded-joint atomid atom-name))
         (input-stub-joints (input-stub-joints joint-template)))
    (kin:put-joint atresidue joint constitution-atoms-index)
    (let ((input-stub0-template (aref input-stub-joints 0))
          (input-stub1-template (aref input-stub-joints 1)))
      (cond
        ((null input-stub0-template)
         ;; Do nothing
         )
        ((null input-stub1-template)
         (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                (input-stub0 (aref (kin:joints atresidue) input-stub0-index)))
           (kin:set-input-stub-joint1 joint input-stub0)))
        (t
         (let* ((input-stub0-index (constitution-atoms-index input-stub0-template))
                (input-stub0 (aref (kin:joints atresidue) input-stub0-index))
                (input-stub1-index (constitution-atoms-index input-stub1-template))
                (input-stub1 (aref (kin:joints atresidue) input-stub1-index)))
           (kin:set-input-stub-joint1 joint input-stub0)
           (kin:set-input-stub-joint2 joint input-stub1))))
      (when parent-joint (kin:joint/add-child parent-joint joint))
      joint)))

(defmethod write-into-joint-tree ((joint-template bonded-joint-template) parent-joint atresidue atmolecule-index atresidue-index)
  (let* ((constitution-atoms-index (constitution-atoms-index joint-template))
         (atom-name (atom-name joint-template))
         (atomid (list atmolecule-index atresidue-index constitution-atoms-index))
         (joint (kin:make-bonded-joint atomid atom-name)))
    (kin:put-joint atresidue joint constitution-atoms-index)
    (when parent-joint (kin:joint/add-child parent-joint joint))
    joint))

;;
;; Build an AtomTreeTemplate recursively using the properties defined
;; for each atom
;;
(defun build-atom-tree-template-recursively (parent root child-indexes in-plug residue constitution-atoms constitution-name topology-name )
  (let ((root-template (new-joint-template-factory parent root child-indexes in-plug constitution-atoms constitution-name topology-name))
        (children (progn
                    (chem:matter-get-property root :children))))
    (loop for child in children
          for sub-child-index from 0
          when (chem:has-content residue child) 
            do (let ((child-template (build-atom-tree-template-recursively root-template
                                                                           child
                                                                           (cons sub-child-index child-indexes)
                                                                           in-plug
                                                                           residue
                                                                           constitution-atoms
                                                                           constitution-name
                                                                           topology-name )))
                 (add-child root-template child-template)))
    root-template))


(defun build-internal-coordinate-joint-template-tree (prepare-topology)
;;;chemdraw-fragment constitution-atoms plugs constitution-name topology-name)
  (with-accessors ((name name) (residue residue) (plugs plugs) (constitution constitution))
      prepare-topology
    (let* ((in-plug (let ((in-plug (find-in-plug plugs)))
                      (if in-plug
                          in-plug
                          (error "There has to be an in-plug in topology ~a" name))))
           (outplugs (find-out-plugs plugs))
           (constitution-atoms (chem:constitution/get-constitution-atoms constitution))
           (root-atom-name (chem:root-atom-name in-plug))
           (root-atom (chem:atom-with-name residue root-atom-name))
           (spanning-loop (chem:spanning-loop/make root-atom))
           (all-spanning-atoms (chem:all-atoms spanning-loop)))
      (chem:set-property root-atom :root-atom in-plug)
      (loop for atom in all-spanning-atoms
            do (chem:set-property atom :out-plug nil)
               (chem:set-property atom :bond1-atom nil)
               (chem:set-property atom :weight 1)
               (chem:set-property atom :children nil))
      ;;
      ;; Now iterate over the plugs and assign atom properties that will 
      ;; define the build order and identify for OutPlugs
      ;; 1) :highestPriority atoms (atoms that bond0 atoms depend on)
      ;; 2) :highPriority atoms (bond0 atoms)
      ;; 3) atoms that need to delay building their children until another atom is defined
      ;;     identify by setting :entityToDelayChildrenFor XXXX where
      ;;           XXXX == :nextResidueBond1
      ;;        or XXXX == the stubPivotAtom that you need to build your children
      (loop for plug in plugs
            do (cond
                 ((typep plug 'chem:out-plug)
                  (let ((bond0-atom (chem:atom-with-name residue (chem:plug/get-b0 plug)))
                        (bond1-atom (when (chem:plug/get-b1 plug)
                                      (chem:atom-with-name residue (chem:plug/get-b1 plug))))
                        (stub-pivot-atom (when (chem:has-stub-pivot-atom plug)
                                           (chem:atom-with-name residue (chem:get-stub-pivot-atom plug)))))
                    (unless (not (and (or bond1-atom stub-pivot-atom) (eq bond1-atom stub-pivot-atom)))
                      (error "The same plug cannot define a bond1Atom[~a] and a stubPivotAtom[~a] at the same time" bond1-atom stub-pivot-atom))
                    (chem:set-property bond0-atom :out-plug plug)
                    (chem:set-property bond0-atom :high-priority t)
                    (cond
                      (bond1-atom
                       (chem:set-property bond1-atom :bond1-atom t)
                       (chem:set-property bond1-atom :highest-priority t)
                       (chem:set-property bond1-atom :build-before bond0-atom)
                       (chem:set-property bond1-atom :entity-to-delay-children-for plug ))
                      (stub-pivot-atom
                       (let ((bond0-depends-on-atom (chem:matter-get-property-or-default bond0-atom :build-before nil)))
                         (if (eq bond0-depends-on-atom stub-pivot-atom)
                             (progn
                               ;; (% "Setting atom[%s] :entity-to-delay-children-for to-> %s" (description bond0-atom) stub-pivot-atom)
                               (chem:set-property bond0-atom :entity-to-delay-children-for stub-pivot-atom))
                             (progn
                               ;; (% "Setting stubPivotAtom[%s] :highestPriority=t and :buildBefore to-> %s" (description stubPivotAtom) (description bond0-atom))
                               (chem:set-property stub-pivot-atom :build-before bond0-atom)
                               (chem:set-property stub-pivot-atom :highest-priority t))))))))))
      ;; Now iterate over the atoms from farthest to the nearest to the root
      ;; and add yourself to your backSpanning atoms :children
      (dolist (atom all-spanning-atoms)
        (let ((parent-atom (chem:get-back-span spanning-loop atom)))
          (when parent-atom
            (let ((parent-children (chem:matter-get-property parent-atom :children)))
              (setf parent-children (cons atom parent-children))
              (chem:set-property parent-atom :children parent-children)))))
      ;; Now iterate over the atoms from farthest to the nearest and assign a branch weight
      ;; to each atom
      ;; 1) Get the largest weight of this atoms children
      ;; 2) If the current atom is a bond0Atom then add 1000 to its weight
      ;; 3) otherwise if the current atom is a bond1Atom then add 10000 to its weight
      ;; 4) otherwise add 1 to its weight
      ;;   This is to try and ensure that branches that involve Plugs that make two bonds are
      ;;   built first and that the bond1Atom of these plugs gets built before the bond0Atom.
      ;;   This is so that the bond1Atom will define the dihedral to the root atom of the next residue
      ;;   that is connected through the bond0Atom.
      ;;   This will help to isolate residues from each other in terms of internal dihedral coordinates
      ;;   so that bis-amino acids are built from the amino acid atoms of the previous residue rather
      ;;   than building from the inside of the previous bis-amino acid. 
      ;;   Here is a concrete example - If we had a pro4(SS) residue X modeled on another pro4(SS) residue Y and
      ;;   the dihedral for X's root were defined relative to the CB of Y then if we built X onto a pro4(SR) 
      ;;   residue Z then that dihedral would be wrong.  If instead we define X's root atom dihedral in terms
      ;;   of the CGC of Y then it won't matter what stereoisomer of Z we put on.
      ;;
      (dolist (atom all-spanning-atoms)
        (let ((atom-weight (cond
                             ((chem:matter-get-property-or-default atom :highest-priority nil) 10000)
                             ((chem:matter-get-property-or-default atom :high-priority nil)     1000)
                             (t                                                                    1))))
          (let ((max-child-weight 0))
            (loop for child in (chem:matter-get-property atom :children)
                  do (setf max-child-weight (max max-child-weight (chem:matter-get-property child :weight))))
            (setf atom-weight (+ atom-weight max-child-weight))
            (chem:set-property atom :weight atom-weight))))
      ;;
      ;; Now sort the children of each atom in descending order of their weight
      ;; and if there are ties (like methylene or methyl hydrogens) sort them
      ;; in increasing alphabetical order by their name
      ;;
      (dolist (atom all-spanning-atoms)
        (let* ((children (chem:matter-get-property atom :children))
               (sorted-children (sort children (lambda (a b)
                                                 (let ((wa (chem:matter-get-property a :weight))
                                                       (wb (chem:matter-get-property b :weight)))
                                                   (cond
                                                     ((> wa wb) t)
                                                     ((< wa wb) nil)
                                                     (t (string< (string (chem:atom-name a))
                                                                 (string (chem:atom-name b))))))))))
          (chem:set-property atom :children sorted-children)))
      ;;
      ;; Ok, now every atom has :children defined and they are sorted in order of how they
      ;; should be built - now we are ready to build the InternalCoordinateTable for the residue
      ;;
      (let ((tree-template (build-atom-tree-template-recursively nil
                                                                 root-atom
                                                                 (list 0)
                                                                 in-plug
                                                                 residue
                                                                 constitution-atoms
                                                                 name
                                                                 name)))
        #+(or)
        (let ((build-order (with-output-to-string (sout)
                                   (dump-build-order-recursively nil tree-template constitution-atoms sout))))
                (format t "Build order:~%~a~%" build-order))
        tree-template))))


(defparameter *max-index* 0)
#+(or)
(defun walk-joint-templates (node)
  (setf *max-index* (max *max-index* (kin:id node)))
  (loop for child in (kin:children node)
        do (walk-joint-templates child)))

#+(or)
(defun save-joint-templates (node vec)
  (format t "Setting ~a -> ~a~%" (kin:id node) (kin:joint/name node))
  (setf (aref vec (kin:id node)) (kin:joint/name node))
  (loop for child in (kin:children node)
        do (save-joint-templates child vec)))

#+(or)
(defun sort-joint-templates (tree)
  (let ((*max-index* 0))
    (walk-joint-templates tree)
    (let ((vec (make-array (1+ *max-index*))))
      (save-joint-templates tree vec)
      vec)))

      
  
