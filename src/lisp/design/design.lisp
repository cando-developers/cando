
(in-package :design)

(defvar *complex-plugs* (make-hash-table))
(defvar *caps*)

(defclass complex-plug ()
  ((name :initarg :name :accessor name)
   (smarts :initarg :smarts :accessor smarts)
   (compiled-smarts :initarg :compiled-smarts :accessor compiled-smarts)))

(defun complex-plug (name smarts)
  "Define a complex plug that uses a SMARTS string to recognize bonding atoms"
  (setf (gethash name *complex-plugs*) (make-instance 'complex-plug
                                                      :name name
                                                      :smarts smarts
                                                      :compiled-smarts (chem:compile-smarts smarts))))

(defun complex-plug-or-nil (name complex-plugs)
  (gethash name complex-plugs))

(defclass design ()
  ((chemdraw :initarg :chemdraw :accessor chemdraw)
   (cap-name-map :initarg :cap-name-map :accessor cap-name-map)
   (complex-plugs :initarg :complex-plugs :accessor complex-plugs)))

(defun configure (chemdraw)
  (error "You must provide a function (defun design:configure (chemdraw) ...) that returns a design:design object"))

(defun make-design (chemdraw plugs)
  (let ((complex-plugs (make-hash-table))
        (caps nil))
    (loop for plug-spec in plugs
          do (destructuring-bind (plug-name &key smarts cap)
                 plug-spec
               (when smarts
                 (setf (gethash plug-name complex-plugs) (make-instance 'complex-plug
                                                                        :name plug-name
                                                                        :smarts smarts
                                                                        :compiled-smarts (chem:compile-smarts smarts))))
               (when cap
                 (push (cons plug-name cap) caps))))
    (make-instance 'design
                   :chemdraw chemdraw
                   :cap-name-map (alexandria:alist-hash-table caps)
                   :complex-plugs complex-plugs)))

(defclass grammar ()
  ((topologys :initarg :topologys :accessor topologys)
   (cap-name-map :initarg :cap-name-map :accessor cap-name-map)
   (trainers :initarg :trainers :accessor trainers)))


#|
(with-design (design ((:+default :smarts "[C](-[N:0])-[C:1]=O" :cap :pro)
                      (:-default :smarts "[C](-[N:0])-[C:1]=O" :cap :phe)
                      (:-side :cap :pro4)
                      (:+side :cap :4py)))
  (design.load:setup design))
|#

(cando:make-class-save-load design)

(defun save-design (design file-name)
  (let ((pn (pathname file-name)))
    (ensure-directories-exist pn)
    (cando:save-cando design file-name)))

(defun load-design (file-name)
  (let ((design (cando:load-cando file-name)))
    (loop for topology in (topologys design)
          do (cando:register-topology topology))
    design))


(defun deftop (residue in-plug-info &rest out-plugs-info)
  (let ((topology (cando:make-simple-topology-from-residue residue)))
    (when in-plug-info
      (destructuring-bind (in-plug-name &optional (atom-name in-plug-name) (bond-order :single-bond))
          (if (consp in-plug-info)
              in-plug-info
              (list in-plug-info))
        (chem:add-plug topology in-plug-name (chem:make-in-plug in-plug-name nil atom-name bond-order))))
    (loop for out-plug-info in out-plugs-info
          collect (destructuring-bind (out-plug-name &optional (atom-name out-plug-name) (bond-order :single-bond))
                      (if (consp out-plug-info)
                          out-plug-info
                          (list out-plug-info))
                    (chem:add-plug topology out-plug-name (chem:make-out-plug out-plug-name nil nil atom-name bond-order))))
    (cando:register-topology topology)
    topology))

(defun lookup-topology (name)
  "Lookup a named object and if it is a chem:topology return it.
Otherwise signal an error."
  (cando:lookup-topology name))

(defun lookup-maybe-part (name)
  "Lookup a named object and return it.
This is for looking up parts but if the thing returned is not a part then return nil."
  (let ((maybe-part (cando:lookup-entity name nil)))
    (if (typep maybe-part 'part)
        maybe-part
        nil)))


(defvar *alpha-carbon* (chem:compile-smarts "[C&H1](-N)(-C=O)"))

(defun find-ends (molecule)
  (let (ends)
    (cando:do-residues (residue molecule)
      (let (start end)
        (cando:do-atoms (atom residue)
          (when (eq (chem:get-name atom) :start) (setf start atom))
          (when (eq (chem:get-name atom) :end) (setf end atom)))
        (when (and start end) (push (list start end) ends))))
    ends))

(defun bonds-between (start end)
  (let* ((spanning-tree (chem:make-spanning-loop start)))
    (loop for ok = (chem:advance-loop-and-process spanning-tree)
          for atom = (chem:get-atom spanning-tree)
          until (eq atom end))
    (chem:get-back-count spanning-tree end)))

(defun add-guidance-restraints (molecule)
  (chem:clear-restraints molecule)
  (let ((ends (find-ends molecule)))
    (loop for (start end) in ends
          for bonds-between = (bonds-between start end)
          do (chem:add-restraint molecule (core:make-cxx-object 'chem:restraint-distance
                                                                :a start
                                                                :b end
                                                                :r0 (float bonds-between 1d0)
                                                                :k 1000.0)))
    (format *debug-io* "Second loop ends: ~a~%" ends)
    (loop for cur on ends
          for (start-from end-from) = (car cur)
          do (format *debug-io* "cur: ~a~%" cur)
          do (format *debug-io* "start-from: ~a   end-from: ~a~%" start-from end-from)
             (format *debug-io* "    (cdr cur) -> ~a~%" (cdr cur))
          do (loop for (start-to end-to) in (cdr cur)
                   for start-bonds-between = (bonds-between start-from start-to)
                   for end-bonds-between = (bonds-between end-from end-to)
                   do (chem:add-restraint molecule (core:make-cxx-object 'chem:restraint-distance
                                                                         :a start-from
                                                                         :b start-to
                                                                         :r0 (float start-bonds-between 1d0)
                                                                         :k 1000.0))
                      (chem:add-restraint molecule (core:make-cxx-object 'chem:restraint-distance
                                                                         :a end-from
                                                                         :b end-to
                                                                         :r0 (float end-bonds-between 1d0)
                                                                         :k 1000.0))))))

(defun add-spiroligomer-guidance-restraints (molecule root-residue-name)
  "Add a stiff restraint across the ends of the spiroligomer
- from the first alpha carbon to the last"
  (let* ((root-residue (chem:content-with-name molecule root-residue-name))
         (start-atom (chem:content-with-name root-residue :+default))
         (spanning-tree (chem:make-spanning-loop start-atom))
         (alpha-carbons (loop for ok = (chem:advance-loop-and-process spanning-tree)
                              for atom = (chem:get-atom spanning-tree)
                              until (null ok)
                              when (and ok (chem:matches *alpha-carbon* atom))
                                collect atom)))
    (chem:add-restraint molecule (core:make-cxx-object 'chem:restraint-distance
                                                       :a (first alpha-carbons)
                                                       :b (car (last alpha-carbons))
                                                       :r0 (* 6.0 (1- (length alpha-carbons)))
                                                       :k 1000.0))))


(defun find-starting-geometry (aggregate &rest arguments)
  (unless arguments
    (setf arguments '(:cg-tolerance 5.0 :max-tn-steps 100 :tn-tolerance 0.5)))
  (let ((mol (chem:content-at aggregate 0)))
    #+(or)(let ((root-residue-name (chem:get-name (chem:|Oligomer_O::rootMonomer| oligomer))))
            (chem:clear-restraints mol)
            (add-spiroligomer-guidance-restraints mol root-residue-name))
    (add-guidance-restraints mol)
    (apply #'energy:minimize aggregate arguments)
    (chem:clear-restraints mol)
    (apply #'energy:minimize aggregate arguments)))

(defclass subtree ()
  ((name :initarg :name :accessor name)
   (code :initarg :code :accessor code)))

(defun make-subtree (name code)
  (make-instance 'subtree :name name :code code))

(defclass part ()
  ((name :initform :name :accessor name)
   (tree :initform :tree :accessor tree)))

(defun make-part (name tree)
  (make-instance 'part :name name :tree tree))

(defun define-part (name tree)
  "Create a part and register it by name with leap.core:register-variable."
  (let ((part (make-part name tree)))
    (cando:register-entity name part)))

(defun my-add-monomers (oligomer names)
  (let ((monomer (chem:make-monomer names)))
    (chem:|Oligomer_O::addMonomer| oligomer monomer)
    (values (list monomer) nil)))

(defun translate-part (oligomer names labels)
  (let (maybe-part)
    (cond
      ((and (symbolp names)
            (setf maybe-part (lookup-maybe-part names)))
       (values (interpret-subtree oligomer maybe-part labels) nil))
      ((symbolp names)
       (my-add-monomers oligomer (list names)))
      ((consp names)
       (if (and (symbolp (first names)) (string= (first names) :cycle))
           (values (gethash (second names) labels) (list (third names) (fourth names)))
           (my-add-monomers oligomer names)))
      (t (error "Handle ~a" names)))))
            
(defun interpret-part (oligomer part-info labels)
  (destructuring-bind (names &key label)
      (cond
        ((and (consp part-info) (symbolp (first part-info)) (string= (first part-info) :cycle))
         (list part-info))
        ((consp part-info)
         part-info)
        ((symbolp part-info)
         (list part-info)))
    (multiple-value-bind (new-parts ringp)
        (translate-part oligomer names labels)
      (when label
        (loop for new-part in new-parts
              do (push new-part (gethash label labels))))
      (values new-parts ringp))))

(defun do-coupling (oligomer coupling-name ring-info previous-parts next-parts)
  (unless next-parts
    (error "There are no next parts to add to ~a" previous-parts))
  (let ((in-plug-name (chem:in-plug-name coupling-name))
        (out-plug-name (chem:out-plug-name coupling-name)))
    (let ((previous-monomer (loop for part in previous-parts
                                  for name = (chem:current-stereoisomer-name part)
                                  for topology = (chem:current-topology part) ; lookup-topology name)
#|                                  do (progn
                                       (format *debug-io* "do-coupling       part -> ~s~%" part)
                                       (format *debug-io* "  (chem:current-stereoisomer-name part) -> ~s~%" (chem:current-stereoisomer-name part))
                                       (format *debug-io* "      (type-of part) -> ~s~%" (type-of part))
                                       (format *debug-io* "  name -> ~s~%" name)
                                       (format *debug-io* "  (type-of name) -> ~s~%" (type-of name))
                                       (format *debug-io* "              topology -> ~s~%" topology)
                                       (format *debug-io* "         out-plug-name -> ~s~%" out-plug-name))
|#
                                  when (chem:has-plug-named topology out-plug-name)
                                    collect part))
          (next-monomer (loop for part in next-parts
                              for topology = (chem:current-topology part) ; lookup-topology (chem:monomer-name part))
                              when (chem:has-plug-named topology in-plug-name)
                                collect part)))
      (when (= (length previous-monomer) 0)
        (error "There is no monomer found with the out-plug-name ~a in the parts ~a" out-plug-name previous-parts))
      (unless (= (length previous-monomer) 1)
        (error "There is more than one monomer(~a) with the out-plug-name ~a" previous-monomer out-plug-name))
      (unless (= (length next-monomer) 1)
        (format t "(length next-monomer) -> ~a~%" (length next-monomer))
        (if (< (length next-monomer) 1)
            (error "There is no monomer with the in-plug-name ~a to connect to ~a" in-plug-name previous-parts)
            (error "There is more than one monomer(~a) with the in-plug-name ~a" next-monomer in-plug-name)))
      (if (null ring-info)
          (chem:|Oligomer_O::couple| oligomer
                       (first previous-monomer)
                       coupling-name
                       (first next-monomer))
          (let ((plug1name (or (first ring-info) out-plug-name))
                (plug2name (or (second ring-info) in-plug-name)))
            (chem:ring-couple-with-plug-names oligomer
                                              (first previous-monomer)
                                              plug1name
                                              (first next-monomer)
                                              plug2name))))))
  

(defun interpret-subtree (oligomer subtree labels)
  (let* ((root-monomer-info (pop subtree))
         (previous-parts (interpret-part oligomer root-monomer-info labels))
         (accumulated-parts previous-parts))
    (loop
      (when (null subtree) (return accumulated-parts))
      (let ((info (pop subtree)))
        (if (consp info)
            (let* ((coupling (first info)) ; interpret a branch
                   (node-info (cdr info)))
              (multiple-value-bind (new-parts ringp)
                  (interpret-part oligomer node-info labels)
                (setf accumulated-parts (append accumulated-parts new-parts))
                (do-coupling oligomer coupling ringp previous-parts new-parts)))
            (let* ((coupling info)      ; Interpret a chain
                   (node-info (pop subtree)))
              (multiple-value-bind (new-parts ringp)
                  (interpret-part oligomer node-info labels)
                (setf accumulated-parts (append accumulated-parts new-parts))
                (do-coupling oligomer coupling ringp previous-parts new-parts)
                (setf previous-parts new-parts))))))))

(defun make-oligomer (tree)
  "Make an oligomer from a description in the **tree**.
The tree is a nested list of lists that look like 
(component coupling component coupling component ... ).
It starts with the name of a topology or a design:part.
The component can be a name of a single topology, it can be a design:part,
or it can 

Examples:
(make-oligomer '(:ccap :default :ala :default :ala :default :ser :default :ncap))
(make-oligomer '((:lego3 :label :first)
                 :default :lego3
                 :default :lego3
                 :default (cycle :first)))
"
  (let* ((oligomer (core:make-cxx-object 'chem:oligomer))
         (labels (make-hash-table)))
    (interpret-subtree oligomer tree labels)
    oligomer))

(defun classify-topologys (list-of-topologys)
  (let ((origins (make-hash-table))
        (body (make-hash-table))
        (caps (make-hash-table)))
    (loop for topology in list-of-topologys
          if (chem:out-plugs-as-list topology)
            do (let ((in-plug (chem:get-in-plug topology)))
                 (if (typep in-plug 'chem:origin-plug)
                     (progn
                       (format *debug-io* "origin topology: ~s~%" topology)
                       (setf (gethash (chem:get-name topology) origins) topology))
                     (progn
                       (format *debug-io* "  body topology: ~s~%" topology)
                       (setf (gethash (chem:get-name topology) body) topology))))
          else
            do (progn
                 (format *debug-io* "   cap topology: ~s~%" topology)
                 (setf (gethash (chem:get-name topology) caps) topology)))
    (values (alexandria:hash-table-values origins)
            (alexandria:hash-table-values body)
            (alexandria:hash-table-values caps))))

(defun cluster-topologys-using-out-plugs (list-of-topologys)
  "The input is a list-of-topologys with the same in plug but maybe different out-plugs.
Subdivide the list into a list of lists that contain Topologys with different numbers and kinds
of out-plugs."
  (let ((clusters (make-hash-table :test #'equal)))
    (loop for topology in list-of-topologys
          for out-plug-names = (mapcar #'chem:get-name (chem:out-plugs-as-list topology))
          for sorted-out-plug-names = (sort out-plug-names #'string<)
          do (push topology (gethash sorted-out-plug-names clusters)))
    (alexandria:hash-table-values clusters)))

(defun ensure-one-unique-out-plug-name (other-monomer in-plug-name)
  "Ensure that there is one topology with one out-plug with a name that corresponds to in-plug-name on
   a monomer that this one out-plug will be coupled through"
  (let* ((topology (chem:current-topology other-monomer))
         (out-plug-names (chem:all-out-plug-names-that-match-in-plug-name topology in-plug-name)))
    (format *debug-io* "out-plug-names -> ~s~%" out-plug-names)
    (case (length out-plug-names)
      (0
       (error "There are no out-plugs in ~s that match the in-plug-name ~s - existing plugs: ~s" other-monomer in-plug-name (chem:plugs-as-list topology)))
      (1
       (return-from ensure-one-unique-out-plug-name (first out-plug-names)))
      (otherwise
       (error "There are multiple out-plugs ~s in ~s that match the in-plug-name ~s"
              out-plug-names other-monomer in-plug-name)))))

(defun find-unsatisfied-monomer-plug-pairs (oligomer)
  (let* ((monomers (chem:|Oligomer_O::monomersAsList| oligomer))
         (couplings (chem:|Oligomer_O::couplingsAsList| oligomer))
         (monomer-plug-work-list (make-hash-table :test #'equal)))
    (loop for monomer in monomers
          for topology = (chem:current-topology monomer)
          for plugs = (chem:plugs-as-list topology)
          do (loop for plug in plugs
                   for plug-name = (chem:get-name plug)
                   when (or (and (null (typep plug 'chem:origin-plug)) (chem:get-is-in plug))
                            (chem:get-is-out plug))
                     do (setf (gethash (cons monomer plug-name) monomer-plug-work-list) t)))
    (loop for coupling in couplings
          when (typep coupling 'chem:directional-coupling)
            do (setf (gethash (cons (chem:get-source-monomer coupling)
                                    (chem:get-source-plug-name coupling))
                              monomer-plug-work-list) nil)
               (setf (gethash (cons (chem:get-target-monomer coupling)
                                    (chem:get-target-plug-name coupling))
                              monomer-plug-work-list) nil))
    (let (result)
      (maphash (lambda (monomer-plug-name dangling-p)
                 (when dangling-p
                   (push monomer-plug-name result)))
               monomer-plug-work-list)
      result)))

(defun one-round-extend-oligomer-with-caps (oligomer cap-name-map topology-map &key (verbose t))
  (let ((monomers-plugs (find-unsatisfied-monomer-plug-pairs oligomer)))
    (when verbose (format *debug-io* "monomers-plugs: ~%~a~%" monomers-plugs))
    (if monomers-plugs
        (progn
          (when verbose (format *debug-io* "Extending oligomer with caps~%"))
          (loop for (monomer . plug-name) in monomers-plugs
                for plug = (chem:get-plug-named monomer plug-name)
                for cap = (gethash (chem:get-name plug) cap-name-map)
                do (when verbose (format *debug-io* "Extending monomer ~s  plug ~s~%" monomer plug))
                   (unless cap
                     (error "Could not find cap for monomer ~s plug ~s in cap-map"
                            monomer plug-name))
                   (multiple-value-bind (other-topology found)
                       (gethash cap topology-map)
                     (when verbose (format *debug-io* "other-topology ~s~%" other-topology))
                     (unless found
                       (error "Could not find topology for cap ~s~%" cap))
                     (let ((other-monomer (chem:make-monomer (list (chem:get-name other-topology)))))
                       (when verbose (format *debug-io* "Adding new monomer ~s~%" other-monomer))
                       (chem:|Oligomer_O::addMonomer| oligomer other-monomer)
                       (etypecase plug
                         (chem:out-plug
                          (let ((other-monomer-in-plug-name (chem:in-plug-name (chem:coupling-name plug-name))))
                            (unless (chem:get-plug-named other-monomer other-monomer-in-plug-name)
                              (error "While trying to couple monomer ~s with out-plug named ~s we could not find a corresponding in-plug named ~s in the cap monomer ~s"
                                     monomer plug-name other-monomer-in-plug-name other-monomer))
                            (when verbose (format *debug-io* "Coupling to out-coupling ~s ~s ~s ~s~%" monomer plug-name other-monomer other-monomer-in-plug-name))
                            (chem:|Oligomer_O::couple| oligomer monomer plug-name other-monomer other-monomer-in-plug-name)))
                         (chem:in-plug
                          (when verbose (format *debug-io* "Coupling monomer ~s with in coupling~%" other-monomer))
                          (let ((other-monomer-out-plug-name (ensure-one-unique-out-plug-name other-monomer plug-name)))
                            (when verbose (format *debug-io* "Coupling to in-coupling ~s ~s ~s ~s~%" other-monomer other-monomer-out-plug-name monomer plug-name))
                            (chem:|Oligomer_O::couple| oligomer other-monomer other-monomer-out-plug-name monomer plug-name)))
                         (chem:origin-plug
                          #| do nothing |#)))))
          t)
        nil)))

#+(or)
(defun build-superposable-conformation-collection (conformation monomer-id)
  "Return (values superposable-conformation-collection aggregate)."
  (let* ((fold-tree (kin:get-fold-tree conformation))
         (monomer-node (kin:lookup-monomer-id fold-tree monomer-id))
         (joint-ht (make-hash-table))
         joints)
    (kin:walk-joints monomer-node (lambda (index joint)
                                    (declare (ignore index))
                                    (setf (gethash joint joint-ht) t)
                                    (push joint joints)))
    (let ((root-joint nil))
      (kin:walk-joints monomer-node
                       (lambda (index joint)
                         (declare (ignore index))
                         (when (null (gethash (kin:get-parent joint) joint-ht))
                           (setf root-joint joint))))
      (when (typep root-joint 'kin:root-bonded-joint)
        (let* ((ancestor1 (kin:get-parent root-joint))
               (ancestor2 (kin:get-parent ancestor1))
               (ancestor3 (kin:get-parent ancestor2)))
          (push ancestor1 joints)
          (push ancestor2 joints)
          (push ancestor3 joints)))
      ;; joints contains the joints that we are interested in superposing
      ;;  now convert that to 
      joints)))


(defun build-one-training-oligomer (focus-topologys-in-list cap-name-map topology-map &key test)
  "Build a training oligomer with a list of topologys in a focus-residue and then repeatedly
add cap monomers until no more cap monomers are needed."
  (let ((oligomer (core:make-cxx-object 'chem:oligomer))
        (focus-residue (chem:make-monomer (mapcar #'chem:get-name focus-topologys-in-list))))
    (format *debug-io* "Focus monomer: ~a~%" focus-residue)
    (chem:|Oligomer_O::addMonomer| oligomer focus-residue)
    ;; Now repeatedly cap the focus monomer until it's finished.
    (unless test
      (loop for step from 0
            while (one-round-extend-oligomer-with-caps oligomer cap-name-map topology-map)
            when (> step 10)
              do (error "Too many rounds of extending oligomer")))
    oligomer))

(defun build-training-oligomers (design)
  "Create oligomers that describe training molecules"
  (let* ((list-of-topologys (topologys design))
         (cap-name-map (cap-name-map design))
         (topology-map (let ((ht (make-hash-table)))
                         (loop for top in list-of-topologys 
                               do (setf (gethash (chem:get-name top) ht) top))
                         ht)))
    (multiple-value-bind (origins body caps)
        (classify-topologys list-of-topologys)
      ;; Build training oligomers for origins
      ;; First make the focus monomer and add it to the oligomer
      (let ((clustered-origins (cluster-topologys-using-out-plugs origins))
            (clustered-bodys (cluster-topologys-using-out-plugs body)))
        (format *debug-io* "clustered-origins -> ~s~%" clustered-origins)
        (format *debug-io* "clustered-bodys -> ~s~%" clustered-bodys)
        (let ((all-oligomers (append
                              (loop for focus-topologys in (append clustered-origins clustered-bodys)
                                    collect (build-one-training-oligomer focus-topologys
                                                                         cap-name-map
                                                                         topology-map))
                              (loop for focus-cap-topology in caps
                                    collect (build-one-training-oligomer (list focus-cap-topology)
                                                                         cap-name-map
                                                                         topology-map)))))
          ;; Eliminate oligomers that cover the same sequence space
          (let ((unique-oligomers (make-hash-table :test 'equal)))
            (loop for oligomer in all-oligomers
                  for canonical-sequence = (canonical-sequence oligomer)
                  do (setf (gethash canonical-sequence unique-oligomers) oligomer))
            (alexandria:hash-table-values unique-oligomers)))))))

#+(or)
(defun monomer-node-context (monomer-node)
  (let* ((parent-node (kin:parent monomer-node))
         (parent (if parent-node
                     (kin:stereoisomer-name parent-node)
                     nil))
         (coupling (chem:coupling-name (kin:parent-plug-name monomer-node))))
    (declare (ignore parent))
    (list coupling (kin:stereoisomer-name monomer-node))))

#+(or)
(defun get-conformation (monomer-node conformations)
  (let* ((context (monomer-node-context monomer-node))
         (entry (gethash context conformations))
         (conformation-index (kin:conformation-index monomer-node))
         (internals (elt entry conformation-index)))
    (format *debug-io* "monomer-node ~s context: ~s ~s~%" monomer-node context internals)
    internals))

#+(or)
(defun set-conformation (monomer-node coordinates)
  (let ((internals (get-conformation monomer-node coordinates)))
    (loop for index = 0 then (+ index 5)
          until (>= index (length internals))
          do (let* ((atom-index (elt internals index))
                    (atom-name (elt internals (+ index 1)))
                    (dist (elt internals (+ index 2)))
                    (theta (elt internals (+ index 3)))
                    (phi (elt internals (+ index 4)))
                    (joint (kin:joint-at monomer-node atom-index)))
               (format *debug-io* "joint: ~a  ~a  name: ~a index: ~a dist: ~a ~a ~a~%" joint (kin:atom-id joint)
                       atom-name atom-index dist theta phi)
               (when (typep joint 'kin:bonded-joint)
                 (kin:set-distance joint dist)
                 (kin:set-theta joint theta)
                 (kin:set-phi joint phi))))))

#|

(defparameter *build-trainer*
  (format nil "~{~s~%~}"
          '((defparameter *design* (load ":%DESIGN-INPUT-FILE%"))
            (def
|#
