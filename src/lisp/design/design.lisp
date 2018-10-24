
(in-package :design)

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
    (cando:register-topology (chem:get-name topology) topology)
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

(defun recursively-build-molecule (prev-monomer
                                   prev-topology
                                   prev-residue
                                   monomer-out-couplings
                                   molecule
                                   monomers-to-residues)
  (loop for out-coupling in (gethash prev-monomer monomer-out-couplings)
        for coupling-name = (chem:get-name out-coupling)
        for next-monomer = (chem:get-out-monomer out-coupling)
        for next-topology = (let ((next-top (chem:current-topology next-monomer)))
                              (unless (typep next-top 'chem:topology)
                                (error "Unexpected object - expected a chem:topology - got ~s" next-top))
                              next-top)
        for next-residue = (progn
                             (format t "About to build residue for next-topology ~a~%" next-topology)
                             (chem:build-residue next-topology))
        do (setf (gethash next-monomer monomers-to-residues) next-residue)
           (chem:add-matter molecule next-residue)
           (chem:connect-residues prev-topology 
                                  prev-residue
                                  (chem:get-out-plug-name out-coupling)
                                  next-topology
                                  next-residue
                                  (chem:get-in-plug-name out-coupling))
           (recursively-build-molecule next-monomer
                                       next-topology
                                       next-residue
                                       monomer-out-couplings
                                       molecule
                                       monomers-to-residues)))

(defun build-molecule (oligomer)
  (let ((root-monomer (chem:root-monomer oligomer))
        (monomer-out-couplings (make-hash-table))
        (monomers-to-residues (make-hash-table))
        (ring-couplings nil))
    (format t "Gathering couplings~%")
    (loop for coupling in (chem:couplings-as-list oligomer)
          if (typep coupling 'chem:directional-coupling)
            do (push coupling (gethash (chem:get-in-monomer coupling) monomer-out-couplings))
          else
            do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
    (format t "Making molecule~%")
    (let* ((molecule (chem:make-molecule :mol))
           (root-topology (chem:current-topology root-monomer))
           (stereoisomer-name (chem:current-stereoisomer-name root-monomer))
           (root-residue (progn
                           (format t "About to build first residue ~a~%" root-topology)
                           (chem:build-residue-for-monomer-name root-topology stereoisomer-name))))
      (unless (typep root-topology 'chem:topology)
        (error "Unexpected object - expected a chem:topology - got ~s" root-topology))
      (setf (gethash root-monomer monomers-to-residues) root-residue)
      (chem:add-matter molecule root-residue)
      (recursively-build-molecule root-monomer
                                  root-topology
                                  root-residue
                                  monomer-out-couplings
                                  molecule
                                  monomers-to-residues)
      ;; Now close the rings
      (loop for ring-coupling in ring-couplings
            for monomer1 = (chem:get-monomer1 ring-coupling)
            for topology1 = (chem:current-topology monomer1)
            for residue1 = (gethash monomer1 monomers-to-residues)
            for monomer2 = (chem:get-monomer2 ring-coupling)
            for topology2 = (chem:current-topology monomer2)
            for residue2 = (gethash monomer2 monomers-to-residues)
            for plug1name = (chem:get-plug1 ring-coupling)
            for plug2name = (chem:get-plug2 ring-coupling)
            do
            #+(or)(progn
                    (format t "ring-coupling: ~a~%" ring-coupling)
                    (format t "    topology1: ~a~%" topology1)
                    (format t "     residue1: ~a~%" residue1)
                    (format t "    plug1name: ~a~%" plug1name)
                    (format t "    topology2: ~a~%" topology2)
                    (format t "     residue2: ~a~%" residue2)
                    (format t "    plug2name: ~a~%" plug2name))
                  (chem:connect-residues topology1
                                         residue1
                                         plug1name
                                         topology2
                                         residue2
                                         plug2name))
      molecule)))

(defvar *alpha-carbon* (smarts:make-chem-info :smarts "[C&H1](-N)(-C=O)"))

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
  (let* ((spanning-tree (chem:make-spanning-loop start))
         (bond-count 0))
    (loop for ok = (chem:advance-loop-and-process spanning-tree)
          for atom = (chem:get-atom spanning-tree)
          until (eq atom end))
    (chem:get-back-count end)))

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
    (format t "Second loop ends: ~a~%" ends)
    (loop for cur on ends
          for (start-from end-from) = (car cur)
          do (format t "cur: ~a~%" cur)
          do (format t "start-from: ~a   end-from: ~a~%" start-from end-from)
             (format t "    (cdr cur) -> ~a~%" (cdr cur))
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

(defun build-all-molecules (oligomer &optional (number-of-sequences (chem:number-of-sequences oligomer)))
  "Return a list of all built molecules for the oligomer.
The number of entries can be limited by passing a lower value for the optional argument **number-of-sequences**
than the (chem:number-of-sequences oligomer)."
  (loop for index from 0 below number-of-sequences
        for aggregate = (chem:make-aggregate (intern (format nil "seq-~a" index) :keyword))
        for molecule = (progn
                         (chem:goto-sequence oligomer index)
                         (build-molecule oligomer))
        do (chem:add-matter aggregate molecule)
           (cando:jostle aggregate)
        collect aggregate))

(defun find-starting-geometry (aggregate &rest arguments)
  (unless arguments
    (setf arguments '(:cg-tolerance 5.0 :max-tn-steps 100 :tn-tolerance 0.5)))
  (let ((mol (chem:content-at aggregate 0)))
    #+(or)(let ((root-residue-name (chem:get-name (chem:root-monomer oligomer))))
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
    (chem:add-monomer oligomer monomer)
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
  (let ((in-plug-name (chem:in-plug-name coupling-name))
        (out-plug-name (chem:out-plug-name coupling-name)))
    (let ((previous-monomer (loop for part in previous-parts
                                  for name = (chem:current-stereoisomer-name part)
                                  for topology = (chem:current-topology part) ; lookup-topology name)
                                  do (format t "do-coupling       part -> ~s~%" part)
                                     (format t "  (chem:current-stereoisomer-name part) -> ~s~%" (chem:current-stereoisomer-name part))
                                     (format t "      (type-of part) -> ~s~%" (type-of part))
                                     (format t "  name -> ~s~%" name)
                                     (format t "  (type-of name) -> ~s~%" (type-of name))
                                     (format t "              topology -> ~s~%" topology)
                                     (format t "         out-plug-name -> ~s~%" out-plug-name)
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
        (error "There is more than one monomer(~a) with the in-plug-name ~a" next-monomer in-plug-name))
      (if (null ring-info)
          (chem:couple oligomer
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
(make-oligomer '(:ccap :ala :default :ala :default :ser :default :ncap))
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
                       (format t "origin topology: ~s~%" topology)
                       (setf (gethash (chem:get-name topology) origins) topology))
                     (progn
                       (format t "  body topology: ~s~%" topology)
                       (setf (gethash (chem:get-name topology) body) topology))))
          else
            do (progn
                 (format t "   cap topology: ~s~%" topology)
                 (setf (gethash (chem:get-name topology) caps) topology)))
    (values origins body caps)))


(defun find-unsatisfied-monomer-plug-pairs (oligomer)
  (let ((monomers (chem:monomers-as-list oligomer)))
    (loop for monomer in monomers
          for monomer-name = (chem:current-stereoisomer-name monomer)
          for topology = (chem:current-topology monomer)
          do (format t "topology: ~s~%" topology)
          for plugs = (chem:plugs-as-list topology)
          do (format t "plugs: ~s~%" plugs)
          append (loop for plug in plugs
                       for plug-name = (chem:get-name plug)
                       when (and (or (typep plug 'chem:in-plug) (typep plug 'chem:out-plug))
                                 (not (chem:has-coupling-with-plug-name monomer plug-name)))
                         collect (cons monomer plug)))))


(defun one-round-extend-oligomer-with-caps (oligomer cap-name-map topology-map &key (verbose t))
  (let ((monomers-plugs (find-unsatisfied-monomer-plug-pairs oligomer)))
    (when verbose (format t "monomers-plugs: ~%~a~%" monomers-plugs))
    (if monomers-plugs
        (progn
          (when verbose (format t "Extending oligomer with caps~%"))
          (loop for (monomer . plug) in monomers-plugs
                for plug-name = (chem:get-name plug)
                for cap = (gethash (chem:get-name plug) cap-name-map)
                do (when verbose (format t "Extending monomer ~s  plug ~s~%" monomer plug))
                   (unless cap
                     (error "Could not find cap for monomer ~s plug ~s in cap-map"
                            monomer plug))
                   (multiple-value-bind (other-topology found)
                       (gethash cap topology-map)
                     (when verbose (format t "other-topology ~s~%" other-topology))
                     (unless found
                       (error "Could not find topology for cap ~s~%" cap))
                     (let ((other-monomer (chem:make-monomer (list (chem:get-name other-topology)))))
                       (when verbose (format t "Adding new monomer ~s~%" other-monomer))
                       (chem:add-monomer oligomer other-monomer)
                       (etypecase plug
                         (chem:out-plug
                          (when verbose (format t "Coupling monomer ~s with out coupling~%" other-monomer))
                          (chem:couple oligomer monomer (chem:coupling-name plug-name) other-monomer))
                         (chem:in-plug
                          (when verbose (format t "Coupling monomer ~s with in coupling~%" other-monomer))
                          (chem:couple oligomer other-monomer (chem:coupling-name plug-name) monomer))))))
          t)
        nil)))

#|
(format t "monomer ~s   plug ~s   cap ~s~%" monomer plug other-topology))))))

(let ((out-monomer (chem:make-monomer 
(chem:couple oligomer

|#

(defun build-superposable-conformation-collection (conformation monomer-id)
  "Return (values superposable-conformation-collection aggregate)."
  (let* ((fold-tree (kin:get-fold-tree conformation))
         (monomer-node (kin:lookup-monomer-id fold-tree monomer-id))
         (joint-ht (make-hash-table))
         joints)
    (kin:walk-joints monomer-node (lambda (index joint)
                                    (setf (gethash joint joint-ht) t)
                                    (push joint joints)))
    (let ((root-joint nil))
      (kin:walk-joints monomer-node
                       (lambda (index joint)
                         (when (null (gethash (kin:get-parent joint) joint-ht))
                           (setf root-node joint))))
      (when (typep root-joint 'kin:root-bonded-atom)
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
  (let ((oligomer (core:make-cxx-object 'chem:oligomer))
        (focus-monomer (chem:make-monomer (mapcar #'chem:get-name focus-topologys-in-list))))
    (format t "Focus monomer: ~a~%" focus-monomer)
    (chem:add-monomer oligomer focus-monomer)
    ;; Now repeatedly cap the focus monomer until it's finished.
    (unless test
      (loop while (one-round-extend-oligomer-with-caps oligomer cap-name-map topology-map)))
    oligomer))

(defun build-training-oligomers (list-of-topologys cap-name-map)
  "Create oligomers that describe training molecules"
  (let ((topology-map (let ((ht (make-hash-table)))
                        (loop for top in list-of-topologys 
                              do (setf (gethash (chem:get-name top) ht) top))
                        ht)))
    (multiple-value-bind (origins body caps)
        (classify-topologys list-of-topologys)
      ;; Build training oligomers for origins
      ;; First make the focus monomer and add it to the oligomer
      (append
       (loop for focus-topologys in (list origins body)
             collect (build-one-training-oligomer (alexandria:hash-table-values focus-topologys)
                                                  cap-name-map
                                                  topology-map))
       (loop for focus-cap-topology being the hash-values in caps
             collect (build-one-training-oligomer (list focus-cap-topology)
                                                  cap-name-map
                                                  topology-map))))))


