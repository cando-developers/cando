(in-package :topology)

(defvar *complex-plugs* (make-hash-table))
(defvar *caps* (make-hash-table))

(defun connect-residues (prev-topology
                         prev-residue
                         out-plug-name
                         next-topology
                         next-residue
                         in-plug-name)
  (let* ((out-plug (topology:plug-named prev-topology out-plug-name))
         (in-plug (topology:plug-named next-topology in-plug-name)))
    (loop for bond-index below (length (topology:atom-names out-plug))
          do (let* ((out-plug-atom-name (elt (topology:atom-names out-plug) bond-index))
                    (in-plug-atom-name (elt (topology:atom-names in-plug) bond-index))
                    (out-atom (chem:atom-with-name prev-residue out-plug-atom-name t))
                    (in-atom (chem:atom-with-name next-residue in-plug-atom-name t))
                    (bo-in0 :single-bond)
                    (bo-out0 :single-bond))
               (declare (ignore bo-out0))
               (chem:bond-to in-atom out-atom bo-in0)))))

(defun out-plug-name-p (name)
  (when (and name (symbolp name))
    (char= #\+ (elt (symbol-name name) 0))))

(defun in-plug-name-p (name)
  (when (and name (symbolp name))
    (char= #\- (elt (symbol-name name) 0))))

(defun is-coupling-name-p (name)
  (when (and name (symbolp name))
    (let ((first-char (schar (symbol-name name) 0)))
      (and (char/= #\- first-char)
           (char/= #\+ first-char)
           )
      )))

(defun coupling-name (plug-name)
  (if (and plug-name
           (symbolp plug-name)
           (not (is-coupling-name-p plug-name)))
      (intern (subseq (symbol-name plug-name) 1) :keyword)
      plug-name))

(defun in-plug-name (name)
  (unless (is-coupling-name-p name)
    (error "Not coupling name ~a" name))
  (intern (concatenate 'base-string "-" (symbol-name name)) :keyword))

(defun out-plug-name (name)
  (unless (is-coupling-name-p name)
    (error "Not coupling name ~a" name))
  (intern (concatenate 'base-string "+" (symbol-name name)) :keyword))


(defun validate-couple (oligomer source-mon coupling-or-source-plug-name target-mon target-plug-name)
  (let (source-plug-name)
    (cond
      ((consp coupling-or-source-plug-name)
       (setf source-plug-name (car coupling-or-source-plug-name))
       (when target-plug-name
         (error "coupling-or-source-plug-name was a pair and so target-plug-name must be nil - but it is ~a" target-plug-name))
       (setf target-plug-name (cdr coupling-or-source-plug-name)))
      ((symbolp coupling-or-source-plug-name)
       (when (and (is-coupling-name-p coupling-or-source-plug-name)
                  target-plug-name)
         (error "coupling-or-source-plug-name ~s was a coupling name and so target-plug-name (~s) must be NIL" coupling-or-source-plug-name target-plug-name))
       (if (is-coupling-name-p coupling-or-source-plug-name)
           (setf source-plug-name (out-plug-name coupling-or-source-plug-name)
                 target-plug-name (in-plug-name coupling-or-source-plug-name))
           (progn
             (unless (out-plug-name-p coupling-or-source-plug-name)
               (error "coupling-or-source-plug-name must be a source-plug-name - it is a ~a" coupling-or-source-plug-name))
             (setf source-plug-name coupling-or-source-plug-name)
             (unless (in-plug-name-p target-plug-name)
               (error "target-plug-name must be a in-plug-name - it is a ~a" target-plug-name)))))
      (t (error "coupling-or-source-plug-name must be a cons or a symbol - it is ~a" coupling-or-source-plug-name)))
    (let ((found-in (position source-mon (monomers oligomer)))
          (found-out (position target-mon (monomers oligomer))))
      (unless found-in (error "Could not find ~a in oligomer" source-mon))
      (unless found-out (error "Could not find ~a in oligomer" target-mon)))
    (values source-plug-name target-plug-name)))

(defun couple (oligomer source-mon coupling-or-source-plug-name target-mon &optional target-plug-name)
  "Couple two monomers together in the oligomer"
  (multiple-value-bind (source-plug-name target-plug-name)
      (validate-couple oligomer source-mon coupling-or-source-plug-name target-mon target-plug-name)
    (let ((coupling (make-instance 'directional-coupling
                                   :name (coupling-name coupling-or-source-plug-name)
                                   :source-monomer source-mon
                                   :target-monomer target-mon
                                   :source-plug-name source-plug-name
                                   :target-plug-name target-plug-name)))
      (setf (gethash source-plug-name (couplings source-mon)) coupling
            (gethash target-plug-name (couplings target-mon)) coupling)
      (vector-push-extend coupling (couplings oligomer)))))

(defun ring-couple-with-plug-names (oligomer mon1 plug1-name mon2 plug2-name)
  (let ((found1 (position mon1 (monomers oligomer)))
        (found2 (position mon2 (monomers oligomer))))
    (unless found1 (error "Could not find ~a in oligomer" mon1))
    (unless found2 (error "Could not find ~a in oligomer" mon2)))
  (let ((coupling (make-instance 'ring-coupling
                                 :plug1 plug1-name
                                 :plug2 plug2-name
                                 :monomer1 mon1
                                 :monomer2 mon2)))
    (setf (gethash plug1-name (couplings mon1)) coupling
          (gethash plug2-name (couplings mon2)) coupling)
    (vector-push-extend coupling (couplings oligomer))))



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




(defun do-define-cap (name topology-name)
  (setf (gethash name *caps*) topology-name))

(defmacro define-cap (name topology-name)
  `(do-define-cap ',name ',topology-name))

(defclass foldamer ()
  ((topology-list :initarg :topology-list :accessor topology-list)
   (cap-name-map :initarg :cap-name-map :accessor cap-name-map)))

(defun make-foldamer (topology-list &optional (cap-name-map *caps*))
  (make-instance 'foldamer :topology-list topology-list
                 :cap-name-map cap-name-map))

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

(cando:make-class-save-load foldamer)

(defun save-foldamer (foldamer file-name)
  (let ((pn (pathname file-name)))
    (ensure-directories-exist pn)
    (cando:save-cando foldamer file-name)))

(defun load-foldamer (file-name)
  (let ((foldamer (cando:load-cando file-name)))
    (loop for topology in (topologys foldamer)
          do (cando:register-topology topology))
    foldamer))

(defun lookup-maybe-part (name parts)
  "Lookup a named object and return it.
This is for looking up parts but if the thing returned is not a part then return nil."
  (gethash name parts))

(defclass subtree ()
  ((name :initarg :name :accessor name)
   (code :initarg :code :accessor code)))

(defun make-subtree (name code)
  (make-instance 'subtree :name name :code code))

(defclass part ()
  ((name :initarg :name :accessor name)
   (tree :initarg :tree :accessor tree)))

(defun make-part (name tree)
  (make-instance 'part :name name :tree tree))

(defparameter *parts* (make-hash-table))

(defun define-part (name tree)
  "Create a part and register it by name with leap.core:register-variable."
  (let ((part (make-part name tree)))
    (setf (gethash name *parts*) part)))

(defun my-add-monomers (oligomer names)
  (let ((monomer (make-instance 'monomer :monomers names)))
    (vector-push-extend monomer (monomers oligomer))
    (values (list monomer) nil)))

(defun translate-part (oligomer names labels &key (parts *parts*))
  (let (maybe-part)
    (cond
      ((and (symbolp names)
            (setf maybe-part (lookup-maybe-part names parts)))
       (values (interpret-subtree oligomer (tree maybe-part) labels) nil))
      ((symbolp names)
       (my-add-monomers oligomer (list names)))
      ((consp names)
       (if (and (symbolp (first names)) (string= (first names) :cycle))
           (values (gethash (second names) labels) (list (third names) (fourth names)))
           (my-add-monomers oligomer names)))
      (t (error "Handle ~a" names)))))

(defun interpret-part (oligomer part-info labels &key (parts *parts*))
  (destructuring-bind (names &key label)
      (cond
        ((and (consp part-info) (symbolp (first part-info)) (string= (first part-info) :cycle))
         (list part-info))
        ((consp part-info)
         part-info)
        ((symbolp part-info)
         (list part-info)))
    (multiple-value-bind (new-parts ringp)
        (translate-part oligomer names labels :parts parts)
      (when label
        (loop for new-part in new-parts
              do (push new-part (gethash label labels))))
      (values new-parts ringp))))

(defun parts-with-plugs (parts plug-name)
  (loop for part in parts
        when (monomer-topologys-all-have-plug-named part plug-name)
          collect part))

(defun do-coupling (oligomer coupling-name ring-info previous-parts next-parts)
  (unless next-parts
    (error "There are no next parts to add to ~a" previous-parts))
  (let ((in-plug-name (chem:in-plug-name coupling-name))
        (out-plug-name (chem:out-plug-name coupling-name)))
    (let ((previous-monomer (parts-with-plugs previous-parts out-plug-name))
          (next-monomer (parts-with-plugs next-parts in-plug-name)))
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
          (couple oligomer
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

(defun interpret-subtree (oligomer subtree labels &key (parts *parts*))
  (let* ((root-monomer-info (pop subtree))
         (previous-parts (interpret-part oligomer root-monomer-info labels :parts parts))
         (accumulated-parts previous-parts))
    (loop
      (when (null subtree) (return accumulated-parts))
      (let ((info (pop subtree)))
        (if (consp info)
            (progn
              #+(or)(format *debug-io* "interpret-subtree info: ~s~%" info)
              (let* ((coupling (first info)) ; interpret a branch
                     (node-info (cdr info)))
                (multiple-value-bind (new-parts ringp)
                    (interpret-part oligomer node-info labels :parts parts)
                  (setf accumulated-parts (append accumulated-parts new-parts))
                  (do-coupling oligomer coupling ringp previous-parts new-parts))))
            (progn
              #+(or)(format *debug-io* "interpret-subtree info: ~s~%" info)
              (let* ((coupling info)    ; Interpret a chain
                     (node-info (pop subtree)))
                (multiple-value-bind (new-parts ringp)
                    (interpret-part oligomer node-info labels :parts parts)
                  (setf accumulated-parts (append accumulated-parts new-parts))
                  (do-coupling oligomer coupling ringp previous-parts new-parts)
                  (setf previous-parts new-parts)))))))))

(defun make-oligomer-space (foldamer tree &key (parts *parts*))
  "Make an oligomer-space from a description in the **tree**.
The tree is a nested list of lists that look like 
(component coupling component coupling component ... ).
It starts with the name of a topology or a design:part.
The component can be a name of a single topology, it can be a design:part,
or it can 

Examples:
(make-oligomer-space '(:ccap :default :ala :default :ala :default :ser :default :ncap))
(make-oligomer-space '((:lego3 :label :first)
                 :default :lego3
                 :default :lego3
                 :default (cycle :first)))
"
  (let* ((oligomer-space (make-instance 'oligomer-space
                                        :foldamer foldamer))
         (labels (make-hash-table)))
    (interpret-subtree oligomer-space tree labels :parts parts)
    oligomer-space))

(defun classify-topologys (topology-hash-table)
  (let ((origins (make-hash-table))
        (body (make-hash-table))
        (caps (make-hash-table)))
    (maphash (lambda (name topology)
               (declare (ignore name))
               (if (out-plugs-as-list topology)
                   (let ((in-plug (find-in-plug topology)))
                     (if (null in-plug)
                         (progn
                           (format *debug-io* "starting topology: ~s~%" topology)
                           (setf (gethash (name topology) origins) topology))
                         (progn
                           (format *debug-io* "  body topology: ~s~%" topology)
                           (setf (gethash (name topology) body) topology))))
                   (progn
                     (format *debug-io* "   end topology: ~s~%" topology)
                     (setf (gethash (name topology) caps) topology))))
             topology-hash-table)
    (values (alexandria:hash-table-values origins)
            (alexandria:hash-table-values body)
            (alexandria:hash-table-values caps))))

(defun cluster-topologys-using-out-plugs (list-of-topologys)
  "The input is a list-of-topologys with the same in plug but maybe different out-plugs.
Subdivide the list into a list of lists that contain Topologys with different numbers and kinds
of out-plugs."
  (let ((clusters (make-hash-table :test #'equal)))
    (loop for topology in list-of-topologys
          for out-plug-names = (mapcar #'name (out-plugs-as-list topology))
          for sorted-out-plug-names = (sort out-plug-names #'string<)
          do (push topology (gethash sorted-out-plug-names clusters)))
    (alexandria:hash-table-values clusters)))

(defun ensure-one-unique-out-plug-name (other-monomer in-plug-name)
  "Ensure that there is one topology with one out-plug with a name that corresponds to in-plug-name on
   a monomer that this one out-plug will be coupled through"
  (let* ((topology (monomer-topology other-monomer))
         (out-plug-names (all-out-plug-names-that-match-in-plug-name topology in-plug-name)))
    (format *debug-io* "out-plug-names -> ~s~%" out-plug-names)
    (case (length out-plug-names)
      (0
       (error "There are no out-plugs in ~s that match the in-plug-name ~s - existing plugs: ~s" other-monomer in-plug-name (chem:topology/plugs-as-list topology)))
      (1
       (return-from ensure-one-unique-out-plug-name (first out-plug-names)))
      (otherwise
       (error "There are multiple out-plugs ~s in ~s that match the in-plug-name ~s"
              out-plug-names other-monomer in-plug-name)))))

(defun find-unsatisfied-monomer-plug-pairs (oligomer-space)
  (let ((monomer-plug-work-list (make-hash-table :test #'equal)))
    (loop for monomer-index below (length (monomers oligomer-space))
          for monomer = (elt (monomers oligomer-space) monomer-index)
          for topology = (monomer-topology monomer)
          for plugs = (plugs-as-list topology)
          do (loop for plug in plugs
                   for plug-name = (name plug)
                   when (or (typep plug 'in-plug)
                            (typep plug 'out-plug))
                     do (setf (gethash (cons monomer plug-name) monomer-plug-work-list) t)))
    (loop for coupling-index below (length (couplings oligomer-space))
          for coupling = (elt (couplings oligomer-space) coupling-index)
          when (typep coupling 'directional-coupling)
            do (setf (gethash (cons (source-monomer coupling)
                                    (source-plug-name coupling))
                              monomer-plug-work-list) nil)
               (setf (gethash (cons (target-monomer coupling)
                                    (target-plug-name coupling))
                              monomer-plug-work-list) nil))
    (let (result)
      (maphash (lambda (monomer-plug-name dangling-p)
                 (when dangling-p
                   (push monomer-plug-name result)))
               monomer-plug-work-list)
      result)))

(defun one-round-extend-oligomer-space-with-caps (oligomer-space cap-name-map topology-map &key (verbose t))
  (let ((monomers-plugs (find-unsatisfied-monomer-plug-pairs oligomer-space)))
    (when verbose (format *debug-io* "monomers-plugs: ~%~a~%" monomers-plugs))
    (if monomers-plugs
        (progn
          (when verbose (format *debug-io* "Extending oligomer-space with caps~%"))
          (loop for (monomer . plug-name) in monomers-plugs
                for topology = (monomer-topology monomer)
                for plug = (plug-named topology plug-name)
                for cap = (gethash (name plug) cap-name-map)
                do (when verbose (format *debug-io* "Extending monomer ~s  plug ~s~%" monomer plug))
                   (unless cap
                     (error "Could not find cap for monomer ~s plug ~s in cap-map"
                            monomer plug-name))
                   (multiple-value-bind (other-topology found)
                       (gethash cap topology-map)
                     (when verbose (format *debug-io* "other-topology ~s~%" other-topology))
                     (unless found
                       (error "Could not find topology for cap ~s~%" cap))
                     (let ((other-monomer (make-instance 'monomer
                                                         :monomers (list (name other-topology)))))
                       (when verbose (format *debug-io* "Adding new monomer ~s~%" other-monomer))
                       (add-monomer oligomer-space other-monomer)
                       (etypecase plug
                         (out-plug
                          (let* ((other-monomer-in-plug-name (in-plug-name (coupling-name plug-name)))
                                 (other-monomer-topology (monomer-topology other-monomer)))
                            (unless (plug-named other-monomer-topology other-monomer-in-plug-name)
                              (error "While trying to couple monomer ~s with out-plug named ~s we could not find a corresponding in-plug named ~s in the cap monomer ~s"
                                     monomer plug-name other-monomer-in-plug-name other-monomer))
                            (when verbose (format *debug-io* "Coupling to out-coupling ~s ~s ~s ~s~%" monomer plug-name other-monomer other-monomer-in-plug-name))
                            (couple oligomer-space monomer plug-name other-monomer other-monomer-in-plug-name)))
                         (in-plug
                          (when verbose (format *debug-io* "Coupling monomer ~s with in coupling~%" other-monomer))
                          (let ((other-monomer-out-plug-name (ensure-one-unique-out-plug-name other-monomer plug-name)))
                            (when verbose (format *debug-io* "Coupling to in-coupling ~s ~s ~s ~s~%" other-monomer other-monomer-out-plug-name monomer plug-name))
                            (couple oligomer-space other-monomer other-monomer-out-plug-name monomer plug-name)))
                         ))))
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


(defun build-one-training-oligomer-space (focus-topologys-in-list cap-name-map topology-map &key test)
  "Build a training oligomer-space with a list of topologys in a focus-residue and then repeatedly
add cap monomers until no more cap monomers are needed."
  (let ((oligomer-space (make-instance 'oligomer-space))
        (focus-monomer (make-instance 'monomer :monomers (mapcar #'name focus-topologys-in-list))))
    (format *debug-io* "focus-monomer: ~a~%" focus-monomer)
    (add-monomer oligomer-space focus-monomer)
    ;; Now repeatedly cap the focus monomer until it's finished.
    (unless test
      (loop for step from 0
            while (one-round-extend-oligomer-space-with-caps oligomer-space cap-name-map topology-map)
            when (> step 10)
              do (error "Too many rounds of extending oligomer-space")))
    oligomer-space))

(defun build-training-oligomer-spaces (foldamer)
  "Create oligomer-spaces that describe training molecules"
  (let* ((list-of-topology-names (topology-list foldamer))
         (topology-map (let ((ht (make-hash-table)))
                         (loop for top-name in list-of-topology-names
                               for top = (chem:find-topology top-name t)
                               do (setf (gethash (name top) ht) top))
                         ht)))
    (multiple-value-bind (origins body caps)
        (classify-topologys topology-map)
      ;; Build training oligomer-spaces for origins
      ;; First make the focus monomer and add it to the oligomer-space
      (let ((clustered-origins (cluster-topologys-using-out-plugs origins))
            (clustered-bodys (cluster-topologys-using-out-plugs body))
            (cap-name-map (cap-name-map foldamer)))
        (format *debug-io* "clustered-origins -> ~s~%" clustered-origins)
        (format *debug-io* "clustered-bodys -> ~s~%" clustered-bodys)
        (let ((all-oligomer-spaces (append
                              (loop for focus-topologys in (append clustered-origins clustered-bodys)
                                    collect (build-one-training-oligomer-space focus-topologys
                                                                         cap-name-map
                                                                         topology-map))
                              (loop for focus-cap-topology in caps
                                    collect (build-one-training-oligomer-space (list focus-cap-topology)
                                                                         cap-name-map
                                                                         topology-map)))))
          ;; Eliminate oligomer-spaces that cover the same sequence space
          (let ((unique-oligomer-spaces (make-hash-table :test 'equal)))
            (loop for oligomer-space in all-oligomer-spaces
                  for canonical-sequence = (canonical-sequence oligomer-space)
                  do (setf (gethash canonical-sequence unique-oligomer-spaces) oligomer-space))
            (alexandria:hash-table-values unique-oligomer-spaces)))))))

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


#||
(define-foldamer-schema :spiro
  '(bnz pbe pro4 apro4 cbala bala ace cgly)
  '((:side :-side pro4 (:-dkp pro4 :-dkp cgly) :+dkp pro4 (:+side (:cap bnz)) :+dkp (:cap pros))
    (:side :-side apro4 (:-amide cbala) (:+amide ace) :+dkp pro4 (:+side (:cap bnz)) :+dkp (:cap pros))
    (:side :-side apro4 (:-amide cbala) (:+amide ace) :+dkp pro))
    (:pro4 (:-dkp pro4 (:-dkp (:cap cgly)) (:+side (:cap bnz))) (:+side (:cap bnz)) :+dkp (:cap pros))
    (:pro4 (:-dkp apro4 (:-amide (:cap cbala)) (:+side (:cap bnz)) (:+amide (:cap ace))) :+dkp (:cap pros))
    (:apro4 (:-amide cbala) (:+amide (:cap ace)) (:+side (:cap bnz)) :+dkp pro)
    (:cbala (:+amide (:cap ace)))
    (:bala (:-amide (:cap cbala)) :+amide (:cap ace))
    (:ace :-amide (:cap cbala))
(:cgly :+dkp pro))
||#
