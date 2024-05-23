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
    (unless out-plug
      (error "Could not find out-plug in ~s named ~s - available plugs: ~s" prev-topology out-plug-name (alexandria:hash-table-keys (topology:plugs prev-topology))))
    (unless in-plug
      (error "Could not find in-plug in ~s named ~s - available plugs: ~s" next-topology in-plug-name (alexandria:hash-table-keys (topology:plugs next-topology))))
    (unless (= (length (topology:plug-bonds out-plug)) (length (topology:plug-bonds in-plug)))
      (error "There is a mismatch between the number of plug-bonds in ~s and ~s" out-plug in-plug))
    (loop for bond-index below (length (topology:plug-bonds out-plug))
          for out-plug-bond = (elt (topology:plug-bonds out-plug) bond-index)
          for in-plug-bond = (elt (topology:plug-bonds in-plug) bond-index)
          for out-plug-atom-name = (topology:atom-name out-plug-bond)
          for in-plug-atom-name = (topology:atom-name in-plug-bond)
          for out-atom = (chem:atom-with-name prev-residue out-plug-atom-name t)
          for in-atom = (chem:atom-with-name next-residue in-plug-atom-name t)
          do (chem:bond-to in-atom out-atom :single-bond))))

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

(cando.serialize:make-class-save-load foldamer)

(defun save-foldamer (foldamer file-name)
  (let ((pn (pathname file-name)))
    (ensure-directories-exist pn)
    (cando.serialize:save-cando foldamer file-name)))

(defun load-foldamer (file-name)
  (let ((foldamer (cando.serialize:load-cando file-name)))
    (loop for topology in (topologys foldamer)
          do (chem:register-topology topology))
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

(defun my-add-monomers (oligomer names &key id)
  (when (consp (car names))
    (error "Illegal names - must be simple list - instead got ~s" names))
  (when (and (= (length names) 1) (null (car names)))
    (error "Illegal monomer names ~s" names))
  (let ((monomer (make-instance 'monomer :monomers names :id id)))
    (vector-push-extend monomer (monomers oligomer))
    (values (list monomer) nil)))

(defun translate-part (oligomer names labels part-info &key (parts *parts*))
  (let (maybe-part)
    (cond
      ((and (symbolp names)
            (setf maybe-part (lookup-maybe-part names parts)))
       (values (interpret-rooted-tree oligomer (tree maybe-part) labels) nil))
      ((symbolp names)
       (when (null names)
         (error "In translate-part illegal name ~s part-info ~s" names part-info))
       (my-add-monomers oligomer (list names) :id part-info))
      ((consp names)
       (if (and (symbolp (first names)) (string= (first names) :cycle))
           (let ((ringa (third names))
                 (ringb (fourth names)))
             (unless (and ringa ringb)
               (error "You must provide two additional arguments for ~s" names))
             (values (let ((res (gethash (second names) labels)))
                       (unless (atom res)
                         (error "The value in labels must be a single atom"))
                       res)
                       (list ringa ringb)))
           (my-add-monomers oligomer names :id part-info)))
      (t (error "Handle translate-part for ~s" names)))))

(defun interpret-part (oligomer part-info labels &key (parts *parts*))
  (destructuring-bind (names &key label)
      (cond
        ;; Handle (:cycle ...)
        ((and (consp part-info) (symbolp (first part-info)) (string= (string (first part-info)) "CYCLE"))
         (list part-info))
        ;; Handle (group-name ... plist)
        ((and (consp part-info) (gethash (first part-info) *topology-groups*))
         (let ((group-names (gethash (first part-info) *topology-groups*)))
           (list* group-names (cdr part-info))))
        ((consp part-info)
         part-info)
        ((symbolp part-info)
         (let ((group-names (gethash part-info *topology-groups*)))
           (if (null group-names)
               (list part-info)
               (list group-names))))
        (t (error "Handle part-info ~s" part-info)))
    (multiple-value-bind (new-parts ringp)
        (translate-part oligomer names labels part-info :parts parts)
      (when label
        (loop for new-part in new-parts
              do (let ((res (gethash label labels)))
                   (when res
                     (error "There is already a monomer ~s with the label ~s" res label))
                   (setf (gethash label labels) new-part))))
      (values new-parts ringp))))

(defun parts-with-plugs (parts plug-name)
  (loop for part in parts
        when (monomer-topologys-all-have-plug-named part plug-name)
          collect part))

(defun do-coupling (oligomer coupling-name ring-info previous-parts next-parts)
  (unless next-parts
    (error "There are no next parts to add to ~a" previous-parts))
  (let ((in-plug-name (in-plug-name coupling-name))
        (out-plug-name (out-plug-name coupling-name)))
    (let ((previous-monomer (parts-with-plugs previous-parts out-plug-name))
          (next-monomer (parts-with-plugs next-parts in-plug-name)))
      (when (= (length previous-monomer) 0)
        (error "There is no monomer found with the out-plug-name ~a in the parts ~a" out-plug-name previous-parts))
      #+(or)(unless (= (length previous-monomer) 1)
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
            (ring-couple-with-plug-names oligomer
                                         (first previous-monomer)
                                         plug1name
                                         (first next-monomer)
                                         plug2name))))))

(defun interpret-subtree (oligomer subtree labels previous-parts &key (parts *parts*))
  (let ((accumulated-parts previous-parts))
    (loop
      (when (null subtree) (return accumulated-parts))
      (let ((info (pop subtree)))
        (format t "interpret-subtree info: ~s~%" info)
        (cond
          ((and (consp info) (consp (car info)) (string= (string (car (car info))) "RING"))
           (let ((ring-info (cdr (car info))))
             (unless (= (length ring-info) 3)
               (error "A RING specifier must have 3 arguments (RING prev-plug-name next-plug-name next-label) but only has ~s" ring-info))
             (let* ((previous-plug-name (first ring-info))
                    (next-plug-name (second ring-info))
                    (label (third ring-info))
                    (next-monomer (gethash label labels))
                    (previous-monomer (parts-with-plugs previous-parts previous-plug-name))
                    )
               (unless previous-monomer
                 (error "Could not find previous monomer ~s with plug-name ~s" previous-parts previous-plug-name))
               (ring-couple-with-plug-names oligomer
                                            (first previous-monomer)
                                            previous-plug-name
                                            next-monomer ; (first next-monomer)
                                            next-plug-name)
               )))
          ((and (consp info) (consp (car info)))
           (error "Add support to handle ~s - currently only RING is allowed as a special edge"))
          ((consp info)
           (interpret-subtree oligomer info labels accumulated-parts :parts parts)
           #+(or)(let* ((coupling (first info)) ; interpret a branch
                        (node-info (cadr info))) ;; CHECK
                   (multiple-value-bind (new-parts ringp)
                       (interpret-part oligomer node-info labels :parts parts)
                     (setf accumulated-parts (append accumulated-parts new-parts))
                     (do-coupling oligomer coupling ringp previous-parts new-parts))))
          (t
           (let* ((coupling info)       ; Interpret a chain
                  (node-info (pop subtree)))
             (multiple-value-bind (new-parts ringp)
                 (interpret-part oligomer node-info labels :parts parts)
               (setf accumulated-parts (append new-parts accumulated-parts))
               (do-coupling oligomer coupling ringp previous-parts new-parts)
               (setf previous-parts new-parts)))))))))

(defun interpret-rooted-tree (oligomer subtree labels &key (parts *parts*))
  (let* ((root-monomer-info (pop subtree))
         (previous-parts (interpret-part oligomer root-monomer-info labels :parts parts))
         (accumulated-parts previous-parts))
    (interpret-subtree oligomer subtree labels accumulated-parts :parts parts)
    ))

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
       (error "There are no out-plugs in ~s that match the in-plug-name ~s - existing plugs: ~s" other-monomer in-plug-name (plugs-as-list topology)))
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
                     (when (null (name other-topology))
                       (error "Illegal monomer with name ~s" (name other-topology)))
                     (let ((other-monomer (make-instance 'monomer
                                                         :monomers (list (name other-topology)))))
                       (break "one-round-extend Assign the id to the monomer")
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


(defun build-one-training-oligomer-space (focus-topologys-in-list cap-name-map topology-map &key test)
  "Build a training oligomer-space with a list of topologys in a focus-residue and then repeatedly
add cap monomers until no more cap monomers are needed."
  (let ((oligomer-space (make-instance 'oligomer-space))
        (focus-monomer (let ((names (mapcar #'name focus-topologys-in-list)))
                         (when (null (car names))
                           (error "Illegal names for monomer ~s" names))
                         (break "build-one-training-oligomer-space Assign the id to the monomer using names")
                         (make-instance 'monomer :monomers names))))
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

