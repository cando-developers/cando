
(in-package :design.oligomer)

(defun deftop (parts-database residue in-plug-info &rest out-plugs-info)
  (let ((topology (cando:make-topology-from-residue residue)))
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
    (setf (gethash (chem:get-name topology) (topologies parts-database)) topology)
    (leap.core:register-variable (chem:get-name topology) topology)
    topology))

#+(or)
(defmacro define-topology (parts-database residue in-plug &rest out-plugs)
  (let ((topology (gensym "topology")))
    `(let ((,topology (cando:make-topology-from-residue ,residue)))
       ,(when in-plug
          (destructuring-bind (in-plug-name &optional (atom-name in-plug-name) (bond-order :single-bond))
              in-plug
            `(chem:add-plug ,topology ,in-plug-name (chem:make-in-plug ,in-plug-name nil ,atom-name ,bond-order))))
       ,@(loop for out-plug in out-plugs
               collect (destructuring-bind (out-plug-name &optional (atom-name out-plug-name) (bond-order :single-bond))
                           out-plug
                         `(chem:add-plug ,topology ,out-plug-name (chem:make-out-plug ,out-plug-name nil nil ,atom-name ,bond-order))))
       (setf (gethash (chem:get-name ,topology) (topologies ,parts-database)) ,topology)
       ,topology)))

(defun recursively-build-molecule (prev-monomer
                                   prev-topology
                                   prev-residue
                                   monomer-out-couplings
                                   molecule
                                   topologies
                                   monomers-to-residues)
  (loop for out-coupling in (gethash prev-monomer monomer-out-couplings)
        for coupling-name = (chem:get-name out-coupling)
        for next-monomer = (chem:get-out-monomer out-coupling)
        for next-topology = (gethash (chem:monomer-name next-monomer) topologies)
        for next-residue = (chem:build-residue next-topology)
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
                                       topologies
                                       monomers-to-residues)))

(defun build-molecule (oligomer topologies)
  (let ((root-monomer (chem:root-monomer oligomer))
        (monomer-out-couplings (make-hash-table))
        (monomers-to-residues (make-hash-table))
        (ring-couplings nil))
    (loop for coupling in (chem:couplings-as-list oligomer)
          if (typep coupling 'chem:directional-coupling)
            do (push coupling (gethash (chem:get-in-monomer coupling) monomer-out-couplings))
          else
            do (pushnew coupling ring-couplings)) ; Only add ring coupling when unique
    (let* ((molecule (chem:make-molecule :mol))
           (root-topology (gethash (chem:monomer-name root-monomer) topologies))
           (root-residue (chem:build-residue root-topology)))
      (setf (gethash root-monomer monomers-to-residues) root-residue)
      (chem:add-matter molecule root-residue)
      (recursively-build-molecule root-monomer
                                  root-topology
                                  root-residue
                                  monomer-out-couplings
                                  molecule
                                  topologies
                                  monomers-to-residues)
      ;; Now close the rings
      (loop for ring-coupling in ring-couplings
            for monomer1 = (chem:get-monomer1 ring-coupling)
            for topology1 = (gethash (chem:monomer-name monomer1) topologies)
            for residue1 = (gethash monomer1 monomers-to-residues)
            for monomer2 = (chem:get-monomer2 ring-coupling)
            for topology2 = (gethash (chem:monomer-name monomer2) topologies)
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



(defun build-all-molecules (oligomer part-database &optional (number-of-sequences (chem:number-of-sequences oligomer)))
  (loop for index from 0 below number-of-sequences
        for aggregate = (chem:make-aggregate (intern (format nil "seq-~a" index) :keyword))
        for molecule = (progn
                         (chem:goto-sequence oligomer index)
                         (build-molecule oligomer (topologies part-database)))
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

(defclass parts-database ()
  ((topologies :initform (make-hash-table) :accessor topologies)
   (parts :initform (make-hash-table) :accessor parts)))

(defun make-parts-database ()
  (make-instance 'parts-database))

(defun define-part (database name subtree)
    (setf (gethash name (parts database)) subtree))

(defun my-add-monomers (oligomer names parts-database)
  (let ((monomer (chem:make-monomer nil)))
    (loop for name in names
          for topology = (gethash name (topologies parts-database))
          for isoname = (chem:make-isoname name topology 0)
          do (chem:add-isoname monomer isoname))
    (chem:add-monomer oligomer monomer)
    (values (list monomer) nil)))

(defun translate-part (oligomer names parts-database labels)
  (cond
    ((and (symbolp names) (gethash names (parts parts-database)))
     (values (interpret-subtree oligomer maybe-part parts-database labels) nil))
    ((symbolp names)
     (my-add-monomers oligomer (list names) parts-database))
    ((consp names)
     (if (and (symbolp (first names)) (string= (first names) :cycle))
         (values (gethash (second names) labels) (list (third names) (fourth names)))
         (my-add-monomers oligomer names parts-database)))
    (t (error "Handle ~a" names))))
            
(defun interpret-part (oligomer part-info parts-database labels)
  (destructuring-bind (names &key label)
      (cond
        ((and (consp part-info) (symbolp (first part-info)) (string= (first part-info) :cycle))
         (list part-info))
        ((consp part-info)
         part-info)
        ((symbolp part-info)
         (list part-info)))
    (multiple-value-bind (new-parts ringp)
        (translate-part oligomer names parts-database labels)
      (when label
        (loop for new-part in new-parts
              do (push new-part (gethash label labels))))
      (values new-parts ringp))))

(defun do-coupling (oligomer coupling-name ring-info previous-parts next-parts parts-database)
  (let ((in-plug-name (chem:in-plug-name coupling-name))
        (out-plug-name (chem:out-plug-name coupling-name)))
    (let ((previous-monomer (loop for part in previous-parts
                                  for name = (chem:monomer-name part)
                                  for topology = (gethash name (topologies parts-database))
                                  do (format t "do-coupling       part -> ~s~%" part)
                                     (format t "  (chem:get-one-monomer part) -> ~s~%" (chem:get-one-monomer part))
                                     (format t "      (type-of part) -> ~s~%" (type-of part))
                                     (format t "  name -> ~s~%" name)
                                     (format t "  (type-of name) -> ~s~%" (type-of name))
                                     (format t "              topology -> ~s~%" topology)
                                     (format t "         out-plug-name -> ~s~%" out-plug-name)
                                  when (chem:has-plug-named topology out-plug-name)
                                    collect part))
          (next-monomer (loop for part in next-parts
                              for topology = (gethash (chem:monomer-name part) (topologies parts-database))
                              when (chem:has-plug-named topology in-plug-name)
                                collect part)))
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
  
(defun make-oligomer (parts-database tree)
  (let* ((oligomer (core:make-cxx-object 'chem:oligomer))
         (labels (make-hash-table)))
    (interpret-subtree oligomer tree parts-database labels)
    oligomer))

(defun interpret-subtree (oligomer subtree parts-database labels)
  (let* ((root-monomer-info (pop subtree))
         (previous-parts (interpret-part oligomer root-monomer-info parts-database labels))
         (accumulated-parts previous-parts))
    (loop
      (when (null subtree) (return accumulated-parts))
      (let ((info (pop subtree)))
        (if (consp info)
            (let* ((coupling (first info)) ; interpret a branch
                   (node-info (cdr info)))
              (multiple-value-bind (new-parts ringp)
                  (interpret-part oligomer node-info parts-database labels)
                (setf accumulated-parts (append accumulated-parts new-parts))
                (do-coupling oligomer coupling ringp previous-parts new-parts parts-database)))
            (let* ((coupling info)      ; Interpret a chain
                   (node-info (pop subtree)))
              (multiple-value-bind (new-parts ringp)
                  (interpret-part oligomer node-info parts-database labels)
                (setf accumulated-parts (append accumulated-parts new-parts))
                (do-coupling oligomer coupling ringp previous-parts new-parts parts-database)
                (setf previous-parts new-parts))))))))
    

