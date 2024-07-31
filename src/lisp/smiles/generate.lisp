(in-package :smiles)

(defclass smiles-node ()
  ((atm :initarg :atm :reader atm)
   (children :initarg :children :accessor children)
   (begin-ring-other-atoms :initform nil :accessor begin-ring-other-atoms)
   (begin-rings :initform nil :accessor begin-rings)
   (end-rings :initarg :end-rings :initform nil :accessor end-rings)))

(defmethod print-object ((object smiles-node) stream)
    (print-unreadable-object (object stream :type t)
      (format stream "~s" (atm object))))

(defclass spanning-tree-builder ()
  ((root-node :accessor root-node)
   (ranks :initarg :ranks :reader ranks)
   (ring-ids :initform (make-array 10 :initial-element nil) :accessor ring-ids)
   (atoms-to-nodes :initform (make-hash-table) :reader atoms-to-nodes)
   (ring-bonds :initform nil :accessor ring-bonds)))

(defclass depth-first-search-builder ()
  ((ranks :initarg :ranks :reader ranks)
   (root-atom :accessor root-atom)
   (seen :initform (make-hash-table) :reader seen)
   (seen-smiles :initform (make-hash-table) :reader seen-smiles)
   (closure-bonds :initform (make-hash-table :test 'eql) :reader closure-bonds)
   (unclosed-closures :initform (make-hash-table :test 'eql) :reader unclosed-closures)
   (camgen-nodes :initform (make-array 16 :adjustable t :fill-pointer 0) :reader camgen-nodes)
   (outgoing-edges :initform (make-array 16 :adjustable t :fill-pointer 0) :reader outgoing-edges)
   (available-closures :initform '("1" "2" "3" "4" "5" "6" "7" "8" "0" "%10")
                       :accessor available-closures)))

(defclass camgen-node ()
  ((atm :initarg :atm :reader atm)
   (outgoing-edges :initform (make-array 4 :adjustable t :fill-pointer 0) :reader outgoing-edges)))

(defmethod print-object ((obj camgen-node) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~s" (chem:atom/atom-name (atm obj)))))

(defclass outgoing-edge ()
  ((atm :initarg :atm :reader atm)
   (other-atom :initarg :other-atom :reader other-atom)
   (other-node-index :accessor other-node-index)
   (bond-order :initarg :bond-order :reader bond-order)
   (bond-index :initarg :bond-index :accessor bond-index)
   ))

(defmethod print-object ((obj outgoing-edge) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~s ~a" (atm obj) (bond-order obj) (other-atom obj))))

(defun next-available-ring-id (ring-ids)
  (let ((next-pos (position nil ring-ids :start 1)))
    (if next-pos
        (progn
          (setf (aref ring-ids next-pos) t)
          next-pos)
        (error "Ran out of ring-ids"))))

#+(or)
(defmacro dbg (fmt &rest args)
  `(format t ,fmt ,@args))

;;#+(or)
(defmacro dbg (fmt &rest args)
  nil)

(defun calculate-other-node-index (dfs)
  (loop for outgoing-edge across (outgoing-edges dfs)
        for index = (position (other-atom outgoing-edge) (camgen-nodes dfs) :key #'atm)
        do (setf (other-node-index outgoing-edge) index)))

(defun atom-smarts (atm)
  (let ((element (chem:atom/get-element atm))
        (aromaticp (gethash atm chem:*current-aromaticity-information*)))
    (if (member element '(:B :C :N :O :S :P))
        (if aromaticp
            (string-downcase element)
            (string element))
        (if (member element '(:F :|Cl| :|Br| :I))
            (string element)
            (format nil "[~a]" element)))))


(defun get-closure-label (bond-smarts closure)
  (format nil "~a" closure))

(defun bond-smarts (outgoing-edge)
  (case (bond-order outgoing-edge)
    (:single-bond "")
    (:double-bond
     (if (and (gethash (atm outgoing-edge) chem:*current-aromaticity-information*)
              (gethash (other-atom outgoing-edge) chem:*current-aromaticity-information*))
         ""
         "="))
    (:triple-bond "#")))

(defun depth-first-search-ring-finder-impl (prev-atm atm dfs)
  (setf (gethash atm (seen dfs)) t)
  (let ((node (make-instance 'camgen-node
                             :atm atm))
        unsorted-children)
    (vector-push-extend node (camgen-nodes dfs))
    (chem:map-bonds
     nil
     (lambda (a1 a2 order bond)
       (declare (ignore a1 a2))
       (let* ((other-atom (chem:bond/get-other-atom bond atm))
              (other-element (chem:atom/get-element other-atom)))
         (when (and (not (eq prev-atm other-atom))
                    (not (chem:element-is-hydrogen other-element)))
           (let ((outgoing-edge (let ((oe-pos (position-if (lambda (oe) (and (eq (atm oe) other-atom)
                                                                             (eq (other-atom oe) atm)))
                                                           (outgoing-edges dfs))))
                                  (if oe-pos
                                      (aref (outgoing-edges dfs) oe-pos)
                                      (let ((oe (make-instance 'outgoing-edge
                                                               :atm atm
                                                               :other-atom other-atom
                                                               :bond-order order
                                                               :bond-index (length (outgoing-edges dfs)))))
                                        (vector-push-extend oe (outgoing-edges dfs))
                                        oe)))))
             (push (cons (gethash other-atom (ranks dfs)) outgoing-edge)
                   unsorted-children)))))
     atm)
    (let ((sorted-children (sort unsorted-children #'< :key #'car)))
      (loop for outgoing-edge in (mapcar #'cdr sorted-children)
            for idx = (bond-index outgoing-edge)
            do (cond
                 ((not (gethash (other-atom outgoing-edge) (seen dfs)))
                  (let ((idx (position outgoing-edge (outgoing-edges dfs))))
                    (setf (bond-index outgoing-edge) idx))
                  #+(or)(vector-push-extend outgoing-edge (outgoing-edges node))
                  (when (not (gethash (other-atom outgoing-edge) (seen dfs)))
                    (dbg "About to dfs ~s ~s~%" atm (other-atom outgoing-edge))
                    (depth-first-search-ring-finder-impl atm (other-atom outgoing-edge) dfs)))
                 (t (setf (gethash idx (closure-bonds dfs)) t)))
            do (vector-push-extend outgoing-edge (outgoing-edges node))
            ))))


(defun depth-first-smiles (dfs)
  (let ((smiles-terms (make-array 16 :adjustable t :fill-pointer 0))
        (stack (list (cons :0extend (list 0 -1))))) ; stack = [(0, (start_index, -1))]
    (loop
      while stack
      do (block outer-while
           (let* ((pop (pop stack))
                  bond-smarts
                  (action (car pop))
                  (data (cdr pop)))
             (case action
               (:0extend
                ;;# Add an atom.

                ;;# The 'while 1:' emulates a goto for the special case
                ;;# where the atom is connected to only one other atom.  I
                ;;# don't need to use the stack machinery for that case, and
                ;;# can speed up this function by about 10%.
                (loop named while
                      do (let* ((num-neighbors 0)
                                (atom-idx (first data))
                                (prev-bond-idx (second data))
                                (camgen-node (aref (camgen-nodes dfs) atom-idx))
                                (atm (atm camgen-node)))
                           (vector-push-extend (atom-smarts atm) smiles-terms) ; smiles_terms.append(cangen_nodes[atom_idx].atom_smarts)
                           (let ((outgoing-edges (outgoing-edges camgen-node))) ; outgoing_edges = cangen_nodes[atom_idx].outgoing_edges
                             (dbg "atm ~s~%" atm)
                             (loop for outgoing-edge across outgoing-edges
                                   for bond-idx = (bond-index outgoing-edge) ;; bond_idx = outgoing_edge.bond_index
                                                                             ;;# Is this a ring closure bond?
                                   do (dbg "    ~s   bond-idx ~d~%" outgoing-edge bond-idx)
                                   do (block continue
                                        (if (gethash bond-idx (closure-bonds dfs))
                                            ;; # Have we already seen it before?
                                            (progn
                                              (dbg "unclosed-closures ~s~%" (alexandria:hash-table-alist (unclosed-closures dfs)))
                                              (if (not (gethash bond-idx (unclosed-closures dfs)))
                                                  ;;# This is new. Add as a ring closure.
                                                  (let ((closure (pop (available-closures dfs))))
                                                    (dbg "Starting closure ~s~%" closure)
                                                    (vector-push-extend (get-closure-label (bond-smarts outgoing-edge) closure) smiles-terms)
                                                    (setf (gethash bond-idx (unclosed-closures dfs)) closure))
                                                  (let ((closure (gethash bond-idx (unclosed-closures dfs))))
                                                    (dbg "Ending closure ~s~%" closure)
                                                    (vector-push-extend (get-closure-label (bond-smarts outgoing-edge) closure) smiles-terms)
                                                    (push closure (available-closures dfs))
                                                    (remhash bond-idx (unclosed-closures dfs)))))
                                            ;; # This is a new outgoing bond.
                                            (progn
                                              (when (= bond-idx prev-bond-idx)
                                                ;; Don't go backwards along the bond I just came in on
                                                (return-from continue nil))
                                              ;; # This is the first bond. There's a good chance that
                                              ;; # it's the only bond.
                                              (if (= 0 num-neighbors)
                                                  (setf data (list (other-node-index outgoing-edge) bond-idx)
                                                        bond-smarts (bond-smarts outgoing-edge))
                                                  ;;# There are multiple bonds. Can't shortcut.
                                                  (progn
                                                    (when (= num-neighbors 1)
                                                      ;;# Capture the information for the first bond
                                                      ;;# This direction doesn't need the (branch) characters.
                                                      (push (cons :0extend data) stack)
                                                      (push (cons :1no-paren-bond bond-smarts) stack))
                                                    (push (cons :3close nil) stack)
                                                    (push (cons :0extend (list (other-node-index outgoing-edge) bond-idx)) stack)
                                                    (push (cons :4open (bond-smarts outgoing-edge)) stack)))
                                              (incf num-neighbors))))))
                           ;;# If there's only one item then goto action==0 again.
                           (when (/= num-neighbors 1)
                             (dbg "Returning from while stack = ~s~%" stack)
                             (return-from while nil))
                           (vector-push-extend bond-smarts smiles-terms))))
               (:1no-paren-bond;; Process a bond which does not need '()'s
                (vector-push-extend data smiles-terms))
               (:3close
                (vector-push-extend ")" smiles-terms))
               (:4open
                (dbg "cmd 4 data -> ~s~%" data)
                (vector-push-extend (format nil "(~a" data) smiles-terms)) ; 'data' is bond_smarts
               (otherwise
                (error "Bad code"))))))
    smiles-terms))


(defun generate (mol)
  (check-type mol chem:molecule)
  (let* ((chem:*current-rings* (chem:identify-rings mol)))
    (chem:with-aromaticity-information (mol :am1bcc)
      (let* ((ranks (cangen mol))
             (root-atm (block find-root
                         (maphash (lambda (atm rank)
                                    (when (= rank 1)
                                      (return-from find-root atm)))
                                  ranks)))
             (dfs (make-instance 'depth-first-search-builder
                                 :ranks ranks)))
        (depth-first-search-ring-finder-impl nil root-atm dfs)
        (setf (root-atom dfs) root-atm)
        (calculate-other-node-index dfs)
        (format nil "~a" (apply 'concatenate 'string (coerce (depth-first-smiles dfs) 'list)))
        ))))

