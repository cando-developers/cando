(in-package :topology.protein)

;;;
;;;ARG  -180 -180    10     1  2  2  1  0.249730    62.5   176.9   176.6    85.7       6.9    11.1    10.5     9.9
;;;     s   es   e                      s       es      es      es      es      e
;;;.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
;;;.         1         2         3         4         5         6         7         8         9        10        11
;;;

(defclass rotamer-values ()
  ((name :initarg :name :accessor name)
   (probability :initarg :probability :accessor probability)
   (phi :initarg :phi :accessor phi)
   (psi :initarg :psi :accessor psi)
   (chi1 :initarg :chi1 :accessor chi1)
   (chi2 :initarg :chi2 :accessor chi2)
   (chi3 :initarg :chi3 :accessor chi3)
   (chi4 :initarg :chi4 :accessor chi4)))

(defun parse-rotamer (line)
  (let ((name (intern (subseq line 0 3) :keyword))
        (phi-str (subseq line 5 9))
        (psi-str (subseq line 10 14))
        (probability-str (subseq line 37 45))
        (chi1-str (subseq line 46 53))
        (chi2-str (subseq line 54 61))
        (chi3-str (subseq line 62 69))
        (chi4-str (subseq line 70 77)))
    (let ((phi (parse-integer phi-str :junk-allowed t))
          (psi (parse-integer psi-str :junk-allowed t))
          (probability (fortran:parse-vecreal probability-str))
          (chi1 (fortran:parse-vecreal chi1-str))
          (chi2 (fortran:parse-vecreal chi2-str))
          (chi3 (fortran:parse-vecreal chi3-str))
          (chi4 (fortran:parse-vecreal chi4-str)))
      (make-instance 'rotamer-values
                     :name name
                     :phi phi
                     :psi psi
                     :probability probability
                     :chi1 chi1
                     :chi2 chi2
                     :chi3 chi3
                     :chi4 chi4))))


(defun read-rotamer-file (filename)
  (let ((rotamers (make-hash-table)))
    (with-open-file (fin filename :direction :input)
      (loop named reader
            for line = (read-line fin nil :eof)
            do (cond
                 ((eq line :eof)
                  (return-from reader rotamers))
                 ((char= (elt line 0) #\#)
                  nil)
                 (t
                  (let* ((rotamer-values (parse-rotamer line))
                         (aa-name (name rotamer-values))
                         (residue-rotamers (alexandria:ensure-gethash aa-name rotamers (make-hash-table :test 'equal))))
                    (let* ((backbone-key (cons (phi rotamer-values) (psi rotamer-values)))
                           (rotamers (alexandria:ensure-gethash backbone-key residue-rotamers (make-array 16 :adjustable t :fill-pointer 0))))
                      (vector-push-extend rotamer-values rotamers)))))))
    rotamers))

(defun find-and-clone (name vector &optional key dict)
  (unless vector
    (break "Could not find ~a in ~a" key dict))
  (let* ((pos (position name vector :key #'topology:name))
         (internal (elt vector pos))
         (new-internal (topology:copy-internal internal)))
    (setf (elt vector pos) new-internal)
    new-internal))

(defun adjust-phi-psi (first-fragment phi-deg psi-deg)
  (error "deprecated")
  #+(or)(format t "adjust-phi-psi ~a ~a~%" phi-deg psi-deg)
  #+(or)
  (let* ((new-fragment (topology:copy-fragment-internals first-fragment))
         (phi-c (find-and-clone :c (topology:internals new-fragment)))
         (phi-ha (find-and-clone :ha (topology:internals new-fragment)))
         (oofi (topology:out-of-focus-internals new-fragment))
         (side (gethash :-side oofi))
         (phi-cb (find-and-clone :cb side :-side oofi))
         (out-aa (gethash :-amide oofi))
         (psi-n (find-and-clone :n out-aa :-amide oofi))
         (psi-o (find-and-clone :o (topology:internals new-fragment)))
         (delta-phi (topology:degree-difference (topology:rad-to-deg (topology:dihedral-rad phi-c)) phi-deg))
         (delta-psi (topology:degree-difference (topology:rad-to-deg (topology:dihedral-rad psi-n)) psi-deg)))
    #+(or)
    (progn
      (format t "---- target phi deg = ~a~%" phi-deg)
      (format t "---- target psi deg = ~a~%" psi-deg)
      (format t "before phi = ~a  psi = ~a~%" phi-c psi-n)
      (format t "delta-phi deg = ~a~%" delta-phi)
      (format t "delta-psi deg = ~a~%" delta-psi))
    (topology:radians-incf (topology:dihedral-rad phi-c) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad phi-cb) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad phi-ha) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad psi-n) (topology:deg-to-rad delta-psi))
    (topology:radians-incf (topology:dihedral-rad psi-o) (topology:deg-to-rad delta-psi))
    #+(or)(format t "after  phi = ~a  psi = ~a~%" phi-c psi-n)
    new-fragment))


#|
(defun adjust-rotamer (first-fragment dihedrals topology)
 (let ((new-fragment (topology:copy-fragment-internals first-fragment)
  (

|#

(defun get-dihedral-atom-names (to-topology dihedral-name name-to-parent)
  (error "We should use the constituion residue-properties rather than atom :dihedral properties")
  (loop for constitution-atom across (topology:constitution-atoms (topology:constitution to-topology))
        for ca-name = (topology:atom-name constitution-atom)
        for properties = (topology:properties constitution-atom)
        for dihedral-property = (getf properties :dihedral)
        when (eq dihedral-name dihedral-property)
          do (return (let* ((parent-joint-template (gethash ca-name name-to-parent)))
                       (unless parent-joint-template
                         (break "nil parent-joint-template for ~a ~a" ca-name name-to-parent))
                       (when parent-joint-template
                         (let* ((children (copy-seq (topology:children parent-joint-template)))
                                (sorted-children (sort children
                                                       (lambda (n1 n2) 
                                                         (declare (ignore n2))
                                                         (eq n1 ca-name))
                                                       :key #'topology:atom-name)))
                           (mapcar #'topology:atom-name sorted-children )))))))

(defun adjust-chi (first-fragment phi-deg psi-deg)
  #+(or)(format t "adjust-phi-psi ~a ~a~%" phi-deg psi-deg)
  (error "deprecated")
  #+(or)
  (let* ((new-fragment (topology:copy-fragment-internals first-fragment))
         (phi-c (find-and-clone :c (topology:internals new-fragment)))
         (phi-ha (find-and-clone :ha (topology:internals new-fragment)))
         (oofi (topology:out-of-focus-internals new-fragment))
         (side (gethash :-side oofi))
         (phi-cb (find-and-clone :cb side :-side oofi))
         (out-aa (gethash :-amide oofi))
         (psi-n (find-and-clone :n out-aa :-amide oofi))
         (psi-o (find-and-clone :o (topology:internals new-fragment)))
         (delta-phi (topology:degree-difference (topology:rad-to-deg (topology:dihedral-rad phi-c)) phi-deg))
         (delta-psi (topology:degree-difference (topology:rad-to-deg (topology:dihedral-rad psi-n)) psi-deg)))
    #+(or)
    (progn
      (format t "---- target phi deg = ~a~%" phi-deg)
      (format t "---- target psi deg = ~a~%" psi-deg)
      (format t "before phi = ~a  psi = ~a~%" phi-c psi-n)
      (format t "delta-phi deg = ~a~%" delta-phi)
      (format t "delta-psi deg = ~a~%" delta-psi))
    (topology:radians-incf (topology:dihedral-rad phi-c) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad phi-cb) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad phi-ha) (topology:deg-to-rad delta-phi))
    (topology:radians-incf (topology:dihedral-rad psi-n) (topology:deg-to-rad delta-psi))
    (topology:radians-incf (topology:dihedral-rad psi-o) (topology:deg-to-rad delta-psi))
    #+(or)(format t "after  phi = ~a  psi = ~a~%" phi-c psi-n)
    new-fragment))

(defun maybe-adjust-dihedral (to-topology dihedral-name name-to-parent dup-fragment-internals target-dihedral-deg)
  (let ((atom-names (get-dihedral-atom-names to-topology dihedral-name name-to-parent)))
    (when atom-names
      (let* ((internals (mapcar (lambda (name) (find-and-clone name (topology:internals dup-fragment-internals))) atom-names))
             (current-dihedral (topology:rad-to-deg (topology:dihedral-rad (elt internals 0))))
             (delta-deg (topology:angle-sub target-dihedral-deg current-dihedral))
             (new-dihedral (topology:degrees-add current-dihedral delta-deg)))
        #+(or)(format t "maybe-adjust-dihedral current-dihedral ~a ~a ~a delta-deg ~a -> new-dihedral ~a target: ~a~%"
                (topology:name to-topology)
                dihedral-name
                current-dihedral
                delta-deg
                new-dihedral
                target-dihedral-deg)
        (when (> (abs (topology:angle-sub new-dihedral target-dihedral-deg)) 0.1)
          (error "Bad dihedral calculation ~a vs ~a" new-dihedral target-dihedral-deg))
        (setf (topology:dihedral-rad (elt internals 0)) new-dihedral)))))

(defun lookup-rotamers (rotamers name)
  (let* ((key (cond
                ((member (string name) '("LYN") :test #'string=) :lys)
                ((member (string name) '("HIE" "HID") :test #'string=) :his)
                ((member (string name) '("GLH") :test #'string=) :glu)
                ((member (string name) '("ASH") :test #'string=) :asp)
                (t (intern (string name) :keyword))))
         (rots (gethash key rotamers)))
    (unless rots
      (error "Could not find (key: ~s) in rotamers keys ~s" key (alexandria:hash-table-keys rotamers)))
    rots))

(defun adjust-backbone-dependent-rotamers (from-monomer-context from-oligomer from-focus-monomer
                                           to-monomer-context to-oligomer to-focus-monomer
                                           rotamers protein-conformations-map)
  (let* ((from-name (topology:oligomer-monomer-name-for-monomer from-oligomer from-focus-monomer))
         (to-name (topology:oligomer-monomer-name-for-monomer to-oligomer to-focus-monomer))
         (from-fragments (gethash from-monomer-context (topology:monomer-context-to-context-rotamers protein-conformations-map)))
         (to-fragments (gethash to-monomer-context (topology:monomer-context-to-context-rotamers protein-conformations-map)))
         (new-from-fragment-vector (make-array 256 :adjustable t :fill-pointer 0))
         (new-to-fragment-vector (make-array 256 :adjustable t :fill-pointer 0))
         (from-fragment-template (elt (topology:fragments from-fragments) 0))
         (to-fragment-template (elt (topology:fragments to-fragments) 0))
         (to-topology (chem:find-topology to-name))
         (name-to-parent (let ((name-to-parent (make-hash-table)))
                           (topology:walk-joint-template (topology:joint-template to-topology)
                                                         (lambda (joint-template)
                                                           (setf (gethash (topology:atom-name joint-template)
                                                                          name-to-parent)
                                                                 (topology:parent joint-template))))
                           name-to-parent))
         (max-backbone-indices (* 36 36))
         (indices (make-array max-backbone-indices))
         (backbone-index -1))
    (cond
      ((string= (string to-name) "ALA")
       (warn "What do we do with ALA???????"))
      (t (let* ((to-rotamers (lookup-rotamers rotamers to-name)))
           (format t "adjust-backbone-dependent-rotamers ~a -> ~a :topology ~a~%" from-name to-name to-topology)
           (loop for iphi from -180 below 180 by 10
                 for phi-rad = (float iphi 1.0s0)
                 do (loop for ipsi from -180 below 180 by 10
                          for psi-rad = (float ipsi 1.0s0)
                          for new-fragment = (adjust-phi-psi from-fragment-template phi-rad psi-rad)
                          ;; from-index is the from monomer-context index
                          for from-index = (vector-push-extend new-fragment new-from-fragment-vector)
                          for side-chain-rotamers = (gethash (cons iphi ipsi) to-rotamers)
                          for backbone-index-vector = (make-array 4 :adjustable t :fill-pointer 0)
                          do (incf backbone-index)
                          do (loop for side-chain-rotamer across side-chain-rotamers
                                   for dup-fragment-internals = (topology:copy-fragment-internals to-fragment-template)
                                   do (setf (topology:probability dup-fragment-internals) (probability side-chain-rotamer))
                                   do (progn
                                        (maybe-adjust-dihedral to-topology :chi1 name-to-parent dup-fragment-internals (chi1 side-chain-rotamer))
                                        (maybe-adjust-dihedral to-topology :chi2 name-to-parent dup-fragment-internals (chi2 side-chain-rotamer))
                                        (maybe-adjust-dihedral to-topology :chi3 name-to-parent dup-fragment-internals (chi3 side-chain-rotamer))
                                        (maybe-adjust-dihedral to-topology :chi4 name-to-parent dup-fragment-internals (chi4 side-chain-rotamer)))
                                   do (let ((index (vector-push-extend dup-fragment-internals new-to-fragment-vector)))
                                        (vector-push-extend index backbone-index-vector)))
                          unless (< backbone-index max-backbone-indices)
                            do (error "backbone-index >= max-backbone-indices ~a" max-backbone-indices)
                          do (setf (elt indices backbone-index) backbone-index-vector))))))
    (values new-from-fragment-vector new-to-fragment-vector indices)))

(defun apply-backbone-rotamers (foldamer rotamers protein-conformations-map)
  "Loop through every backbone residue and generate context-rotamers for every pair of phi/psi angles at 10 degree increments"
  (let ((context-pairs (foldamer:all-matching-monomer-contexts foldamer))
        (monomer-context-to-context-rotamers (make-hash-table))
        (rotamer-context-connections (make-rotamer-context-connections)))
    (loop for monomer-contexts in context-pairs
          for from-monomer-context = (car monomer-contexts)
          for to-monomer-context = (cdr monomer-contexts)
          do (multiple-value-bind (from-oligomer from-focus-monomer)
                 (foldamer:find-oligomer-for-monomer-context foldamer from-monomer-context)
               (multiple-value-bind (to-oligomer to-focus-monomer)
                   (foldamer:find-oligomer-for-monomer-context foldamer to-monomer-context)
                 (let ((from-name (topology:oligomer-monomer-name-for-monomer from-oligomer from-focus-monomer)))
                   (cond
                     ((and (string= (string from-name) "AA")
                           (string= (string (topology:id to-focus-monomer)) "AASIDE"))
                      (multiple-value-bind (from-fragments to-fragments indices)
                          (adjust-backbone-dependent-rotamers from-monomer-context from-oligomer from-focus-monomer
                                                              to-monomer-context to-oligomer to-focus-monomer
                                                              rotamers protein-conformations-map)
                        ;; set the monomer-context-to-context-rotamers for the;
                        ;; from-monomer-context -> from-fragments
                        ;; to-monomer-context -> to-fragments
                        ;; and the rotamer-context-connections object
                        ;; (cons from-monomer-context to-monomer-context) -> indices
                        (setf (gethash from-monomer-context monomer-context-to-context-rotamers) from-fragments
                              (gethash to-monomer-context monomer-context-to-context-rotamers) to-fragments
                              )
                        (set-rotamer-context-connections rotamer-context-connections from-monomer-context to-monomer-context indices)
                        ))
                     ((string= (string from-name) "AA")
                      (let ((to-name (topology:oligomer-monomer-name-for-monomer to-oligomer to-focus-monomer)))
                        (format t "Handle in next pass ~a -> ~a~%" from-name to-name)))
                     ((string= (string from-name) "AP")
                      (let ((to-name (topology:oligomer-monomer-name-for-monomer to-oligomer to-focus-monomer)))
                        (format t "Ignore in next pass ~a -> if PRO ~a~%" from-name to-name)))
                     ((string= (string from-name) "ACE")
                      (let ((to-name (topology:oligomer-monomer-name-for-monomer to-oligomer to-focus-monomer)))
                        #+(or)(format t "Handle in next pass ~a -> ~a~%" from-name to-name)))
                     (t
                      (break "check from-focus-monomer ~a to-focus-monomer ~a" from-focus-monomer to-focus-monomer))
                     ))
                 )))
    (make-instance 'topology:connected-rotamers-map
                   :monomer-context-to-context-rotamers monomer-context-to-context-rotamers
                   :rotamer-context-connections rotamer-context-connections)
    ))


(defun extract-conformations-apply-rotamers (&key rotamers (path foldamer::*conformations-path*) foldamer (verbose t))
  (format t "Extracting conformations for the path: ~a~%" path)
  (let ((protein-conformations-map (foldamer:extract-context-rotamers-map "./")))
    (format t "Matching fragment conformations~%")
    #+(or)(let ((matched (protein:optimize-context-rotamers-map protein-conformations-map protein verbose)))
      (loop while (> (protein-report-conformations matched) 0)
            do (format t "Eliminating missing matches~%")
            do (setf matched (protein-eliminate-missing-rotamer-context-connections matched protein)))
      (format t "Saving ~a~%" path)
      (save-protein-conformations-map matched path)
      (values matched protein-conformations-map)
            )
    protein-conformations-map))

(defun residue-name-to-graph (residue-name foldamer)
  (cond
    ((string= "GLY" (string residue-name))
     (list :amide :gly))
    ((string= "PRO" (string residue-name))
     `(:amide :ap (:side :pro)))
    ((string= "LYS" (string residue-name))
     `(:amide :aa (:side :LYN)))
    ((string= "CYS" (string residue-name))
     `(:amide :aa (:side :LYN)))
    ((string= "NMET" (string residue-name))
     (let ((side-chain-name ':MET))
       `(:naa (:side ,side-chain-name))))
    (t
     (let ((side-chain-name (intern (string residue-name) :keyword)))
       `(:amide :aa (:side ,side-chain-name))))))

(defun amber-protein-to-oligomer-space (molecule foldamer)
  (let ((tree (loop for seq-index below (chem:content-size molecule)
                    for res = (chem:content-at molecule seq-index)
                    for res-name = (chem:get-name res)
                    for part = (residue-name-to-graph res-name foldamer)
                    append part
                    )))
    (topology:make-oligomer-space foldamer tree)))
