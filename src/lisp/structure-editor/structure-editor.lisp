(in-package :structure-editor)


(defun summarize-stereochemistry (aggregate)
  (cando:do-molecules (mol aggregate)
    (format t "Molecule: ~a~%" mol)
    (cando:do-residues (res mol)
      (cando:do-atoms (atm res)
        (when (eq (chem:get-stereochemistry-type atm) :chiral)
          (format t "~a  ~a ~a~%" atm (chem:get-stereochemistry-type atm) (chem:get-configuration atm))
          (format t "~a~%" (chem:bonds-as-list atm)))))
    (chem:map-bonds nil
                    (lambda (a1 a2 o bb)
                      (declare (ignore bb))
                      (when (member o '(:single-wedge-begin :single-wedge-end :single-hash-begin :single-hash-end))
                        (format t "~a -> ~a  ~a~%" a1 a2 o)))
                    mol)))

(defun json-bond-order (edge)
  (let ((bond-order (gethash "bondOrder" edge))
        (stereo (gethash "stereo" edge 0)))
    (case bond-order
      (1 (case stereo
           (0 :single-bond)
           (1 :single-wedge-begin)
           (3 :single-hash-begin)
           (otherwise (error "Illegal bond-order stereo (~a/~a) combination" bond-order stereo))))
      (2 :double-bond)
      (3 :triple-bond)
      (otherwise (warn "Unrecognized bond order ~a" bond-order)))))

(defun calculate-element (jisotope)
  (cond
    ((null jisotope)
     :du)
    ((stringp jisotope)
     (intern jisotope :keyword))
    (t (warn "Could not calculate element using ~s" jisotope))))

(defun calculate-name (type symbol element index)
  (cond
    ((string= type "Kekule.Atom")
     (if element
         (intern (symbol-name (gensym (symbol-name element))) :keyword)
         (error "Kekule.Atom doesn't have a valid element ~s~%" element)))
    ((string= type "Kekule.Pseudoatom")
     (intern symbol :keyword))
    ((string= type "Kekule.SubGroup")
     (error "Handle Kekule.SubGroup symbol: ~s  element: ~s  index: ~d"
            symbol element index))
    (t (error "Cannot calculate name ~a" type))))


(defun set-properties-from-annotation (json matter)
  (loop with *package* = (find-package :keyword)
        for cur = (read-from-string (format nil "(~a)" (gethash "annotation" json ""))) then (cddr cur)
        for label = (first cur)
        for value = (second cur)
        while label do (chem:set-property matter label value)))


(defun parse-kekule-atom-pseudo-atom (jnode type atoms residue coord-x coord-y)
  (let* ((symbol (gethash "symbol" jnode))
         (jcoord2d (gethash "coord2D" jnode))
         (jisotope (gethash "isotopeId" jnode))
         (element (calculate-element jisotope))
         (atom (chem:make-atom (calculate-name type symbol element (length atoms)) element))
         (xp (if jcoord2d (float (gethash "x" jcoord2d)) 0.0))
         (yp (if jcoord2d (float (gethash "y" jcoord2d)) 0.0)))
    (chem:set-position atom (geom:vec (- (+ xp coord-x)) (+ yp coord-y) 0.0)) ; flip X coordinate for stereochem
    (chem:set-charge atom (gethash "charge" jnode))
    (set-properties-from-annotation jnode atom)
    (vector-push-extend atom atoms)
    (chem:add-matter residue atom)))

(defun atom-reference (atoms reference)
  (etypecase reference
    (integer
     (aref atoms reference))
    (cons
     (if (null (cdr reference))
         (aref atoms (car reference))
         (let ((subarray (aref atoms (car reference))))
           (atom-reference subarray (cdr reference)))))))
    
(defun parse-kekule-bonds (json atoms)
  (loop for jedge across (gethash "connectors" json #())
        for connected-objs = (gethash "connectedObjs" jedge)
        for from = (elt connected-objs 0)
        for to = (elt connected-objs 1)
        for from-atom = (atom-reference atoms from)
        for to-atom = (atom-reference atoms to)
        do (unless from-atom
             (error "Could not identify from-atom at ~a in ~a" from atoms))
        do (unless to-atom
             (error "Could not identify to-atom at ~a in ~a" to atoms))
        do (when (and from-atom to-atom)
             (chem:bond-to from-atom to-atom (json-bond-order jedge)))))

(defun parse-kekule-json-subgroups-and-atoms (json atoms residue coord-x coord-y molecule)
  (loop for jnode across (gethash "nodes" json #())
        for type = (gethash "__type__" jnode)
        do (progn
             (cond
               ((string= type "Kekule.SubGroup")
                (let ((subgroup-atoms (make-array 16 :adjustable t :fill-pointer 0)))
                  (vector-push-extend subgroup-atoms atoms)
                  (parse-kekule-json-group jnode subgroup-atoms molecule)))
               (t
                (parse-kekule-atom-pseudo-atom jnode type atoms residue coord-x coord-y))))))

(defun parse-kekule-json-group (json atoms molecule)
  (let* ((residue (chem:make-residue nil))
         (name-string (gethash "id" json))
         (name (intern name-string :keyword))
         (coord-2d (gethash "coord2D" json))
         (coord-x (if coord-2d
                      (gethash "x" coord-2d)
                      0.0))
         (coord-y (if coord-2d
                      (gethash "y" coord-2d)
                      0.0))
         (ctab (gethash "ctab" json)))
    (set-properties-from-annotation json molecule)
    (parse-kekule-json-subgroups-and-atoms ctab atoms residue coord-x coord-y molecule)
    (parse-kekule-bonds ctab atoms)
    (chem:add-matter molecule residue)))

(defun parse-kekule-json-molecule (json)
  (let ((molecule (chem:make-molecule nil))
        (atoms (make-array 256 :adjustable t :fill-pointer 0)))
    (parse-kekule-json-group json atoms molecule)
    molecule))


(defun parse-kekule-json (json)
  (loop with agg = (chem:make-aggregate nil)
        for jmolecule across (gethash "items" (gethash "children" (gethash "root" json)))
        do (chem:add-matter agg (parse-kekule-json-molecule jmolecule))
        finally (return agg)))

