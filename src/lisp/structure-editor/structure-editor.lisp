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
                    (lambda (a1 a2 o)
                      (when (member o '(:single-wedge-begin :single-wedge-end :single-hash-begin :single-hash-end))
                        (format t "~a -> ~a  ~a~%" a1 a2 o)))
                    mol)))

(defun json-bond-order (jbo stereo)
  (case jbo
    (1 (cond
         ((null stereo) :single-bond)
         ((= stereo 1) :single-wedge-begin)
         ((= stereo 3) :single-hash-begin)
         (t (error "Illegal bond-order stereo (~a/~a) combination" jbo stereo))))
    (2 :double-bond)
    (3 :triple-bond)
    (otherwise (warn "Unrecognized bond order ~a" jbo))
    ))

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

(defun property-list-from-annotation (annotation)
  (let ((*package* (find-package :keyword)))
    (with-input-from-string (sin (format nil "(~a)" annotation))
      (let ((plist (read sin)))
        plist))))

(defun parse-kekule-atom-pseudo-atom (jnode type atoms residue coord-x coord-y)
  (let* ((symbol (gethash "symbol" jnode))
         (jid (gethash "id" jnode))
         (jcoord2d (gethash "coord2D" jnode))
         (jcharge (gethash "charge" jnode))
         (charge jcharge)
         (jparity (gethash "parity" jnode))
         (jisotope (gethash "isotopeId" jnode))
         (jannotation (gethash "annotation" jnode))
         (element (calculate-element jisotope))
         (id (intern jid :keyword))
         (name (calculate-name type symbol element (length atoms)))
         (atom (chem:make-atom name element))
         (xp (if jcoord2d (float (gethash "x" jcoord2d)) 0.0))
         (yp (if jcoord2d (float (gethash "y" jcoord2d)) 0.0))
         )
    (chem:set-position atom (geom:vec (* -1.0 (+ xp coord-x)) (+ yp coord-y) 0.0)) ; flip X coordinate for stereochem
    (chem:set-charge atom charge)
    (when (stringp jannotation)
      (let ((property-list (property-list-from-annotation jannotation)))
        (loop for cur = property-list then (cddr cur)
              for label = (first cur)
              for value = (second cur)
              while label do (chem:set-property atom label value))))
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
    
(defun parse-kekule-bonds (jconnectors atoms)
  (loop for jedge in jconnectors
        for jbond-type = (gethash "bondType" jedge)
        for jbond-order = (gethash "bondOrder" jedge)
        for connected-objs = (gethash "connectedObjs" jedge)
        for stereo = (gethash "stereo" jedge)
        for from = (first connected-objs)
        for to = (second connected-objs)
        for from-atom = (atom-reference atoms from)
        for to-atom = (atom-reference atoms to)
        for bond-order = (json-bond-order jbond-order stereo)
        do (unless from-atom
             (error "Could not identify from-atom at ~a in ~a" from atoms))
        do (unless to-atom
             (error "Could not identify to-atom at ~a in ~a" to atoms))
        do (when (and from-atom to-atom)
             (chem:bond-to from-atom to-atom bond-order))))

(defun parse-kekule-json-subgroups-and-atoms (jnodes jconnectors atoms residue coord-x coord-y molecule)
  (loop for jnode in jnodes
        for type = (gethash "__type__" jnode)
        do (progn
             (cond
               ((string= type "Kekule.SubGroup")
                (let ((subgroup-atoms (make-array 16 :adjustable t :fill-pointer 0)))
                  (vector-push-extend subgroup-atoms atoms)
                  (parse-kekule-json-group jnode subgroup-atoms molecule)))
               (t
                (parse-kekule-atom-pseudo-atom jnode type atoms residue coord-x coord-y)
                ))))
  (parse-kekule-bonds jconnectors atoms))

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
         (annotation (gethash "annotation" json))
         (ctab (gethash "ctab" json))
         (jnodes (gethash "nodes" json))
         (jconnectors (gethash "connectors" ctab)))
    (when (stringp annotation)
      (let ((property-list (property-list-from-annotation annotation)))
        (loop for cur = property-list then (cddr cur)
              for label = (first cur)
              for value = (second cur)
              while cur do (chem:set-property molecule label value))))
    (parse-kekule-json-subgroups-and-atoms jnodes jconnectors atoms residue coord-x coord-y molecule)
    (chem:add-matter molecule residue)))

(defun parse-kekule-json-molecule (json)
  (let ((molecule (chem:make-molecule nil))
        (atoms (make-array 256 :adjustable t :fill-pointer 0)))
    (parse-kekule-json-group json atoms molecule)
    molecule))


(defun parse-kekule-json (json)
  (let* ((root (gethash "root" json))
         (children (gethash "children" root))
         (jmolecules (gethash "items" children)))
    (let ((agg (chem:make-aggregate nil)))
      (loop for jmolecule in jmolecules
            for mol = (parse-kekule-json-molecule jmolecule)
            do (chem:add-matter agg mol)
            )
      agg)))


