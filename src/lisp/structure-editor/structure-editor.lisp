(in-package :structure-editor)

(defun json-lookup (json key)
  (unless (eq (car json) :obj)
    (error "The json must start with :obj but it doesn't json: ~s" json))
  (loop for entry in (cdr json)
        when (string= (car entry) key)
          do (return (cdr entry))))

(defun describe-json (json)
  (loop for entry in json
        do (format t "entry: ~s~%" entry)))

(defun json-bond-order (jbo)
  (case jbo
    (1 :single-bond)
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
         (intern (format nil "~a~a" element (1+ index)) :keyword)
         (error "Kekule.Atom doesn't have a valid element ~s~%" element)))
    ((string= type "Kekule.Pseudoatom")
     (intern symbol :keyword))
    (t (error "Cannot calculate name ~a" type))))

(defun parse-kekule-json-molecule (json)
  (let* ((molecule (chem:make-molecule nil))
         (residue (chem:make-residue nil))
         (name-string (json-lookup json "id"))
         (name (intern name-string :keyword))
         (coord-2d (json-lookup json "coord2D"))
         (coord-x (if coord-2d
                      (json-lookup coord-2d "x")
                      0.0))
         (coord-y (if coord-2d
                      (json-lookup coord-2d "y")
                      0.0))
         (annotation (json-lookup json "annotation"))
         (ctab (json-lookup json "ctab"))
         (jnodes (json-lookup ctab "nodes"))
         (jconnectors (json-lookup ctab "connectors"))
         (atoms (make-array (length jnodes))))
    (chem:add-matter molecule residue)
    (when annotation
      (chem:set-property molecule :annotation annotation))
    (loop for jnode in jnodes
          for type = (json-lookup jnode "__type__")
          for symbol = (json-lookup jnode "symbol")
          for jid = (json-lookup jnode "id")
          for index from 0
          for jcoord2d = (json-lookup jnode "coord2D")
          for jcharge = (json-lookup jnode "charge")
          for charge = jcharge
          for jparity = (json-lookup jnode "parity")
          for jisotope = (json-lookup jnode "isotopeId")
          for jannotation = (json-lookup jnode "annotation")
          for element = (calculate-element jisotope)
          for id = (intern jid :keyword)
          for name = (calculate-name type symbol element index)
          for atom = (chem:make-atom name element)
          for xp = (if jcoord2d
                       (float (json-lookup jcoord2d "x"))
                       0.0)
          for yp = (if jcoord2d
                       (float (json-lookup jcoord2d "y"))
                       0.0)
          do (chem:set-position atom (geom:vec (+ xp coord-x) (+ yp coord-y) 0.0))
          do (chem:set-charge atom charge)
          do (when jannotation
               (chem:set-property atom :annotation jannotation))
          do (setf (aref atoms index) atom)
          do (chem:add-matter residue atom))
    (loop for jedge in jconnectors
          for jbond-type = (json-lookup jedge "bondType")
          for jbond-order = (json-lookup jedge "bondOrder")
          for connected-objs = (json-lookup jedge "connectedObjs")
          for stereo = (json-lookup jedge "stereo")
          for from = (first connected-objs)
          for to = (second connected-objs)
          for from-atom = (aref atoms from)
          for to-atom = (aref atoms to)
          for bond-order = (json-bond-order jbond-order)
          for bond = (chem:bond-to from-atom to-atom bond-order)
          do (when stereo
               (chem:bond-set-property bond :stereo (list from-atom to-atom stereo))))
    molecule))

(defun parse-kekule-json (json)
  (let* ((root (json-lookup json "root"))
         (children (json-lookup root "children"))
         (jmolecules (json-lookup children "items")))
    (let ((agg (chem:make-aggregate nil)))
      (loop for jmolecule in jmolecules
            for mol = (parse-kekule-json-molecule jmolecule)
            do (chem:add-matter agg mol)
            )
      agg)))

