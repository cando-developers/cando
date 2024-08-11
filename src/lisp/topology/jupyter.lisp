(in-package :topology-jupyter)

(defun sketch-svg (topology &rest args)
  "Generate an svg sketch of the topology - send args to sketch2d"
  (let ((mol (topology:build-one-molecule-for-topology topology)))
    (sketch2d:svg (apply 'sketch2d:sketch2d mol args))))

(defclass node ()
  ((label :initarg :label :reader label)
   (uid :initarg :uid :reader uid)
   (monomer :initarg :monomer :reader monomer)
   ))

(defclass edge ()
  ((label :initarg :label :reader label)
   (class :initarg :class :reader class)
   (coupling :initarg :coupling :reader coupling)
   (uid :initarg :uid :reader uid)
   (weight :initarg :weight :reader weight)
   (width :initarg :width :reader width)
   (source :initarg :source :reader source)
   (target :initarg :target :reader target)
   ))

(defgeneric make-node-label (monomer oligomer-or-space))


(defmethod make-node-label (monomer (oligomer-space topology:oligomer-space))
  (string-downcase
   (cond
     ((and (consp (topology:id monomer)) (symbolp (car (topology:id monomer))))
      (format nil "~a" (first (topology:id monomer))))
     ((consp (topology:id monomer))
      (format nil "~{~a~%~}" (first (topology:id monomer))))
     ((symbolp (topology:id monomer))
      (format nil "~a" (topology:id monomer)))
     (t (break "check monomer to generate name") (format nil "~a" monomer))))
  )

(defmethod make-node-label (monomer (oligomer topology:oligomer))
  (string-downcase (topology:oligomer-monomer-name oligomer monomer)))


(defun make-edge-label (coupling)
  (typecase coupling
    (topology:directional-coupling
     (format nil "~a" (string-downcase (topology:name coupling))))
    (topology:ring-coupling
     (format nil "RING~%~a~%~a" (string-downcase (topology:plug1 coupling)) (string-downcase (topology:plug2 coupling))))
    (t (break "Handle coupling ~s" coupling))))

(defgeneric make-element (element))
(defmethod make-element ((element node))
  (cytoscape:make-element :group "nodes" :data (j:make-object "id" (string (uid element))
                                                              "label" (label element)
                                                              )))

(defmethod make-element ((element edge))
  (cytoscape:make-element :group "edges" :data (j:make-object "id" (string (uid element))
                                                              "weight" (weight element)
                                                              "width" (width element)
                                                              "label" (label element)
                                                              "source" (string (uid (source element)))
                                                              "target" (string (uid (target element)))
                                                              )))

(defun make-weight (coupling)
  (typecase coupling
    (topology:directional-coupling
     "100")
    (topology:ring-coupling
     "1")
    (t (break "Handle coupling ~s" coupling))))

(defun make-width (coupling)
  (typecase coupling
    (topology:directional-coupling
     "6")
    (topology:ring-coupling
     "2")
    (t (break "Handle coupling ~s" coupling))))

(defun generate-elements (oligomer-or-space)
  (let ((monomer-to-node (make-hash-table))
        (rev-elements nil))
    (loop for monomer across (topology:monomers oligomer-or-space)
          for label = (make-node-label monomer oligomer-or-space)
          for uid = (gensym)
          for monomer-node = (make-instance 'node
                                            :uid uid
                                            :label label
                                            :monomer monomer)
          do (push monomer-node rev-elements)
          do (setf (gethash monomer monomer-to-node) monomer-node))
    (loop for coupling across (topology:couplings oligomer-or-space)
          for source-node = (gethash (topology:source-monomer coupling) monomer-to-node)
          for target-node = (gethash (topology:target-monomer coupling) monomer-to-node)
          for edge = (make-instance 'edge
                                    :label (make-edge-label coupling)
                                    :uid (gensym)
                                    :weight (make-weight coupling)
                                    :width (make-width coupling)
                                    :coupling coupling
                                    :source source-node
                                    :target target-node
                                    )
          for uid = (gensym)
          do (push edge rev-elements)
          )
    (loop for element in (nreverse rev-elements)
          collect (make-element element))))

(defun build-cytoscape (elements)
  (cytoscape:make-cytoscape-widget
   :box-selection-enabled nil
   :auto-unselectify t
   :layout (make-instance 'jw:layout :height "640px")
   :graph-layouts (list (cytoscape:make-fcose-layout :directed t :padding 10 :roots "#a"))
   :elements elements
   :graph-style "
 node { content: data(label);
        font-size: 12;
        text-valign: center;
        text-halign: center;
        text-wrap: wrap;
      }
 edge { content: data(label);
        curve-style: bezier;
        target-arrow-shape: triangle;
        width: data(width);
        line-color: #0dd;
        target-arrow-color: #d00;
        font-size: 10;
        text-wrap: wrap;
      }
 .highlighted {
         backgrond-color: #61bffc;
         line-color: #61bffc;
         target-arrow-color: #61bffc;
         transition-property: background-color, line-color, target-arrow-color;
         transition-duration: 0.5s;
 }"
   ))


(defmethod cando-widgets::show-on-pane (pane-instance (object topology:oligomer-space) &rest rest &key &allow-other-keys)
  (let ((elements (generate-elements object)))
    (build-cytoscape elements)))

(defmethod cando-widgets::show-on-pane (pane-instance (object topology:oligomer) &rest rest &key &allow-other-keys)
  (let ((elements (generate-elements object)))
    (build-cytoscape elements)))

(defmethod cando-widgets::show-on-pane (pane-instance (object topology:topology) &rest rest &key &allow-other-keys)
  (apply 'cando-widgets::show-on-pane pane-instance (topology-jupyter:sketch-svg object) rest))

(defmethod cando-widgets::show-on-pane (pane-instance (object symbol) &rest rest &key &allow-other-keys)
  (let ((top (chem:find-topology object)))
    (apply 'cando-widgets::show-on-pane pane-instance (topology-jupyter:sketch-svg top) rest)))


(defmethod cando-widgets::show-on-pane (pane-instance (object topology:oligomer-shape) &rest rest &key &allow-other-keys)
  (let* ((assembler (topology:make-assembler (list object)))
         (coords (topology:make-coordinates-for-assembler assembler)))
    (topology:fill-internals-from-oligomer-shape assembler object)
    (topology:build-external-coordinates assembler :coords coords)
    (cando-widgets::show-on-pane pane-instance (cando:mol (topology:aggregate* assembler coords) 0))))
