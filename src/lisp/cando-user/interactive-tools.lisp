
(in-package :cando-user)

(defparameter *energy-function* nil)
(defparameter *picked* nil)
(defparameter *pick-aggregate* nil)

(defun identify-atom (atom-ht aggregate)
  (let ((xpos (gethash "x" atom-ht))
        (ypos (gethash "y" atom-ht))
        (zpos (gethash "z" atom-ht)))
    (cando:do-atoms (atm aggregate)
      (let ((deltax (abs (- xpos (geom:vx (chem:get-position atm)))))
            (deltay (abs (- ypos (geom:vy (chem:get-position atm)))))
            (deltaz (abs (- zpos (geom:vz (chem:get-position atm))))))
        (when (and (< deltax 0.01)
                   (< deltay 0.01)
                   (< deltaz 0.01))
          (return-from identify-atom atm))))
    nil))

(defun log-interactions (picked-atoms energy-function)
  (cond
    ((and (= (length picked-atoms) 4)
          (first picked-atoms)
          (second picked-atoms)
          (third picked-atoms)
          (fourth picked-atoms))
     (destructuring-bind (a1 a2 a3 a4)
         picked-atoms
       (let ((term (chem:get-dihedral-component energy-function))
             (atom-table (chem:atom-table energy-function))
             (dihedral-terms (chem:lookup-dihedral-terms atom-table a1 a2 a3 a4)))
         (format t "dihedral-terms ~s~%" dihedral-terms))

    ))))

(defun stage-pick (inst data)
  (declare (ignore inst))
  (let ((control-key (gethash "ctrl_key" data)))
    (if control-key
      (let ((atm-ht (gethash "atom" data)))
        (when atm-ht
          (let ((atm (identify-atom atm-ht *pick-aggregate*)))
            (when (> (length *picked*) 3)
              (rplacd (nthcdr 2 *picked*) nil))
            (push atm *picked*)
            (format t "Added to selection: ~s~%" *picked*))
          (log-interactions *picked* *energy-function*)))
      (progn
        (format t "Clearing selection~%")
        (setf *picked* nil)))))

(defun force-field-viewer (aggregate)
  (let* ((energy-function (chem:make-energy-function :matter aggregate :assign-types t))
         (view (cando-user:show aggregate))
         (stage (find-if (lambda (wid) (typep wid 'ngl:stage)) (jupyter-widgets:widget-children view))))
    (setf *energy-function* energy-function)
    (setf *pick-aggregate* aggregate)
    (ngl:on-stage-pick stage #'stage-pick)
    view))
    
