;;
;; Recognize aromaticity in Matter
;;
;;

;(print "Loading aromaticity")
;(include "CANDO/tools/common.l")


;(make-package :aromaticity)
(in-package :aromaticity)
;(warn "(use-package 'core)")
;(warn "(use-package 'u)")

;(warn "(export 'identify-aromatic-rings
;        'is-aromatic
;        )
;")


;;
;;
;; Describe ring membership in a string
;;
(defun ring-membership-as-string (a)
  (let ((member (format t  "~a~a~a~a~a~%"
                   (if (chem:in-ring-size a 3) "3" "")
                   (if (chem:in-ring-size a 4) "4" "")
                   (if (chem:in-ring-size a 5) "5" "")
                   (if (chem:in-ring-size a 6) "6" "")
                   (if (chem:in-ring-size a 7) "7" ""))))
    (if (/= member "") (format t "ring~a~%" member) "")))


(defstruct aromatic-info aromaticity-type rule-name)


;; Return the aromaticity type of the atom followed by the rule-name or (values nil nil)
(defun is-aromatic (atom &optional (aromaticity-info chem:*current-aromaticity-information*))
  (gethash atom aromaticity-info))

(defun is-aromatic-type (atom aromatic-type &optional (aromaticity-info chem:*current-aromaticity-information*))
  (loop named test
        for info in (gethash atom aromaticity-info)
        when (eq (aromatic-info-aromaticity-type info) aromatic-type)
          do (return-from is-aromatic-type t))
  (null aromatic-type))

(defun is-aromatic-rule (atom aromatic-rule &optional (aromaticity-info chem:*current-aromaticity-information*))
  (loop named test
        for info in (gethash atom aromaticity-info)
        when (eq (aromatic-info-rule-name info) aromatic-rule)
          do (return-from is-aromatic-rule t))
  (null aromatic-rule))

;; Set the aromatic type of an atom
(defun set-aromaticity-type (aromaticity-info atom arom rule aromatic-atoms)
  (unless (is-aromatic-type atom arom)
    (pushnew (make-aromatic-info :aromaticity-type arom :rule-name rule)
             (gethash atom aromaticity-info)
             :test #'equalp))
  (pushnew atom aromatic-atoms)
  aromatic-atoms)

(defun aromaticity-information-as-string (atom &optional (aromaticity-info chem:*current-aromaticity-information*))
  (format nil "~a" (gethash atom aromaticity-info)))

(defun remove-hits-from-list (atoms-in-rings hits)
  (loop for atm in atoms-in-rings
        unless (member atm hits)
          collect atm))

(progn
  (defparameter *artests* (list (cons :ar6 (lambda (atom)
                                             (let ((is-aromatic-type (is-aromatic-type atom :ar6)))
                                               is-aromatic-type)))
                                (cons :abx #'chem:am1bcc-x)
                                (cons :aby #'chem:am1bcc-y)))

;; Apply aromaticity rule1 of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
  (defparameter *rule1* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<abx>:5]=[<abx>:6]1" :tests *artests*))))


(defun exhaustively-apply-aromatic-rule (aromaticity-info molecule-graph rule aromaticity-type rule-name)
  (when (chem:verbose 1)
    (format t "Applying aromatic rule ~a~%" rule-name))
  (let ((aromatic-atoms nil))
    (loop
      ;; Keep applying the rule until we don't find any new atoms
      (let ((hits (chem:boost-graph-vf2 rule molecule-graph))
            (new-aromatic-atoms aromatic-atoms))
        (loop for hit-vector in hits
              do (loop for index from 0 below (length hit-vector)
                       for atom = (elt hit-vector index)
                       when atom
                         do (setf new-aromatic-atoms (set-aromaticity-type aromaticity-info atom aromaticity-type rule-name new-aromatic-atoms))))
        (when (eq new-aromatic-atoms aromatic-atoms)
          (return-from exhaustively-apply-aromatic-rule nil))
        (setf aromatic-atoms new-aromatic-atoms)))))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defparameter *rule2* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<ar6>]:[<ar6>]1" :tests *artests*)))

(defparameter *rule3* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<ar6>]:[<ar6>]:[<ar6>]:[<ar6>]1" :tests *artests*)))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defparameter *rule4* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[C&+:1]1-[<abx>:2]=[<abx>:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*)))


(defparameter *ar67test*
  (list* (cons :ar67
               (lambda (atom)
                 (or (is-aromatic-type atom :ar6)
                     (is-aromatic-type atom :ar7))))
         *artests*))

(defparameter *rule5* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[<aby>:1]1-[<abx>&!<ar67>:2]=[<abx>&!<ar67>:3]-[<abx>&!<ar67>:4]=[<abx>&!<ar67>:5]1" :tests *ar67test*)))

(defparameter *rule6* #+nosmarts nil #-nosmarts(chem:make-chem-info-graph (chem:compile-smarts "[<aby>:1]1-[<abx>&!<ar67>:2]=[<abx>&!<ar67>:3]-[<abx>&!<ar67>:4]=[<abx>&!<ar67>:5]1" :tests *ar67test*)))

(defun ring-bonds-aromatic-p (ring)
  (let ((ring-set (make-hash-table :test 'eq))
        (ok t))
    (dolist (a ring) (setf (gethash a ring-set) t))
    (dolist (a ring)
      (let* ((ring-neighbors
               (remove-if-not (lambda (b) (gethash b ring-set))
                              (chem:atom/bonded-atoms-as-list a))))
        (when (/= (length ring-neighbors) 2)
          (setf ok nil))
        (dolist (b ring-neighbors)
          (unless (eq :aromatic-bond (chem:bond-order-to a b))
            (setf ok nil)))))
    ok))

(defun ring-aromaticity-type (ring)
  (case (length ring)
    (5 :ar5)
    (6 :ar6)
    (7 :ar7)
    (t :ar6)))

(defun mark-aromatic-rings-from-bonds (aromaticity-info rings)
  (let ((aromatic-atoms nil))
    (dolist (ring rings)
      (when (ring-bonds-aromatic-p ring)
        (let ((arom (ring-aromaticity-type ring)))
          (dolist (atom ring)
            (setf aromatic-atoms
                  (set-aromaticity-type aromaticity-info atom arom :rosetta-bonds aromatic-atoms))))))
    aromatic-atoms))


(defun mark-aromatic-atoms-from-properties (aromaticity-info matter)
  (let ((aromatic-atoms nil))
    (chem:do-atoms (atom matter)
      (let ((ar-type (getf (chem:properties atom) :aromatic)))
        (when ar-type
          (setf aromatic-atoms
                (set-aromaticity-type aromaticity-info atom ar-type :property-list aromatic-atoms))))
      aromatic-atoms)))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defun is-ring-aromatic (aromaticity-info ring)
  (dolist (atom ring)
    (when (not (is-aromatic atom aromaticity-info)) (return-from is-ring-aromatic nil))
    )
  t)

(defun select-aromatic-rings (aromaticity-info all-rings)
  (let (aromatic-rings)
    (dolist (ring all-rings)
      (when (is-ring-aromatic aromaticity-info ring)
        (setq aromatic-rings (cons ring aromatic-rings))))
    aromatic-rings))


(defun all-ring-atoms (matter)
  (let* ((all-rings (chem:identify-rings matter))
         (atoms-in-rings (loop for ring in all-rings
                               append ring)))
    (values atoms-in-rings all-rings)))


(defparameter *save-aromatic-info* nil)

(defgeneric chem:identify-aromatic-rings (matter force-field-name))

(defmethod chem:identify-aromatic-rings (matter (force-field-name t))
  (make-hash-table))

(defmethod chem:identify-aromatic-rings (matter (force-field-name (eql :am1bcc)))
  "Use the atoms-in-rings - found using RingFinder_O::identifyRings(matter) and run all of the
aromaticity tests on the atoms in the rings.   Assign aromaticity flags of aromatic atoms by
associating the atom with its aromaticity info in a hash-table and return the hash-table."
  (unless (boundp 'chem:*current-rings*)
    (error "chem:*current-rings* must be bound to a list of rings in the molecule identified using chem:identify-rings"))
  (let* ((aromaticity-info (make-hash-table))
         (chem:*current-aromaticity-information* aromaticity-info))
    (mark-aromatic-rings-from-bonds aromaticity-info chem:*current-rings*)
    (mark-aromatic-atoms-from-properties aromaticity-info matter)
    (chem:do-molecules (molecule matter)
      (let ((molecule-graph (chem:make-molecule-graph-from-molecule molecule)))
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule1* :ar6 :rule1)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule2* :ar6 :rule2)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule3* :ar6 :rule3)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule4* :ar7 :rule4)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule5* :ar5 :rule5)))
    (setf *save-aromatic-info* aromaticity-info)
    aromaticity-info))


(defmethod chem:identify-aromatic-rings (matter (force-field-name (eql :mdl)))
  "Use the atoms-in-rings - found using RingFinder_O::identifyRings(matter) and run all of the
aromaticity tests on the atoms in the rings.   Assign aromaticity flags of aromatic atoms by
associating the atom with its aromaticity info in a hash-table and return the hash-table."
  (unless (boundp 'chem:*current-rings*)
    (error "chem:*current-rings* must be bound to a list of rings in the molecule identified using chem:identify-rings"))
  (let* ((aromaticity-info (make-hash-table))
         (chem:*current-aromaticity-information* aromaticity-info))
    (mark-aromatic-rings-from-bonds aromaticity-info chem:*current-rings*)
    (mark-aromatic-atoms-from-properties aromaticity-info matter)
    (chem:do-molecules (molecule matter)
      (let ((molecule-graph (chem:make-molecule-graph-from-molecule molecule)))
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule1* :ar6 :rule1)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule2* :ar6 :rule2)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule3* :ar6 :rule3)))
    (setf *save-aromatic-info* aromaticity-info)
  aromaticity-info))

(defmethod chem:identify-aromatic-rings (matter (force-field-name (eql :default)))
  (chem:identify-aromatic-rings matter :mdl))

(defmethod chem:identify-aromatic-rings (matter (force-field-name (eql :smirnoff)))
  (chem:identify-aromatic-rings matter :mdl))




