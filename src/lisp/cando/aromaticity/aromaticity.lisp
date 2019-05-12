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


(defstruct aromatic-info aromaticity-type step-assigned)

;; Set the aromatic type of an atom
(defun set-aromaticity-type (aromaticity-info atom arom rule)
  ;;  (set-is-aromatic a t)
  (pushnew (make-aromatic-info :aromaticity-type arom :step-assigned rule)
           (gethash atom aromaticity-info)
           :test #'equalp))

;; Return the aromaticity type of the atom followed by the step-assigned or (values nil nil)
(defun is-aromatic (atom &optional (aromaticity-info chem:*current-aromaticity-information*))
  (gethash atom aromaticity-info))

(defun aromaticity-information-as-string (atom &optional (aromaticity-info chem:*current-aromaticity-information*))
  (format nil "~a" (gethash atom aromaticity-info)))

(defun remove-hits-from-list (atoms-in-rings hits)
  (loop for atm in atoms-in-rings
        unless (member atm hits)
          collect atm))

(progn
  (defparameter *artests* (list (cons :ar6 (lambda (atom)
                                             (loop named test
                                                   for info in (gethash atom chem:*current-aromaticity-information*)
                                                   when (eq (aromatic-info-aromaticity-type info) :ar6)
                                                     do (return-from test t))
                                             nil))
                                (cons :abx #'chem:am1bcc-x)
                                (cons :aby #'chem:am1bcc-y)))

;; Apply aromaticity rule1 of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
  (defparameter *rule1* (chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<abx>:5]=[<abx>:6]1" :tests *artests*))))


(defun exhaustively-apply-aromatic-rule (aromaticity-info molecule-graph rule aromaticity-type rule-name)
  (let ((hits (chem:boost-graph-vf2 rule molecule-graph)))
    (loop for hit-vector in hits
          do (loop for index from 0 below (length hit-vector)
                   for atom = (elt hit-vector index)
                   when atom
                     do (set-aromaticity-type aromaticity-info atom aromaticity-type rule-name)))))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defparameter *rule2* (chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<ar6>]-[<ar6>]1" :tests *artests*)))

(defparameter *rule3* (chem:make-chem-info-graph (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<ar6>]-[<ar6>]-[<ar6>]-[<ar6>]1" :tests *artests*)))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defparameter *rule4* (chem:make-chem-info-graph (chem:compile-smarts "[C&+:1]1-[<abx>:2]=[<abx>:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*)))


(defparameter *ar67test*
  (list* (cons :ar67
               (lambda (atom)
                 (loop named test
                       for info in (gethash atom chem:*current-aromaticity-information*)
                       when (or (eq (aromatic-info-aromaticity-type info) :ar6)
                                (eq (aromatic-info-aromaticity-type info) :ar7))
                         do (return-from test t))
                 nil))
         *artests*))

(defparameter *rule5* (chem:make-chem-info-graph (chem:compile-smarts "[<aby>:1]1-[<abx>&!<ar67>:2]=[<abx>&!<ar67>:3]-[<abx>&!<ar67>:4]=[<abx>&!<ar67>:5]1" :tests *ar67test*)))

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

(defun chem:identify-aromatic-rings (matter)
  "Use the atoms-in-rings - found using RingFinder_O::identifyRings(matter) and run all of the
aromaticity tests on the atoms in the rings.   Assign aromaticity flags of aromatic atoms by
associating the atom with its aromaticity info in a hash-table and return the hash-table."
  (let* ((aromaticity-info (make-hash-table))
         (chem:*current-aromaticity-information* aromaticity-info))
    (cando:do-molecules (molecule matter)
      (let ((molecule-graph (chem:make-molecule-graph molecule)))
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule1* :ar6 'rule1)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule2* :ar6 'rule2)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule3* :ar6 'rule3)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule4* :ar7 'rule4)
        (exhaustively-apply-aromatic-rule aromaticity-info molecule-graph *rule5* :ar5 'rule5)))
    (setf *save-aromatic-info* aromaticity-info)
    aromaticity-info))

;;
;; Identify all rings, isolate atoms in rings and apply all of the aromaticity rules
;; from Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
;; Return the rings for bond type assignment
;; Make this work using the hash-table returned by chem:identify-aromatic-rings
#+(or)
(defun identify-aromatic-rings (matter)
  (unless (or (typep matter 'chem:molecule) (typep matter 'chem:aggregate))
    (error "identify-aromatic-rings only accepts aggregate or molecule - you passed ~s" matter))
  (cando:do-molecules (molecule matter)
    (multiple-value-bind (atoms-in-rings all-rings)
        (all-ring-atoms molecule)
      (chem:identify-aromatic-rings atoms-in-rings)
      (select-aromatic-rings all-rings))))


(defmacro with-aromaticity-information ((matter) &body body)
  "Provide a dynamic environment where the aromaticity:is-aromatic function works.
This macro first calculates the aromaticity information of the matter and then
evaluates the body in that dynamic environment."
  `(let ((chem:*current-aromaticity-information* (chem:identify-aromatic-rings ,matter)))
     ,@body))


