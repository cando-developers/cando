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
  (push (make-aromatic-info :aromaticity-type arom :step-assigned rule)
        (gethash atom aromaticity-info)))


;; Return the aromaticity type of the atom followed by the step-assigned or (values nil nil)
(defun is-aromatic (aromaticity-info atom)
  (gethash atom aromaticity-info))

(defun one-atom-aromatic-rule (aromaticity-info atm rule aromaticity-type rule-name)
  ;; (print (% "test1 on %s" (description a)))
  (let ((match (chem:matches rule atm))
        hits)
    (when match
      (loop for match-atom across (chem:tags-as-vector match)
            when match-atom
              do (progn
                   (set-aromaticity-type aromaticity-info match-atom aromaticity-type rule-name)
                   (push match-atom hits))))
    ;; At this point if we look at bonds between the atoms in hits we could
    ;; set the bond types to 'aromatic' if that was the right thing to do.
    hits))

(defun remove-hits-from-list (atoms-in-rings hits)
  (loop for atm in atoms-in-rings
        unless (member atm hits)
          collect atm))

(defparameter *artests* (list (cons :ar6 (lambda (atom)
                                           (loop named test
                                                 for info in (gethash atom chem:*current-aromaticity-information*)
                                                 when (eq (aromatic-info-aromaticity-type info) :ar6)
                                                   do (return-from test t))
                                           nil))
                              (cons :abx #'chem:am1bcc-x)
                              (cons :aby #'chem:am1bcc-y)))

(defparameter *rule1* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<abx>:5]=[<abx>:6]1" :tests *artests*))


;; Apply aromaticity rule1 of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule (aromaticity-info atoms-in-rings rule aromaticity-type rule-name)
  (let ((chem:*current-aromaticity-information* aromaticity-info))
    (loop
      (let ((aromatic-rule-hits (loop for atm in atoms-in-rings
                                      append (one-atom-aromatic-rule aromaticity-info atm rule aromaticity-type rule-name))))
        (unless aromatic-rule-hits
          (return-from aromatic-rule atoms-in-rings))
        (setf atoms-in-rings (remove-hits-from-list atoms-in-rings aromatic-rule-hits))))))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defparameter *rule2* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<ar6>]-[<ar6>]1" :tests *artests*))

(defparameter *rule3* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<ar6>]-[<ar6>]-[<ar6>]-[<ar6>]1" :tests *artests*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defparameter *rule4a* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[C&+:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*))
(defparameter *rule4b* (chem:compile-smarts "[C&+:1]1-[<abx>:2]=[<abx>:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*))


(defparameter *ar67test* (list* (cons :ar67 (lambda (atom)
                                              (loop named test
                                                    for info in (gethash atom chem:*current-aromaticity-information*)
                                                    when (or (eq (aromatic-info-aromaticity-type info) :ar6)
                                                             (eq (aromatic-info-aromaticity-type info) :ar7))
                                                      do (return-from test t))
                                             nil))
                                *artests*))

(defparameter *rule5* (chem:compile-smarts "[<aby>:1]1-[<abx>&!<ar67>:2]=[<abx>&!<ar67>:3]-[<abx>&!<ar67>:4]=[<abx>&!<ar67>:5]1" :tests *ar67test*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defun is-ring-aromatic (aromaticity-info ring)
  (dolist (atom ring)
    (when (not (is-aromatic aromaticity-info atom)) (return-from is-ring-aromatic nil))
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



(defun chem:identify-aromatic-rings (atoms-in-rings)
  "Use the atoms-in-rings - found using RingFinder_O::identifyRings(matter) and run all of the
aromaticity tests on the atoms in the rings.   Assign aromaticity flags of aromatic atoms by
associating the atom with its aromaticity info in a hash-table and return the hash-table."
  (format t "atoms in rings ~a~%" atoms-in-rings)
  (let ((aromaticity-info (make-hash-table)))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule1* :ar6 'rule1))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule2* :ar6 'rule2))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule3* :ar6 'rule3))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule4a* :ar7 'rule4a))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule4b* :ar7 'rule4b))
    (setf atoms-in-rings (aromatic-rule aromaticity-info atoms-in-rings *rule5* :ar5 'rule5))
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


