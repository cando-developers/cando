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


;; Set the aromatic type of an atom
(defun set-aromaticity-type (a arom rule)
;;  (set-is-aromatic a t)
  (chem:set-property a 'aromaticity arom)
  (chem:set-property a 'step-assigned rule))


;; Return the aromaticity type of the atom followed by the step-assigned or (values nil nil)
(defun is-aromatic (a)
  (values (chem:matter-get-property-or-default a 'aromaticity ())
          (chem:matter-get-property-or-default a 'step-assigned ())))

(defun one-atom-aromatic-rule (atm rule aromaticity-type rule-name)
  ;; (print (% "test1 on %s" (description a)))
  (let ((match (chem:matches rule atm))
        hits)
    (when match
      (loop for match-atom across (chem:tags-as-vector match)
            when match-atom
              do (progn
                   (set-aromaticity-type match-atom aromaticity-type rule-name)
                   (push match-atom hits))))
    hits))

(defun remove-hits-from-list (atoms-in-rings hits)
  (loop for atm in atoms-in-rings
        unless (member atm hits)
          collect atm))

(defparameter *artests* (list (cons :ar6 (lambda (a) (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar6)))
                              (cons :abx #'chem:am1bcc-x)
                              (cons :aby #'chem:am1bcc-y)))

(defparameter *rule1* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<abx>:5]=[<abx>:6]1" :tests *artests*))


;; Apply aromaticity rule1 of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule (atoms-in-rings rule aromaticity-type rule-name)
  (loop
    (let ((aromatic-rule-hits (loop for atm in atoms-in-rings
                                     append (one-atom-aromatic-rule atm rule aromaticity-type rule-name))))
      (unless aromatic-rule-hits
        (return-from aromatic-rule atoms-in-rings))
      (setf atoms-in-rings (remove-hits-from-list atoms-in-rings aromatic-rule-hits)))))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defparameter *rule2* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<abx>:3]=[<abx>:4]-[<ar6>]-[<ar6>]1" :tests *artests*))

(defparameter *rule3* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[<ar6>]-[<ar6>]-[<ar6>]-[<ar6>]1" :tests *artests*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defparameter *rule4a* (chem:compile-smarts "[<abx>:1]1=[<abx>:2]-[C&+:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*))
(defparameter *rule4b* (chem:compile-smarts "[C&+:1]1-[<abx>:2]=[<abx>:3]-[<abx>:4]=[<abx>:5]-[<abx>:6]=[<abx>:7]1" :tests *artests*))


(defparameter *ar67test* (list* (cons :ar67 (lambda (a) (or (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar6)
                                                            (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar7))))
                                *artests*))

(defparameter *rule5* (chem:compile-smarts "[<aby>:1]1-[<abx>:2&!<ar67>]=[<abx>:3&!<ar67>]-[<abx>:4&!<ar67>]=[<abx>:5&!<ar67>]1" :tests *ar67test*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry

(defun is-ring-aromatic (ring)
  (dolist (atom ring)
    (when (not (is-aromatic atom)) (return-from is-ring-aromatic nil))
    )
  t)
             
(defun select-aromatic-rings (all-rings)
  (let (aromatic-rings)
    (dolist (ring all-rings)
      (when (is-ring-aromatic ring)
        (setq aromatic-rings (cons ring aromatic-rings))))
    aromatic-rings))

(defun generate-atoms-in-rings-list (mo)
  (let ((atoms-in-rings nil))
    (chem:map-atoms
     nil
     (lambda (a)
       (when (chem:is-in-ring a)
         (push a atoms-in-rings)))
     mol)
    atoms-in-rings))


(defun all-ring-atoms (matter)
  (let* ((all-rings (chem:identify-rings matter))
         (atoms-in-rings (loop for ring in all-rings
                               append ring)))
    (values atoms-in-rings all-rings)))

;;
;; Identify all rings, isolate atoms in rings and apply all of the aromaticity rules
;; from Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
;; Return the rings for bond type assignment
(defun identify-aromatic-rings (mol)
  (multiple-value-bind (atoms-in-rings all-rings)
      (all-ring-atoms mol)
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule1* 'ar6 'rule1))
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule2* 'ar6 'rule2))
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule3* 'ar6 'rule3))
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule4a* 'ar7 'rule4a))
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule4b* 'ar7 'rule4b))
    (setf atoms-in-rings (aromatic-rule atoms-in-rings *rule5* 'ar5 'rule5))
    (select-aromatic-rings all-rings)))


