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

(defparameter *artests* (list (cons :ar6 (lambda (a) (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar6)))
                              (cons :abx #'chem:am1bcc-x)
                              (cons :aby #'chem:am1bcc-y)))

(defparameter *rule1a* (chem:compile-smarts "[<abx>]1-[<abx>]=[<abx>]-[<abx>]=[<abx>]-[<abx>]=[<abx>]1" :tests *artests*))

(defun one-atom-aromatic-rule1 (a)
;        (print (% "test1 on %s" (description a)))
  (when (chem:matches *rule1a* a)
    (set-aromaticity-type a 'ar6 'rule1)))

;; Apply aromaticity rule1 of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule1 (atoms-in-rings)
  (loop for a in atoms-in-rings
;        (print (% "test1 on %s" (description a)))
        do (one-atom-aromatic-rule1 a)))
  
(defparameter *rule2a* (chem:compile-smarts "[<abx>]1=[<abx>]-[<abx>]=[<abx>]-[<ar6>]-[<ar6>]-[<abx>]1" :tests *artests*))
(defparameter *rule2b* (chem:compile-smarts "[<abx>]1-[<abx>]=[<abx>]-[<ar6>]-[<ar6>]-[<abx>]=[<abx>]1" :tests *artests*))


;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule2 (atoms-in-rings)
  (let ((atoms-left atoms-in-rings)
        (failed-atoms ())
        (rule-pass 0)
        (found-some t))
    (loop while found-some
       do (setq found-some nil)
       do (loop for a in atoms-left
             do (if (or (chem:matches *rule2a* a)
                        (chem:matches *rule2b* a))
                    (progn
                      (setq found-some t)
                      (set-aromaticity-type a 'ar6 'rule2))
                    (setq failed-atoms (cons a failed-atoms))))
       do (setq atoms-left failed-atoms)
       do (setq rule-pass (+ rule-pass 1)))))


(defparameter *rule3a* (chem:compile-smarts "[<abx>]1=[<abx>]-[<ar6>]-[<ar6>]-[<ar6>]-[<ar6>]-[<abx>]1" :tests *artests*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule3 (atoms-in-rings)
  (let ((atoms-left atoms-in-rings)
	(failed-atoms ())
	(rule-pass 0)
	(found-some t))
    (loop while found-some
       do (setq found-some nil)
       do (loop for a in atoms-left
	     do (if (chem:matches *rule3a* a)
		    (progn
		      (setq found-some t)
		      (set-aromaticity-type a 'ar6 'rule3)
		      )
		    (setq failed-atoms (cons a failed-atoms))))
       do (setq atoms-left failed-atoms)
       do (setq rule-pass (+ rule-pass 1)))))

(defparameter *rule4a* (chem:compile-smarts "[<abx>]1-[<abx>]=[<abx>]-[C&+]-[<abx>]=[<abx>]-[<abx>]=[<abx>]1" :tests *artests*))
(defparameter *rule4b* (chem:compile-smarts "[<abx>]1=[<abx>]-[C&+]-[<abx>]=[<abx>]-[<abx>]=[<abx>]-[<abx>]1" :tests *artests*))
(defparameter *rule4c* (chem:compile-smarts "[<abx>]1-[C&+]-[<abx>]=[<abx>]-[<abx>]=[<abx>]-[<abx>]=[<abx>]1" :tests *artests*))
(defparameter *rule4d* (chem:compile-smarts "[C&+]1-[<abx>]=[<abx>]-[<abx>]=[<abx>]-[<abx>]=[<abx>]-[C]1" :tests *artests*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule4 (atoms-in-rings)
    (loop for a in atoms-in-rings
         do (when (or (chem:matches *rule4a* a)
                   (chem:matches *rule4b* a)
                   (chem:matches *rule4c* a)
                   (chem:matches *rule4d* a)
                   )
           (set-aromaticity-type a 'ar7 'rule4)
           )))

(defparameter *ar67test* (list* (cons :ar67 (lambda (a) (or (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar6)
                                                              (eq (chem:matter-get-property-or-default a 'aromaticity nil) 'ar7))))
                                *artests*))

(defparameter *rule5a* (chem:compile-smarts "[<aby>]1-[<abx>&!<ar67>]=[<abx>&!<ar67>]-[<abx>&!<ar67>]=[<abx>&!<ar67>]-[<aby>]1" :tests *ar67test*))
(defparameter *rule5b* (chem:compile-smarts "[<abx>&!<ar67>]1-[<aby>]-[<abx>&!<ar67>]=[<abx>&!<ar67>]-[<abx>&!<ar67>]=[<abx>]1" :tests *ar67test*))
(defparameter *rule5c* (chem:compile-smarts "[<abx>&!<ar67>]1=[<abx>&!<ar67>]-[<aby>]-[<abx>&!<ar67>]=[<abx>&!<ar67>]-[<abx>]1" :tests *ar67test*))

;; Apply aromaticity rule of Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun aromatic-rule5 (atoms-in-rings)
    (loop for a in atoms-in-rings
        do (when (or (chem:matches *rule5a* a)
                   (chem:matches *rule5b* a)
                   (chem:matches *rule5c* a))
           (set-aromaticity-type a 'ar5 'rule5)
           )))


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

(defun generate-atoms-in-rings-list (mol)
  (let ((atoms-in-rings nil))
    (chem:map-atoms
     nil
     (lambda (a)
       (when (chem:is-in-ring a)
         (push a atoms-in-rings)))
     mol)
    atoms-in-rings))



;;
;; Identify all rings, isolate atoms in rings and apply all of the aromaticity rules
;; from Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
;; Return the rings for bond type assignment
(defun identify-aromatic-rings (mol)
  (let ((all-rings (chem:identify-rings mol))
        #+(or)(atoms-in-rings (mapcar (lambda (x) (isInRing x)) (all-atoms-as-cons mol nil)))
        (atoms-in-rings (generate-atoms-in-rings-list mol)))
    (format t "rings ~a~%" atoms-in-rings)
    (aromatic-rule1 atoms-in-rings)
    (format t "rings ~a~%" atoms-in-rings)
    (aromatic-rule2 atoms-in-rings)
    (format t "rings ~a~%" atoms-in-rings)
    (aromatic-rule3 atoms-in-rings)
    (format t "rings ~a~%" atoms-in-rings)
    (aromatic-rule4 atoms-in-rings)
    (format t "rings ~a~%" atoms-in-rings)
    (aromatic-rule5 atoms-in-rings)
    (format t "rings ~a~%" atoms-in-rings)
    (select-aromatic-rings all-rings)))


