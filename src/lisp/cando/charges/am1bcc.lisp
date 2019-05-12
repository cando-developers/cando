;(include "CANDO/pattern/aromaticity.l")
;(include "CANDO/db/build-molecule.l")

(in-package :charges)

(defvar *rr* 1)

(defvar *am1-bcc-raw-data*
  '((110111 0.0000)
    (110112 0.0042) (110113 -0.0753) (110114 -0.0500) (110115 -0.0269)
    (110116 0.0073) (110117 -0.0943) (110121 0.1582) (110122 0.0374) (110123 -0.0193)
    (110124 0.1252) (110125 0.1825) (110131 0.0718) (110141 0.2579)
    (110142 0.4078) (110151 0.1821) (110152 0.4263) (110153 0.4608) (110161 0.1540)
    (110171 0.0713) (110172 0.0734) (110173 0.1274) (110174 0.3009)
    (110191 0.0393) (120112 0.0000) (120113 -0.0866) (120114 -0.0192) (120115 0.0333)
    (120116 0.0118) (120117 -0.1098) (120121 0.0820) (120122 0.0026)
    (120123 -0.0271) (120124 0.1406) (120125 0.4897) (120131 0.0441) (120141 0.2986)
    (120142 0.4154) (120151 0.2406) (120152 0.4307) (120153 0.4509)
    (120161 0.2128) (120171 0.0625) (120172 0.0844) (120173 0.1257) (170151 0.3273)
    (170152 0.5416) (170153 0.5169) (170171 0.1402) (170172 0.1660)
    (170173 0.2312) (170174 0.3528) (170191 0.1369) (210121 0.0000) (210122 -0.0302)
    (210123 -0.1185) (210124 0.0867) (210131 -0.0517) (120174 0.2728)
    (120191 0.0000 ) (130113 0.0000) (130114 0.0693) (130115 0.0302) (130116 0.0839)
    (130117 -0.0297) (130121 0.1522) (130122 0.0069) (130123 0.0388)
    (130124 0.2193) (130125 0.2406) (130131 0.1336) (130141 0.3529) (130142 0.5230)
    (130151 0.2242) (130152 0.4296) (130153 0.5077) (130171 0.1055)
    (130172 0.1123) (130173 0.1253) (130174 0.2752) (130191 0.1292) (140114 0.0000)
    (140115 -0.0308) (140116 0.0206) (140117 -0.0895) (140122 0.0670)
    (140123 -0.0432) (140124 0.1394) (140125 0.0070) (140131 0.0901) (140141 0.2828)
    (140142 0.4820) (140151 0.1826) (140152 0.4319) (140153 0.4840)
    (140171 0.1115) (140172 0.0486) (140173 0.1106) (140174 0.2601) (140191 0.0928)
    (150115 0.0000) (150116 0.0040) (150117 -0.1070) (230153 0.3418)
    (230161 0.1766) (230171 0.0561) (230172 0.0642) (230173 0.1802) (230191 -0.0497)
    (240124 0.0000) (240125 0.1043) (240131 -0.0985) (240141 0.1287)
    (240142 0.1498) (240151 0.0315) (240152 0.1731) (150121 0.0558) (150122 0.0062)
    (150123 -0.0543) (150124 0.0520) (150125 0.0798) (150131 0.0286)
    (150141 0.3987) (150142 0.3545) (150151 0.2251) (150152 0.4586) (150153 0.4018)
    (150161 0.2343) (150171 0.0021) (150172 0.0899) (150173 0.2338)
    (150174 0.2886) (150191 0.0575) (160116 0.0000) (160117 -0.1075) (160121 0.0816)
    (160122 -0.0074) (160123 -0.0452) (160124 0.1394) (160125 0.1877)
    (160131 0.0451) (160141 0.3270) (160142 0.3895) (160151 0.2269) (160152 0.4285)
    (160153 0.4445) (160161 0.1938) (160171 0.0349) (160172 0.0784)
    (160173 0.1401) (160174 0.2859) (160191 0.0000) (170117 0.0000) (170121 0.1929)
    (170122 0.0799) (170123 0.0635) (170124 0.2190) (170125 0.2100)
    (170131 0.1703) (170141 0.4597) (170142 0.4543) (420171 -0.1976) (420172 -0.3077)
    (420173 -0.2425) (420191 -0.3756) (510151 0.0000) (510152 0.2576)
    (510153 0.2971) (510161 -0.0936) (510171 -0.0547) (510172 -0.0260) (510173 0.1105)
    (510191 -0.1718) (520152 0.0000) (210141 0.2527) (210142 0.2947)
    (210151 0.0757) (210152 0.3531) (210153 0.3359) (210161 0.1037) (210171 -0.0204)
    (210172 -0.0369) (210173 0.0876) (210174 0.1898) (210191 -0.2048)
    (220122 0.0000) (220123 -0.0327) (220124 0.1326) (220125 0.0236) (220131 -0.0059)
    (220141 0.1682) (220142 0.3046) (220151 0.1593) (220152 0.3657)
    (220153 0.3322) (220171 -0.0155) (220172 -0.0015) (220173 0.0838) (220174 0.2067)
    (220191 -0.0865) (230123 0.0000) (230124 0.1794 ) (230131 0.0202)
    (230141 0.2538) (230142 0.4175) (230151 0.2024) (240141 0.0535) (240142 0.4287)
    (240151 0.1684) (240152 0.3764) (240153 0.4265) (250225 0.0000)
    (310241 0.1315) (310242 0.2707) (310251 0.1583) (310252 0.2718) (310253 0.2792)
    (410241 0.0000) (410251 -0.5630) (510252 0.3295) (150215 0.0000)
    (150325 -0.1168) (210631 0.0465) (330631 0.1317) (240631 -0.0817) (250631 -0.0543)
    (250651 0.1473) (160716 0.0000) (160717 -0.0653) (160723 -0.0243)
    (160731 -0.1518) (160741 0.2240) (240153 0.2714) (240161 0.1150) (240171 -0.0554)
    (240172 -0.0603) (240173 0.0192) (240174 0.1171) (240191 -0.2444)
    (250191 -0.1634) (310131 0.0000) (310141 0.2229) (310142 0.2336) (310151 0.0442)
    (310152 0.2152) (310153 0.2093) (310161 0.0083) (310171 0.0375)
    (310172 0.0014) (310173 0.0699) (310174 0.1712) (310191 -0.2010) (410141 0.0000)
    (410142 0.2762) (410151 -0.1072) (410161 0.0592) (410171 -0.2769)
    (410172 -0.2065) (410173 -0.1309) (410174 0.0181) (410191 -0.1593) (420142 0.0000)
    (420151 -0.3626) (420152 0.0675) (160751 0.2699) (170717 0.0000)
    (170723 0.0885) (170724 0.2630) (170731 0.1294) (170741 0.2464) (170751 0.3144)
    (230723 0.0000) (230724 0.1698) (230731 0.0677) (230741 0.1806)
    (230751 0.1565) (240724 0.0000) (240731 -0.0723) (240741 0.0456) (240751 0.0419)
    (310741 0.0710) (310751 0.0983) (410751 0.0587) (510751 0.0000)
    (160816 0.0000) (160817 -0.0653) (160823 -0.0243) (160831 -0.1518) (160841 0.2240)
    (160851 0.2699) (170817 0.0000) (170823 0.0885) (520172 -0.2350)
    (520191 -0.4040) (530171 -0.2041) (530172 -0.2761) (530173 -0.1970) (530174 -0.1490)
    (530191 -0.4588) (610172 -0.0555) (610191 -0.0258) (120212 0.0000)
    (120215 -0.0631) (130223 0.0547) (130224 0.2877) (130225 0.1624) (130241 0.2496)
    (130242 0.6836) (140231 0.2391) (140232 0.1890) (140233 0.2755)
    (140251 0.2834) (140252 0.4447) (140253 0.5617) (150215 0.0000) (150224 0.2204)
    (150231 0.2335) (150242 0.7150) (150251 0.3429) (150252 0.3052)
    (230224 0.1731) (240224 0.0000) (240225 0.0040) (240231 -0.0467) (170824 0.2630)
    (170831 0.1294) (170841 0.2464) (170851 0.3144) (230824 0.1698)
    (240824 0.0000) (240841 0.0456) (110931 0.1615) (110951 0.0833) (120931 0.1844)
    (120951 0.1879) (130931 0.3046) (130951 0.1969) (140931 0.2653)
    (140951 0.1713) (150931 0.2549) (150951 0.3281) (160931 0.2797) (160951 0.2541)
    (170931 0.3370) (170951 0.3476) (230931 -0.1500) (310941 0.2354)
    (310942 0.3195) (310951 0.1530) (310952 0.3916) (310953 0.3228) (510953 0.5218)
    ))

(defvar *am1-bcc-parameters* (make-hash-table :test #'eql))

;(loop for pair in *am1-bcc-table*
;                  (setf (gethash (car pair) *am1-bcc-parameters*) (second pair)))

(defun set-am1-bcc-type (a id)
  (chem:set-property a 'am1-bcc-type id))

(defun get-am1-bcc-type (a)
  (chem:matter-get-property-or-default a 'am1-bcc-type ()))



;;
;; Return a text description of an atom
(defun describe-atom (a)
  (if a
      (multiple-value-bind (aromatic step-assigned) (aromaticity:is-aromatic a)
        (format t "~a type[~a] ~a ~a ~a~%" (chem:get-name a)
           (get-am1-bcc-type a)
           (aromaticity:aromaticity-information-as-string a)
           aromatic
           step-assigned)
        )
    (debug))
  )




;;
;; Set the a property in the bond that registers the Am1Bcc bond type
(defun set-bond-type (b id)
;  (print (% "Setting bond type[%d] for bond: %s" id (repr b)))
  (chem:bond-set-property b 'am1-bcc-type id))

(defun has-bond-type (b)
  (chem:bond-has-property b 'am1-bcc-type))

(defun get-bond-type (b)
  (chem:bond-get-property b 'am1-bcc-type ()))




(defparameter *carboxyl-nitro-pattern* (chem:compile-smarts "[$([C:1](=[O:2])-[O&-;X1:3]),$([N;+:1](=[O:2])-[O&-;X1:3])]"))

;;
;; Assign Am1Bcc bond types for all bonds
;; For the rules see Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun assign-am1-bcc-bond-types (mol)
                                        ; Identify the aromatic bonds
                                        ; any pair of atoms in the same ring that has a single or
                                        ; double bond is an aromatic bond
  (let ((all-rings (chem:identify-rings mol)))
    (loop for ring in all-rings
          for ring-bonds = (chem:ring-bonds ring)
          do (loop for bond in ring-bonds
                   for a1 = (chem:get-atom1 bond)
                   for a2 = (chem:get-atom2 bond)
                   when (and (aromaticity:is-aromatic a1)
                             (aromaticity:is-aromatic a2))
                     do (cond
                          ((eq (chem:get-order bond) :single-bond)
                           (set-bond-type bond 7))
                          ((eq (chem:get-order bond) :double-bond)
                           (set-bond-type bond 8))
                          (t nil)))))
                                        ; First look for delocalized bonds - carboxyl groups, nitro groups
  ;;    (loop for a in (all-atoms-as-cons mol nil)
  (chem:map-atoms
   nil
   (lambda (a)
     (let (match)
       (when (setf match (chem:matches *carboxyl-nitro-pattern* a))
         (let* ((a1 (chem:tag match "1"))
                (a2 (chem:tag match "2"))
                (a3 (chem:tag match "3"))
                (double-bond (chem:bond-to a1 a2))
                (single-bond (chem:bond-to a1 a3)))
           (set-bond-type double-bond 9)
           (set-bond-type single-bond 9)))))
   mol)
;   (loop for bond in (all-bonds-as-cons mol nil)
  (chem:map-bonds
   nil
   (lambda (a1 a2 order)
     (let ((bond (chem:get-bond-to a1 a2)))
       (when (not (has-bond-type bond))
         (cond
           ((eq (chem:get-order bond) :single-bond)
            (cond
              ((or (and (= (chem:get-ionization a1) -1) (= (chem:get-ionization a2) +1))
                   (and (= (chem:get-ionization a1) 1) (= (chem:get-ionization a2) -1)))
               (set-bond-type bond 06)) ; +- polarized bond
              ((or (and (= (chem:get-ionization a1) 0) (= (chem:get-ionization a2) -1))
                   (and (= (chem:get-ionization a1) -1) (= (chem:get-ionization a2) 0)))
               (set-bond-type bond 09))     ; 0- polarized bond
              (t (set-bond-type bond 1)))) ; simple single bond
           ((eq (chem:get-order bond) :double-bond)
            (set-bond-type bond 02))          ; simple double bond
           (t (set-bond-type bond 03))))))  ; triple bond
   mol)
 #+(or)  (chem:map-bonds
   nil
   (lambda (a1 a2 order)
     (when (not (has-bond-type bond))
       (let ((a1 (chem:get-atom1 atom1))
             (a2 (chem:get-atom2 atom2)))
         (cond
           ((= (chem:get-order bond) :single-bond)
            (cond
              ((or (and (= (chem:get-ionization a1) -1) (= (chem:get-ionization a2) +1))
                   (and (= (chem:get-ionization a1) 1) (= (chem:get-ionization a2) -1)))
               (set-bond-type bond 06)) ; +- polarized bond
              ((or (and (= (chem:get-ionization a1) 0) (= (chem:get-ionization a2) -1))
                   (and (= (chem:get-ionization a1) -1) (= (chem:get-ionization a2) 0)))
               (set-bond-type bond 09))     ; 0- polarized bond
              (t (set-bond-type bond 1)))) ; simple single bond
           ((= (chem:get-order bond) :double-bond)
            (set-bond-type bond 02))          ; simple double bond
           (t (set-bond-type bond 03))))))  ; triple bond
   mol)
   )
  


;; CSMARTS rules for carbon atom types
(defvar *cx3=c* (chem:compile-smarts "[C&X3]=C"))
(defvar *cx3=norp* (chem:compile-smarts "[C&X3]=[N,P]"))
(defvar *cx3=oors* (chem:compile-smarts "[C&X3]=[O,S]"))
(defvar *cx3-ox2_cx3-nx2* (chem:compile-smarts "[$([C&X3]-[O&X2]),$([C&X3]~[N&X2])]"))

#|
(defun make-chem-info (&key smarts)
  (let ((ci (core:make-cxx-object 'chem:chem-info)))
    (chem:compile-smarts ci smarts)
    ci))
|#
(defun apply-carbon-atom-types (a)
  (let ((num-bonds (chem:number-of-bonds a)))
    (cond
     ((= num-bonds 4) (set-am1-bcc-type a 11))
     ((/= num-bonds 3) (set-am1-bcc-type a 15))
     ((not (aromaticity:is-aromatic a))
      (cond
       ((chem:matches *cx3=c* a)    (set-am1-bcc-type a 12)) 
       ((chem:matches *cx3=norp* a) (set-am1-bcc-type a 13))
       ((chem:matches *cx3=oors* a)  (set-am1-bcc-type a 14))
       (t
        (error "Could not assign type"))))
     (t
      (cond
       ((chem:matches *cx3-ox2_cx3-nx2* a) (set-am1-bcc-type a 17))
       (t (set-am1-bcc-type a 16)))))))



(defvar *match-lactone-lactam* (chem:compile-smarts "[$(O=[C]1-[O,N]~[*]~[*]1),$(O=[C]1-[O,N]~[*]~[*]~[*]1),$(O=[C]1-[O,N]~[*]~[*]~[*]~[*]1),$(O=[C]1-[O,N]~[*]~[*]~[*]~[*]~[*]1)]"))

(defvar *match-ester-acid* (chem:compile-smarts "O=C-[O]"))

(defun apply-oxygen-atom-types (a)
  (cond
   ((= (chem:number-of-bonds a) 2)         (set-am1-bcc-type a 31))
   ((chem:matches *match-lactone-lactam* a) (set-am1-bcc-type a 33))
   ((chem:matches *match-ester-acid* a)     (set-am1-bcc-type a 32))
   (t (set-am1-bcc-type a 31))))

(defvar *n-highly-deloc* (chem:compile-smarts "N1(-[C]=[C])-[C]=[C]1"))
(defvar *ndeloc* (chem:compile-smarts "N-[C]=[O]"))


;; Figure out what these properties deloc and hdeloc really mean
(defun is-highly-deloc (a) (chem:matches *n-highly-deloc* a))
(defun is-deloc (a) (chem:matches *ndeloc* a))


(defun apply-nitrogen-atom-types (a)
  (let ((num-bonds (chem:number-of-bonds a)))
    (cond
     ((= num-bonds 4) (set-am1-bcc-type a 21))
     ((= num-bonds 3)
      (cond
       ((is-highly-deloc a) (set-am1-bcc-type a 23))
       ((is-deloc a)       (set-am1-bcc-type a 22))
       (t              (set-am1-bcc-type a 21))))
     ((/= num-bonds 2) (set-am1-bcc-type a 25))
     (t
      (cond
       ((>= (chem:get-ionization a) 0) ; not anion
        (if (chem:has-bond-with-order a :double-bond)
            (set-am1-bcc-type a 24)
          (set-am1-bcc-type a 25)))
       ((is-deloc a) (set-am1-bcc-type a 22))
       (t (set-am1-bcc-type a 21)))))))


(defun apply-phosphorous-atom-types (a)
  (cond
   ((= (chem:number-of-bonds a) 4) (set-am1-bcc-type a 42))
   ((/= (chem:number-of-bonds a) 3) (set-am1-bcc-type a 41))
   ((chem:has-bond-with-order a :double-bond) (set-am1-bcc-type a 42))
   (t (set-am1-bcc-type a 41))))


;; Calculate all atom types
;; See Jakalian, Jack, and Bayly • Vol. 23, No. 16 • Journal of Computational Chemistry
(defun apply-atom-types (mol)
  (chem:map-atoms
   nil
   (lambda (a)
     (let ((ae (chem:get-element a)))
       (cond
         ((eq ae :C) (apply-carbon-atom-types a))
         ((eq ae :N) (apply-nitrogen-atom-types a))
         ((eq ae :O) (apply-oxygen-atom-types a))
         ((eq ae :S)
          (let ((num-bonds (chem:number-of-bonds a)))
            (cond
              ((= num-bonds 4) (set-am1-bcc-type a 53))
              ((= num-bonds 3) (set-am1-bcc-type a 52))
              (t  (set-am1-bcc-type a 51)))))
         ((eq ae :P) (apply-phosphorous-atom-types a))
         ((eq ae :Si) (set-am1-bcc-type a 61))
         ((eq ae :H)  (set-am1-bcc-type a 91))
         ((eq ae :F)  (set-am1-bcc-type a 71))
         ((eq ae :Cl) (set-am1-bcc-type a 72))
         ((eq ae :Br) (set-am1-bcc-type a 73))
         ((eq ae :I)  (set-am1-bcc-type a 74))
         (t (warn "Unhandled element in am1-bcc charge calculation [~s]" ae)
            (set-am1-bcc-type a nil))
         )))
   mol))



(defun lookup-am1-bcc-key (a1 bond a2)
  (let ((type1 (get-am1-bcc-type a1))
        (type2 (get-am1-bcc-type a2))
        (bond-type (get-bond-type bond))
        (sign 1.0))
    (if (and type1 type2)
        (list (+ (* 10000 type1) (* 100 bond-type) type2) sign)
        #+(or) (progn
          (when ( > type1 type2 )
            (let (sw)
              (setq sw a1) (setq a1 a2) (setq a2 sw)
              (setq type1 (get-am1-bcc-type a1))
              (setq type2 (get-am1-bcc-type a2))
              (setq sign -1.0)))
          (list (+ (* 10000 type1) (* 100 bond-type) type2) sign))
        (list nil sign))))

(defun lookup-am1-bcc-correction (a1 bond a2)
  (let* ((key-sign (lookup-am1-bcc-key a1 bond a2))
         (bcc-correction 0.0)
         (key (car key-sign))
         (sign (cadr key-sign)))
   ;; (format t "key ~a~%" key)
    (if key
        (progn
          (loop for pair in *am1-bcc-raw-data*
                do (when (= (car pair) key)
                     (setf bcc-correction (second pair))))
          (setf *am1-bcc-parameters* bcc-correction))
        0.0 #| There was no correction - so return 0.0 |#)))

;;
;; Given a molecule (mol) assign the Am1Bcc atom and bond types
;; these type values are put into the 'am1-bcc-type property for each
;; atom and bond
(defun calculate-bcc-corrections (mol)
  "Return an ObjectMap from Atoms to Am1Bcc charge corrections"
  (aromaticity:with-aromaticity-information (mol)
    (let ((correction-map (make-hash-table)))
      (assign-am1-bcc-bond-types mol)
      (apply-atom-types mol)
      (chem:map-atoms
       nil
       (lambda (a)
         ;;       (let ((atom-correction 0.0))
         (loop for b in (chem:bonds-as-list a)
               do (let* ((previous-atom-correction 0.0)
                         (previous-atom-correction-other 0.0)
                         (atom-correction 0.0)
                         (atom-correction-other 0.0)
                         (other (chem:get-other-atom b a))
                         (bond-correction (lookup-am1-bcc-correction a b other)))
                    (if (gethash a correction-map)
                        (setf previous-atom-correction (gethash a correction-map)))
                    (if (gethash other correction-map)
                        (setf previous-atom-correction-other (gethash other correction-map))) 
                    (setq atom-correction (+ previous-atom-correction bond-correction))
                    (setq atom-correction-other (- previous-atom-correction-other bond-correction)) 
                    (setf (gethash a correction-map) atom-correction)
                    (setf (gethash other correction-map) atom-correction-other))))
       mol)
      ;;    (loop for a in (all-atoms-as-cons mol nil)
      ;;         do (let ((atom-correction 0.0))
      ;;           (loop for b in (chem:bonds-as-list a)
      ;;                do (let* ((other (chem:get-other-atom b a))
      ;;                          (bond-correction (lookup-am1-bcc-correction a b other)))
      ;;                  (setq atom-correction (+ atom-correction bond-correction))))
      ;;           do (extend correction-map a atom-correction)))
      correction-map
      )))




;;; This takes a hashtable of atom->am1charges and another of atoms->bcc-corrections
;;;   and adds the am1charge to the bcc-correction and creates another hash-table
;;;   of atom->(+ am1charge bcc-correction)
(defun combine-am1-bcc-charges (atoms-to-am1-charges atoms-to-bcc-corrections)
  "Combine Am1 charges with Bcc corrections and return an ObjectMap of atoms to am1bcc charges"
  (let ((atoms-to-am1-bcc-charges (make-hash-table)))
    (maphash (lambda (atom am1charge)
                (let ((bcc-correction (gethash atom atoms-to-bcc-corrections)))
                  (setf (gethash atom atoms-to-am1-bcc-charges) (+ am1charge bcc-correction))))
                atoms-to-am1-charges)
    atoms-to-am1-bcc-charges))

(defun calculate-am1-bcc-charges (aggregate &key (maxcyc 9999))
  "Calculate Am1-Bcc charges and add the results to the aggregate."
  (let ((bcc (calculate-bcc-corrections aggregate))
        (order (charges:write-sqm-calculation (open "/tmp/sqm-input.txt" :direction :output) aggregate
                                              :maxcyc maxcyc))
        (output-file-name "/tmp/sqm-output.out"))
    (ext:vfork-execvp (list (namestring (translate-logical-pathname #P"amber:bin;sqm"))
                            "-O"
                            "-i" "/tmp/sqm-input.txt"
                            "-o" output-file-name))
    (unless (probe-file "/tmp/sqm-output.out")
      (error "The AM1 charge calculation using the external sqm program did not generate output ~a - did sqm run?" output-file-name))
    (let* ((am1 (charges:read-am1-charges "/tmp/sqm-output.out" order))
           (am1-bcc (combine-am1-bcc-charges am1 bcc)))
      (maphash (lambda (atom charge)
                 (chem:set-charge atom charge))
               am1-bcc))))

(defun calculate-mopac-am1-bcc-charges (aggregate &key (version :2016))
  "Calculate Am1-Bcc charges (Am1 charges are calculated in mopac) and add the results to the aggregate."
  (let ((bcc (calculate-bcc-corrections aggregate))
        (order (charges:write-mopac-calculation (open "/tmp/mopac.mop" :direction :output) aggregate)))
    (ext:vfork-execvp (list "/opt/mopac/MOPAC2016.exe"
                            "/tmp/mopac.mop"))
    (unless (probe-file "/tmp/mopac.out")
      (error "The AM1 charge calculation using the mopac program did not generate mopac.out - did mopac run?"))
    (let* ((am1 (charges:read-mopac-am1-charges "/tmp/mopac.out" order version))
           (am1-bcc (combine-am1-bcc-charges am1 bcc)))
      (maphash (lambda (atom charge)
                 (chem:set-charge atom charge))
               am1-bcc))))

#||
(defun describe-fragment-and-residue (frag res am1-charges bcc-corrections am1bcc-charges)
  (loop for a in (contents-as-cons res)
     do (print (% "        atom: %s" (describe-atom a))))
  (print "Describing bonds")
  (chem:map-atoms nil (lambda (a)
                        
                        (format t "Atom ~a ----~%" (chem:get-name a)) ;(print (% "Atom %s -----" (chem:get-name a)))
                        (loop for b in (chem:bonds-as-list a)
                           do (let* ((other (chem:get-other-atom b a))
                                     (key-sign (lookup-am1-bcc-key a b other)))
                                (format t " am1-bcc bond key: ~a~%" key-sign))) ; (print (% " am1-bcc bond key: %s" (repr key-sign)))))
                        (let ((am1-charge (lookup am1-charges a))
                              (bcc-correction (lookup bcc-corrections))
                              (charge (lookup am1bcc-charges a)))
                          (print ( % "   am1[%8.3f] bcc-correction[%8.3f] am1bccCharge[%8.3f]" am1-charge bcc-correction charge))))
                  res))
||#



#||
(when (is-top-level-script)
  (let ((cd (chem-draw :fileName "aspirin.cdxml")))
    ;;
    ;; Put the residues into molecules and
    ;; start classifying the atoms
    ;;
    (loop for frag in (chem:get-fragments cd)
        do  (let* ((mol (moleule))
                (res (get-entire-residues frag)))
           (let* ((name (chem:matter-get-propterty frag :name)))
             (print "-------------------")
             (print (% "Fragment: %s" frag))
             (print (% "    name: %s" name))
             (print (% "    properties====="))
             (print (% "%s" (get-propterties frag))))
           (chem:add-matter mol res)
           (build-molecule mol)
           (let ((bcc-corrections (calculate-am1-bcc-corrections mol))))
           (describe-fragment-and-residue frag res)
           ))))
||#
