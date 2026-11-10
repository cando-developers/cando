#+tests(in-package #:clasp-tests)

(defparameter agg (cando:load-mol2 "sys:extensions;cando;src;lisp;regression-tests;data;hexapeptide-smirnoff.mol2"))
(chem:setf-force-field-name (cando:mol agg 0)  :smirnoff)

(leap:load-smirnoff-params (probe-file "sys:extensions;cando;src;lisp;regression-tests;data;force-field.offxml"))


(core:set-simd-width 1)

(defparameter ef (chem:make-energy-function :matter agg))


(defparameter pos (chem:make-nvector (chem:get-nvector-size ef)))
(defparameter force (chem:make-nvector (chem:get-nvector-size ef)))

(chem:load-coordinates-into-vector ef pos)
(format t "pos = ~s~%" pos)
#+(or)(progn
;;  (gctools:wait-for-user-signal "signal")
  (time (dotimes (i 1) (defparameter energy (chem:evaluate-energy-force ef pos :calc-force t :force force)))))
(format t "About to calculate energy~%")
(defparameter energy (chem:evaluate-energy-force ef pos :calc-force t :force force))
(defparameter expected-energy 110.45436709880707d0)
(format t "         energy = ~f~%" energy)
(format t "expected-energy = ~f~%" expected-energy)

(defparameter delta (abs (- energy expected-energy)))

(format t "delta = ~f~%" delta)
;; These two tests are not yet reliable.
#+tests(test-true energy-delta (< delta 0.2))

(defparameter correct-energy-stretch 12.748283740690873)
(defparameter correct-energy-angle 83.66812723256757)
(defparameter correct-energy-dihedral 19.549603494795992)
(defparameter correct-energy-nonbond -5.538841840435541d0)
(let ((components (chem:all-components ef))
      (ht (make-hash-table)))
  (format t "Components ~s~%" components)
  (loop for comp in components
        for energy = (chem:energy-component-evaluate-energy ef comp pos)
        for class-name = (class-name (class-of comp))
        do (setf (gethash class-name ht) energy)
        do (format t "~s -> ~f~%" (class-name (class-of comp)) energy))
  #+tests(test-true energy-stretch  (< (abs (print (- correct-energy-stretch (gethash 'chem:energy-stretch ht)))) 0.01))
  #+tests(test-true energy-angle    (< (abs (print (- correct-energy-angle (gethash 'chem:energy-angle ht)))) 0.01))
  #+tests(test-true energy-dihedral (< (abs (print (- correct-energy-dihedral (gethash 'chem:energy-dihedral ht)))) 0.01))
  #+tests(test-true energy-nonbond  (< (abs (print (- correct-energy-nonbond (gethash 'chem:energy-nonbond ht)))) 0.01))
  )


(defun hand-calculate-stretch-energy (stretch-component pos)
  (let ((ebtotal 0.0))
    (chem:walk-stretch-terms
     stretch-component
     (lambda (idx a1 a2 i3x1 i3x2 kb r0)
       (let* ((pos1 (geom:vec-array pos i3x1))
              (pos2 (geom:vec-array pos i3x2))
              (dist (geom:vlength (geom:v- pos1 pos2)))
              (eb (* kb (expt (- dist r0) 2))))
         (incf ebtotal eb)
         (format t "idx ~3d  a1 = ~10s a2 = ~10s kb=~f r0=~f eb=~f~%" idx (chem:get-name a1) (chem:get-name a2) kb r0 eb))))
    (format t "ebtotal = ~f~%" ebtotal)
    ebtotal))

#+tests(let ((hand-stretch-energy (hand-calculate-stretch-energy (chem:get-stretch-component ef) pos)))
         (test-true energy-stretch-terms-match (< (abs (- hand-stretch-energy correct-energy-stretch)) 0.1)))


(defparameter expected-force (coerce
                              #(0.731319310655266d0 -0.7843311796634809d0 0.11428318967122131d0
                                -0.34291397274597024d0 0.3276675700150138d0 -0.32990218672574967d0
                                -0.7110660751472764d0 0.6592797769984873d0 0.4034257747458261d0
                                0.6315708875902306d0 -0.454396848652765d0 0.1278687132577171d0
                                -0.026085905116969857d0 -0.10100959824016205d0 -0.2746502892395881d0
                                0.4059262114094644d0 -0.29610550184605605d0 0.053861437311081506d0
                                -0.9567373033867423d0 -0.2651436300796668d0 -0.31254257619849257d0
                                -0.7439880744349008d0 -0.08287205229295314d0 -0.4481724660319921d0
                                0.8060150098147107d0 0.3450456953402482d0 0.3584798419169815d0
                                0.03948980746933861d0 0.09860133198526544d0 0.19604315201816114d0
                                -0.5634399839851341d0 0.28022883026556d0 -0.12994819225705068d0
                                0.6409832024283713d0 -0.5687773607303566d0 -0.8980260426520714d0
                                -0.713421211564184d0 -0.31687749710773194d0 0.2662266623983198d0
                                -0.49290197442539885d0 0.36761951422076233d0 0.7816662277554582d0
                                -0.16073989446191095d0 0.4319870599231008d0 0.45570497932117293d0
                                0.15083685578018446d0 -0.3122485004991003d0 -0.2606601022344863d0
                                0.4365321693331242d0 -0.1947140325176802d0 -0.7115330294162182d0
                                0.4034004799020559d0 0.06441437196574963d0 -0.8098411084063514d0
                                0.13554466607241777d0 0.07415253106130781d0 0.039742702094795884d0
                                0.19653687533440367d0 0.2123236479738253d0 -0.20829397089059842d0
                                -0.6014614483076015d0 0.3302195146639502d0 0.21775149852190426d0
                                -0.07361895416343245d0 -0.4782162555468972d0 0.5924793130095125d0
                                -0.13751781686270625d0 -0.011136987795237863d0 0.1874776898688073d0
                                0.5070870617260752d0 -0.08970240106698535d0 0.05380567478744182d0
                                -0.4823071620960435d0 0.5445806740856343d0 -0.5089601029182287d0
                                0.17238158892887556d0 -0.0736436775136875d0 0.21107475625716465d0
                                -0.1421739674939752d0 -0.5071506944268407d0 -0.4407616141530343d0
                                0.3608545814821641d0 -0.11667894492417563d0 0.3155008000155455d0
                                0.5275612133491934d0 -0.09186762069473964d0 0.607566785409247d0
                                0.14243774611747095d0 0.4088063065731248d0 -0.14201175711974917d0
                                -0.6331397502489926d0 0.14964879158389288d0 0.25065251789692444d0
                                53.70074041697893d0 -12.828397087479281d0 35.65679140245699d0
                                -0.32043306127042337d0 -0.10354896052204143d0 0.034031022109602294d0
                                -0.04884730245004373d0 -0.015337940554351105d0 -0.05743713248965603d0
                                0.665762856999269d0 -0.5069683295534729d0 0.06573491856447754d0
                                -16.407326997387656d0 -30.21932276734742d0 40.15748108104006d0
                                40.98949135347185d0 24.699676017303148d0 -22.633620198781188d0
                                -1.1414589829763677d0 0.4468671535394476d0 0.0013078143954289162d0
                                -57.484854753494865d0 -5.237974044187037d0 -4.969474443276542d0
                                -0.6813647555169753d0 -0.14114596461455697d0 0.5702190768990647d0
                                -20.73064096756879d0 24.60883055707393d0 -47.82323867742666d0
                                0.3838071300554763d0 0.20526511215917154d0 -0.4437636225168631d0
                                0.2397514445786697d0 0.1566922921365812d0 -0.18433642437022746d0
                                0.6950221241816816d0 0.23467537676404593d0 -0.1610399407006362d0
                                -1.0826245915482442d0 1.2320855251416125d0 0.677655993177801d0
                                -0.44599734505265926d0 0.49947063883085896d0 -0.11666577108309496d0
                                0.15542240961679799d0 -0.3333357932517944d0 -0.14005714027631766d0
                                0.42924257844381d0 -0.8788763352914799d0 -0.3453331520575605d0
                                -0.30745835648659287d0 0.17770784094393183d0 -0.2877437237684371d0
                                -0.052632810379684064d0 0.26684871804379573d0 0.011289360405689328d0
                                0.08866446654085117d0 -0.08340840235338687d0 0.04940569151470878d0
                                0.09466895626189925d0 -0.07769712408212211d0 0.2864498098820425d0
                                -0.01637603851325932d0 0.257998936686512d0 0.040817281957742736d0
                                0.06963829126401855d0 -0.025358537018870503d0 0.0033217749176716847d0
                                0.6143212050329835d0 0.13749490619807736d0 0.020859028382874203d0
                                0.1306899393464552d0 0.24114585593997515d0 0.0476619044113666d0
                                -0.22894175410104325d0 -0.18573900637826207d0 0.2234762138723116d0
                                0.1301423550157204d0 -0.21333232397931096d0 0.15489064502339583d0
                                0.4491742896745319d0 0.13520188728200389d0 0.07763242576119632d0
                                0.042679163208733426d0 -0.006245875974112192d0 -0.07511336248090159d0
                                -0.09237613965744684d0 0.11150346510564789d0 -0.13141177547696062d0
                                0.009290860042796875d0 -0.06634418348260611d0 0.11346828489428795d0
                                -0.49066840692151137d0 0.07430670361705963d0 0.7034596448634842d0
                                0.026328745598186168d0 -0.009845742577265237d0 -0.04555989368471433d0
                                0.22004318595871433d0 -0.07543966126180118d0 -0.04057604028891251d0
                                0.2696796654587654d0 -0.002604136535056298d0 -0.4345343184142471d0
                                -0.016201423428105932d0 -0.10027583360144697d0 0.013785688700951038d0
                                -0.10995079895003d0 -0.15986613169716155d0 -0.14031827801100233d0
                                -0.0010338700485759689d0 -0.005371956767050726d0 0.06033185786703074d0
                                -0.09821368246926804d0 -0.009662144937015768d0 -0.07589270949798878d0
                                -0.10066657368825285d0 0.15231063421847313d0 -0.28345644468713627d0
                                -0.33284042132474767d0 -0.04316273427569363d0 0.017128041705433316d0
                                0.1523904805406699d0 0.1572295157854951d0 -0.0497765941358859d0
                                -0.003598988368529852d0 0.012946571099886749d0 0.0153633453387422d0
                                0.13172046021317615d0 -0.2468776550411409d0 -0.395924736868541d0
                                -0.05262651137484208d0 -0.22564132542377371d0 -0.1400416610756559d0
                                0.14197302815520346d0 -0.12168117058297656d0 -0.10307993393558534d0
                                0.9451104540789164d0 -0.08770612184211582d0 0.047526759076065905d0
                                -0.25704342780207834d0 -0.16735891923718188d0 -0.06364899328428969d0
                                0.07081856316434666d0 -0.13686871821804997d0 -0.0850680088378832d0
                                0.046073543873994344d0 -0.06185109700046438d0 0.049214921101621456d0
                                -0.06260520634426621d0 -0.003514409099446031d0 0.13436866936254913d0
                                -0.11012121498133326d0 -0.2791097404847669d0 0.0544878559042824d0
                                -0.09674906900236632d0 -0.2174380282624657d0 0.1688519379584673d0
                                0.22395943143249356d0 -0.08002437339734723d0 -0.09523690223317229d0
                                0.03303621398728089d0 -0.12702968596135703d0 0.2114429784359123d0
                                0.008791197383839608d0 -0.3592198046880002d0 -0.07373597393385999d0
                                0.10825446755820213d0 0.08562152603138089d0 -0.030677854442546174d0)
                              (list 'vector (geom:vecreal-type))))
(format t "force = ~s~%" force)

(defparameter expected-force-mag (chem:nvector-magnitude expected-force))
(defparameter force-mag (chem:nvector-magnitude force))

(format t "         force-mag = ~f~%" force-mag)
(format t "expected-force-mag = ~f~%" expected-force-mag)
#+tests(test-true force-mag (< (abs (- expected-force-mag force-mag)) 1.0))

;;; Calculate the angle between expected-force and force - they should be parallel
;;; so the angle should be near zero
(let* ((dotdot (chem:nvector-dot expected-force force))
       (dd1 (/ dotdot expected-force-mag))
       (cos-angle (/ dd1 force-mag)))
  (format t "dotdot = ~f~%" dotdot)
  (format t "dd1 = ~f~%" dd1)
  (format t "cos-angle = ~f~%" cos-angle)

  (defparameter force-acos (acos (min 1.0 (max cos-angle 0.0))))
  (format t "force-acos = ~f (angle between force vectors in radians)~%" force-acos)

  #+tests(test-true force-acos (< force-acos 0.01)))


(defparameter minimizer (chem:make-minimizer ef))
#+(or)(time (dotimes (i 10)
        (chem:save-coordinates-from-vector ef pos)
        (chem:minimize minimizer)))


(defparameter spiro (cando:load-mol2 (merge-pathnames "data/struct-0000.mol2"
                                                      *load-truename*)))

(chem:setf-force-field-name (cando:mol spiro 0) :smirnoff)

(defparameter ef2 (chem:make-energy-function :matter spiro))

(defparameter mask (chem:matter/make-mask spiro))

(chem:matter/walk-atoms spiro (lambda (index atm res mol)
                                (when (eq :pro4ss (chem:get-name res))
                                  (setf (aref mask index) 1))))

(defparameter pos2 (chem:make-nvector (chem:get-nvector-size ef2)))
(chem:load-coordinates-into-vector ef2 pos2)


(defun energy-components (name en pos mask)
  (let ((components (chem:all-components en))
        (ht (make-hash-table)))
    (format t "~s energy total = ~f~%" name (chem:evaluate-energy en pos :active-atom-mask mask))
    (loop for comp in components
          for energy = (chem:energy-component-evaluate-energy
                        en comp pos :active-atom-mask mask)
          for class-name = (class-name (class-of comp))
          do (setf (gethash class-name ht) energy)
          do (format t "~s ~s -> ~f~%" name (class-name (class-of comp)) energy))
    ht
    )
  )

(defun energy-multiple-dihedral (num en pos &optional mask)
  (let ((components (chem:all-components en))
        (ht (make-hash-table)))
    (loop for comp in components
          for class-name = (class-name (class-of comp))
          do (setf (gethash class-name ht) comp))
    (let ((dihedral (gethash 'chem:energy-dihedral ht)))
      (loop for idx below num
            for energy = (chem:energy-component-evaluate-energy
                          en dihedral pos :active-atom-mask mask)))))

#+(or)
(progn
  (format t "Running simd1~%")
  (core:set-simd-width 1)
  (energy-components "simd1" ef2 pos2 mask)
  (format t "Running simd8~%")
  (core:set-simd-width 8)
  (energy-components "simd8" ef2 pos2 mask)
  )

#+(or)
(progn
  (core:set-simd-width 1)
  (format t "no mask simd1 timing~%")
  (time (energy-multiple-dihedral 100000 ef2 pos2 mask)))

#+(or)
(progn
  (core:set-simd-width 8)
  (format t "no mask simd8 timing~%")
  (time (energy-multiple-dihedral 100000 ef2 pos2 mask)))


(let* ((butane (cando:load-mol2 (merge-pathnames "data/butane.mol2" *load-pathname*))))
  (chem:setf-force-field-name (cando:mol butane 0)  :smirnoff)
  (leap:load-smirnoff-params (probe-file "sys:extensions;cando;src;lisp;regression-tests;data;force-field.offxml"))
  (let* ((butane-ef (chem:make-energy-function :matter butane))
         (dh (chem:energy-function/get-dihedral-restraint-component butane-ef))
         (a1 (chem:matter/first-atom-with-name butane :C1))
         (a2 (chem:matter/first-atom-with-name butane :C2))
         (a3 (chem:matter/first-atom-with-name butane :C3))
         (a4 (chem:matter/first-atom-with-name butane :C4))
         (ri (chem:add-dihedral-restraint dh butane-ef
                                          100000.0 (* -180.0 0.0174533)
                                          a1 a2 a3 a4 :proper
                                          ))
         (minimizer (chem:make-minimizer butane-ef))
         )
    (format t "dh = ~s~%" dh)
    #+tests
    (flet ((rejected-p (atoms kind)
             (let ((before (chem:energy-component/number-of-terms dh)))
               (and (handler-case
                        (progn
                          (apply #'chem:add-dihedral-restraint dh butane-ef 1d0 0d0
                                 (append atoms (list kind)))
                          nil)
                      (error () t))
                    (= before (chem:energy-component/number-of-terms dh))))))
      (test-true dihedral-restraint-proper-kind
                 (eq :proper (nth-value 6 (chem:energy-dihedral-restraint/get-dihedral-restraint dh ri))))
      (test-true dihedral-restraint-rejects-broken-chain
                 (rejected-p (list a1 a3 a2 a4) :proper))
      (test-true dihedral-restraint-rejects-wrong-improper-center
                 (rejected-p (list a1 a2 a3 a4) :improper))
      (test-true dihedral-restraint-rejects-unknown-kind
                 (rejected-p (list a1 a2 a3 a4) :unknown))
      (test-true dihedral-restraint-rejects-repeated-atoms
                 (rejected-p (list a1 a2 a1 a2) :proper))
      (let ((neighbor nil))
        (chem:map-atoms nil
                       (lambda (atom)
                         (when (and (not (member atom (list a1 a2 a3)))
                                    (chem:is-bonded-to a2 atom))
                           (setf neighbor atom)))
                       butane)
        (assert neighbor)
        ;; Atom 3 is central; weight zero avoids altering the following scan.
        (let ((ii (chem:add-dihedral-restraint dh butane-ef 0d0 0d0 a1 a3 a2 neighbor :improper)))
          (chem:energy-dihedral-restraint/update-dihedral-restraint dh ii 0d0 0.1d0)
          (test-true dihedral-restraint-improper-kind-survives-update
                     (eq :improper (nth-value 6 (chem:energy-dihedral-restraint/get-dihedral-restraint dh ii)))))))
    (with-open-file (sdf-stream "/tmp/butane.sdf" :direction :output :if-exists :supersede)
      (loop for angle from -180 by 20 below 180
            do (format t "angle = ~d~%" angle)
            do (chem:update-dihedral-restraint dh ri 100000.0 (* angle 0.0174533))
            do (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
                 (chem:minimize minimizer))
            do (sdf:write-sdf-stream butane sdf-stream))))
  (defparameter *butane* butane)
  )


;;;(ql:quickload :open-force-field/load)

(let* ((agg (chem:sample-aggregate))
       (ef (chem:make-energy-function :matter agg :use-excluded-atoms nil))
       (coords (chem:matter/extract-coordinates agg))
       (components (chem:make-energy-components))
       (energy (chem:scoring-function/evaluate-all ef coords :energy-components components))
       (expected-energy 669.7740554130356d0)
       (expected-components '((CHEM::ENERGY-STRETCH . 95.7978)
                              (CHEM::ENERGY-ANGLE . 301.018)
                              (CHEM::ENERGY-DIHEDRAL . 230.37)
                              (CHEM::ENERGY-NONBOND-TOTAL . -55.4358)
                              (CHEM::ENERGY-NONBOND14 . 75.18347128717448d0)
                              (CHEM::ENERGY-CHIRAL-RESTRAINT . 0)
                              (CHEM::ENERGY-ANCHOR-RESTRAINT . 0)
                              (CHEM::ENERGY-DIHEDRAL-RESTRAINT . 22.840712650543228d0)
                              (CHEM::MONOMER-CORRECTIONS . 0)))
       (nonbond-component (chem:get-nonbond-component ef))
       (dihedral-component (chem:get-dihedral-component ef))
       (angle-component (chem:get-angle-component ef))
       (stretch-component (chem:get-stretch-component ef))
       )
  (format t "components = ~s~%" components)
  (format t "nonbond description ~s~%" (chem:energy-component/implementation-details nonbond-component))
  (format t "dihedral description ~s~%" (chem:energy-component/implementation-details dihedral-component))
  (format t "Total energy = ~f  expected energy: ~f~%" energy expected-energy)
  #+tests(test-true sample-agg-energy (< (abs (- energy expected-energy)) 0.001))
  (loop for (comp . en) in (chem:energy-components/components components)
        for expected-en = (cdr (assoc comp expected-components))
        for index from 0
        for name = (intern (format nil "energy-rosetta-nonbond-~d" index))
        for new-name = (intern (format nil "new-energy-rosetta-nonbond-~d" index))
        if expected-en
          do (progn
               (format t "~s (~a)-> energy: ~f  expected: ~f ~%" comp
                       (if (< (abs (- expected-en en)) 0.001)
                           "matched"
                           "MISMATCHED")
                       en expected-en)
               #+tests(%test name `(< (abs (- ,expected-en ,en)) 0.001) (lambda () (not (not (< (abs (- expected-en en)) 0.001)))) '(t)))
        else
          do (progn
               (warn "There is a new energy component ~s that isn't in the regression test" comp)
               #+tests(test-true new-name t)))
  (format t "About to run gradient and hessian finite-difference tests~%")
  (chem:run-test-calls nonbond-component *standard-output* coords)
  (chem:run-test-calls dihedral-component *standard-output* coords)
  (chem:run-test-calls angle-component *standard-output* coords)
  (chem:run-test-calls stretch-component *standard-output* coords)
  #|
  #+tests(nonbond-fd-tests (= 0 (chem:run-test-calls nonbond-component *standard-output* coords)))
  #+tests(dihedral-fd-tests (= 0 (chem:run-test-calls dihedral-component *standard-output* coords)))
  #+tests(angle-fd-tests (= 0 (chem:run-test-calls angle-component *standard-output* coords)))
  #+tests(stretch-fd-tests (= 0 (chem:run-test-calls stretch-component *standard-output* coords)))
|#
  )


;;; ======================================================================
;;; Regression tests for the restraint energy kernels:
;;;   anchor restraint, flat-bottom anchor restraint, chiral restraint.
;;;
;;; The flat-bottom and chiral restraints are ONE-SIDED, so each is tested
;;; in BOTH domains:
;;;   flat-bottom : atom inside R0 (energy 0) vs outside R0 (harmonic)
;;;   chiral      : the active cubic domain (E>0) vs the clamped domain (E=0)
;;;
;;; The restraint parameters (center / R0 / K / CO) are chosen so the
;;; expected energy is exact and independent of the butane geometry:
;;;   - placing the anchor center a known OFFSET from the atom fixes the
;;;     distance exactly,
;;;   - the normalized chiral measure Q satisfies |Q| <= 1, so |CO| > 1
;;;     makes the sign of (CO+Q) - hence which chiral domain is active -
;;;     depend only on CO and the sign of K.
;;; ======================================================================

(defun %fresh-butane-ef ()
  "Load a fresh butane aggregate with smirnoff parameters.
   Returns (values energy-function aggregate)."
  (let ((butane (cando:load-mol2 (merge-pathnames "data/butane.mol2" *load-truename*))))
    (chem:setf-force-field-name (cando:mol butane 0) :smirnoff)
    (values (chem:make-energy-function :matter butane) butane)))

;;; ---- Anchor restraint:  E = ka * |r - center|^2 ----------------------
(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (let* ((coords (chem:make-nvector (chem:get-nvector-size ef)))
         (c1 (chem:matter/first-atom-with-name butane :C1))
         (anchor (chem:make-energy-anchor-restraint ef))
         (ka 100.0d0)
         (offset 3.0d0))
    (chem:load-coordinates-into-vector ef coords)
    ;; center placed exactly OFFSET away from the atom along x -> distance = offset
    (let* ((center (geom:v+ (chem:get-position c1) (geom:vec offset 0.0d0 0.0d0)))
           (expected (* ka offset offset)))
      (chem:add-anchor-restraint-term anchor ef c1 center ka)
      (let ((e (chem:energy-component-evaluate-energy ef anchor coords)))
        (format t "anchor energy = ~f  expected = ~f~%" e expected)
        #+tests(test-true anchor-restraint-energy (< (abs (- e expected)) 1d-6))
        ;; analytic gradient/hessian vs finite-difference
        (format t "anchor FD test:~%")
        #+tests(test-true anchor-restraint-fd
                          (= 0 (chem:run-test-calls anchor *standard-output* coords)))))))

;;; ---- Flat-bottom anchor restraint: domain 1 (inside R0) -> energy 0 ---
(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (let* ((coords (chem:make-nvector (chem:get-nvector-size ef)))
         (c1 (chem:matter/first-atom-with-name butane :C1))
         (fb (chem:make-energy-flat-bottom-anchor-restraint ef))
         (ka 100.0d0)
         (r0 1.0d0))
    (chem:load-coordinates-into-vector ef coords)
    ;; atom is 0.2 from the center, R0 = 1.0  ->  inside flat region -> E = 0
    (let ((center (geom:v+ (chem:get-position c1) (geom:vec 0.2d0 0.0d0 0.0d0))))
      (chem:add-flat-bottom-anchor-restraint-term fb ef c1 center r0 ka)
      (let ((e (chem:energy-component-evaluate-energy ef fb coords)))
        (format t "flat-bottom INSIDE energy = ~f  expected = 0~%" e)
        #+tests(test-true flat-bottom-inside-domain (< (abs e) 1d-9))))))

;;; ---- Flat-bottom anchor restraint: domain 2 (outside R0) -> harmonic -
(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (let* ((coords (chem:make-nvector (chem:get-nvector-size ef)))
         (c1 (chem:matter/first-atom-with-name butane :C1))
         (fb (chem:make-energy-flat-bottom-anchor-restraint ef))
         (ka 100.0d0)
         (r0 1.0d0)
         (offset 5.0d0))
    (chem:load-coordinates-into-vector ef coords)
    ;; atom is 5.0 from the center, R0 = 1.0 -> outside -> E = ka*(d-R0)^2
    (let* ((center (geom:v+ (chem:get-position c1) (geom:vec offset 0.0d0 0.0d0)))
           (expected (* ka (expt (- offset r0) 2))))
      (chem:add-flat-bottom-anchor-restraint-term fb ef c1 center r0 ka)
      (let ((e (chem:energy-component-evaluate-energy ef fb coords)))
        (format t "flat-bottom OUTSIDE energy = ~f  expected = ~f~%" e expected)
        #+tests(test-true flat-bottom-outside-domain (< (abs (- e expected)) 1d-6))
        ;; analytic gradient/hessian vs finite-difference (harmonic branch)
        (format t "flat-bottom FD test:~%")
        #+tests(test-true flat-bottom-restraint-fd
                          (= 0 (chem:run-test-calls fb *standard-output* coords)))))))

;;; ---- Chiral restraint: one-sided  E = K*(CO+Q)^3 when positive, else 0.
(defun %chiral-energy (k co)
  "Build a fresh butane energy function with a single chiral restraint term
   over C1-C2-C3(center)-C4 and return (values energy component coords)."
  (multiple-value-bind (ef butane) (%fresh-butane-ef)
    (let ((coords (chem:make-nvector (chem:get-nvector-size ef)))
          (cc (chem:make-energy-chiral-restraint ef))
          (a1 (chem:matter/first-atom-with-name butane :C1))
          (a2 (chem:matter/first-atom-with-name butane :C2))
          (a3 (chem:matter/first-atom-with-name butane :C3))
          (a4 (chem:matter/first-atom-with-name butane :C4)))
      (chem:load-coordinates-into-vector ef coords)
      (chem:add-chiral-restraint-term cc ef a1 a2 a3 a4 k co)
      (values (chem:energy-component-evaluate-energy ef cc coords) cc coords))))

(let ((w 100.0d0))
  ;; CO = +2 : (CO+Q) > 0  ->  K>0 active (E>0),  K<0 clamped to 0
  ;; CO = -2 : (CO+Q) < 0  ->  K<0 active (E>0),  K>0 clamped to 0
  (let ((e-pos-active   (%chiral-energy    w     2.0d0))
        (e-pos-clamped  (%chiral-energy (- w)    2.0d0))
        (e-neg-active   (%chiral-energy (- w)  (- 2.0d0)))
        (e-neg-clamped  (%chiral-energy    w   (- 2.0d0))))
    (format t "chiral CO=+2 K>0 (active)  = ~f~%" e-pos-active)
    (format t "chiral CO=+2 K<0 (clamped) = ~f~%" e-pos-clamped)
    (format t "chiral CO=-2 K<0 (active)  = ~f~%" e-neg-active)
    (format t "chiral CO=-2 K>0 (clamped) = ~f~%" e-neg-clamped)
    ;; Active domain: CO+Q in [1,3] -> E = w*(CO+Q)^3 in [w, 27w], strictly > 0.
    #+tests(test-true chiral-active-domain-pos-k (> e-pos-active 0.0d0))
    #+tests(test-true chiral-active-domain-neg-k (> e-neg-active 0.0d0))
    ;; Clamped domain: the one-sided guard yields EXACTLY 0 (never negative,
    ;; which is what the old two-sided cubic would have produced here).
    #+tests(test-true chiral-clamped-domain-pos-k (< (abs e-pos-clamped) 1d-9))
    #+tests(test-true chiral-clamped-domain-neg-k (< (abs e-neg-clamped) 1d-9))
    ;; The active energy must not depend on the sign convention used to reach
    ;; it: swapping the sign of both K and CO maps Q-domain to Q-domain with the
    ;; same |CO|, so the two active energies bracket the same cubic family.
    #+tests(test-true chiral-active-finite (and (< e-pos-active 1d6)
                                                (< e-neg-active 1d6)))))

;;; analytic gradient/hessian vs finite-difference on an active chiral term
(multiple-value-bind (e cc coords) (%chiral-energy 100.0d0 2.0d0)
  (declare (ignore e))
  (format t "chiral FD test:~%")
  #+tests(test-true chiral-restraint-fd
                    (= 0 (chem:run-test-calls cc *standard-output* coords))))


;;; ======================================================================
;;; SMIRNOFF routing: a 180-degree (sp/linear) angle parameter must be
;;; assigned to the singularity-free linear-bend term (energy-linear-angle,
;;; E = kt*(1+cos theta)) rather than the standard harmonic energy-angle term.
;;;
;;; Acetylene (H-C#C-H) has two H-C#C angles, both matching SMIRNOFF a16
;;; ([*:1]~[#6X2:2]~[*:3], angle = 180 deg).  The reduced test force field has
;;; no triple-bond stretch / sp torsion parameters, so those are skipped
;;; (chem:skip-missing-parameters) while the (covered) angle is still assigned.
;;; Geometry is slightly trans-bent (~175 deg) so the finite-difference check
;;; is not evaluated exactly at the degenerate 180-degree point.
;;; ======================================================================

(leap:load-smirnoff-params (probe-file "sys:extensions;cando;src;lisp;regression-tests;data;force-field.offxml"))

(let* ((acetylene (cando:load-mol2 (merge-pathnames "data/acetylene.mol2" *load-truename*))))
  (chem:setf-force-field-name (cando:mol acetylene 0) :smirnoff)
  (let ((ef (handler-bind ((chem:parameterization-error
                             (lambda (c)
                               ;; skip the sp bond/torsion params the reduced FF lacks
                               (let ((r (find-restart 'chem:skip-missing-parameters c)))
                                 (when r (invoke-restart r))))))
              (chem:make-energy-function :matter acetylene))))
    (let* ((coords (chem:make-nvector (chem:get-nvector-size ef)))
           (angle-comp (chem:get-angle-component ef))
           ;; ensureComponent: returns the linear-angle component built during assignment
           (linear-comp (chem:make-energy-linear-angle ef))
           (n-harmonic (if angle-comp (chem:number-of-terms angle-comp) 0))
           (n-linear (chem:number-of-terms linear-comp)))
      (chem:load-coordinates-into-vector ef coords)
      (format t "acetylene: harmonic-angle terms = ~d  linear-angle terms = ~d~%"
              n-harmonic n-linear)
      ;; Both H-C#C angles (SMIRNOFF a16, theta0 = 180) route to the linear term ...
      #+tests(test-true linear-angle-routing-count (= n-linear 2))
      ;; ... and none stayed in the singular harmonic angle term.
      #+tests(test-true linear-angle-routing-harmonic-empty (= n-harmonic 0))
      (let ((e (chem:energy-component-evaluate-energy ef linear-comp coords)))
        (format t "acetylene linear-angle energy = ~f~%" e)
        ;; E = kt*(1+cos theta) >= 0, finite at/near linearity
        #+tests(test-true linear-angle-routing-energy (>= e 0.0d0)))
      (format t "linear-angle FD test:~%")
      #+tests(test-true linear-angle-routing-fd
                        (= 0 (chem:run-test-calls linear-comp *standard-output* coords))))))


;;; ======================================================================
;;; Hessian preconditioner for the truncated-Newton minimizer.
;;;
;;; chem:setup-hessian-preconditioner assembles each *bonded* component's local
;;; analytic Hessian into a sparse symmetric matrix (the preconditioner that
;;; truncated Newton LDLt-factorizes).  The nonbond / nonbond14 components
;;; contribute nothing to it.  This test verifies the assembled matrix equals a
;;; finite-difference Hessian of the bonded gradient -- catching scatter / global
;;; index bugs in the per-component setupHessianPreconditioner code, independent
;;; of the kernels themselves.
;;;
;;; NOTE: the default dihedral kernel (dihedral_fast_floored) supplies a
;;; Gauss-Newton (outer-product-only) Hessian -- it drops the
;;; dE/dphi * d2phi/dq^2 curvature term -- so dihedral-coupled entries do NOT
;;; match the exact FD Hessian (they differ by a few percent).  That is
;;; intentional (a preconditioner only needs to approximate H; truncated Newton
;;; still converges -- see the TN test below).  We therefore relax the tolerance
;;; on dihedral entries and keep it tight on the exact (stretch/angle/...) ones.
;;;
;;; H[i][j] = d(grad_i)/dq_j = -(F_i(q + h e_j) - F_i(q - h e_j)) / (2h),
;;; where F = -grad is the force returned by evaluate-energy-force.
;;; ======================================================================

(defun %nv-ref (v idx)
  "Read scalar element IDX of NVector V."
  (let* ((base (* 3 (floor idx 3)))
         (comp (mod idx 3))
         (vv (geom:vec-array v base)))
    (ecase comp
      (0 (geom:get-x vv))
      (1 (geom:get-y vv))
      (2 (geom:get-z vv)))))

(defun %nv-set! (v idx val)
  "Set scalar element IDX of NVector V to VAL."
  (let* ((base (* 3 (floor idx 3)))
         (comp (mod idx 3))
         (vv (geom:vec-array v base))
         (x (geom:get-x vv)) (y (geom:get-y vv)) (z (geom:get-z vv)))
    (ecase comp (0 (setf x val)) (1 (setf y val)) (2 (setf z val)))
    (geom:vec-set (geom:vec x y z) v base)))

(defun %rms-force (ef pos)
  "RMS of the gradient (force) of energy-function EF at coordinates POS."
  (let* ((nn (chem:get-nvector-size ef))
         (f  (chem:make-nvector nn)))
    (chem:evaluate-energy-force ef pos :calc-force t :force f)
    (/ (chem:nvector-magnitude f) (sqrt (float nn 1d0)))))

(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (declare (ignore butane))
  ;; Restrict the force to the bonded terms that actually populate the
  ;; preconditioner, so the finite-difference Hessian matches it exactly.
  (loop for comp in (chem:all-components ef)
        when (member (class-name (class-of comp))
                     '(chem:energy-nonbond chem:energy-nonbond14))
          do (chem:disable comp))
  (let* ((n (chem:get-nvector-size ef))
         (pos (chem:make-nvector n))
         (fp (chem:make-nvector n))
         (fm (chem:make-nvector n))
         (h 1d-5)
         (hfd (make-array (list n n) :initial-element 0d0)))
    (chem:load-coordinates-into-vector ef pos)
    ;; central-difference Hessian of the (bonded) gradient
    (dotimes (j n)
      (let ((orig (%nv-ref pos j)))
        (%nv-set! pos j (+ orig h))
        (chem:evaluate-energy-force ef pos :calc-force t :force fp)
        (%nv-set! pos j (- orig h))
        (chem:evaluate-energy-force ef pos :calc-force t :force fm)
        (%nv-set! pos j orig)
        (dotimes (i n)
          (setf (aref hfd i j)
                (- (/ (- (%nv-ref fp i) (%nv-ref fm i)) (* 2d0 h)))))))
    ;; assemble the full bonded preconditioner
    (let ((m (chem:make-sparse-large-square-matrix n 'chem:symmetric-diagonal-lower))
          (mismatches 0)
          (max-err 0d0)
          (dih-pairs (make-hash-table :test 'equal)))
      (chem:setup-hessian-preconditioner ef pos m nil)
      ;; The default dihedral kernel (dihedral_fast_floored) contributes a
      ;; Gauss-Newton (outer-product-only) Hessian, so its (row,col) entries
      ;; intentionally differ from the exact FD Hessian by the dropped
      ;; dE/dphi * d2phi/dq^2 curvature term.  Identify those entries by
      ;; assembling a dihedral-only preconditioner, and relax the tolerance there
      ;; (keep it tight everywhere else, where the Hessian is exact).
      (let ((dih (chem:get-dihedral-component ef)))
        (loop for comp in (chem:all-components ef)
              unless (eq comp dih) do (chem:disable comp))
        (let ((mdih (chem:make-sparse-large-square-matrix n 'chem:symmetric-diagonal-lower)))
          (chem:setup-hessian-preconditioner ef pos mdih nil)
          (chem:walk-matrix
           mdih
           (lambda (col row val)
             (when (> (abs val) 1d-12)
               (setf (gethash (cons row col) dih-pairs) t))))))
      (chem:walk-matrix
       m
       (lambda (col row val)                        ; val = stored preconditioner entry (col,row)
         (let* ((ref       (aref hfd row col))
                (err       (abs (- val ref)))
                (dihedralp (gethash (cons row col) dih-pairs))
                ;; exact entries: tight (0.1% + 1e-2).  Gauss-Newton dihedral
                ;; entries: relaxed to 10% (still catches sign/scatter bugs, which
                ;; are O(100%), while tolerating the dropped curvature term ~few%).
                (tol       (if dihedralp
                               (+ 1d-2 (* 1d-1 (abs ref)))
                               (+ 1d-2 (* 1d-3 (abs ref))))))
           (when (> err max-err) (setf max-err err))
           (when (> err tol)
             (incf mismatches)
             (format t "  precon mismatch (~d,~d)~a: precon=~f  fd=~f  err=~f~%"
                     row col (if dihedralp " [dihedral/GN]" "") val ref err)))))
      (format t "preconditioner vs FD Hessian: max-err = ~f  mismatches = ~d (dihedral entries are Gauss-Newton)~%"
              max-err mismatches)
      #+tests(test-true preconditioner-matches-fd-hessian (= mismatches 0)))))


;;; ======================================================================
;;; Truncated-Newton convergence test.
;;;
;;; Starting near a relaxed minimum, verify that truncated Newton reaches
;;; a tight gradient tolerance. Explicitly disable the default L-BFGS path.
;;; ======================================================================

(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (declare (ignore butane))
  (let ((n (chem:get-nvector-size ef))
        (tol 1d-4))
    (let* ((m0 (chem:make-minimizer ef))
           (refpos (progn
                     (chem:disable-print-intermediate-results m0)
                     (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
                       (chem:minimize m0)))))
      (dotimes (k (min 9 n))
        (%nv-set! refpos k (+ (%nv-ref refpos k) 0.15d0)))
      (let ((mtn (chem:make-minimizer ef)))
        (chem:set-use-lbfgs mtn nil)
        (chem:set-maximum-number-of-steepest-descent-steps mtn 0)
        (chem:set-maximum-number-of-conjugate-gradient-steps mtn 0)
        (chem:set-maximum-number-of-truncated-newton-steps mtn 10000)
        (chem:set-truncated-newton-tolerance mtn tol)
        (chem:disable-print-intermediate-results mtn)
        (let* ((tn-final (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
                           (chem:minimize mtn :coords refpos)))
               (tn-rms (%rms-force ef tn-final)))
          (format t "truncated-Newton: ~d steps, rms-force = ~e~%"
                  (chem:get-iteration mtn) tn-rms)
          #+tests(test-true truncated-newton-converged (< tn-rms (* 10d0 tol))))))))


;;; ======================================================================
;;; Cartesian active-coordinate masks.
;;;
;;; A mask has one bit per Cartesian coordinate.  A bonded term crossing the
;;; active/fixed boundary must remain in the objective, while minimization must
;;; change only coordinates whose bits are set.  Exercise a single movable x
;;; coordinate so this also catches accidental atom-index (i3x/3) masking.
;;; ======================================================================

(multiple-value-bind (ef butane) (%fresh-butane-ef)
  (declare (ignore butane))
  (let* ((n (chem:get-nvector-size ef))
         (pos (chem:make-nvector n))
         (mask (make-array n :element-type 'bit :initial-element 0))
         (stretch (chem:get-stretch-component ef))
         (expected 0d0))
    (chem:load-coordinates-into-vector ef pos)
    ;; Only x of the atom at i3x=0 is movable.
    (setf (aref mask 0) 1)
    (chem:walk-stretch-terms
     stretch
     (lambda (index atom1 atom2 i3x1 i3x2 kb r0)
       (declare (ignore index atom1 atom2))
       (when (or (= i3x1 0) (= i3x2 0))
         (let* ((p1 (geom:vec-array pos i3x1))
                (p2 (geom:vec-array pos i3x2))
                (distance (geom:vlength (geom:v- p1 p2))))
           (incf expected (* kb (expt (- distance r0) 2)))))))
    (let ((masked (chem:energy-component-evaluate-energy
                   ef stretch pos :active-atom-mask mask)))
      #+tests(test-true active-coordinate-mask-keeps-cross-boundary-terms
                        (< (abs (- masked expected)) 1d-8)))

    ;; Displace the one active coordinate and run projected TN.  Every other
    ;; coordinate, including y and z of the same atom, must remain bit-identical.
    (%nv-set! pos 0 (+ (%nv-ref pos 0) 0.15d0))
    (let* ((before (copy-seq pos))
           (minimizer (chem:make-minimizer ef)))
      (chem:set-maximum-number-of-steepest-descent-steps minimizer 0)
      (chem:set-maximum-number-of-conjugate-gradient-steps minimizer 0)
      (chem:set-maximum-number-of-truncated-newton-steps minimizer 10000)
      (chem:set-truncated-newton-tolerance minimizer 1d-4)
      (chem:disable-print-intermediate-results minimizer)
      (let* ((after
               (ext:with-float-traps-masked
                   (:underflow :overflow :invalid :inexact :divide-by-zero)
                 (chem:minimize minimizer :coords pos :active-atom-mask mask)))
             (energy (chem:evaluate-energy ef after)))
        #+tests(test-true projected-truncated-newton-moves-active-coordinate
                          (> (abs (- (%nv-ref after 0) (%nv-ref before 0))) 1d-10))
        #+tests(test-true projected-truncated-newton-preserves-fixed-coordinates
                          (loop for coordinate from 1 below n
                                always (= (%nv-ref after coordinate)
                                          (%nv-ref before coordinate))))
        #+tests(test-true projected-truncated-newton-energy-is-finite
                          (and (not (ext:float-nan-p energy))
                               (not (ext:float-infinity-p energy))))))))


;;; BEGIN L-BFGS REGRESSION TESTS
;;;
;;; Known harmonic minima provide an independent numerical reference.  The
;;; coupled stretch test additionally leaves a nonzero force on a fixed atom;
;;; merely masking the initial force cannot preserve that constrained minimum.

(defun %lbfgs-harmonic-fixture (&key negative)
  (multiple-value-bind (ef butane) (%fresh-butane-ef)
    (let* ((size (chem:get-nvector-size ef))
           (coords (chem:make-nvector size))
           (weights (make-array size :initial-element 0d0))
           (indices (make-hash-table :test 'eq))
           (anchor (chem:make-energy-anchor-restraint ef)))
      (chem:walk-stretch-terms
       (chem:get-stretch-component ef)
       (lambda (index atom1 atom2 i1 i2 kb r0)
         (declare (ignore index kb r0))
         (setf (gethash atom1 indices) i1 (gethash atom2 indices) i2)))
      (dolist (component (chem:all-components ef)) (chem:disable component))
      (chem:enable anchor)
      (chem:load-coordinates-into-vector ef coords)
      (let ((target (copy-seq coords)))
        (chem:matter/walk-atoms
         butane
         (lambda (index atom residue molecule)
           (declare (ignore residue molecule))
           (let ((weight (if negative -1d0 (float (1+ (mod index 7)) 1d0)))
                 (base (gethash atom indices)))
             (assert base)
             (chem:add-anchor-restraint-term anchor ef atom (chem:get-position atom) weight)
             (dotimes (axis 3)
               (let ((ii (+ base axis)))
                 (setf (aref weights ii) weight)
                 (%nv-set! coords ii (+ (%nv-ref coords ii)
                                       (nth axis '(0.7d0 -0.4d0 0.25d0)))))))))
        (values ef coords target weights)))))

(defun %lbfgs-analytic-anchor-rms (coords target weights)
  (sqrt (/ (loop for ii below (length coords)
                 sum (expt (* 2d0 (aref weights ii)
                              (- (%nv-ref coords ii) (%nv-ref target ii))) 2))
           (length coords))))

(defun %lbfgs-analytic-anchor-energy (coords target weights)
  (loop for ii below (length coords)
        sum (* (aref weights ii)
               (expt (- (%nv-ref coords ii) (%nv-ref target ii)) 2))))

(defun %lbfgs-analytic-anchor-slope (coords target weights direction)
  (loop for ii below (length coords)
        sum (* 2d0 (aref weights ii)
               (- (%nv-ref coords ii) (%nv-ref target ii))
               (%nv-ref direction ii))))

(defun %lbfgs-test-minimizer (ef &key (enabled t) (sd-steps 0) (limit 10000)
                                    (tolerance 1d-7))
  (let ((minimizer (chem:make-minimizer ef)))
    (chem:set-use-lbfgs minimizer enabled)
    (chem:set-maximum-number-of-lbfgs-steps minimizer limit)
    (chem:set-maximum-number-of-steepest-descent-steps minimizer sd-steps)
    (chem:set-truncated-newton-tolerance minimizer tolerance)
    (chem:set-conjugate-gradient-tolerance minimizer tolerance)
    (chem:disable-print-intermediate-results minimizer)
    minimizer))

(defun %lbfgs-test-run (minimizer coords &key mask callback)
  (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
    (chem:minimize minimizer :coords coords :active-atom-mask mask :callback callback)))

(defun %lbfgs-test-constructed-minimizer (ef steps)
  ;; Ordinary record fields require a complete initializer. Supply those
  ;; from a prototype, but pass each tested L-BFGS keyword exactly once.
  (let* ((prototype (%lbfgs-test-minimizer ef :enabled nil :tolerance 1d-12))
         (arguments (loop for (key . value) in (core:encode prototype)
                          unless (member key '(:use-lbfgs :number-of-lbfgs-steps))
                            append (list key value))))
    (apply #'core:make-cxx-object 'chem:minimizer
           :use-lbfgs t :number-of-lbfgs-steps steps arguments)))

(defun run-lbfgs-regression-tests ()
  (multiple-value-bind (ef start target weights) (%lbfgs-harmonic-fixture)
    (let ((default (chem:make-minimizer ef)))
      #+tests(test-true lbfgs-enabled-by-default (chem:get-use-lbfgs default))
      #+tests(test-true lbfgs-default-step-limit
                        (= 10000 (chem:get-maximum-number-of-lbfgs-steps default))))
    (let* ((configured (%lbfgs-test-constructed-minimizer ef 7))
           (fields (core:encode configured))
           (restored (core:decode (chem:make-minimizer ef) fields)))
      #+tests(test-true lbfgs-constructor-step-limit
                        (and (chem:get-use-lbfgs configured)
                             (= 7 (chem:get-maximum-number-of-lbfgs-steps configured))))
      #+tests(test-true lbfgs-step-limit-record-roundtrip
                        (and (equal (assoc :number-of-lbfgs-steps fields)
                                    '(:number-of-lbfgs-steps . 7))
                             (chem:get-use-lbfgs restored)
                             (= 7 (chem:get-maximum-number-of-lbfgs-steps restored)))))
    (let* ((minimizer (%lbfgs-test-constructed-minimizer ef 0))
           (initial-iteration (chem:get-iteration minimizer))
           (final (%lbfgs-test-run minimizer start)))
      #+tests(test-true lbfgs-constructor-zero-skips
                        (and (every #'= final start)
                             (= initial-iteration (chem:get-iteration minimizer)))))
    (let* ((minimizer (%lbfgs-test-constructed-minimizer ef 1))
           (initial-iteration (chem:get-iteration minimizer))
           (condition (handler-case (progn (%lbfgs-test-run minimizer start) nil)
                        (chem:minimizer-exceeded-max-steps (condition) condition))))
      #+tests(test-true lbfgs-constructor-limit-is-enforced
                        (and condition (= (1+ initial-iteration) (chem:get-iteration minimizer)))))
    (let* ((minimizer (%lbfgs-test-minimizer ef))
           (stages nil)
           (previous-energy nil) (previous-slope nil)
           (previous-step nil) (previous-direction nil)
           (wolfe-ok t) (checked-steps 0) (history-direction-seen nil))
      ;; Observe the chosen step without changing any callback arguments.
      ;; Both energies and directional derivatives below come from the known
      ;; quadratic, independently of the energy function and line search.
      (chem:set-step-callback
       minimizer
       (lambda (stage coords force step direction)
         (assert (eq stage 'chem:lbfgs))
         (setf previous-energy (%lbfgs-analytic-anchor-energy coords target weights)
               previous-slope (%lbfgs-analytic-anchor-slope coords target weights direction)
               previous-step step previous-direction (copy-seq direction))
         (let ((dot (loop for ii below (length coords)
                          sum (* (%nv-ref force ii) (%nv-ref direction ii))))
               (force2 (loop for ii below (length coords) sum (expt (%nv-ref force ii) 2)))
               (direction2 (loop for ii below (length coords) sum (expt (%nv-ref direction ii) 2))))
           (when (and (> force2 1d-12) (> direction2 1d-12)
                      (< (* dot dot) (* (- 1d0 1d-8) force2 direction2)))
             (setf history-direction-seen t)))))
      (let* ((final (%lbfgs-test-run
                     minimizer start
                     :callback
                     (lambda (stage energy-function coords mask)
                       (declare (ignore energy-function mask))
                       (push stage stages)
                       (let ((energy (%lbfgs-analytic-anchor-energy coords target weights))
                             (slope (%lbfgs-analytic-anchor-slope coords target weights previous-direction)))
                         (incf checked-steps)
                         (setf wolfe-ok
                               (and wolfe-ok (> previous-step 0d0) (< previous-slope 0d0)
                                    (<= energy (+ previous-energy
                                                  (* 1d-4 previous-step previous-slope)
                                                  (* 1d-12 (max 1d0 (abs previous-energy)))))
                                    (<= (abs slope)
                                        (+ (* 0.9d0 (abs previous-slope))
                                           (* 1d-12 (max 1d0 (abs previous-slope)))))))))))
           (energy (chem:evaluate-energy ef final))
           (rms (%lbfgs-analytic-anchor-rms final target weights)))
        (format t "L-BFGS harmonic: energy=~e analytic-rms=~e steps=~d~%"
                energy rms (chem:get-iteration minimizer))
        #+tests(test-true lbfgs-opt-in-enabled (chem:get-use-lbfgs minimizer))
        #+tests(test-true lbfgs-harmonic-energy (< (abs energy) 1d-12))
        #+tests(test-true lbfgs-harmonic-analytic-force (< rms 1.01d-7))
        #+tests(test-true lbfgs-accepted-steps-satisfy-strong-wolfe
                          (and (> checked-steps 0) wolfe-ok))
        #+tests(test-true lbfgs-observation-callbacks-preserve-history history-direction-seen)
        #+tests(test-true lbfgs-bypasses-cg-and-tn
                          (and stages (every (lambda (stage) (eq stage 'chem:lbfgs)) stages)))))
    ;; Keeping the opt-in disabled must retain the legacy CG route.
    (let* ((minimizer (%lbfgs-test-minimizer ef :enabled nil))
           (stages nil))
      (chem:set-maximum-number-of-truncated-newton-steps minimizer 0)
      (let ((final (%lbfgs-test-run
                    minimizer start
                    :callback (lambda (stage energy-function coords mask)
                                (declare (ignore energy-function coords mask))
                                (push stage stages)))))
        #+tests(test-true lbfgs-disabled-legacy-converges
                          (< (%lbfgs-analytic-anchor-rms final target weights) 1.01d-7))
        #+tests(test-true lbfgs-disabled-legacy-callbacks
                          (and (member 'chem:conjugate-gradient stages)
                               (not (member 'chem:lbfgs stages))))))
    ;; At the exact minimum no curvature update or line search is necessary.
    (let* ((minimizer (%lbfgs-test-minimizer ef))
           (initial-iteration (chem:get-iteration minimizer))
           (final (%lbfgs-test-run minimizer target)))
      #+tests(test-true lbfgs-zero-force-start
                        (and (= initial-iteration (chem:get-iteration minimizer))
                             (every #'= final target)
                             (zerop (chem:get-rmsforce minimizer)))))
    (let* ((minimizer (%lbfgs-test-minimizer ef))
           (initial-iteration (chem:get-iteration minimizer))
           (mask (make-array (length start) :element-type 'bit :initial-element 0))
           (final (%lbfgs-test-run minimizer start :mask mask)))
      #+tests(test-true lbfgs-no-active-coordinates
                        (and (every #'= final start)
                             (= initial-iteration (chem:get-iteration minimizer))
                             (zerop (chem:get-rmsforce minimizer)))))
    ;; Unequal spring constants cannot all relax in one first-order step.
    (let* ((minimizer (%lbfgs-test-minimizer ef :limit 1 :tolerance 1d-12))
           (initial-iteration (chem:get-iteration minimizer))
           (condition (handler-case (progn (%lbfgs-test-run minimizer start) nil)
                        (chem:minimizer-exceeded-max-steps (condition) condition))))
      #+tests(test-true lbfgs-set-step-limit
                        (= 1 (chem:get-maximum-number-of-lbfgs-steps minimizer)))
      #+tests(test-true lbfgs-reports-step-limit
                        (and condition (= (1+ initial-iteration) (chem:get-iteration minimizer)))))
    ;; The optional SD prelude still executes before L-BFGS.
    (let ((minimizer (%lbfgs-test-minimizer ef :sd-steps 10000))
          (stages nil))
      (chem:set-steepest-descent-tolerance minimizer 1d0)
      (let ((final (%lbfgs-test-run
                    minimizer start
                    :callback (lambda (stage energy-function coords mask)
                                (declare (ignore energy-function coords mask))
                                (push stage stages)))))
        #+tests(test-true lbfgs-preserves-steepest-descent-prelude
                          (and (member 'chem:steepest-descent stages)
                               (member 'chem:lbfgs stages)
                               (not (member 'chem:conjugate-gradient stages))
                               (not (member 'chem:truncated-newton stages))
                               (< (%lbfgs-analytic-anchor-rms final target weights) 1.01d-7))))))
  ;; Changing the objective in a pre-step callback invalidates the selected
  ;; trial even when the callback leaves every coordinate and force untouched.
  (multiple-value-bind (ef start target weights) (%lbfgs-harmonic-fixture)
    (declare (ignore target weights))
    (let ((minimizer (%lbfgs-test-minimizer ef))
          (anchor (chem:make-energy-anchor-restraint ef))
          (changed nil))
      (chem:set-step-callback
       minimizer
       (lambda (stage coords force step direction)
         (declare (ignore stage coords force step direction))
         (unless changed (chem:disable anchor) (setf changed t))))
      (let ((final (%lbfgs-test-run minimizer start)))
        #+tests(test-true lbfgs-callback-objective-change-discards-selected-trial
                          (and changed (every #'= final start)
                               (zerop (chem:evaluate-energy ef final))
                               (zerop (%rms-force ef final)))))))
  ;; E = 4 (x-1)^2 + (x-3)^2 for x>0, with the second atom fixed at x=0.
  ;; The constrained minimum is x=7/5, E=16/5; the fixed atom's force is 16/5.
  (multiple-value-bind (ef butane) (%fresh-butane-ef)
    (let* ((c1 (chem:matter/first-atom-with-name butane :C1))
           (c2 (chem:matter/first-atom-with-name butane :C2))
           (stretch (chem:get-stretch-component ef))
           (anchor (chem:make-energy-anchor-restraint ef))
           (size (chem:get-nvector-size ef))
           (start (chem:make-nvector size))
           (mask (make-array size :element-type 'bit :initial-element 0))
           (c1-index nil) (c2-index nil) (found nil))
      (dolist (component (chem:all-components ef)) (chem:disable component))
      (chem:enable stretch)
      (chem:enable anchor)
      (chem:walk-stretch-terms
       stretch
       (lambda (index atom1 atom2 i1 i2 kb r0)
         (declare (ignore kb r0))
         (chem:modify-stretch-term-kb stretch index 0.0)
         (when (or (and (eq atom1 c1) (eq atom2 c2))
                   (and (eq atom1 c2) (eq atom2 c1)))
           (setf found t c1-index (if (eq atom1 c1) i1 i2)
                 c2-index (if (eq atom1 c2) i1 i2))
           (chem:modify-stretch-term-kb stretch index 4.0)
           (chem:modify-stretch-term-r0 stretch index 1.0))))
      (assert found)
      (chem:add-anchor-restraint-term anchor ef c1 (geom:vec 3d0 0d0 0d0) 1d0)
      (chem:load-coordinates-into-vector ef start)
      (geom:vec-set (geom:vec 2d0 0d0 0d0) start c1-index)
      (geom:vec-set (geom:vec 0d0 -0d0 0d0) start c2-index)
      (setf (aref mask c1-index) 1)
      (let* ((before (copy-seq start))
             (minimizer (%lbfgs-test-minimizer ef))
             (final (%lbfgs-test-run minimizer start :mask mask))
             (force (chem:make-nvector size))
             (energy (chem:evaluate-energy-force ef final :calc-force t :force force)))
        (format t "L-BFGS coupled mask: x=~e energy=~e active-force=~e fixed-force=~e~%"
                (%nv-ref final c1-index) energy (%nv-ref force c1-index) (%nv-ref force c2-index))
        #+tests(test-true lbfgs-coupled-constrained-minimum
                          (and (< (abs (- (%nv-ref final c1-index) 1.4d0)) 1d-7)
                               (< (abs (- energy 3.2d0)) 1d-10)
                               (< (abs (%nv-ref force c1-index)) 1.01d-7)))
        #+tests(test-true lbfgs-coupled-fixed-force-is-nonzero
                          (> (abs (%nv-ref force c2-index)) 3d0))
        #+tests(test-true lbfgs-preserves-every-fixed-coordinate
                          (loop for ii below size
                                always (or (= ii c1-index)
                                           (= (ext:double-float-to-bits (%nv-ref before ii))
                                              (ext:double-float-to-bits (%nv-ref final ii)))))))))
  ;; A negative harmonic potential has no strong-Wolfe step: its directional
  ;; derivative becomes more negative forever.  The failed search must leave
  ;; the last accepted coordinates recoverable, without accepting a trial.
  (multiple-value-bind (ef start target weights) (%lbfgs-harmonic-fixture :negative t)
    (declare (ignore target weights))
    (let* ((minimizer (%lbfgs-test-minimizer ef))
           (initial-iteration (chem:get-iteration minimizer))
           (condition (handler-case (progn (%lbfgs-test-run minimizer start) nil)
                        (chem:minimizer-stuck (condition) condition))))
      #+tests(test-true lbfgs-line-search-failure-is-reported condition)
      #+tests(test-true lbfgs-failed-search-retains-accepted-coordinates
                        (and condition
                             (every #'= start (chem:coordinates condition))
                             (= initial-iteration (chem:get-iteration minimizer))))))
  ;; At this translated origin, a unit displacement rounds back to the
  ;; starting double.  The initial line search must expand to a representable
  ;; step, rather than repeatedly shrinking an already unrepresentable one.
  (multiple-value-bind (ef butane) (%fresh-butane-ef)
    (let* ((atom (chem:matter/first-atom-with-name butane :C1))
           (anchor (chem:make-energy-anchor-restraint ef))
           (size (chem:get-nvector-size ef))
           (start (chem:make-nvector size))
           (mask (make-array size :element-type 'bit :initial-element 0))
           (target-x 1d16)
           (coordinate nil))
      (chem:walk-stretch-terms
       (chem:get-stretch-component ef)
       (lambda (index atom1 atom2 i1 i2 kb r0)
         (declare (ignore index kb r0))
         (cond ((eq atom atom1) (setf coordinate i1))
               ((eq atom atom2) (setf coordinate i2)))))
      (assert coordinate)
      (dolist (component (chem:all-components ef)) (chem:disable component))
      (chem:enable anchor)
      (chem:add-anchor-restraint-term anchor ef atom (geom:vec target-x 0d0 0d0) 1d0)
      (chem:load-coordinates-into-vector ef start)
      (geom:vec-set (geom:vec (+ target-x 8d0) 0d0 0d0) start coordinate)
      (setf (aref mask coordinate) 1)
      #+tests(test-true lbfgs-translated-fixture-has-unrepresentable-unit-step
                        (and (= (%nv-ref start coordinate) (- (%nv-ref start coordinate) 1d0))
                             (= (chem:evaluate-energy ef start) 64d0)))
      (let* ((minimizer (%lbfgs-test-minimizer ef))
             (final (%lbfgs-test-run minimizer start :mask mask)))
        #+tests(test-true lbfgs-expands-unrepresentable-initial-step
                          (and (= (%nv-ref final coordinate) target-x)
                               (zerop (chem:evaluate-energy ef final))
                               (zerop (%rms-force ef final)))))))
  ;; Finally exercise the complete nonlinear molecular objective, rather than
  ;; only a quadratic.  The force is reevaluated independently after minimize.
  (multiple-value-bind (ef butane) (%fresh-butane-ef)
    (declare (ignore butane))
    (let* ((start (chem:make-nvector (chem:get-nvector-size ef)))
           (minimizer (%lbfgs-test-minimizer ef :tolerance 1d-4)))
      (chem:load-coordinates-into-vector ef start)
      (dotimes (ii 9) (%nv-set! start ii (+ (%nv-ref start ii) 0.15d0)))
      (let* ((initial-energy (chem:evaluate-energy ef start))
             (final (%lbfgs-test-run minimizer start))
             (energy (chem:evaluate-energy ef final))
             (rms (%rms-force ef final)))
        (format t "L-BFGS butane: energy ~e -> ~e rms-force=~e steps=~d~%"
                initial-energy energy rms (chem:get-iteration minimizer))
        #+tests(test-true lbfgs-physical-butane-lowers-energy (< energy initial-energy))
        #+tests(test-true lbfgs-physical-butane-converges (< rms 1.01d-4)))))
  (values))

(run-lbfgs-regression-tests)

;;; END L-BFGS REGRESSION TESTS

;;; ===========================================================================
;;; Dihedral kernel agreement: dihedral_fast vs dihedral_fast_floored
;;;
;;; chem:compare-dihedral-kernels (dihedralKernelCompare.cc) compiles BOTH
;;; kernels into the binary and compares them term-by-term, so this needs only
;;; ONE build -- no recompiling to swap CANDO_DIHEDRAL_KERNEL.  It returns
;;; (max |dEnergy|, max |dForce|) over every dihedral term.  Use a fresh energy
;;; function + freshly loaded coordinates so the geometry is the (non-degenerate)
;;; input structure and is independent of any earlier in-file minimization.
;;;
;;; A sign error in the floored kernel's Blondel-Karplus dphi/dq would flip the
;;; force and blow max|dForce| up to ~force magnitude, so the force gate doubles
;;; as the sign check.
;;; ===========================================================================
(let* ((cmp-ef   (chem:make-energy-function :matter agg))
       (cmp-pos  (chem:make-nvector (chem:get-nvector-size cmp-ef))))
  (chem:load-coordinates-into-vector cmp-ef cmp-pos)
  (let ((cmp-dih (chem:get-dihedral-component cmp-ef)))
    (multiple-value-bind (max-de max-df)
        (chem:compare-dihedral-kernels cmp-dih cmp-pos)
      (format t "dihedral kernel agreement: max|dEnergy| = ~e   max|dForce| = ~e~%"
              max-de max-df)
      #+tests(test-true dihedral-kernels-energy-agree (< max-de 1d-6))
      #+tests(test-true dihedral-kernels-force-agree  (< max-df 1d-6)))))


;;; ===========================================================================
;;; Near-linear (degenerate) dihedral: dihedral_fast blows up, floored stays finite
;;;
;;; chem:dihedral-kernel-near-linear builds one dihedral whose atoms 1-2-3 are
;;; exactly collinear (eps = 0), so the AD kernel's 1/(t1^2+t2^2) factor is 1/0
;;; (non-finite) while the floored Blondel-Karplus kernel produces a finite,
;;; bounded force.  It returns (max|force_fast|, max|force_floored|).  Wrap the
;;; call so the fast kernel's divide-by-zero doesn't trap.
;;; ===========================================================================
(multiple-value-bind (nl-fast nl-floored)
    (ext:with-float-traps-masked (:underflow :overflow :invalid :inexact :divide-by-zero)
      (chem:dihedral-kernel-near-linear 0d0))
  (format t "near-linear dihedral: max|force| fast = ~e   floored = ~e~%"
          nl-fast nl-floored)
  ;; The floored kernel must stay finite (the whole point of the floor) ...
  #+tests(test-true near-linear-floored-finite
                    (< nl-floored most-positive-double-float))
  ;; ... while the AD kernel diverges (inf/nan) at the exact degeneracy.
  #+tests(test-true near-linear-fast-diverges
                    (not (< nl-fast most-positive-double-float))))
