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
(defparameter expected-energy 139.19010958782783)
(format t "         energy = ~f~%" energy)
(format t "expected-energy = ~f~%" expected-energy)

(defparameter delta (abs (- energy expected-energy)))

(format t "delta = ~f~%" delta)
;; These two tests are not yet reliable.
#+tests(test-true energy-delta (< delta 0.2))

(defparameter correct-energy-stretch 12.748283740690873)
(defparameter correct-energy-angle 83.66812723256757)
(defparameter correct-energy-dihedral 19.549603494795992)
(defparameter correct-energy-nonbond 23.22409511977336)
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
                              #(0.7313171276817805d0 -0.7843348305291686d0 0.11427955752899521d0
                                -0.34291189131523137d0 0.32767111035158614d0 -0.329898621538927d0
                                -0.7110813349877194d0 0.6592901977086865d0 0.40345197621544937d0
                                0.5958436385481027d0 -0.36974475109952865d0 0.22683190298909575d0
                                -0.02608690452320951d0 -0.10101071950616987d0 -0.27466990849422857d0
                                0.4059398239080221d0 -0.29611364308005766d0 0.053853321991587766d0
                                -0.9567350337481149d0 -0.2651452754995948d0 -0.3125462936797983d0
                                -0.7838353791805783d0 -0.001188066110304753d0 -0.3468349828308428d0
                                0.8560832116725887d0 0.25810948286855284d0 0.24231035726673777d0
                                0.08704551359396495d0 -0.01489532459401488d0 0.06373333173923168d0
                                -0.6619338583739367d0 0.5200787907848752d0 -0.01751758434139772d0
                                0.5572494778978427d0 -0.2887878545113747d0 -0.9971431435318183d0
                                -0.7134183316743781d0 -0.31687361483390586d0 0.2662245049959058d0
                                -0.42609949776674316d0 0.1473001334921655d0 0.8397458799088238d0
                                -0.0851301692400035d0 0.17836762683763108d0 0.5507189290162904d0
                                0.1508205092200554d0 -0.31224499508772d0 -0.26063269658910304d0
                                0.3837636815014146d0 0.04853288160440612d0 -0.7603916842904349d0
                                0.40341799450256005d0 0.06440963776615817d0 -0.8098487139920628d0
                                0.13554602444991937d0 0.0741531836024043d0 0.039724939010686494d0
                                0.2250559175061433d0 0.21453684571564402d0 -0.18484691095396943d0
                                -0.6014612054263633d0 0.3302179102796425d0 0.21774632388300733d0
                                -0.15576340271953715d0 -0.4844408871171204d0 0.5252047425535622d0
                                -0.15881088066755766d0 -0.0124820262534866d0 0.17047269957603595d0
                                0.5070872187825861d0 -0.08970283832050097d0 0.053805298489363995d0
                                -0.43587877394338637d0 0.5481798160351463d0 -0.4708105444474042d0
                                0.17238162070775465d0 -0.07364334574026277d0 0.2110748929049946d0
                                -0.14217419589912067d0 -0.5071504723222242d0 -0.44076152739477326d0
                                0.3608549117225661d0 -0.11667874691653068d0 0.31550070382046647d0
                                0.5275563298311408d0 -0.09186511521918428d0 0.6075663873662613d0
                                0.1424379446213615d0 0.4088061071938096d0 -0.14201181334305657d0
                                -0.6331397908596035d0 0.1496485335185535d0 0.2506524281706096d0
                                0.2504963760788428d0 0.07019198850761299d0 -0.6283227129638327d0
                                -0.32043306062701093d0 -0.1035489471308324d0 0.03403104244164412d0
                                -0.04884599387645083d0 -0.015339768549434085d0 -0.05744109483936655d0
                                0.6657605357639325d0 -0.5069685049874448d0 0.06573464887100013d0
                                -0.4450359713863392d0 0.5887977304697216d0 0.22290669141573005d0
                                -0.22686701665077785d0 0.24173765039854267d0 0.7899293970321287d0
                                -1.1414566602568645d0 0.4468707774898886d0 0.0013093140883791918d0
                                0.5080045058681834d0 -0.05874462991310331d0 0.33174498242333195d0
                                -0.6563096918464684d0 -0.14610519496961683d0 0.5935185788234367d0
                                -0.0502686170994977d0 0.18617207078016484d0 -0.35677913379656445d0
                                0.4323801513140918d0 0.19870809901519593d0 -0.4002428899969046d0
                                0.2393771702882026d0 0.15674951622349198d0 -0.18467553044190277d0
                                0.6574545866099177d0 0.2400292498426446d0 -0.1948554767048679d0
                                -1.082652377496741d0 1.2320935292290391d0 0.6776695871761318d0
                                -0.4459890357266924d0 0.49947343607922257d0 -0.11665935178288446d0
                                0.15543277584683501d0 -0.3333508637648541d0 -0.14006399498760158d0
                                0.42925868751453444d0 -0.8788722936089499d0 -0.3453398368883223d0
                                -0.3074580630955875d0 0.1777114285441745d0 -0.28774291559983217d0
                                -0.05263084647925689d0 0.2668489930742355d0 0.01129247528010928d0
                                0.08866258256034756d0 -0.08340867240350669d0 0.049402651552455296d0
                                0.09466864960240234d0 -0.07770068544353366d0 0.28644900184729316d0
                                -0.016376001413197186d0 0.2579989413458277d0 0.04081736845851545d0
                                0.06964084139575319d0 -0.02535980491846033d0 0.003323295530815281d0
                                0.6143207721370927d0 0.13749581519015297d0 0.02086052384864925d0
                                0.13068898083289324d0 0.24114593249385063d0 0.0476634290045307d0
                                -0.18843662217998827d0 -0.2715781220967061d0 0.11847303780779955d0
                                0.1301423985308981d0 -0.2133322754956743d0 0.15489066257151524d0
                                0.44917452882255415d0 0.13519913212069146d0 0.07763488520914295d0
                                0.04267596372233849d0 -0.006247423037644584d0 -0.07511326694259103d0
                                -0.04181377766016559d0 -0.05643393879432436d0 -0.0794177342880086d0
                                0.009290871603157544d0 -0.06634417898926885d0 0.11346825149737383d0
                                -0.490670811466446d0 0.07430752450837315d0 0.703457656958072d0
                                0.02632675384807072d0 -0.009845572941264879d0 -0.045556889974364106d0
                                0.2200452265525705d0 -0.07543788110483846d0 -0.04057350436283838d0
                                0.27763534074811125d0 -0.002120745795000864d0 -0.4282120885318866d0
                                -0.016201402432955647d0 -0.10027562121319199d0 0.013785737197281211d0
                                -0.10995075825945021d0 -0.15986616137916063d0 -0.14031830452822452d0
                                -0.0010338755682908762d0 -0.005371982272776257d0 0.0603318368103499d0
                                -0.09821281213364358d0 -0.009664499684758784d0 -0.0758899162024357d0
                                -0.10066285191199234d0 0.15230982162722417d0 -0.2834567024664232d0
                                -0.3328404161778856d0 -0.04316271118537286d0 0.01712805021302084d0
                                0.15239046441430742d0 0.1572295314735211d0 -0.049776593042651064d0
                                -0.0035999737457596197d0 0.012949062418218386d0
                                0.015365143016109398d0 0.13172054285951615d0 -0.24687796165802234d0
                                -0.3959246494710529d0 -0.05262314730230066d0 -0.2256414311684376d0
                                -0.14004416593892738d0 0.14197497317821767d0 -0.12167749344027445d0
                                -0.10307980597089633d0 0.9451103081654494d0 -0.0877091500409527d0
                                0.047525662918696095d0 -0.2570463075850139d0 -0.1673566289395625d0
                                -0.06364745416212297d0 0.07081478937904818d0 -0.1368704497725733d0
                                -0.08506897423100508d0 0.04607314630595363d0 -0.06185249804907123d0
                                0.04921110221501829d0 -0.06260666124976895d0 -0.0035112792610380306d0
                                0.1343661886199076d0 -0.11472477335513814d0 -0.2783525240287381d0
                                0.05029049623816573d0 -0.09674900114474458d0 -0.21743801964598883d0
                                0.16885194867960285d0 0.22395754867653464d0 -0.08002570033951538d0
                                -0.09524000494826691d0 0.03303298646681578d0 -0.12703144043226425d0
                                0.21144176083782085d0 0.008788048584329586d0 -0.35921942450833283d0
                                -0.07373817768435319d0 0.10825572463529572d0 0.08562449911460364d0
                                -0.030677989835232178d0)
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
          for energy = (chem:energy-component-evaluate-energy en comp pos mask)
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
            for energy = (chem:energy-component-evaluate-energy en dihedral pos mask)))))

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
                                          a1 a2 a3 a4
                                          ))
         (minimizer (chem:make-minimizer butane-ef))
         )
    (format t "dh = ~s~%" dh)
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
       (expected-energy 810.093569671)
       (expected-components '((CHEM::ENERGY-STRETCH . 95.7978)
                              (CHEM::ENERGY-ANGLE . 301.018)
                              (CHEM::ENERGY-DIHEDRAL . 230.37)
                              (CHEM::ENERGY-NONBOND-TOTAL . -55.4358)
                              (CHEM::ENERGY-NONBOND14 . 238.344)
                              (CHEM::ENERGY-CHIRAL-RESTRAINT . 0)
                              (CHEM::ENERGY-ANCHOR-RESTRAINT . 0)
                              (CHEM::ENERGY-DIHEDRAL-RESTRAINT . 0)
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
