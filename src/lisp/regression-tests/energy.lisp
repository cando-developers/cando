#+tests(in-package #:clasp-tests)

(defparameter agg (cando:load-mol2 "sys:extensions;cando;src;lisp;regression-tests;data;hexapeptide.mol2"))
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
(defparameter expected-energy 178.83531388394212)
(format t "         energy = ~f~%" energy)
(format t "expected-energy = ~f~%" expected-energy)

(defparameter delta (abs (- energy expected-energy)))

(format t "delta = ~f~%" delta)
;; These two tests are not yet reliable.
#+tests(test-true energy-delta (< delta 3.0))

(let ((components (chem:all-components ef))
      (ht (make-hash-table)))
  (format t "Components ~s~%" components)
  (loop for comp in components
        for energy = (chem:energy-component-evaluate-energy ef comp pos)
        for class-name = (class-name (class-of comp))
        do (setf (gethash class-name ht) energy)
        do (format t "~s -> ~f~%" (class-name (class-of comp)) energy))
  #+tests(test-true energy-stretch  (< (abs (print (- 14.873171 (gethash 'chem:energy-stretch ht)))) 0.01))
  #+tests(test-true energy-angle    (< (abs (print (- 86.5912   (gethash 'chem:energy-angle ht)))) 0.01))
  #+tests(test-true energy-dihedral (< (abs (print (- 49.70431  (gethash 'chem:energy-dihedral ht)))) 0.01))
;;  (test-true energy-nonbond  (< (abs (print (- 30.878897 (gethash 'chem:energy-nonbond ht)))) 0.01)) ;; old
  #+tests(test-true energy-nonbond  (< (abs (print (- 27.67031857189474 (gethash 'chem:energy-nonbond ht)))) 0.01))
  )

(defparameter expected-force (coerce
                              #(2.648438613813617d0 1.4691883842230746d0 2.0984212884230526d0
                                -1.85037235384908d0 -0.740820573680203d0 -2.4285662607173113d0
                                -8.697339926007915d0 15.030147773191908d0 12.311608823205955d0
                                -0.044330493133313914d0 -9.151267913775813d0 -5.4576267149918944d0
                                4.744442001870059d0 -7.946210794362269d0 -21.801174477134083d0
                                6.461062082875487d0 -7.489652040747462d0 0.8173248283669575d0
                                7.749823842111675d0 6.212960464519702d0 10.829454493380583d0
                                -24.76344010290134d0 -4.103864095882832d0 -29.94523115988545d0
                                26.29342148574872d0 7.691758798183855d0 17.212164925484707d0
                                0.05826614843678057d0 -1.4812934303497438d0 0.2170538128987216d0
                                -3.4110262864426555d0 8.163442747473553d0 2.0248051780343834d0
                                -30.5591045166907d0 8.296602670532263d0 -23.853497849184222d0
                                -0.5816743315393653d0 -8.452191242349228d0 2.88158568327548d0
                                32.67331171186668d0 -3.036707899042751d0 15.864584364775231d0
                                1.9730493898646968d0 3.0842775379157197d0 1.652399192575691d0
                                -16.870093748990357d0 2.7915120743828594d0 10.147050785051876d0
                                0.6829101245827727d0 -4.6665608947817665d0 -2.1640893750909287d0
                                28.191593479532738d0 -5.464074203126487d0 -6.045341499842981d0
                                1.9954354791386755d0 -3.1717263459928375d0 -8.52119143643655d0
                                -28.679785177798237d0 16.71963876568821d0 -20.1898265144166d0
                                -10.442888920059799d0 1.142213447441778d0 4.1270874045157155d0
                                22.780929716763232d0 -21.60239312661955d0 7.359959117140701d0
                                -0.7319202534072868d0 -1.416450969323416d0 0.9231121964458304d0
                                13.182647543923759d0 -23.896689473320457d0 14.452255827932953d0
                                -3.837712650943934d0 3.2180939418059014d0 5.085749959169358d0
                                -9.352513706191706d0 21.512060696424594d0 -7.5252866411319435d0
                                4.647368686199324d0 -13.99340568130923d0 -7.793553886910733d0
                                -13.451807826686249d0 20.09245294914637d0 10.14805361754992d0
                                2.4315097643592596d0 13.626072613887576d0 -0.7737048005380218d0
                                3.9329046230177522d0 -9.883530072366437d0 -13.356873618936177d0
                                -6.616723308016074d0 15.751632686117823d0 -2.396972944448554d0
                                21.21653997941251d0 -19.488523383280217d0 -14.511961191749329d0
                                -4.434417441858745d0 -1.9139401813107708d0 4.149220257995272d0
                                -13.072918769550654d0 -12.77874017601126d0 -12.456758167880833d0
                                -1.4240894015533114d0 3.3576599877562803d0 12.424119321676685d0
                                -7.850525622480408d0 7.104251752142027d0 3.1326413708777507d0
                                2.7704383945082003d0 1.6891125035204355d0 0.929025699137854d0
                                1.3965553734994807d0 5.195141436552376d0 18.88930637610872d0
                                7.507524228121186d0 -3.8624122526458167d0 -0.6920605647145714d0
                                -20.508404376221023d0 -9.262994388600992d0 -25.05836145115043d0
                                -1.9387121833940943d0 3.644474305633089d0 1.1454409515528532d0
                                28.018881238250167d0 13.671927683723277d0 15.26433185432932d0
                                -0.8316599831450806d0 -1.0114375590716511d0 -0.37949306235601066d0
                                -2.523796068385369d0 -9.898423950434626d0 6.644580243522978d0
                                10.015068010559359d0 17.74987737576017d0 14.679595646900676d0
                                -4.745750093665304d0 6.900214858931701d0 -1.5461712238717467d0
                                -13.425987099214073d0 -23.18688778764417d0 -6.9801774669459045d0
                                -1.1216180099317952d0 -0.02261301366683364d0 -11.492459312010071d0
                                -0.9824265468439362d0 -2.830781547293762d0 -1.7624924672250872d0
                                -0.4094269692770521d0 -0.889934474613457d0 -0.33912638751108726d0
                                0.31601318063416467d0 -0.09924772911222196d0 0.7610025704838531d0
                                1.347817725193501d0 2.277646626530442d0 2.422548829148709d0
                                -0.048667669531854786d0 -0.3038298235546239d0 0.3348036930865084d0
                                -1.9005435696034418d0 2.6828141808733696d0 8.917936404039411d0
                                -7.029174316049497d0 -6.1749066892309425d0 1.1465255479598953d0
                                0.9289848189016017d0 3.393158298009854d0 -7.6117749606449046d0
                                -3.319390088764841d0 -1.1739502713167413d0 8.16170323021066d0
                                1.3626957779038242d0 -0.752690706165138d0 2.808381672528942d0
                                2.9710888278680323d0 2.207191413666899d0 -1.3034230098905177d0
                                0.2732978327314467d0 1.404045440727069d0 -3.1770808828898733d0
                                -2.2372765799613465d0 -1.7625240533659727d0 8.669768764818059d0
                                -0.007625886055557182d0 1.6071477793391586d0 2.4649137266261056d0
                                -14.76921221044264d0 4.7482037484929975d0 4.259424334939795d0
                                -0.3708311320639189d0 1.2376054099208762d0 -1.6267604677936995d0
                                6.072348024696663d0 -0.19294937650646368d0 -1.450081129624229d0
                                5.09851704762864d0 5.152263651816458d0 7.795564287251896d0
                                1.1960210256537243d0 -0.37342076992137657d0 -3.2192120410263434d0
                                0.32403129606853465d0 -2.775445056009548d0 -6.153027691464538d0
                                1.339652995450124d0 -2.0661933248750963d0 3.7913008921507756d0
                                -1.8048286855482067d0 0.01710186002389057d0 2.4541068211537915d0
                                -1.7679671626027935d0 0.34104405802371734d0 0.49408949343826186d0
                                -1.080780327857666d0 1.7444446539939904d0 -2.8531300984868917d0
                                -0.08755170350237362d0 1.164084203626301d0 4.403389422001439d0
                                14.101579862492331d0 0.9120318096434773d0 7.077583378125778d0
                                -1.8252068272398179d0 0.1575272464003339d0 -1.7416886580725888d0
                                -1.1469369803622786d0 0.521496806007318d0 1.6456753223089557d0
                                -0.22352746025726944d0 0.094968173324905d0 2.416755116983889d0
                                -2.10218205868649d0 0.47590141533826724d0 -13.791568751341906d0
                                0.5588939470196018d0 -0.6470352523161138d0 -3.066545007392963d0
                                -4.217178589529473d0 -3.5246208554916967d0 -3.1692591035400923d0
                                -0.8487903832580903d0 -5.886033621298784d0 0.4588603560796763d0
                                4.915651077563303d0 -1.2428493752557082d0 0.854927897728811d0
                                -1.670361565106614d0 -3.8427961251853775d0 8.770624599260925d0
                                -1.02292716046444d0 4.488355450112775d0 -2.3240795270426955d0
                                1.032428095369244d0 -1.1993510355111068d0 -1.448654068795384d0
                                0.9342715123458717d0 -3.937896274120281d0 1.9414453837550123d0
                                1.9131484961158332d0 0.022772234687114895d0 -0.831298210645949d0
                                4.582863062974887d0 8.834747895327332d0 3.7773171193227157d0)
                              (list 'vector (geom:vecreal-type))))
(format t "force = ~s~%" force)

(defparameter expected-force-mag (chem:nvector-magnitude expected-force))
(defparameter force-mag (chem:nvector-magnitude force))

(format t "         force-mag = ~f~%" force-mag)
(format t "expected-force-mag = ~f~%" expected-force-mag)
#+tests(test-true force-mag (< (abs (- expected-force-mag force-mag)) 1.0))

(defparameter force-acos (/ (chem:nvector-dot expected-force force) expected-force-mag force-mag))
(format t "force-acos = ~f~%" force-acos)

#+tests(test-true force-acos (< force-acos 0.01))


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


(let* ((butane (cando:load-mol2 "data/butane.mol2")))
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


(ql:quickload :open-force-field/load)


(let* ((agg (chem:sample-aggregate))
       (ef (chem:make-energy-function :matter agg :use-excluded-atoms nil))
       (coords (chem:matter/extract-coordinates agg))
       (components (chem:make-energy-components))
       (energy (chem:scoring-function/evaluate-all ef coords :energy-components components))
       (expected-energy 834.5838566176413)
       (expected-components '((CHEM:ENERGY-ANCHOR-RESTRAINT . 0.0d0)
                              (CHEM:ENERGY-CHIRAL-RESTRAINT . 0.0d0)
                              (CHEM:ENERGY-DIHEDRAL-RESTRAINT . 0.0d0)
                              (CHEM:ENERGY-NONBOND-TOTAL . -55.43579102222518d0)
                              (CHEM:ENERGY-NONBOND14 . 238.34374324574765d0)
                              (CHEM:ENERGY-DIHEDRAL . 219.3372469384417d0)
                              (CHEM:ENERGY-ANGLE . 301.01757487803803d0)
                              (CHEM:ENERGY-STRETCH . 131.32108257763906d0)))
       (nonbond-component (chem:get-nonbond-component ef))
       (dihedral-component (chem:get-dihedral-component ef))
       (angle-component (chem:get-angle-component ef))
       (stretch-component (chem:get-stretch-component ef))
       )
  (format t "nonbond description ~s~%" (chem:energy-component/description nonbond-component))
  (format t "dihedral description ~s~%" (chem:energy-component/description dihedral-component))
  (format t "Total energy = ~f  expected energy: ~f~%" energy expected-energy)
  #+tests(test-true energy (< (abs (- energy expected-energy)) 0.00001))
  (loop for (comp . en) in (chem:energy-components/components components)
        for expected-en = (cdr (assoc comp expected-components))
        if expected-en
          do (progn
               (format t "~s (~a)-> energy: ~f  expected: ~f ~%" comp
                       (if (< (abs (- expected-en en)) 0.001)
                           "matched"
                           "MISMATCHED")
                       en expected-en)
               #+tests(test-true energy-component (< (abs (- expected-en en)) 0.001)))
        else
          do (progn
               (warn "There is a new energy component ~s that isn't in the regression test" comp)
               #+tests(test-true new-energy-component t)))
  (format t "About to run gradient and hessian finite-difference tests~%")
  (chem:run-test-calls nonbond-component *standard-output* coords)
  (chem:run-test-calls dihedral-component *standard-output* coords)
  (chem:run-test-calls angle-component *standard-output* coords)
  (chem:run-test-calls stretch-component *standard-output* coords)
  #+tests(nonbond-fd-tests (= 0 (chem:run-test-calls nonbond-component *standard-output* coords)))
  #+tests(dihedral-fd-tests (= 0 (chem:run-test-calls dihedral-component *standard-output* coords)))
  #+tests(angle-fd-tests (= 0 (chem:run-test-calls angle-component *standard-output* coords)))
  #+tests(stretch-fd-tests (= 0 (chem:run-test-calls stretch-component *standard-output* coords)))

  )
