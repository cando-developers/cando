(in-package #:clasp-tests)

(defparameter agg (chem:load-mol2 "sys:extensions;cando;src;lisp;regression-tests;data;hexapeptide.mol2"))
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
  (time (dotimes (i 1) (defparameter energy (chem:evaluate-energy-force ef pos t force)))))
(format t "About to calculate energy~%")
(defparameter energy (chem:evaluate-energy-force ef pos t force))
(defparameter expected-energy 20937032.6233387)
(format t "         energy = ~f~%" energy)
(format t "expected-energy = ~f~%" expected-energy)

(defparameter delta (abs (- energy expected-energy)))

(format t "delta = ~f~%" delta)
;; These two tests are not yet reliable.
(test-true energy-delta (< delta 3.0))

(let ((components (chem:all-components ef))
      (ht (make-hash-table)))
  (format t "Components ~s~%" components)
  (loop for comp in components
        for energy = (chem:energy-component-evaluate-energy ef comp pos)
        for class-name = (class-name (class-of comp))
        do (setf (gethash class-name ht) energy)
        do (format t "~s -> ~f~%" (class-name (class-of comp)) energy))
  (test-true energy-stretch  (< (abs (print (- 14.873171 (gethash 'chem:energy-stretch ht)))) 0.01))
  (test-true energy-angle    (< (abs (print (- 86.5912   (gethash 'chem:energy-angle ht)))) 0.01))
  (test-true energy-dihedral (< (abs (print (- 49.70431  (gethash 'chem:energy-dihedral ht)))) 0.01))
  (test-true energy-nonbond  (< (abs (print (- 30.878897 (gethash 'chem:energy-nonbond ht)))) 0.01))
  )

(defparameter expected-force (coerce
                              #(2.6174526 1.4443848 2.0994837 -1.8517512 -0.75102925 -2.432416
                                -8.5939865 15.114094 12.309942 -0.047213633 -9.139131 -5.459774
                                4.710642 -7.9739184 -21.801796 6.423041 -7.5212755 0.8190722 7.756392
                                6.2105813 10.829898 -24.751427 -4.121009 -29.946684 26.274252
                                7.7273483 17.21085 0.055465523 -1.477793 0.21809128 -3.4077888
                                8.153169 2.0224102 -30.561218 8.284508 -23.848635 -0.58561605
                                -8.473848 2.8910348 32.6781 -3.019022 15.858077 1.97434 3.0862103
                                1.6496463 -16.856556 2.8574927 10.113729 0.67898744 -4.673321
                                -2.160513 28.180765 -5.485939 -6.028311 1.9882052 -3.1959887
                                -8.5077095 -28.67401 16.732014 -20.196583 -10.418779 1.1545113
                                4.122681 22.763374 -21.61034 7.3822007 -0.73201305 -1.4176271
                                0.9215267 15.21612 -28.872467 13.385801 -3.8336127 3.2060728
                                5.0825415 -9.04366 19.65255 -12.886541 6.726975 -17.260607 -3.466814
                                -13.445648 20.107267 10.146154 2.4285238 13.6339655 -0.77199274
                                1.9192798 -6.590973 -17.69256 -6.8651495 17.641161 2.9584932
                                21.212105 -19.502378 -14.507019 -4.437197 -1.9174879 4.146387
                                -13.066817 -12.779705 -12.458004 -3.5705705 8.2768955 13.5090475
                                -7.840398 7.105832 3.1228426 2.766078 1.684678 0.9286907 1.4099352
                                5.206618 18.887606 7.509185 -3.8592947 -0.69349605 -20.522873
                                -9.241659 -25.045328 -1.9408325 3.6344519 1.1450552 28.018919
                                13.662445 15.258868 -0.82603025 -1.017411 -0.3827567 -2.533032
                                -9.897736 6.642954 10.044313 17.745338 14.678083 -4.7409983 6.9025574
                                -1.5460207 -13.431607 -23.181091 -6.9765716 -1.1340516 -0.025261158
                                -11.490736 -0.98068196 -2.8314056 -1.7598569 -0.40958136 -0.8879694
                                -0.34210607 0.31659028 -0.09658031 0.76061964 1.3496021 2.280601
                                2.4254735 -0.048426416 -0.30820397 0.33929715 -1.9025369 2.6829236
                                8.917992 -7.0347095 -6.177149 1.1471264 0.928406 3.395796 -7.6135273
                                -3.312034 -1.1857884 8.16532 1.3643141 -0.7515218 2.8079681 2.972939
                                2.209898 -1.3026403 0.27428436 1.4050107 -3.1786232 -2.23991
                                -1.7728083 8.6719675 -0.007715225 1.6070875 2.4650762 -14.763551
                                4.74731 4.2538605 -0.37204972 1.2359741 -1.6259735 6.0726347
                                -0.19272625 -1.4523464 5.1062555 5.156993 7.786087 1.1951859
                                -0.371384 -3.2192175 0.30030414 -2.7865481 -6.1512914 1.3151727
                                -2.0746167 3.7924454 -1.8045332 0.016275978 2.4537725 -1.7647833
                                0.33926266 0.4947312 -1.0805506 1.7444935 -2.8531337 -0.08640733
                                1.1643641 4.4030867 14.096523 0.91154593 7.077344 -1.825268
                                0.15667284 -1.7416267 -1.1449207 0.5241355 1.644149 -0.22344476
                                0.095424354 2.4192533 -2.1054282 0.47404295 -13.791035 0.55793864
                                -0.6467289 -3.0642304 -4.2187953 -3.525773 -3.168699 -0.8448597
                                -5.879744 0.4610304 4.91517 -1.2414533 0.85410076 -1.6704061
                                -3.841198 8.771013 -1.0210749 4.4894695 -2.324779 1.031044 -1.2010121
                                -1.4488937 0.9349262 -3.9381697 1.9427706 1.9104602 0.022106707
                                -0.8322981 4.580284 8.833524 3.774918)
                              (list 'vector (geom:vecreal-type))))
(format t "force = ~s~%" force)

(defparameter expected-force-mag (chem:nvector-magnitude expected-force))
(defparameter force-mag (chem:nvector-magnitude force))

(format t "         force-mag = ~f~%" force-mag)
(format t "expected-force-mag = ~f~%" expected-force-mag)
(test-true force-mag (< (abs (- expected-force-mag force-mag)) 1.0))

(defparameter force-acos (/ (chem:nvector-dot expected-force force) expected-force-mag force-mag))
(format t "force-acos = ~f~%" force-acos)

(test-true force-acos (< force-acos 0.01))


(defparameter min (chem:make-minimizer ef))
#+(or)(time (dotimes (i 10)
        (chem:save-coordinates-from-vector ef pos)
        (chem:minimize min)))


(defparameter spiro (cando:load-mol2 "data/struct-0000.mol2"))

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
    (format t "~s energy total = ~f~%" name (chem:evaluate-energy en pos mask))
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
(progn
  (format t "Running simd1~%")
  (core:set-simd-width 1)
  (energy-components "simd1" ef2 pos2 mask)
  (format t "Running simd8~%")
  (core:set-simd-width 8)
  (energy-components "simd8" ef2 pos2 mask)
  )

(progn
  (core:set-simd-width 1)
  (format t "no mask simd1 timing~%")
  (time (energy-multiple-dihedral 100000 ef2 pos2 mask)))
(progn
  (core:set-simd-width 8)
  (format t "no mask simd8 timing~%")
  (time (energy-multiple-dihedral 100000 ef2 pos2 mask)))

