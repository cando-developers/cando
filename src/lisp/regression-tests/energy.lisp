(in-package #:clasp-tests)

(defparameter agg (chem:load-mol2 "sys:extensions;cando;src;lisp;regression-tests;data;hexapeptide.mol2"))
(chem:setf-force-field-name (cando:mol agg 0)  :smirnoff)

(leap:load-smirnoff-params (probe-file "sys:extensions;cando;src;lisp;regression-tests;data;force-field.offxml"))


(defparameter ef (chem:make-energy-function :matter agg))


(defparameter pos (chem:make-nvector (chem:get-nvector-size ef)))
(defparameter force (chem:make-nvector (chem:get-nvector-size ef)))

(chem:load-coordinates-into-vector ef pos)
(format t "pos = ~s~%" pos)
(progn
;;  (gctools:wait-for-user-signal "signal")
  (time (dotimes (i 10000) (defparameter energy (chem:evaluate-energy-force ef pos t force)))))
(format t "Energy = ~f~%" energy)

(defparameter expected-force (coerce
  #(-86.73449263636415d0 -73.55316904531266d0 -5.346204607031707d0
    -2.0971564256407627d0 -0.7133730245137889d0 -1.7900555664882147d0
    222116.62005163686d0 310.0771401187096d0 1389.8106853626543d0
    0.010446539771560306d0 -9.179454068068946d0 -5.492142900666708d0
    -45.54054304938287d0 0.24740912282260116d0 18.525367831923102d0
    -221973.9537660818d0 -233.14021239748405d0 -1409.7198162783047d0
    7.768839981611046d0 6.197115563085713d0 10.808476201871994d0
    -24.748262742466512d0 -4.107537948149908d0 -29.952863348830995d0
    26.310019977771134d0 7.689710443045155d0 17.20794817105258d0
    0.06357540006066577d0 -1.4823590682777776d0 0.2134485442155784d0
    -3.4064379451535203d0 8.163551246947106d0 2.0237846876883485d0
    -30.557938661805203d0 8.296541292207154d0 -23.853834334424565d0
    -0.5769948458732149d0 -8.4510353183983d0 2.8801330619393077d0
    32.67405157821935d0 -3.036795759377446d0 15.864450880037262d0
    1.973810581235812d0 3.0842537370233662d0 1.6520665087074962d0
    -16.85969782357271d0 2.7958101967161433d0 10.145181946362218d0
    0.683155792605414d0 -4.666583169058494d0 -2.1641327160224564d0
    28.197609807116354d0 -5.461316151141949d0 -6.045179908592856d0
    2.028427353567891d0 -3.1537965728915167d0 -8.530850849024969d0
    -28.67965929357215d0 16.719612739027646d0 -20.189852933047497d0
    -10.442730397434051d0 1.142215822998993d0 4.127046848695567d0
    22.7810150934944d0 -21.602413278466077d0 7.3599482585371d0
    -0.7317872962717615d0 -1.4164873474920825d0 0.9230707537512389d0
    15.30097955929534d0 -28.833826523055702d0 13.372074007447422d0
    -3.83766620873595d0 3.218078623818504d0 5.0857430978598d0
    -9.11071029256718d0 19.619412821175793d0 -12.879653536299807d0
    6.6617606745130225d0 -17.288531797250535d0 -3.459275697889793d0
    -13.45178651495714d0 20.092446382251513d0 10.148051354082774d0
    2.4315688891648826d0 13.626046139088702d0 -0.7737113680869595d0
    1.9189375493829126d0 -6.588279530196939d0 -17.691246186824305d0
    -6.858096799147599d0 17.644406402676506d0 2.9573079491521916d0
    21.21654357520866d0 -19.48851094893725d0 -14.511954628464188d0
    -4.434397670976444d0 -1.9139447048013534d0 4.149219094065561d0
    -13.07288318444486d0 -12.778759053472934d0 -12.456763055243242d0
    -3.541990936783008d0 8.29493261148671d0 13.50420523502712d0
    -7.850515433800668d0 7.104248266303729d0 3.1326391875002795d0
    2.7704484958548132d0 1.689108118743631d0 0.9290250070221205d0
    1.3968071005937657d0 5.195295486062979d0 18.889258741289865d0
    7.507536204404101d0 -3.8624227658528576d0 -0.6920699452454815d0
    -20.50839794887366d0 -9.262996409808741d0 -25.05836281411205d0
    -1.938700107393006d0 3.644461441018225d0 1.1454396405571714d0
    28.018886618953907d0 13.67192651539256d0 15.264330221271132d0
    -0.8316547580003428d0 -1.0114390806905813d0 -0.3794949921805104d0
    -2.5237916465515364d0 -9.898424560879166d0 6.644580090949788d0
    10.01507146184396d0 17.749876776529252d0 14.67959533352954d0
    -4.745747609287943d0 6.900214280284919d0 -1.5461715289601474d0
    -13.425982330386372d0 -23.186888299089983d0 -6.980176520520492d0
    -1.1216152797120134d0 -0.022613383781221224d0 -11.492459894230459d0
    -0.4691732784416666d0 -2.173320089989665d0 -2.2979874315927207d0
    -1.771654917694326d0 -0.25828182299919d0 -0.3112968035008171d0
    0.34649848051214327d0 -0.09924746675117747d0 0.7274632170828793d0
    1.2536147282248662d0 2.2846444035700193d0 2.455534823736825d0
    -0.03968946744470276d0 -0.31007391826059805d0 0.33136313905549236d0
    -6.163138575955535d0 -1.1770229702885453d0 9.035342845140402d0
    -7.02857409395159d0 -6.175722272952579d0 1.145627709456536d0
    0.9292887507635063d0 3.3929641404619453d0 -7.612192400224856d0
    -3.3186492628830035d0 -1.1740750034914573d0 8.161632369953537d0
    1.363261896009406d0 -0.7526689016241428d0 2.808343139397231d0
    2.971411717440667d0 2.207287905457964d0 -1.3035147220741607d0
    0.2741607950170039d0 1.4042685713659875d0 -3.1775134982878863d0
    -2.2372548685189853d0 -1.7625278084286566d0 8.669766877236954d0
    -0.007595520362896499d0 1.6071454162372742d0 2.4649113603941677d0
    -14.769207169430345d0 4.748206712160729d0 4.25942503555152d0
    -0.37080481764155054d0 1.2376052543894571d0 -1.6267698970571955d0
    6.072363417913433d0 -0.1929491109306484d0 -1.4500846733052872d0
    5.098518766481003d0 5.152263335512919d0 7.795564152523721d0
    1.1960270393930943d0 -0.3734228274320816d0 -3.2192134902565286d0
    0.32405789826726916d0 -2.7754393728183784d0 -6.1530401541909d0
    1.3396888141431604d0 -2.06618830946613d0 3.791300118517204d0
    -1.8048198568989613d0 0.017097958698184777d0 2.454106705036069d0
    -1.767952051532204d0 0.3410371434952514d0 0.4940872271651772d0
    -1.0807481577412639d0 1.7444598408802034d0 -2.8531448315010706d0
    -0.08750753731714339d0 1.1641005255554557d0 4.40339036228301d0
    14.101579876799677d0 0.9120318020268278d0 7.077583376721713d0
    -1.8252046343023105d0 0.15752683612682494d0 -1.7416893634834592d0
    -1.1469351780422732d0 0.5214958845569455d0 1.6456753494156284d0
    -0.22352623014516293d0 0.09496773138764744d0 2.416755228167835d0
    -2.102181923716297d0 0.47590150823211397d0 -13.791568790644904d0
    0.5588946952512743d0 -0.6470353723165534d0 -3.0665451290227637d0
    -4.217177355821328d0 -3.5246215685984343d0 -3.169259424186668d0
    -0.8487893265923566d0 -5.886034226339872d0 0.4588601430349549d0
    4.915651925454647d0 -1.242849783448273d0 0.8549277186766596d0
    -1.6703613939415465d0 -3.842796160116949d0 8.770624579865544d0
    -1.0229262282001228d0 4.488355766877113d0 -2.324079162780976d0
    1.0324288899683762d0 -1.1993512928112695d0 -1.4486544219371218d0
    0.9342722908961351d0 -3.9378964732606065d0 1.9414453903015403d0
    1.9131488481776606d0 0.0227721793821479d0 -0.8312982008110424d0
    4.5828630652258d0 8.83474789498359d0 3.777317119445569d0) (list 'vector (geom:vecreal-type))))

(defparameter expected-force-mag (chem:nvector-magnitude expected-force))
(defparameter force-mag (chem:nvector-magnitude force))

(defparameter delta (abs (- energy 20937035.650304314d0)))

(test-true energy-delta (< delta 0.0000000001))

(test-true force-mag (< (abs (- expected-force-mag force-mag)) 0.000001))

(defparameter force-acos (/ (chem:nvector-dot expected-force force) expected-force-mag force-mag))

(test-true force-acos (< force-acos 0.0001))

(format t "force-acos = ~f~%" force-acos)

(defparameter min (chem:make-minimizer ef))
(time (dotimes (i 10)
        (chem:save-coordinates-from-vector ef pos)
        (chem:minimize min)))
(format t "minimizer done~%")
