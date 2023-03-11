(load-smirnoff-params "~/Development/openff-sage/inputs-and-results/optimizations/vdw-v1/forcefield/force-field.offxml")

(ql:quickload :topology)
;;;(ql:quickload :design)

(topology:define-topology pro4 ((CAC C :plugs (-dkp.0))
                                (= OA)
                                - (CA :C :stereochemistry-type :chiral) (- HA)
                                - CB (- HB1) (- HB2)
                                - (CG :C :stereochemistry-type :chiral)
                                (- (NG N :plugs (+side.0 +dkp.0)))
                                (- (CGC C :plugs (+dkp.1)) = OG)
                                - CD (- HD1) (- HD2) -
                                (NA N :plugs (-dkp.1)) - CA ))

(topology:define-topology apro4 ((CAC C :plugs (-amide))
                                (= OA)
                                - (CA :C :stereochemistry-type :chiral) (- HA)
                                - CB (- HB1) (- HB2)
                                - (CG :C :stereochemistry-type :chiral)
                                (- (NG N :plugs (+side.0 +dkp.0)))
                                (- (CGC C :plugs (+dkp.1)) = OG)
                                - CD (- HD1) (- HD2) -
                                (NA N :plugs (+amide)) - CA ))

(topology::define-topology bala ((C C :plugs (-amide)) (= O) - CA (- HA1) (- HA2) - CB (- HB1) (- HB2) - (N N :plugs (+amide)) - HN))

(topology::define-topology cbala (C (= O) (- OA - HA) - CA (- HA1) (- HA2) - CB (- HB1) (- HB2) - (N N :plugs (+amide)) - HN))

(topology::define-topology cgly ((C C :plugs (+dkp.1))
                               (= O)
                               - CA (- HA1) (- HA2)
                               - (N N :plugs (+dkp.0))
                               - CN (- HN1) (- HN2) - HN3))

(topology::define-topology ace ((C C :plugs (-amide)) (= O) - CA (- H1) (- H2) - H3))

(topology:define-topology pro 
                          ((C C :plugs (-dkp.0))
                                   (= O)
                               - (CA C :stereochemistry-type :chiral) (- HA)
                               - CB (- HB1) (- HB2)
                               - CG (- HG1) (- HG2)
                               - CD (- HD1) (- HD2)
                               - (N N :plugs (-dkp.1)) - CA ))

(topology:define-topology (bnz side)
                          ((CM C :plugs (-side.0))
                               (- HM1) (- HM2)
                               - C1
                               = C2 (- H2)
                               - C3 (- H3)
                               = C4 (- H4)
                               - C5 (- H5)
                               = C6 (- H6)
                               - C1))

(topology:define-topology (pba side)
                          ((CM C :plugs (-side.0))
                               (- HM1) (- HM2)
                               - C1
                               = C2 (- H2)
                               - C3 (- H3)
                               = C4 (- C4M (- H41) (- H42) - H43)
                               - C5 (- H5)
                               = C6 (- H6)
                               - C1))



(define-foldamer :spiroligomer
  '(bnz pbe pro4 apro4 cbala bala ace cgly)
  '((side :-side pro4 (:-dkp pro4 :-dkp cgly) :+dkp pro4 (:+side (:cap bnz)) :+dkp (:cap pros))
    (side :-side apro4 (:-amide cbala) (:+amide ace) :+dkp pro4 (:+side (:cap bnz)) :+dkp (:cap pros))
    (side :-side apro4 (:-amide cbala) (:+amide ace) :+dkp pro))
    (pro4 (:-dkp pro4 (:-dkp (:cap cgly)) (:+side (:cap bnz))) (:+side (:cap bnz)) :+dkp (:cap pros))
    (pro4 (:-dkp apro4 (:-amide (:cap cbala)) (:+side (:cap bnz)) (:+amide (:cap ace))) :+dkp (:cap pros))
    (apro4 (:-amide cbala) (:+amide (:cap ace)) (:+side (:cap bnz)) :+dkp pro)
    (cbala (:+amide (:cap ace)))
    (bala (:-amide (:cap cbala)) :+amide (:cap ace))
    (ace :-amide (:cap cbala))
    (cgly :+dkp pro))

(topology:define-cap :+dkp pros)
(topology:define-cap :-dkp gly)
(topology:define-cap :-side pro4ss)
(topology:define-cap :+side bnz)

(topology:define-part 'myspiro '(apro4ss (:side bnz)
                                 :dkp ((pro4ss pro4rr)) (:side bnz)
                                 :dkp ((pro4ss pro4rr)) (:side bnz)
                                 :dkp ((pro4ss pro4rr)) (:side bnz)
                                 :dkp pros))

(defparameter olig (topology::make-oligomer '(cbala
                                              :amide myspiro
                                              :amide bala
                                              :amide myspiro
                                              :amide bala
                                              :amide myspiro
                                              :amide bala
                                              :amide ace)))

(topology::number-of-sequences olig)

(topology::goto-sequence olig 456)

conf = (topology:make-conformation olig)

agg = (topology::aggregate conf)

(chem:force-field-name (mol agg 0))

(jostle agg)

(starting-geometry agg)

show agg

(topology::copy-atom-positions-into-joints conf)

(topology::update-joint-tree-internal-coordinates conf)

(save-mol2 agg "/tmp/symba.mol2" :use-sybyl-types t )


