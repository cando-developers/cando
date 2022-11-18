(ql:quickload :topology)
;;;(ql:quickload :design)

(topology::define-topology gly ((C C :plugs (+dkp.1))
                               (= O)
                               - CA (- HA1) (- HA2)
                               - (N N :plugs (+dkp.0))
                               - CN (- HN1) (- HN2) - HN3))

(topology:define-topology pro4 ((CAC C :plugs (-dkp.0))
                                (= OA)
                                - (CA :C :stereochemistry-type :chiral) (- HA)
                                - CB (- HB1) (- HB2)
                                - (CG :C :stereochemistry-type :chiral)
                                (- (NG N :plugs (+side.0 +dkp.0)))
                                (- (CGC C :plugs (+dkp.1)) = OG)
                                - CD (- HD1) (- HD2) -
                                (NA N :plugs (-dkp.1)) - CA ))

(topology:define-topology pro ((C C :plugs (-dkp.0))
                               (= O)
                               - (CA C :stereochemistry-type :chiral) (- HA)
                               - CB (- HB1) (- HB2)
                               - CG (- HG1) (- HG2)
                               - CD (- HD1) (- HD2)
                               - (N N :plugs (-dkp.1)) - CA ))

(topology:define-topology bnz ((CM C :plugs (-side.0))
                               (- HM1) (- HM2)
                               - C1
                               = C2 (- H2)
                               - C3 (- H3)
                               = C4 (- H4)
                               - C5 (- H5)
                               = C6 (- H6)
                               - C1))


(topology:define-cap :+dkp pros)
(topology:define-cap :-dkp gly)
(topology:define-cap :-side pro4ss)
(topology:define-cap :+side bnz)

(defparameter foldamer (topology::make-foldamer '(gly pro4ss pro4sr pro4rs pro4rr bnz pros pror)))

(defparameter oligs (topology::build-training-oligomers foldamer))

(topology::build-all-molecules (first oligs))

(chem:find-topology 'pros t)




(defparameter olig (topology::make-oligomer '(:gly :dkp :pro4ss (:side :bnz) :dkp :pro4rr (:side :bnz) :dkp :pro4ss (:side :bnz) :dkp :pro4ss (:side :bnz) :dkp :pro4rr (:side :bnz) :dkp :pro4ss (:side :bnz) :dkp :pros)))


(defparameter mol (topology::build-molecule olig))
(easy-gaff2)
(assign-atom-types mol)
(jostle mol)
(starting-geometry mol)

(save-mol2 mol "/tmp/mol.mol2" :use-sybyl-types t )
