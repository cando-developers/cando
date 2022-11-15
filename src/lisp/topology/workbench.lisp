c(ql:quickload :cando-user)


(defparameter gr (topology::interpret 'pro4 '((CAC C :plugs (-main.0)) (= OA) - (CA :C :stereochemistry-type :chiral) (- HA) - CB (- HB1) (- HB2) - (CG :C :stereochemistry-type :chiral) (- (NG N :plugs (+side.0 +main.0))) (- (CGC C :plugs (+main.1)) = OG) - CD (- HD1) (- HD2) - (NA N :plugs (-main.1)) - CA )))


(defparameter start (topology::interpret 'pro '((CA C :stereochemistry-type :chiral) (- HA) - (C C :plugs (+main.1)) (= O) - CB (- HB1) (- HB2) - CG (- HG1) (- HG2) - CD (- HD1) (- HD2) - (N N :plugs (+main.0)) - CA )))

(defparameter side (topology::interpret 'bnz '((CM C :plugs (-side.0)) (- HM1) (- HM2) - C1 = C2 (- H2) - C3 (- H3) = C4 (- H4) - C5 (- H5) = C6 (- H6) - C1)))

(defparameter end (topology::interpret 'gly '((C C :plugs (-main.0)) - CA (- HA1) (- HA2) - (N N :plugs (-main.1)) - H)))


(defparameter ttt (topology::topologies-from-graph gr))
(topology::topologies-from-graph start)
(topology::topologies-from-graph side)
(topology::topologies-from-graph end)


