(
 (:H
  ;; Sidechain hydrogens
  (:|Haro| "[#1]-[#6;a]")
  (:|Hpol| "[#1]-[#7]")
  (:|Hpol| "[#1]-[#8]")
  (:|Hapo| "[#1]-[#6;!a]")
  (:|HNbb| "[#1]-[#7;X3](-[#6;X4])-[#6;X3](=[#8])~[#6]") ; Backbone
  )
 (:C
  ;; Sidechain carbons
  (:|aroC| "[#6;a]")
  (:COO "[#6;X3](~[#8])~[#8]")
  (:CNH2 "[#6;X3](=[#8])-[#7]")
  (:CH3 "[#6;X4;H3]")
  (:CH2 "[#6;X4;H2]")
  (:CH1 "[#6;X4;H1]")
  (:CH1 "[#6;X4;H0;!a]")
  (:|CAbb| "[#6;X4](-[#7;X3])-[#6;X3](=[#8])~[#7]") ; peptoid backbone
  (:|CObb| "[#6;X3](=[#8])(~[#7;X3])~[#6;X4]")      ; peptoid backbone
  )
 (:N
  ;; Sidechain nitrogens
  (:|Npro| "[#7+](=[#8])([#8-])[#8-]")
  (:|Ntrp| "[#7;a;H1]")
  (:|Nhis| "[#7;a;H0]")
  (:|Nhis| "[#7;X3;H0;!$([#7]=[#8])]")
  (:|Nlys| "[#7+;X4]")
  (:|Nlys| "[#7;X4;H3]")
  (:|Nlys| "[#7;X4;H1;!a;!$([#7]=[#8])]")
  (:|Nbb|  "[#7;X3;H0](-[#6;X4])-[#6;X3](=[#8])~[#6]") ; peptoid backbone
  (:|Nbb|  "[#7;X3;H1](-[#6;X4])-[#6;X3](=[#8])~[#6]") ; peptoid backbone
  (:|Npro| "[#7;X3](~[#8;X1])(~[#6])~[#8;X1]") ; nitro nitrogen groups peptoid012
  )
 (:O
  ;; Sidechain oxygens
  (:|Oaro| "[#8;a]")
  (:ONH2 "[#8;X1]~[#7;X3](~[#6])~[#8;X1]") ; nitro oxygen groups
  (:OOC "[#8-]-[#6](=O)")
  (:OOC "[#8;X1]-[#6](=O)")
  (:OOC "[#8;X1]=[#6](-O)")
  (:|ONH2| "[#8]=[#7+]")
  (:OH "[#8;X2;H1]")
  (:OH "[#8;X2;H0;!R](-[#6])-[#6]") ; acyclic diphenyl ether, pheny-alkyl ethers, alkyl ethers
  (:OH "[#8;X2;H0;R;!a]") ; cyclic alkyl ethers
;;;  (:OH "[#8;X2;H0]")
  (:|OCbb| "[#8;X1]=[#6;X3](~[#7;X3])~[#6;X4]") ; peptoid backbone
  (:|Oaro| "[#8;X2;R]~[#6;a]")
  )
 (:F
  (:F "[#9]"))
 (:|Cl|
   (:|Cl| "[#17]"))
 )
