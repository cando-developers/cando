(
 (:H
  (:HS "[#1](~[#16])")
  (:|Hpol| "[#1](~[#8])")
  (:|Hpol| "[#1](~[#7])")
  (:HNBB "[#1](~[#7;X2])")
  (:|Hapo| "[#1]-[#6;A;!a]")
  (:|Haro| "[#1]-[#6;a]")
  (:|Hapo| "[#1]-[#6;r5](~[#7])~[#7]") ; His HE1
  (:|Hapo| "[#1]-[#6;r5](~[#6])~[#7]") ; His HD2
  (:|Haro| "[#1]-[#6;r5](~[#7])~[#6;H0]~[#6;H0]~[#6]") ; Trp HD1
  (:|HNbb| "[#1]-[#7;X3](-[#6]=[#8])-[#6;X4]") ; Amide nitrogen backbone
  )
 (:C
  (:|aroC| "[#6](~[#7])(~[#7])~[#7]") ; Arg
  (:|aroC| "[#6;r5](~[#7])~[#7]")      ; His
  (:|aroC| "[#6;r5](~[#7])~[#6]")      ; His
  (:|aroC| "[#6;a]") ; Trp
  (:CH0 "[#6;r5;H0]") ; His
  (:CH0 "[#6;a](-[#6;X4;H2])~[#6;a]~[#6;a]~[#6;a]") ; Phe/Tyr uses CH0 for the para (from OH) aromatic carbon????
  (:CH0 "[#6;a]-[#8]") ; Tyr uses CH0 for the ipso (from OH) aromatic carbon????
  (:CH1 "[#6;H1;X4]")
  (:CH2 "[#6;H2;X4]")
  (:CH3 "[#6;H3;X4]")
  (:COO "[#6;!R](~[#8])~[#8]")
  (:|CAbb| "[#6;X4](-[#7])-[#6]=[#8]") ; Backbone alpha carbon
  (:|CObb| "[#6;X3](=[#8])(-[#6])-[#7;X3]") ; Backbone amide carbonyl
  (:CNH2 "[#6;X3](=[#8])(-[#6])-[#7;H2]") ; Backbone amide carbonyl
  )
 (:N
  (:|Narg| "[#7]~[#6](~[#7])~[#7]")
  (:NH2O "[#7](~[#16])")
  (:NHIS "[#7](~[#7])")
  (:NTRP "[#7;a;H1]")
  (:|Nlys| "[#7;H3]-[#6;X4]")
  (:|NtrR| "[#7;!R;H1;X3]")
  (:|Nbb| "[#7;X3;H1](-[#6;X4])-[#6;X3]=[#8]") ; Amide nitrogen backbone
  (:|Nbb| "[#7;X3;H0](-[#6;X4])-[#6;X3]=[#8]") ; Peptoid Amide nitrogen backbone
  (:NH2O "[#7;X3;H2]-[#6;X3](=[#8])-[#6;X4]")
  (:|Nhis| "[#7;r5](~[#6])~[#6]")
  (:|Npro| "[#7;r5](~[#6;X4])~[#6;X4]") ; Pro N
  (:|Ntrp| "[#7;r5;H1](~[#6])~[#6]")
  )
 (:O
  (:OOC "[#8](~[#16])")
  (:OH "[#8](~[#6;H2])") ; Ser
  (:OH "[#8](~[#6;H1])") ; Thr
  (:OH "[#8]-[#6;a]") ; Tyr
  (:|Oaro| "[#8;R](~[#6;H0])")
  (:|OCbb| "[#8]=[#6;H0;X3](-[#6])-[#7;X3]") ; Amide carbonyl backbone
  (:|ONH2| "[#8]=[#6;H0;X3](-[#6])-[#7;H2]") ; Terminal amide carbonyl
  (:OOC "[#8]~[#6;X3](~[#8])~[#6]") ; Sidechain carboxyl
  )
 (:F
  (:F "[#9]"))
 (:S
  (:S "[#16;!R](-[#6;H2;X4])-[#6;H3;X4]")
  (:SH1 "[#16;H1]"))
 (:|Cl|
  (:CL "[#17]"))
)
