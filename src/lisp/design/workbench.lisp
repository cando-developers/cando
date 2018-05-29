

(:= *olig-simple* (design:make-oligomer *parts* '(:ncap :default :legodkp (:fg1>fg :mbenzoic) (:fg2>fg :mbenzoic) (:fg3>fg :acetyl) :default :ccap)))

(design:define-part *parts* :lego3 '(((:legodkp :legohyd))
                                     (:fg1>fg (:mbenzoic :pbenzoic :acetyl))
                                     (:fg2>fg (:mbenzoic :pbenzoic :acetyl))
                                     (:fg3>fg (:mbenzoic :pbenzoic :acetyl))))

(:= *olig1* (design:make-oligomer *parts* '(:ncap
                                             :default (:lego3 :label :first)
                                             :default (:lego3 :label :second)
                                             :default (:lego3 :label :third)
                                             :default :ccap)))

(:= *olig2* (design:make-oligomer *parts* '((:lego3 :label :first)
                                             :default (:lego3 :label :second)
                                             :default (:lego3 :label :third)
                                             :default (cycle :first))))

(:= *olig3* (design:make-oligomer *parts* '((:lego3 :label :first)
                                            :default :lego3
                                            :default :lego3
                                            :default (cycle :first))))







(:= *parts* (design:make-parts-database))
(design::deftop *parts* (chem:content-with-name *lego3dkp* :lego3dkp) :-default :+default :fg1 :fg2 :fg3)
(design:deftop *parts* (chem:content-with-name *lego3hyd* :lego3hyd) :-default :+default :fg1 :fg2 :fg3)
(design:deftop *parts* (chem:content-with-name *lego2hyd* :lego2hyd) :-default :+default :fg1 :fg2)
(design:deftop *parts* (chem:content-with-name (chem:content-with-name *all* :ends) :ncap) nil :+default)
(design:deftop *parts* (chem:content-with-name (chem:content-with-name *all* :ends) :ccap) :-default)
(design:deftop *parts* (chem:content-with-name *pbenzoic* :pbenzoic) :fg)
(design:deftop *parts* (chem:content-with-name *mbenzoic* :mbenzoic) :fg)
(design:deftop *parts* (chem:content-with-name *acetyl* :acetyl) :fg)
(design:deftop *parts* (chem:content-with-name *benzyl* :benzyl) :fg)
*parts*
