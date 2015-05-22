(defvar *a* (Atom :name "C1"))
(chem:setElement *a* :C)
(print (% "element = %s" (chem:getElement *a*)))

(defvar *c* (ChemInfo :smarts "C-[C&+]"))
(print (asXmlString *c*))