(in-package :cando-user)

(defun show (matter)
  (cond
    ((typep matter 'chem:aggregate))
    ((typep matter 'chem:molecule)
     (let ((agg (chem:make-aggregate nil)))
       (chem:add-matter agg matter)
       (setf matter agg)))
    ((typep matter 'chem:residue)
     (let ((agg (chem:make-aggregate nil))
           (mol (chem:make-molecule nil)))
       (chem:add-matter agg mol)
       (chem:add-matter mol matter)
       (setf matter agg)))
    (t (error "You cannot show a ~a" matter)))
  (cando-jupyter:show-aggregate matter))

(defun repr (widget representation &optional (selection "all"))
  (funcall (find-symbol "ADD-REPRESENTATION" :nglv) widget representation :selection selection))



