;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :jsme)

(defparameter *sketch* "")

;;; Old url http://peter-ertl.com/jsme/JSME_2017-02-26/jsme/jsme.nocache.js

(defparameter *jsme-form* "<script type=\"text/javascript\" src=\"jsme/jsme.nocache.js\"></script>
<script type=\"text/javascript\">
function jsmeOnLoad() {
    //arguments: HTML id, width, height (must be string not number!)
    jsmeApplet12345 = new JSApplet.JSME(\"appletContainer12345\", \"640px\", \"400px\", {
                     //optional parameters
                     \"options\" : \"number,query,hydrogens,paste,zoom,atommovebutton,border\"
    });
    var jme = `~A`;
    jsmeApplet12345.readMolFile(jme);
}
function onSubmit() {
    var drawing = jsmeApplet12345.molFile();
    // document.getElementById('jsme_smiles12345').value = drawing;
    var command = '(defparameter jsme::*sketch* \"' + drawing + '\")';
    console.log(\"Executing Command: \" + command);
    var kernel = IPython.notebook.kernel;
    kernel.execute(command);
}
</script>
<table align=\"left\" style=\"border: none;\">
<tr style=\"border: none;\">
<td id=\"appletContainer12345\" style=\"border: none;\"></td>
<td style=\"vertical-align: bottom; border: none;\">
<button onclick=\"onSubmit()\">done !</button>
</td>
</tr>
</table>
")

(esrap:defrule skippable
    (+ (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(esrap:defrule skippable?
    (* (or parser.common-rules:shell-style-comment
           parser.common-rules:whitespace)))

(parser.common-rules:defrule/s number
    (+ (esrap:character-ranges (#\0 #\9)))
  (:lambda (x)
    (parse-integer (esrap:text x))))

(parser.common-rules:defrule/s label
    (+ (esrap:character-ranges (#\A #\Z) (#\a #\z) #\?))
  (:text t))

(defparameter *elements*
  '("He" "Li" "Be" "Ne" "Na" "Mg" "Al" "Si" "Cl" "Ar" "Ca" "Sc" "Ti" "Cr" "Mn" "Fe" "Co"
    "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh"
    "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "Xe" "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu"
    "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Hf" "Ta" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl"
    "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac" "Th" "Pa" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es"
    "Fm" "Md" "No" "Lr"
    "H" "W" "U" "B" "C" "N" "O" "F" "P" "S" "K" "V" "I" "Y"))

(esrap:defrule atom-line
    (and parser.common-rules:float-literal
         parser.common-rules:float-literal
         parser.common-rules:float-literal
         parser.common-rules:whitespace
         label/s
         number/s number/s number/s number/s number/s number/s number/s number/s number/s number/s number/s number
         )
  (:lambda (x)
    x))

(esrap:defrule bond-line
    (and number/s number/s number/s number/s number/s number/s number)
  (:lambda (x)
    x))

(esrap:defrule header-line
    (and number/s number/s number/s number/s number/s number/s number/s number/s number/s "0999 V2000"))

#|
(esrap:parse 'header-line "8  8  0  0  0  0  0  0  0  0999 V2000")
(esrap:parse 'atom-line "  0.0000    4.9000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0")
|#


(defun sketch (&optional (mol *sketch*))
  (unless mol
    (setf mol ""))
  (funcall (find-symbol "HTML" :cl-jupyter-user) (format nil *jsme-form* mol)))

(defun load-sketch (filename)
  (with-open-file (fin filename :direction :input)
    (let ((str (make-string (file-length fin))))
      (read-sequence str fin)
      (setf *sketch* str)
      '*sketch*
      )))

(defun save-sketch (filename &optional (sketch *sketch*))
  (unless sketch
    (error "There is no sketch"))
  (with-open-file (fout filename :direction :output :if-exists :supersede)
    (write sketch :stream fout))
  nil)




(defun build-atom-vector (atoms bonds)
  (let ((atom-vector (make-array (length atoms)))
        (map-numbers (make-hash-table)))
    (loop for atom-info in atoms
          for index from 0
          for xpos = (elt atom-info 0)
          for ypos = (elt atom-info 1)
          for zpos = (elt atom-info 2)
          for raw-map-number = (elt atom-info 14)
          for map-number = (if (= raw-map-number 0) nil raw-map-number)
          for label-string = (string (elt atom-info 4))
          for match-element-string = (progn
                                       (find label-string *elements* :test #'string=))
          for element-string = (if match-element-string match-element-string "C")
          for element = (intern element-string :keyword)
          for label = (if match-element-string nil (intern label-string :keyword))
          for atom-name = (intern (format nil "~a~a" element-string index) :keyword)
          for atom = (chem:make-atom atom-name element)
          do (chem:set-position atom (geom:vec xpos ypos zpos))
          do (when label (chem:set-property atom :label label))
          do (when map-number (setf (gethash atom-name map-numbers) map-number))
          do (setf (elt atom-vector index) atom))
    (loop for bond-info in bonds
          for i1 = (1- (elt bond-info 0))
          for i2 = (1- (elt bond-info 1))
          for order-index = (elt bond-info 2)
          for order = (elt '(:no-bond :single-bond :double-bond :triple-bond) order-index)
          for atom1 = (elt atom-vector i1)
          for atom2 = (elt atom-vector i2)
          do (chem:bond-to atom1 atom2 order))
    (values atom-vector map-numbers)))

(defun parse-groups (mol-string)
  (with-input-from-string (sin mol-string)
    (let* ((smiles (read-line sin))
           (label (read-line sin))
           (blank (read-line sin))
           (header-line (string-trim '(#\space) (read-line sin)))
           (header (esrap:parse 'header-line header-line))
           (num-atoms (first header))
           (num-bonds (second header))
           (atoms (loop for atom-index below num-atoms
                        for line = (string-trim '(#\space) (read-line sin))
                        collect (esrap:parse 'atom-line line)))
           (bonds (loop for bond-index below num-bonds
                        for line = (string-trim '(#\space) (read-line sin))
                        collect (esrap:parse 'bond-line line))))
      (multiple-value-bind (atom-vector map-names-numbers)
          (jsme:build-atom-vector atoms bonds)
        (let ((atom-to-index (make-hash-table))
              groups)
          (loop for index below (length atom-vector)
                for atom = (aref atom-vector index)
                do (setf (gethash atom atom-to-index) index))
          (loop for index below (length atom-vector)
                for atom = (aref atom-vector index)
                when atom
                  do (let* ((spanning-loop (chem:make-spanning-loop atom))
                            (group (loop for span-atom in (chem:all-atoms spanning-loop)
                                         for index = (gethash span-atom atom-to-index)
                                         do (setf (aref atom-vector index) nil)
                                         collect span-atom)))
                       (push group groups)))
          (values groups map-names-numbers))))))
