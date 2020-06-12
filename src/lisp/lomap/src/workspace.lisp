(setf *default-pathname-defaults* (pathname (format nil "~a/Development/clomap/" (ext:getenv "HOME"))))

(asdf:load-asd (pathname (format nil "~a/Development/clomap/lomap.asd" (ext:getenv "HOME"))))

(push :debug-lomap *features*)
(progn
  (asdf:load-system :lomap)
  (format t "Done~%"))

(in-package :lomap)

(defparameter *mols* (sdf:load-sdf-as-list-of-molecules (format nil "~a/Development/fep-benchmark/eg5/ligands.sdf" (ext:getenv "HOME"))))

;;; Smaller graph
(defparameter *mols4* (subseq *mols* 6 10))
(defparameter *mat4* (similarity-matrix *mols4*))
(defparameter *graph* (similarity-graph *mols4* *mat4*))

;;; Larger graph
(defparameter *mols11* (subseq *mols* 6 16))
(defparameter *mat11* (similarity-matrix *mols11*))
(defparameter *graph* (similarity-multigraph *mols11* *mat11*))


;;; All molecules graph
(defparameter *mat* (similarity-matrix *mols*))
(defparameter *graph* (similarity-multigraph *mols* *mat*))

(multigraph-set-positions *graph*
                     '(
                       ("CHEMBL1088740" 400 190)
                       ("CHEMBL1078998" 558 276)
                       ("CHEMBL1084935" 1739 18)
                       ("CHEMBL1084676" 1172 104)
                       ("CHEMBL1093088" 202 276)
                       ("CHEMBL1089056" 347 362)
                       ("CHEMBL1096002" 1583 104)
                       ("CHEMBL1084143" 1093 190)
                       ("CHEMBL1085692" 1305 276)
                       ("CHEMBL1085666" 2106 104)
                       ("CHEMBL1083836" 1886 190)
                       ("CHEMBL1077204" 971 362)
                       ("CHEMBL1096003" 1247 448)
                       ("CHEMBL1083517" 2228 276)
                       ("CHEMBL1086410" 1971 362)
                       ("CHEMBL1077227" 2289 448)
                       ("CHEMBL1084677" 1275 534)
                       ("CHEMBL1085895" 173 362)
                       ("CHEMBL1084431" 824 448)
                       ("CHEMBL1084115" 453 448)
                       ("CHEMBL1084678" 1064 620)
                       ("CHEMBL1082248" 1936 448)
                       ("CHEMBL1086409" 2050 534)
                       ("CHEMBL1078691" 287 534)
                       ("CHEMBL1089393" 499 620)
                       ("CHEMBL1082249" 1500 534)
                       ("CHEMBL1093087" 1205 706)
                       ("CHEMBL1078774" 414 706)
                       ))
 *graph*

(progn
  (time (defparameter *new-graph* (lomap-multigraph *graph* :debug t)))
  (format t "Done~%"))



