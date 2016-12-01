(in-package :leap.topology)

(defclass amber-topology ()
  ((aggregate :initarg :aggregate :accessor aggregate)
   ))

(defun save-amber-parm-format (aggregate topology-pathname coordinate-pathname force-field)
  (let* ((energy-function (chem:make-energy-function aggregate force-field
                                                     :use-excluded-atoms t))
         (nonbonds (chem:get-nonbond-component energy-function))
         (number-excluded-atoms (chem:number-excluded-atoms nonbonds))
         (excluded-atom-list (chem:excluded-atom-list nonbonds)))
    ;; Skip assigning MarkMainChainAtoms and MarkSideChain atoms for now
    ;; see (unitio.c:4889).  This won't mean anything for spiroligomers.
    (with-fortran-file (ftop topology-pathname :direction :output)
      (fortran:debug-on ftop)
      (fortran:debug "-1-")        ;
      (fortran:format 1 "%-80s")
      (fortran:write (core:strftime 81 "%%VERSION  VERSION_STAMP = V0002.000  DATE = %m/%d/%y  %H:%M:%S"))
      (fortran:write "%FLAG TITLE")
      (fortran:write "%FORMAT(20a4)")
      (fortran:write (chem:get-name aggregate))
      (fortran:format 1 "%-80s")
      (fortran:write "%FLAG POINTERS")
      (fortran:write "%FORMAT(10I8)")
      (fortran:debug "-2-")
      (fortran:format 10 +intformat+)
      ;; NTOTAT
      (fortran:write (chem:get-number-of-atoms (chem:atom-table energy-function)))
      ;; pick up at unitio.cc:4969
      )))
