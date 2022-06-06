
(in-package :leap.parmchk2)

(defun prepare-for-parmchk2 (amber-force-field)
  (let* ((directory (sys:mkdtemp "/tmp/parmchk-XXXXXX"))
         (input-parm-file (make-pathname :name "input-parameters" :type "dat" :defaults directory)))
    (with-open-file (fout input-parm-file :direction :output)
      (format fout "Write force-field here~%"))
    input-parm-file))

