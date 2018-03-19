(in-package :charges)


(defun write-sqm-header (fout qm-theory qm-charge &optional (header "Run a qm calculation"))
  "Arguments
fout : A stream to write the header to.
qm-theory : The level of theory ( :AM1 | :PM3 | :PM6)
qm-charge : The overall charge of the molecule"
  (format fout "~a~%" header)
  (format fout " &qmmm~%")
  (format fout "    qm_theory='~a',~%" (string qm-theory))
  (format fout "    qmcharge=~d,~%" qm-charge)
  (format fout " /~%"))

(defun write-sqm-molecule (fout molecule)
  (chem:map-atoms
   nil
   (lambda (atom)
     (let ((position (chem:get-position atom))
           (atom-name (cond
                        ((eq :ZA (chem:get-name atom)) :CA)
                        ((eq :ZB (chem:get-name atom)) :CB)
                        (t (chem:get-name atom)))))
       (format fout "~2d ~a ~12,6f ~12,6f ~12,6f~%"
               (chem:atomic-number-for-element (chem:get-element atom))
               atom-name
               (geom:vx position)
               (geom:vy position)
               (geom:vz position))))
   molecule))

(defun write-sqm-calculation (pathname molecule &key (qm-theory :am1) (qm-charge 0))
  (with-open-file (fout pathname :direction :output)
    (write-sqm-header fout qm-theory qm-charge)
    (write-sqm-molecule fout molecule)))



(defun generate-charge-atoms-array (pathname)
  (let ((atom-number-vec (make-array 50 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
        (element-vec (make-array 50 :fill-pointer 0 :adjustable t))
        (charge-vec (make-array 50 :element-type 'double-float :fill-pointer 0 :adjustable t))
        (line-charge 0))
    (with-open-file (in pathname :if-does-not-exist nil)
      (when in
        (loop for line = (read-line in nil) 
           until (and (> (length line) 22)
                      (string= (subseq line 1 22) "Total Mulliken Charge"))
           do (when (= line-charge 1)
                (multiple-value-bind (atom-number element charge)
                    (parse-mulliken-charge line)
                  (format t  "~a~%" element)
                  (vector-push-extend atom-number atom-number-vec)
                  (vector-push-extend element element-vec)
                  (vector-push-extend charge charge-vec))) 
         do (when (string= line "  Atom    Element       Mulliken Charge")
                (setf line-charge 1))
             )))
    (values atom-number-vec element-vec charge-vec)))

(defun parse-mulliken-charge (line)
  (let ((atom-number (parse-integer line :start 0 :end 6))
        (element (string-trim '(#\space) (subseq line 6 14)))
        (mulliken-charge (fortran::parse-double-float line :start 14 :end 36)))
    (format t " number ~a element ~a charge ~a~%" atom-number element mulliken-charge)
    (values atom-number element mulliken-charge)))
                   

