(in-package :charges)

(defun write-mopac-header (fout &optional (header "Run a mopac calculation"))
  "Arguments
fout : A stream to write the header to."
  (format fout "AM1 XYZ~%")
  (format fout "~a~%" header)
  (format fout "~%"))

(defun write-mopac-molecule (fout molecule)
  (let ((atom-order (chem:map-atoms
                'list
                (lambda (atom)
                  (let ((position (chem:get-position atom)))
                    (format fout "~3a~8,3f~2d~8,3f~2d~8,3f~2d~%"
                            (chem:get-element atom)
                            (geom:get-x position)
                            0
                            (geom:get-y position)
                            0
                            (geom:get-z position)
                            0))
                  atom)
                molecule)))
    atom-order))

(defun write-mopac-calculation (fout molecule)
  (write-mopac-header fout)
  (write-mopac-molecule fout molecule))

(defun generate-mopac-charge-atoms-array (pathname version) 
  (let ((atom-number-vec (make-array 50 :element-type '(signed-byte 32) :fill-pointer 0 :adjustable t))
        (element-vec (make-array 50 :fill-pointer 0 :adjustable t))
        (charge-vec (make-array 50 :element-type (geom:vecreal-type) :fill-pointer 0 :adjustable t))
        (line-charge 0)
        data)
    (with-open-file (in pathname :if-does-not-exist nil)
      (setf data (make-string (file-length in)))
      (read-sequence data in))
    (case version
      (:2016
       (with-input-from-string (in data)
         (when in
           (loop for line = (read-line in nil)
              when (null line) do (return)
              until (and line
                         (> (length line) 50)
                         (string= (subseq line 1 51) "DIPOLE           X         Y         Z       TOTAL"))
             ;; do (format t "line ~a line-charge ~a~%" line line-charge)
              do (when (= line-charge 1)
                   (multiple-value-bind (atom-number element charge)
                       (parse-mopac-charge-2016 line)
                     (vector-push-extend atom-number atom-number-vec)
                     (vector-push-extend element element-vec)
                     (vector-push-extend charge charge-vec)))
              do (when (and (> (length line) 71)
                            (string= (subseq line 2 72)
                                     "ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop"))
                   (setf line-charge 1))))))
      (otherwise
       (with-input-from-string (in data)
         (when in
           (loop for line = (read-line in nil)
              when (null line) do (return)
              until (and line
                         (> (length line) 50)
                         (string= (subseq line 1 51) "DIPOLE           X         Y         Z       TOTAL"))
             ;; do (format t "line ~a line-charge ~a~%" line line-charge)
              do (when (= line-charge 1)
                   (multiple-value-bind (atom-number element charge)
                       (parse-mopac-charge line)
                     (vector-push-extend atom-number atom-number-vec)
                     (vector-push-extend element element-vec)
                     (vector-push-extend charge charge-vec)))
              do (when (and (> (length line) 69)
                            (string= (subseq line 9 70)
                                     "ATOM NO.   TYPE          CHARGE        ATOM  ELECTRON DENSITY"))
                   (setf line-charge 1)))))))
    (values atom-number-vec element-vec charge-vec)))
  
(defun parse-mopac-charge-2016 (line)
  (let ((atom-number (parse-integer line :start 0 :end 5))
        (element (intern (string-trim '(#\space) (subseq line 6 16)) :keyword))
        (charge (fortran::parse-vecreal line :start 16 :end 35)))
   ;;; (format t " number ~a element ~a charge ~a~%" atom-number element charge)
    (values atom-number element charge)))

(defun parse-mopac-charge (line)
  (let ((atom-number (parse-integer line :start 0 :end 12))
        (element (intern (string-trim '(#\space) (subseq line 13 24)) :keyword))
        (charge (fortran::parse-vecreal line :start 25 :end 42)))
   ;;; (format t " number ~a element ~a charge ~a~%" atom-number element charge)
    (values atom-number element charge)))

     
(defun read-mopac-am1-charges (pathname atom-order version)
 "Read the AM1 charges from the mopac output file.
Return a hash-table of atoms -> charge or NIL if the charges couldn't be read"
 (let ((correction-map (make-hash-table))
       (count 0))
   (multiple-value-bind (atom-number-vec element-vec charge-vec)
       (generate-mopac-charge-atoms-array pathname version)
     (if (= (length atom-order) (length atom-number-vec))
         (mapc (lambda (a)
                 (if (eq (chem:get-element a ) (aref element-vec count))
                     (setf (gethash a correction-map) (aref charge-vec count))
                     (error "Mismatch atom name ~a element ~a element ~a" (chem:get-name a) (chem:get-element a) (aref atom-number-vec count)))
                 (incf count))
               atom-order)
         nil)
     correction-map)))
