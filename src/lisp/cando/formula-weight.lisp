(in-package :cando)


(defparameter *most-abundant-isotopes*
  '((:H . 1.007825)
    (:D . 2.0141017)
    (:C . 12.00000)
    (:N . 14.003074)
    (:O . 15.994915)
    (:F . 18.998403)
    (:|Cl| . 34.968853)
    (:|Br| . 78.918338)
    (:S . 31.972071)
    (:B . 11.009305)))

    

(defun formula-string-to-alist (formula)
  (loop with formula-alist = nil
        with rest-formula = formula
        with num = 1
        do (let* ((element-symbol (chem:element-from-atom-name-string rest-formula))
                  (element-len (length (string element-symbol))))
             (setf rest-formula (subseq rest-formula element-len))
             (setf num (cond
                         ((= (length rest-formula) 0) 1) ; end of formula N=1
                         ((digit-char-p (aref rest-formula 0))
                          (let ((next-alpha-pos (position-if #'alpha-char-p rest-formula)))
                            (if (null next-alpha-pos)
                                (prog1
                                    (parse-integer rest-formula) ; last element - read N
                                  (setf rest-formula nil))
                                (prog1
                                    (parse-integer rest-formula :junk-allowed t) ;  - more formula - read N
                                  (setf rest-formula (subseq rest-formula next-alpha-pos))))))
                         (t 1)          ; No N - return 1 more formula
                         ))
             (when (eq element-symbol :ne)
               (error "Could not identify an element in formula ~a" formula))
             (push (cons element-symbol num) formula-alist)
             (when (or (null rest-formula) (string= rest-formula ""))
               (return-from formula-string-to-alist (nreverse formula-alist))))))


(defun formula-string (formula)
  (with-output-to-string (sout)
    (loop for (element . number) in formula
          do (format sout "~a~d" element number))))

(defgeneric formula-weight (formula))

(defmethod formula-weight ((formula cons))
  (loop for (element . number) in formula
        for atomic-weight = (chem:atomic-weight-for-element element)
        for weight = (* atomic-weight number)
        sum weight))

(defmethod formula-weight ((formula-string string))
  (let ((formula (formula-string-to-alist formula-string)))
    (formula-weight formula)))


(defgeneric most-abundant-isotope-weight (formula))

(defmethod most-abundant-isotope-weight ((formula cons))
  (loop for (element . number) in formula
        for most-abundant-isotope = (assoc element *most-abundant-isotopes*)
        sum (if most-abundant-isotope
                (let ((weight (cdr most-abundant-isotope)))
                  (* number weight))
                (error "Could not find most abundant isotope weight for ~a" element))))

(defmethod most-abundant-isotope-weight ((formula-string string))
  (let ((formula (formula-string-to-alist formula-string)))
    (most-abundant-isotope-weight formula)))

(defun chemical-formula (matter)
  "Return the chemical formula as an alist for the matter"
  (let ((element-counts (make-hash-table)))
    (chem:do-atoms (atm matter)
      (let* ((element (chem:atom/get-element atm))
             (count (gethash element element-counts)))
        (if count
            (setf (gethash element element-counts) (1+ count))
            (setf (gethash element element-counts) 1))))
    (let* ((formula (alexandria:hash-table-alist element-counts))
           (sorted-formula (sort formula (lambda (aa bb)
                                           (let ((aa-pos (position aa #(:C :H)))
                                                 (bb-pos (position bb #(:C :H))))
                                             (cond
                                               ((and aa-pos bb-pos (< aa-pos bb-pos)))
                                               ((and (null aa-pos) (null bb-pos) (string< aa bb)))
                                               (t aa-pos))))
                                 :key #'car)))
      sorted-formula)))

(defun chem:formula (matter)
  (chemical-formula matter))

(defgeneric chem:formula-weight (obj))
(defmethod chem:formula-weight ((obj cons))
  (formula-weight obj))
(defmethod chem:formula-weight ((obj string))
  (formula-weight obj))

(defmethod chem:formula-string (formula)
  (formula-string formula))

