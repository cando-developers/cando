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
        
        
