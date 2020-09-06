#|
This code was all copied from the rosettacode.org/wiki/Rosetta_Code

|#

(in-package :numscl)




;; Calculates the Cholesky decomposition matrix L 
;; for a positive-definite, symmetric nxn matrix A.
(defun chol (A)
  (let* ((n (car (array-dimensions A)))
         (L (make-array `(,n ,n) :initial-element 0)))
 
    (do ((k 0 (incf k))) ((> k (- n 1)) nil)
      ;; First, calculate diagonal elements L_kk.
      (setf (aref L k k)
            (sqrt (- (aref A k k)
                     (do* ((j 0 (incf j))
                           (sum (expt (aref L k j) 2) 
                                (incf sum (expt (aref L k j) 2))))
                          ((> j (- k 1)) sum)))))
 
      ;; Then, all elements below a diagonal element, L_ik, i=k+1..n.
      (do ((i (+ k 1) (incf i)))
          ((> i (- n 1)) nil)
 
        (setf (aref L i k)
              (/ (- (aref A i k)
                    (do* ((j 0 (incf j))
                          (sum (* (aref L i j) (aref L k j))
                               (incf sum (* (aref L i j) (aref L k j)))))
                         ((> j (- k 1)) sum)))
                 (aref L k k)))))
 
    ;; Return the calculated matrix L.
    L))


(defun mmul (A B)
  (let* ((m (car (array-dimensions A)))
         (n (cadr (array-dimensions A)))
         (l (cadr (array-dimensions B)))
         (C (make-array `(,m ,l) :initial-element 0)))
    (loop for i from 0 to (- m 1) do
      (loop for k from 0 to (- l 1) do
        (setf (aref C i k)
              (loop for j from 0 to (- n 1)
                    sum (* (aref A i j)
                           (aref B j k))))))
    C))

;; Solve a linear system AX=B where A is symmetric and positive definite, so it can be Cholesky decomposed.
(defun linsys (A B)
  (let* ((n (car  (array-dimensions A)))
         (m (cadr (array-dimensions B)))
         (y (make-array n        :element-type 'long-float :initial-element 0.0L0))
         (X (make-array `(,n ,m) :element-type 'long-float :initial-element 0.0L0))
         (L (chol A))) ; A=LL'
    (loop for col from 0 to (- m 1) do
      ;; Forward substitution: y = L\B
      (loop for k from 0 to (- n 1)
            for numerator = (- (aref B k col)
                               (loop for j from 0 to (- k 1)
                                     sum (* (aref L k j)
                                            (aref y j))))
            for denomenator = (aref L k k)
            do (setf (aref y k) (/ numerator denomenator)))
       ;; Back substitution. x=L'\y
       (loop for k from (- n 1) downto 0
             do (setf (aref X k col)
                      (/ (- (aref y k)
                            (loop for j from (+ k 1) to (- n 1)
                                  sum (* (aref L j k)
                                         (aref X j col))))
                         (aref L k k)))))
    X))
 
;; Solve a linear least squares problem. Ax=b, with A being mxn, with m>n.
;; Solves the linear system A'Ax=A'b.
(defun lsqr (A b)
  (linsys (mmul (mtp A) A)
          (mmul (mtp A) b)))

(defun transpose (m)
  (apply #'mapcar #'list m))

;; Transpose a mxn matrix A to a nxm matrix B=A'.
(defun mtp (A)
  (let* ((m (array-dimension A 0))
         (n (array-dimension A 1))
         (B (make-array `(,n ,m) :initial-element 0)))
    (loop for i from 0 below m do
          (loop for j from 0 below n do
                (setf (aref B j i)
                      (aref A i j))))
    B))


;; Least square fit of a polynomial of order n the x-y-curve.
(defun polyfit (x y n)
  (let* ((m (cadr (array-dimensions x)))
         (A (make-array `(,m ,(+ n 1)) :initial-element 0)))
    (loop for i from 0 to (- m 1) do
          (loop for j from 0 to n do
                (setf (aref A i j)
                      (expt (aref x 0 i) j))))
    (let* ((mtpy (mtp y))
           (rev-result (cl-mathstats::svd-solve-linear-system A mtpy nil)))
      (nreverse rev-result))))


