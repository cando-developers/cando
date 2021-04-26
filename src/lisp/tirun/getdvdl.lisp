(ql:quickload :numscl)

(defvar *skip* 5)

(defclass online-av-var ()
  ((step% :initform 0 :initarg :step% :accessor step%)
   (mean :initform 0.0 :initarg :mean :accessor mean)
   (m2   :initform 0.0 :initarg :m2 :accessor m2)
   (store :initform nil :initarg :store :accessor store)
   (data  :initform (make-array 16 :adjustable t :fill-pointer 0))))

(defmethod accumulate ((self online-av-var) x)
  (with-slots (store mean m2 step% data) self
    (incf step%)
    (let ((delta (- x mean)))
      (incf mean (/ delta step%))
      (incf m2 (* delta (- x mean)))
      (when store
        (vector-push-extend x data)))))

(defmethod get-variance ((self online-av-var))
  "Convenience function to return the variance"
  (with-slots (m2 step%) self
    (/ (m2 self) (1- (step% self)))))

(defmethod get-stat ((self online-av-var))
  "Convenience function to return the mean and standard deviation"
  (with-slots (mean m2 step%) self
    (values mean (sqrt (/ m2 (1- step%))))))



(defun onep (num)
  (= num 1.0))


(defun process-linear (x y)
  (unless (find-if #'zerop x)
    (let ((l (/ (- (* (elt x 0) (elt y 1)) (* (elt x 1) (elt y 0)))
                (- (elt x 0) (elt x 1)))))
      (push 0.0 x)
      (push l y)))
  (unless (find-if #'onep x)
    (let* ((x[-2] (elt x (- (length x) 2)))
           (x[-1] (elt x (- (length x) 1)))
           (y[-2] (elt y (- (length y) 2)))
           (y[-1] (elt y (- (length y) 1)))
           (l #n(( (x[-2] - 1.0) * y[-1] + ((1.0 - x[-1]) * y[-2]) ) / (x[-2] - x[-1]))))
      (rplacd (last x) (list 1.0))
      (rplacd (last y) (list l)))))

(defun process-polyfit (x y)
  (let* ((deg (if (< (length x) 6)
                  (1- (length x))
                  6))
         ;; numscl:polyfit returns the coefficients in reverse order from numpy
         (coefs (nreverse
                 (let ((ax (make-array (list 1 (length x))))
                       (ay (make-array (list 1 (length y)))))
                   (setf (aref ax 0 0) 10)
                   (loop for index below (length x)
                         do (setf (aref ax 0 index) (elt x index)))
                   (loop for index below (length y)
                         do (setf (aref ay 0 index) (elt y index)))
                   (numscl:polyfit ax ay deg)))))
    (when (null (find-if #'zerop x))
      (push 0.0 x)
      (push (elt coefs (1- (length coefs))) y))
    (when (null (find-if #'onep x))
      (rplacd (last x) (list 1.0))
      (rplacd (last y) (list (loop for coef in coefs sum coef))))))


(defun read-one-window (en-filename)
  (let ((en-file (open en-filename :direction :input))
        (ln 0)
        (dvdl (make-instance 'online-av-var)))
    (loop for line = (read-line en-file nil :eof)
          until (eq line :eof)
          do (incf ln)
             (when (and (> ln *skip*)
                        (string= line "L9" :start1 0 :end1 2)
                        (null (search "dV/dlambda" line)))
               (let ((val (elt (sys:split line " ") 5)))
                 (accumulate dvdl (read-from-string val)))))
    (multiple-value-bind (mean std)
        (get-stat dvdl)
      (let* ((window-dir (pathname (probe-file en-filename)))
             (window-name (car (last (pathname-directory window-dir))))
             (window-float (read-from-string window-name)))
        (cons window-float (list mean (/ std (sqrt (step% dvdl))) std))))))
    

(defun read-dvdl-windows (windows)
  (loop for window in windows
        for one-data = (read-one-window window)
        collect one-data))

(defun calculate-dvdl (windows)
  (let* ((extrap :polyfit)
         (data (read-dvdl-windows windows))
         (x (mapcar #'car data))
         (y (mapcar (lambda (entry) (car (cdr entry))) data)))
    (case extrap
      (:linear (process-linear x y))
      (:polyfit (process-polyfit x y)))
    (values x y data)))
           
               
(defun run-dvdl (output-filename windows)
  (multiple-value-bind (x y data)
      (calculate-dvdl windows)
    (let ((fout (open output-filename :direction :output :if-exists :supersede)))
      (format fout "(( :lambdas (~%")
      (loop for a in x
            for b in y
            for val-cur = (member a data :key #'car :test #'=)
            do (if val-cur
                   (let ((v (cdr (car val-cur))))
                     (format fout "( ~a ~{~a ~})~%" a v))
                   (format fout ";;; ~a ~a~%" a b)))
      (format fout "))~%")
      (format fout "( :dg . ~a ))~%" (numscl:trapz (coerce y 'vector)
                                                   (coerce x 'vector)))
      (close fout))))

(defun main ()
  (let ((prog (elt sys:*command-line-arguments* 0)))
    (when (< (length sys:*command-line-arguments*) 3)
      (format t "Usage: ~a output-filename windows...~%" prog)
      (sys:exit 1)))
  (let* ((output-filename (elt sys:*command-line-arguments* 0))
         (windows-vec (subseq sys:*command-line-arguments* 1 (length sys:*command-line-arguments*)))
         (windows (coerce windows-vec 'list)))
    (format t "run-dvdl output: ~s\n" output-filename)
    (run-dvdl output-filename windows)))
    
  
(main)
(format t "Done getdvdl~%")
(ext:quit)

