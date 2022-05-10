(in-package :cif)

(defclass cif-block ()
  ((name :initarg :name :accessor name)
   (lines :initform nil :accessor lines)))

(defclass cif ()
  ((data :initform nil :accessor data)
   (dict :initform (make-hash-table) :accessor dict)))

(defclass data-block ()
  ((name :initarg :name :accessor name)
   (dict :initarg :dict :accessor dict)))

(defclass cif-line ()
  ((lineno :initarg :lineno :accessor lineno)
   (line :initarg :line :accessor line)))

(defmethod print-object ((obj cif-line) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~s" (lineno obj) (line obj))))

(defun parse-cif-block (cur cur-block)
  (let ((rev-lines nil))
    (let ((new-cur (loop named lines
                         for pcur = cur
                         for ccur = cur then (cdr ccur)
                         for entry = (car ccur)
                         when (null ccur)
                           do (return-from lines ccur)
                         do (let ((line (line entry)))
                              (when (string= line "data_" :start1 0 :end1 5)
                                (return-from lines pcur))
                              (push entry rev-lines)))))
      (setf (lines cur-block) (nreverse rev-lines))
      new-cur)))

(defun process-line (line cif fin &key debug)
  (declare (ignore fin))
  (cond
    ((char= (char line 0) #\#)
                                        ; nothing - comment
     )
    ((string= line "data_" :start1 0 :end1 5)
     (let ((str (subseq line 5)))
       (when debug (format t "data_ |~a|~%" str))
       (setf (data cif) str)
       ))
    (t
     (format t "???? ~a~%" line))))

(defun parse-multi-line-value (lines)
  (loop named lines
        for pcur = cur
        for cur = lines then (cdr cur)
        for entry = (car cur)
        when (null cur)
          do (error "expecting closing ; ~a" (car pcur))
        collect (let ((line (line entry)))
                  (when (char= (char line 0) #\;)
                    (return-from lines (values cur value)))
                  line) into value))

(defun parse-values (line values)
  (loop for read-value in (core:split line " ")
        for value in values
        do (vector-push-extend read-value value)))

(defun process-loop (lines dict)
  (multiple-value-bind (cur keys)
      (loop named lines
            for cur = lines then (cdr cur)
            for entry = (car cur)
            unless entry
              do (error "File terminated before loop_ data was complete")
            collect (let ((line (line entry)))
                      (when (not (char= (char line 0) #\_))
                        ;; data is done
                        (return-from lines (values cur keys)))
                      line) into keys)
    (let ((values (loop for key in keys
                        collect (make-array 16 :adjustable t :fill-pointer 0))))
      (let ((after-cur (loop named value-lines
                             for pcur = ncur
                             for ncur = cur then (cdr ncur)
                             for entry = (car ncur)
                             unless ncur
                               do (return-from value-lines ncur)
                             do (let ((line (line entry)))
                                  (when (char= (char line 0) #\_)
                                    (return-from value-lines pcur))
                                  (parse-values line values)))))
        (loop for key in keys
              for value in values
              do (setf (gethash key dict) value))
        after-cur
        ))))

(defun convert-number (number-string)
  (let ((open-paren (position #\( number-string)))
    (cond
      (open-paren
       (let ((just-number (subseq number-string 0 open-paren)))
         (fortran:parse-double-float just-number)))
      (t (fortran:parse-double-float number-string)))))

(defun elements-match (el1 el2 ela elb)
  (or (and (eq el1 ela) (eq el2 elb))
      (and (eq el1 elb) (eq el2 ela))))

(defun calculate-bonds (agg)
  (flet ((bonded-order (a1 a2 el1 el2 dist)
           (cond
             ((elements-match el1 el2 :c :c)
              (cond
                ((< 1.45 dist 1.55) :single-bond)
                ((< 1.29 dist 1.45) :double-bond)
                ((< dist 1.29) :triple-bond)))
             ((elements-match el1 el2 :c :o)
              (cond
                ((< 1.23 dist 1.3)
                 :double-bond)
                ))
             ((elements-match el1 el2 :c :s)
              (cond
                ((< 1.7 dist 1.82)
                 :single-bond)
                ))
             ((elements-match el1 el2 :|Rh| :|Rh|)
              (cond
                ((< 2.3 dist 2.5)
                 :single-bond)
                )))))
    (let ((atoms (chem:map-atoms 'list #'identity agg)))
      (loop named outer-loop
            for atm-cur = atoms then (cdr atm-cur)
            for atom1 = (car atm-cur)
            for element1 = (chem:get-element atom1)
            for pos1 = (chem:get-position atom1)
            when (null (cdr atm-cur))
              do (return-from outer-loop nil)
            do (loop for atom2 in (cdr atm-cur)
                     for element2 = (chem:get-element atom2)
                     for pos2 = (chem:get-position atom2)
                     for delta = (geom:v- pos1 pos2)
                     for dist = (geom:vlength delta)
                     for order = (bonded-order atom1 atom2 element1 element2 dist)
                     when order
                       do (chem:bond-to atom1 atom2 order))))))

(defun generate-aggregate (data)
  (let ((sname (intern (name data) :keyword))
        (dict (dict data)))
    (let* ((_cell_angle_alpha (convert-number (gethash "_cell_angle_alpha" dict)))
           (_cell_angle_beta (convert-number (gethash "_cell_angle_beta" dict)))
           (_cell_angle_gamma (convert-number (gethash "_cell_angle_gamma" dict)))
           (alpha (* _cell_angle_alpha 0.0174533))
           (beta (* _cell_angle_beta 0.0174533))
           (gamma (* _cell_angle_gamma 0.0174533))
           (_cell_length_a (convert-number (gethash "_cell_length_a" dict)))
           (_cell_length_b (convert-number (gethash "_cell_length_b" dict)))
           (_cell_length_c (convert-number (gethash "_cell_length_c" dict))))
      (let ((n2 (/ (- (cos alpha) (* (cos gamma) (cos beta))) (sin gamma)))
            (mat (geom:make-matrix t)))
        (geom:set-from-vector mat (vector _cell_length_a
                                          0.0
                                          0.0
                                          0.0
                                          (* _cell_length_b (cos gamma))
                                          (* _cell_length_b (sin gamma))
                                          0.0
                                          0.0
                                          (* _cell_length_c (cos beta))
                                          (* _cell_length_c n2)
                                          (* _cell_length_c (sqrt (- (expt (sin beta) 2) (* n2 n2))))
                                          0.0
                                          0.0
                                          0.0
                                          0.0
                                          1.0))
        (let ((tmat (geom:transposed3x3 mat))
              (res (chem:make-residue sname))
              (mol (chem:make-molecule sname))
              (agg (chem:make-aggregate sname)))
          (chem:add-matter agg mol)
          (chem:add-matter mol res)
          (loop for atom-name across (gethash "_atom_site_label" dict)
                for element across (gethash "_atom_site_type_symbol" dict)
                for sfx across (gethash "_atom_site_fract_x" dict)
                for sfy across (gethash "_atom_site_fract_y" dict)
                for sfz across (gethash "_atom_site_fract_z" dict)
                for fx = (convert-number sfx)
                for fy = (convert-number sfy)
                for fz = (convert-number sfz)
                for fvec = (geom:vec fx fy fz)
                for vec = (geom:m*v tmat fvec)
                for atm = (chem:make-atom (intern atom-name :keyword) (intern element :keyword))
                do (chem:set-position atm vec)
                do (chem:add-matter res atm)
                )
          (calculate-bonds agg)
          agg)))))
                                        
(defun process-cif-block (cif-block)
  (declare (optimize (debug 3)))
  (let ((dict (make-hash-table :test 'equal)))
    (loop named lines
          for cur = (lines cif-block) then (cdr cur)
          for entry = (car cur)
          when (null cur)
            do (return-from lines dict)
          do (let ((line (line entry)))
               (cond
                 ((string= line "loop_" :start1 0 :end1 5)
                  (setf cur (process-loop (cdr cur) dict)))
                 (t
                  (let* ((first-space (or (position-if (lambda (c) (or (char= c #\space) (char= c #\tab))) line)
                                          (length line)))
                         (non-space (if (< first-space (length line))
                                        (position-if (lambda (c) (not (or (char= c #\space) (char= c #\tab)))) line :start first-space)
                                        nil))
                         (key (subseq line 0 first-space))
                         (value (cond
                                  ((and (null non-space) (not (char= (char (line (car (cdr cur))) 0) #\;)))
                                   (let ((data-line (line (car (cdr cur)))))
                                     (setf cur (cdr cur))
                                     data-line))
                                  ((and non-space (< non-space (length line)))
                                   (subseq line non-space (length line)))
                                  (t
                                   (let ((tvalue nil))
                                     (multiple-value-setq (cur tvalue)
                                       (parse-multi-line-value (cdr (cdr cur))))
                                     tvalue)))))
                    (setf (gethash key dict) value))))))
    (make-instance 'data-block
                   :name (name cif-block)
                   :dict dict)))
  
(defun load-cif-blocks (filename &key debug)
  (with-open-file (fin filename :direction :input)
    (let ((liness (loop named lines
                        for line = (read-line fin nil :eof)
                        for count from 1
                        when (eq line :eof)
                          do (return-from lines result)
                        when (and (> (length line) 0) (char/= (char line 0) #\#)) ; skip comments and blank lines
                          collect (make-instance 'cif-line :lineno count :line line)
                            into result)))
      (let ((cif-blocks (loop named parse-blocks
                              for cur = liness then (cdr cur)
                              for entry = (car cur)
                              when (null cur)
                                do (return-from parse-blocks result)
                              collect (let ((line (line entry))
                                            (lineno (lineno entry)))
                                        (cond
                                          ((string= line "data_" :start1 0 :end1 5)
                                           (let* ((name (subseq line 6))
                                                  (cur-block (make-instance 'cif-block :name name)))
                                             (when debug (format t "!~a~%" line))
                                             (setf cur (parse-cif-block (cdr cur) cur-block))
                                             cur-block))
                                          (t (error "Illegal block start at line ~a: ~a" lineno line))))
                                into result)))
        (loop for one in cif-blocks
              collect (process-cif-block one))
        ))))

(defun load-cif (filename &key debug)
  "Load a CIF file containing a single data block"
  (unless (probe-file filename)
    (error "Could not find ~a" filename))
  (let ((cifs (load-cif-blocks filename :debug debug)))
    (unless (= (length cifs) 1) (error "load-cif expects one block - there are ~a - use load-cifs" (length cifs)))
    (values (generate-aggregate (first cifs)) (first cifs))))

(defun load-cifs (filename &key debug)
  "Load a CIF file containing a multiple data blocks and return them as aggregates in a list"
  (unless (probe-file filename)
    (error "Could not find ~a" filename))
  (let ((cifs (load-cif-blocks filename :debug debug)))
    (values (mapcar #'generate-aggregate cifs) cifs)))
