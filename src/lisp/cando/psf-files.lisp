(in-package :cando)

(export 'load-psf-pdb)

(defmacro log (fmt &rest fmtargs)
  `(format t ,fmt ,@fmtargs)) 

(defun read-psf-header (fin)
  (let ((raw-line (loop for rl = (read-line fin nil :eof)
                     do (when (eq rl :eof) (return :eof))
                       do (when (> (length rl) 0) (return rl)))))
    (format t "read header: ~s~%" raw-line)
    (if (eq raw-line :eof)
        nil
        (let* ((*package* (find-package :keyword))
               (line (subseq raw-line 0 (position #\: raw-line)))
               (bang-pos (position #\! line))
               (value-part (subseq line 0 bang-pos))
               (var-part (subseq line (1+ bang-pos)))
               (values (with-input-from-string (sin value-part)
                         (loop for val = (read sin nil :eof)
                            until (eq val :eof)
                            collect val)))
               (vars (with-input-from-string (sin var-part)
                       (loop for var = (read sin nil :eof)
                          until (eq var :eof)
                          collect var))))
          (mapcar (lambda (var val) (cons var val)) vars values)))))

(defun skip-to-blank-line (fin)
  (loop for line = (read-line fin nil :eof)
     until (eq line :eof)
     until (equal line "")))

(defun read-psf-title (psfin header)
  (let ((num (cdr (assoc :ntitle header))))
    (with-output-to-string (sout)
      (dotimes (n num)
        (print (read-line psfin) sout)))))

(defun read-psf-atoms (psfin header)
  (let ((num (cdr (assoc :natom header)))
        (atoms (make-hash-table :test #'eql))
        (residues (make-hash-table :test #'eql)))
    (flet ((process-one-atom (line)
             (with-input-from-string (sin line)
               (let ((atomid (read sin))
                     (segment-name (read sin))
                     (residue-id (read sin))
                     (residue-name (read sin))
                     (atom-name (read sin))
                     (atom-type (read sin))
                     (charge (read sin))
                     (mass (read sin)))
                 (let ((residue (gethash residue-id residues))
                       (atom (chem:make-atom atom-name
                                             (chem:element-from-atom-name-string-case-insensitive (string atom-name)))))
                   (chem:set-type atom atom-type)
                   (unless residue
                     (setq residue (chem:make-residue residue-name))
                     (setf (gethash residue-id residues) residue))
                   (chem:add-atom residue atom)
                   (setf (gethash atomid atoms) atom))))))
      (format t "About to read ~a lines~%" num)
      (dotimes (ln num)
        (let ((line (read-line psfin)))
          (handler-case (process-one-atom line)
            (error (e) (error "Could not process psf atom entry#~a line <~a> as an atom~%" ln line)))))
      (values atoms residues))))

(defun read-psf-bonds (psfin header atoms)
  (let ((num (cdr (assoc :nbond header))))
    (dotimes (bond num)
      (let ((id1 (read psfin))
            (id2 (read psfin)))
        (let ((atom1 (gethash id1 atoms))
              (atom2 (gethash id2 atoms)))
          (chem:bond-to atom1 atom2 :single-bond))))))

(defun read-psf-molnt (psfin header atoms residues)
  (format t "Reading MOLNT section~%")
  (flet ((create-molecule-hash-table (num)
           (let ((mht (make-hash-table :test #'eql)))
             (loop for mi from 1 to num
                do (let ((mol (chem:make-molecule)))
                     (setf (gethash mi mht) (chem:make-molecule))))
             mht))
         (create-atoms-to-residues ()
           (let ((a2r (make-hash-table :test #'eq)))
             (maphash (lambda (ri res)
                        (chem:map-atoms 'vector (lambda (a)
                                                  (setf (gethash a a2r) res))
                                        res))
                      residues)
             a2r)))
    (let* ((num-molecules (cdr (assoc :molnt header)))
           (num-atoms (hash-table-count atoms))
           (molecules (create-molecule-hash-table num-molecules))
           (atoms-to-residues (create-atoms-to-residues))
           (residues-to-molecules (make-hash-table :test #'eq)))
      (format t "About to read ~a MOLNT entries~%" num-atoms)
      (dotimes (i num-atoms)
        (let* ((atomid (1+ i))
               (molidx (read psfin))
               (residue (gethash (gethash atomid atoms) atoms-to-residues))
               (molecule (gethash molidx molecules)))
          (unless (gethash residue residues-to-molecules)
            (setf (gethash residue residues-to-molecules) (gethash molidx molecules))
            (chem:add-matter molecule residue))))
      molecules)))

(defun build-aggregate (molecules)
  "* Arguments
- atoms :: hash-table of ids to atoms.
- residues :: hash-table of residue ids to residues.
* Description build an Aggregate from the atoms and residues"
  (let ((aggregate (chem:make-aggregate)))
    (maphash (lambda (molid mol)
               (chem:add-molecule aggregate mol))
             molecules)
    aggregate))
  
(defparameter *atoms* nil)
(defparameter *residues* nil)
      
(defun load-psf (psfin)
  (block top
    (let ((first-line (read-line psfin))
          atoms residues)
      (read-line psfin) ;; skip blank
      (let (header
            atom-header
            nbond-header
            last-header
            atoms residues
            molnt-header
            molecules)
        (setq header (read-psf-header psfin))
        (read-psf-title psfin header)
        (setq atom-header (read-psf-header psfin))
        (multiple-value-setq (atoms residues)
          (read-psf-atoms psfin atom-header))
        (setq nbond-header (read-psf-header psfin))
        (read-psf-bonds psfin nbond-header atoms)
        (skip-to-blank-line psfin)      ; NTHETA
        (skip-to-blank-line psfin)
        (skip-to-blank-line psfin)
        (skip-to-blank-line psfin)
        (skip-to-blank-line psfin)
        (skip-to-blank-line psfin)
        (skip-to-blank-line psfin)      ; NACC
        (skip-to-blank-line psfin)      ; NNB
        (skip-to-blank-line psfin)      ; weird bank of zeros
        (skip-to-blank-line psfin)      ; NGRP NST2
        (setq molnt-header (read-psf-header psfin))
        (setf molecules (read-psf-molnt psfin molnt-header atoms residues))
        (setq last-header (read-psf-header psfin))
        (format t "Most recent header: ~a~%" last-header)
        (values atoms (build-aggregate molecules))))))


(defun load-pdb-for-atoms (pdbin atoms)
  (let* ((*package* (find-package :keyword)))
    (loop for line = (read-line pdbin nil :eof)
     until (eq line :eof)
       do (with-input-from-string (sin line)
            (let ((head (read sin)))
              (case head
                (:ATOM
                 (let ((atomid (read sin))
                       (dummy (progn
                                (read sin) ; atom name
                                (read sin) ; residue name
                                (read sin) ; residueid
                                nil))
                       (pos (geom:ov3 (read sin) (read sin) (read sin))))
                   (chem:set-position (gethash atomid atoms) pos)))
                (otherwise #|nothing|#)))))))
                  
(defun load-psf-pdb (psf-filename pdb-filename)
  (with-open-file (psfin psf-filename)
    (multiple-value-bind (atoms aggregate)
        (load-psf psfin)
      (with-open-file (pdbin pdb-filename)
        (load-pdb-for-atoms pdbin atoms))
      aggregate)))


