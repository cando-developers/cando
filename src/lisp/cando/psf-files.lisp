(in-package :cando)

(defmacro log (fmt &rest fmtargs)
  `(format t ,fmt ,@fmtargs)) 

(defun read-psf-section-header (fin)
  (let ((line (read-line fin nil :eof)))
    (if (eq line :eof)
        nil
        (with-input-from-string (sin line)
          (let ((num (read sin))
                (section (let ((*package* (find-package :keyword)))
                           (read sin))))
            (values num section))))))

(defun read-psf-title (psfin num)
  (prog1 
      (with-output-to-string (sout)
        (dotimes (n num)
          (print (read-line psfin) sout)))
    (read-line psfin)))

(defun read-psf-atoms (psfin num)
  (let ((atoms (make-hash-table :test #'eql))
        (residues (make-hash-table :test #'eql)))
    (dotimes (line num)
      (let ((line (read-line psfin)))
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
                  (atom (chem:make-atom (string atom-name)
                                   (chem:element-from-atom-name-string-case-insensitive (string atom-name)))))
              (chem:set-type atom atom-type)
              (unless residue
                (setq residue (chem:make-residue (string residue-name)))
                (setf (gethash residue-id residues) residue))
              (chem:add-atom residue atom)
              (setf (gethash atom-id atoms) atom))))))
    (values atoms residues)))

(defun load-psf-pdb (psf-fin pdb-fin)
  (let ((first-line (read-line psf-fin))
        atoms residues)
    (read-line psf-fin) ;; skip blank
    (loop
       (multiple-value-bind (num section)
           (read-psf-section-header psf-fin)
         (unless num
           (return-from load-psf-pdb (values atoms residues)))
         (case section
           (:!NTITLE
            (log "Reading ntitle~%")
            (let ((title (read-psf-title psf-fin num)))
              (log "Read: ~a~%" title)))
           (:!NATOM
            (log "Reading atoms~%")
            (multiple-value-setq (atoms residues)
              (read-psf-atoms psf-fin num)))
           (otherwise
            (format t "Add support for ~a section~%")))))))


(with-open-file (psfin #P"POPC36.psf")
  (with-open-file (pdbin #P"POPC36.pdb")
    (load-psf-pdb psfin pdbin)))

(apropos "make-residue")
(untrace read-line)
