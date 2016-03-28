;;;
;;;    File: psf-files.lisp
;;;

;; Open Source License
;; Copyright (c) 2016, Christian E. Schafmeister
;; Permission is hereby granted, free of charge, to any person obtaining a copy
;; of this software and associated documentation files (the "Software"), to deal
;; in the Software without restriction, including without limitation the rights
;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;; copies of the Software, and to permit persons to whom the Software is
;; furnished to do so, subject to the following conditions:
;; The above copyright notice and this permission notice shall be included in
;; all copies or substantial portions of the Software.
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
;; THE SOFTWARE.
;;  
;; This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.

;; -^-
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
        (residues (make-hash-table :test #'equal))
        residue)
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
                 (let ((unique-residue-id (format nil "~a-~a" segment-name residue-id))
                       (atom (chem:make-atom atom-name
                                             (chem:element-from-atom-name-string-case-insensitive
                                              (string atom-name)))))
                   (chem:set-type atom atom-type)
                   (chem:set-charge atom charge)
                   (unless (setf residue (gethash unique-residue-id residues))
                     (setq residue (chem:make-residue residue-name))
                     (setf (gethash unique-residue-id residues) residue))
                   (chem:add-atom residue atom)
                   (setf (gethash atomid atoms) atom))))))
      (format t "About to read ~a lines~%" num)
      (dotimes (ln num)
        (let ((line (read-line psfin)))
          (handler-case (process-one-atom line)
            (error (e) (error "Could not process psf atom entry#~a line <~a> as an atom - error: ~a~%" ln line e)))))
      (values atoms residues))))

(defun read-psf-bonds (psfin header atoms)
  (let ((num (cdr (assoc :nbond header))))
    (dotimes (bond num)
      (let ((id1 (read psfin))
            (id2 (read psfin)))
        (let ((atom1 (gethash id1 atoms))
              (atom2 (gethash id2 atoms)))
          (chem:bond-to atom1 atom2 :single-bond))))))


(defun create-atoms-to-residues (residues)
  "* Arguments
- residues :: A hashtable of residue ids to residues.
* Description
Return a hash table that maps atoms to residues that contain them."
  (let ((a2r (make-hash-table :test #'eq)))
    (maphash (lambda (ri res)
               (chem:map-atoms 'vector (lambda (a)
                                         (setf (gethash a a2r) res))
                               res))
             residues)
    a2r))


(defun walk-atom-spanning-tree (atom molecule
                                atoms-seen
                                residues-to-molecules
                                atoms-to-residues)
  "* Arguments
- atom :: An atom.
- molecule :: A molecule.
- atoms-seen :: A hash-table of atoms to T.
- residues-to-molecules :: A hash-table of residues to molecules.
- atoms-to-residues :: A hash-table of atoms to the residues that contain them.
* Description
Walk the spanning tree of the ATOM and add the ATOM to ATOMS-SEEN.
For each atom walked look up the residue in ATOMS-TO-RESIDUES
and if the residue is not in RESIDUES-TO-MOLECULES add it."
  (let ((spanning-loop (chem:make-spanning-loop atom)))
    (loop for a = (chem:next spanning-loop (lambda (a b) t))
       while a
       do (unless (gethash a atoms-seen)
            (setf (gethash a atoms-seen) t)
            (let ((residue (gethash a atoms-to-residues)))
              (or residue (error "The residue could not be found for atom: ~a" a))
              (unless (gethash residue residues-to-molecules)
                (setf (gethash residue residues-to-molecules) molecule)
                (chem:add-matter molecule residue)))))))


(defun build-molecules-from-atom-connectivity (atoms residues)
  "* Arguments
- atoms :: Hash-table of atom ids to atoms.
- residues :: Hash-table of residue ids to residues.
* Description
Return a hash table of molecule-ids to molecules that contains molecules identified by walking the
spanning trees of the atoms and identifying the residues that belong to
those atoms."
  (let ((atoms-seen (make-hash-table :test #'eq))
        (atoms-to-residues (create-atoms-to-residues residues))
        (residues-to-molecules (make-hash-table :test #'eq))
        (molecules (make-hash-table :test #'eql))
        (molecule-id 0))
    (maphash (lambda (id atom)
               (unless (gethash atom atoms-seen)
                 ;; Start a new molecule
                 (let ((mol (chem:make-molecule)))
                   (walk-atom-spanning-tree atom mol
                                            atoms-seen
                                            residues-to-molecules
                                            atoms-to-residues)
                   (setf (gethash (incf molecule-id) molecules) mol))))
             atoms)
    molecules))
             



(defun read-psf-molnt (psfin header atoms residues)
  (format t "Reading MOLNT section~%")
  (flet ((create-molecule-hash-table (num)
           (let ((mht (make-hash-table :test #'eql)))
             (loop for mi from 1 to num
                do (let ((mol (chem:make-molecule)))
                     (setf (gethash mi mht) (chem:make-molecule))))
             mht)))
    (let* ((num-molecules (cdr (assoc :molnt header)))
           (num-atoms (hash-table-count atoms))
           (molecules (create-molecule-hash-table num-molecules))
           (atoms-to-residues (create-atoms-to-residues residues))
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
            molnt-header
            molecules)
        (setq header (read-psf-header psfin))
        (read-psf-title psfin header)
        (setq atom-header (read-psf-header psfin))
        (multiple-value-setq (atoms residues)
          (read-psf-atoms psfin atom-header))
        (setq nbond-header (read-psf-header psfin))
        (read-psf-bonds psfin nbond-header atoms)
        (read-psf-header psfin)         ;NTHETA
        (format t "Expected NTHETA~%")
        (skip-to-blank-line psfin)
        (read-psf-header psfin)         ; NPHI
        (format t "Expected NPHI~%")
        (skip-to-blank-line psfin)
        (read-psf-header psfin)         ; NIMPHI
        (format t "Expected NIMPHI~%")
        (skip-to-blank-line psfin)
        (read-psf-header psfin)         ; NDON
        (format t "Expected NDON~%")
        (skip-to-blank-line psfin)
        (read-psf-header psfin)         ; NACC
        (format t "Expected NACC~%")
        (skip-to-blank-line psfin)
        (read-psf-header psfin)         ; NNB
        (format t "Expected NNB~%")
        (skip-to-blank-line psfin)
        (format t "Skipping bank of zeros~%")
        (skip-to-blank-line psfin)     ; Skip mysterious bank of zeros
        (read-psf-header psfin)        ; NGRP NST2
        (format t "Expected NGRP NST2~%")
        (skip-to-blank-line psfin)
        (let ((header (read-psf-header psfin)))
          (when (eq (car molnt-header) :NCRTERM)
            (setq molnt-header header)
            (format t "Expected MOLNT~%")
            (setq *atoms* atoms)
            (setq *residues* residues)
            (format t "molnt-header = ~a~%" molnt-header)))
        (setf molecules (if molnt-header
                            (read-psf-molnt psfin molnt-header atoms residues)
                            (build-molecules-from-atom-connectivity atoms residues)))
        (values atoms (build-aggregate molecules))))))


(defun load-pdb-for-atoms (pdbin atoms)
  (let* ((*package* (find-package :keyword)))
    (loop for line = (read-line pdbin nil :eof)
     until (eq line :eof)
       do (let ((head (read-from-string line nil nil :start 0 :end (min (length line) 4))))
            (case head
              (:ATOM
               (let ((atomid (read-from-string line t nil :start 4 :end 11))
                     (pos (geom:vec (read-from-string line t nil :start 30 :end 38)
                                    (read-from-string line t nil :start 38 :end 46)
                                    (read-from-string line t nil :start 46 :end 54))))
                 (chem:set-position (gethash atomid atoms) pos)))
              (otherwise #|nothing|#))))))
                  
(defun load-psf-pdb (psf-filename pdb-filename)
  (with-open-file (psfin psf-filename)
    (multiple-value-bind (atoms aggregate)
        (load-psf psfin)
      (with-open-file (pdbin pdb-filename)
        (load-pdb-for-atoms pdbin atoms))
      aggregate)))
