;;;
;;;    File: ring-builder.lisp
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
(in-package :linked-rings)

(defclass linkable-ring ()
  ((atoms :initarg :atoms :accessor atoms)
   (join :initarg :join :initform nil :accessor join)))

(defun rings-share-atoms (fring1 fring2)
  (dolist (r1 (atoms fring1))
    (dolist (r2 (atoms fring2))
      (if (eq r1 r2) (return-from rings-share-atoms t))))
  nil)

;; Link the rings together that share atoms
(defun link-rings (fused-rings)
  (loop for (a . rest) on fused-rings
     do (loop for b in rest
           do (when (rings-share-atoms a b)
                (push b (join a))
                (push a (join b))))))

(defun linkable-rings (rings)
  (mapcar (lambda (r)
            (make-instance 'linkable-ring
                           :atoms r))
          rings))

(defun join-rings (rings)
  (let* ((linkable-rings (linkable-rings rings)))
    (link-rings linkable-rings)
    linkable-rings))

(defun build-ring-hash-table (linkable-rings)
  (let ((ht (make-hash-table)))
    (dolist (lring linkable-rings)
      (setf (gethash lring ht) lring))
    ht))

(defun gather-isolated-rings (rings)
  (let (ring-chains)
    (maphash (lambda (k v)
               (when (null (join k))
                 (remhash k rings)
                 (push (cons k nil) ring-chains)))
             rings)
    ring-chains))

(defun find-start (the-rings)
             (maphash (lambda (k v)
                        (when (= (length (join k)) 1)
                          (return-from find-start k)))
                      the-rings))
(defun find-next (one prev the-rings)
             (let* ((joins (join one)))
               (find-if (lambda (a) (not (eq a prev))) joins)))

(defun gather-one-ring-chain (rings)
  (let (ring-chain
        start cur prev)
    (do* ((next (find-start rings) (find-next cur prev rings)))
         ((null next) ring-chain)
      (push next ring-chain)
      (remhash next rings)
      (setq prev cur)
      (setq cur next))))

(defun gather-ring-chains (rings)
  (let (ring-chains)
    (when (> (hash-table-count rings) 0)
      (push (gather-one-ring-chain rings) ring-chains))
    ring-chains))

(defun ring-chains (linkable-rings)
  (let* ((unlinked-rings (build-ring-hash-table linkable-rings))
         (isolated-rings (gather-isolated-rings unlinked-rings))
         (ring-chains (gather-ring-chains unlinked-rings)))
    (append isolated-rings ring-chains)))

(defun identify-all-ring-chains (mol)
  (let* ((rings (chem:identify-rings mol))
         (linkable-rings (join-rings rings)))
    (ring-chains linkable-rings)))



;;;
;;; Coerce a list of atoms to a set (hash-table) of atoms
;;;
(defun coerce-to-atom-set (arg)
  (cond
    ((and (listp arg) (typep (car arg) 'linkable-ring))
     (let ((new-set (make-hash-table)))
       (mapc (lambda (lr)
               (mapc (lambda (a)
                       (setf (gethash a new-set) a))
                     (atoms lr)))
             arg)
       new-set))
    ((listp arg)
     (let ((new-set (make-hash-table)))
       (mapc (lambda (a) (setf (gethash a new-set) a)) arg)
       new-set))
    ((hash-table-p arg)
     arg)
    (t (error "Add support to coerce ~s to atom-set" arg))))

;;;
;;; Expand a collection of atoms by one bond
;;;
(defun expand-atom-set-by-one-bond (passed-atom-set)
  (let ((new-set (make-hash-table)))
    (maphash (lambda (a v)
               (setf (gethash a new-set) a)
               (mapc (lambda (b) (setf (gethash b new-set) b))
                     (chem:atom/bonded-atoms-as-list a)))
             (coerce-to-atom-set passed-atom-set))
    new-set))

;;;
;;; Build geometry for a subset of atoms in a matter
;;;
(defun build-structure-for-atoms (matter atom-set)
)


(defun expanded-atom-set-from-linked-rings (linked-rings)
  (let ((atom-set (make-hash-table)))
    (loop for lr in linked-rings
         do (loop for atom in (atoms lr)
               do (setf (gethash atom atom-set) atom)))
    (expand-atom-set-by-one-bond atom-set)))

(defun expanded-atom-set-for-ring (linked-ring)
  (let ((atom-set (make-hash-table)))
    (loop for atom in (atoms linked-ring)
          do (setf (gethash atom atom-set) atom))
    (expand-atom-set-by-one-bond atom-set)))


(defun scramble-linked-rings (linked-rings)
  "Scramble the coordinates of a set of linked rings"
  (loop for lr in linked-rings
     for x from 0 to (length linked-rings)
     do (cando:scramble-positions (expanded-atom-set-for-ring lr)
                                  :center (geom:vec 0.0 (* x 5.0) 0.0 )
                                  :width 5.0)))
