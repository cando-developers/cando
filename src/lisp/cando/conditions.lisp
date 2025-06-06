;;;
;;;    File: conditions.lisp
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
(in-package :chem)

(defun linear-atoms-reporter (condition stream)
  (let ((indexes (chem::indexes condition))
        (coordinates (coordinates condition)))
    (let ((positions (loop for index in indexes
                           for xpos = (elt coordinates index)
                           for ypos = (elt coordinates (+ 1 index))
                           for zpos = (elt coordinates (+ 2 index))
                           collect (geom:vec xpos ypos zpos)))) 
      (format stream "A ~a occurred~%Atoms ~a~%Coordinates: ~a~%" (class-name (class-of condition)) (chem:atoms condition) positions)
      (when (slot-boundp condition 'simd-index) (format stream ":simd-index ~d" (simd-index condition))))))
  

(define-condition chem:minimizer-error (error)
  ((message :initform "There was a minimizer-error" :initarg :message :reader message)
   (minimizer :initarg :minimizer :reader minimizer)
   (coordinates :initarg :coordinates :reader coordinates))
  (:report (lambda (condition stream)
             (format stream "~a ~a" (class-name (class-of condition)) (message condition)))))

(define-condition chem::linear-atoms-error (chem:minimizer-error)
  ((atoms :initarg :atoms :reader chem:atoms)
   (coordinates :initarg :coordinates :reader coordinates)
   (simd-index :initarg :simd-index :reader simd-index)
   (indexes :initarg :indexes :reader chem::indexes))
  (:report linear-atoms-reporter))

(define-condition chem:linear-angle-error (linear-atoms-error) ())
(define-condition chem:linear-dihedral-error (linear-atoms-error) ())
(define-condition chem:linear-improper-restraint-error (linear-atoms-error) ())


(define-condition chem:minimizer-exceeded-max-steps (chem:minimizer-error)
  ((number-of-steps :initarg :number-of-steps :accessor number-of-steps))
  (:report (lambda (condition stream)
             (format stream "~a :steps ~a"
                     (class-name (class-of condition))
                     (number-of-steps condition)))))

(define-condition chem:minimizer-exceeded-sd-max-steps (chem:minimizer-exceeded-max-steps) ())
(define-condition chem:minimizer-exceeded-cg-max-steps (chem:minimizer-exceeded-max-steps) ())
(define-condition chem:minimizer-exceeded-tn-max-steps (chem:minimizer-exceeded-max-steps) ())

(defun make-minimizer-exceeded-max-steps (minimizer coordinates number-of-steps)
  (make-condition 'chem:minimizer-exceeded-max-steps
                  :minimizer minimizer
                  :coordinates coordinates
                  :number-of-steps number-of-steps))

(define-condition chem:minimizer-stuck (chem:minimizer-error)
  ((local-steps :initarg :local-steps :reader local-steps))
  (:report (lambda (condition stream)
             (format stream "~a ~a :local-steps ~a" (class-name (class-of condition)) (message condition) (local-steps condition)))))

(defun make-minimizer-stuck (minimizer coordinates number-of-steps)
  (make-condition 'chem:minimizer-stuck
                  :minimizer minimizer
                  :local-steps number-of-steps
                  :coordinates coordinates))



;;; ------------------------------------------------------------
;;;
;;; Parameter warnings
;;;

(defvar *parameter-warnings*)

(define-condition estimated-stretch-term (error)
  ((type1 :initarg :type1 :accessor type1)
   (type2 :initarg :type2 :accessor type2)
   (distance :initarg :distance :accessor angle-rad)
   (force :initarg :force :accessor force)))

(defun warn-estimated-stretch-term (t1 t2 distance force)
  (when (string> (string t1) (string t2))
    (let ((temp t1))
      (setf t1 t2
            t2 temp)))
  (let ((term (make-condition 'estimated-stretch-term
                              :type1 t1
                              :type2 t2
                              :distance distance
                              :force force)))
    (push term *parameter-warnings*)))

(defmethod print-object ((condition estimated-stretch-term) stream)
  (print-unreadable-object (condition stream :type t)
    (format stream "types(~a ~a) distance(~3,2f) k(~5,2f)"
            (type1 condition)
            (type2 condition)
            (distance condition)
            (force condition))))

(define-condition estimated-angle-term (error)
  ((type1 :initarg :type1 :accessor type1)
   (type2 :initarg :type2 :accessor type2)
   (type3 :initarg :type3 :accessor type3)
   (angle-rad :initarg :angle-rad :accessor angle-rad)
   (force :initarg :force :accessor force)))

(defun warn-estimated-angle-term (t1 t2 t3 angle-rad force)
  (when (string> (string t1) (string t3))
    (let ((temp t1))
      (setf t1 t3
            t3 temp)))
  (let ((term (make-condition 'estimated-angle-term
                              :type1 t1
                              :type2 t2
                              :type3 t3
                              :angle-rad angle-rad
                              :force force)))
    (push term *parameter-warnings*)))

(defmethod print-object ((condition estimated-angle-term) stream)
  (print-unreadable-object (condition stream :type t)
    (format stream "types(~a ~a ~a) angle(~5,2f deg) k(~5,2f)"
            (type1 condition)
            (type2 condition)
            (type3 condition)
            (/ (angle-rad condition) 0.0174533)
            (force condition))))

(defun compare-estimated-terms (term1 term2)
;;; Assume only angles for now
  (if (eq (class-of term1)
          (class-of term2))
      (etypecase term1
        (estimated-angle-term
         (or (string< (string (type1 term1)) (string (type1 term2)))
             (string< (string (type2 term1)) (string (type2 term2)))
             (string< (string (type3 term1)) (string (type3 term2)))))
        )
      (string< (string (class-name (class-of term1)))
               (string (class-name (class-of term2))))))

(defun estimated-terms-eq (term1 term2)
  ;;; Assume only angles for now
  (if (eq (class-of term1)
          (class-of term2))
      (etypecase term1
        (estimated-angle-term
         (and (eq (type1 term1) (type1 term2))
              (eq (type2 term1) (type2 term2))
              (eq (type3 term1) (type3 term2))))
        )
      nil))

(defun report-parameter-warnings ()
  (when *parameter-warnings*
    (let ((unique (remove-duplicates *parameter-warnings* :test #'estimated-terms-eq)))
      (loop for warning in unique
            do (format t "~a~%" warning)))))

(define-condition chem:missing-atom-type-error (error)
  ((message :initform "There was a missing atom type" :initarg :message :reader message)
   (atm :initarg :atm :reader atm)
   (atom-types-hash-table :initarg :atom-types-hash-table :reader atom-types-hash-table))
  (:report (lambda (condition stream)
             (format stream "~a ~a" (atm condition) (atom-types-hash-table condition)))))
