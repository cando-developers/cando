;;;
;;;    File: select.lisp
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
(in-package :select)

(defun first-atom-with-name (agg name &key debug)
  (let (atom)
    (chem:map-atoms nil (lambda (a)
                          (when debug
                            (format t "Looking for atom with name ~s current atom: ~s~%" name (chem:get-name a)))
                          (when (eq name (chem:get-name a))
                            (setq atom a)))
                    agg)
    atom))


(defun atom-lookup (matter residue-sequence-number atom-name)
  (chem:map-residues
   nil (lambda (residue)
         (when (= (chem:get-file-sequence-number residue) residue-sequence-number)
           (chem:map-atoms
            nil
            (lambda (a)
              (when (string= (chem:get-name a) atom-name)
                (return-from atom-lookup a)))
            residue)))
   matter)
  (error "Could not locate atom with residue-sequence-number ~s name ~s"
         residue-sequence-number atom-name ))


(defun atoms-with-unique-names (matter name-list)
  (mapcar (lambda (n) (first-atom-with-name matter n))
          name-list))

(defun atoms-with-property (matter prop-name)
  (let (atoms)
    (chem:map-atoms nil
                    (lambda (x)
                      (when (not (eq :eof (getf (chem:get-properties x) prop-name :eof)))
                        (push x atoms)))
                    matter)
    atoms))

(defun atoms-with-property-sorted (matter prop-name &key (order #'string<))
  "Return a list of atoms with the propery sorted by the property"
  (format t "In atoms-with-property-sorted~%")
  (let* ((atoms (atoms-with-property matter prop-name)))
    (format t "atoms-with-property ~a~%" atoms)
    (let* ((atoms-with-props (mapcar (lambda (a) ;associate atom with prop
                                       (cons a (getf (chem:get-properties a) prop-name)))
                                     atoms))
           (sorted (sort atoms-with-props (lambda (x y) ; sort by prop
                                            (let ((px (string (cdr x)))
                                                  (py (string (cdr y))))
                                              (funcall order px py))))))
      (mapcar (lambda (x) (car x)) sorted))))


(defun atoms-inside-sphere (matter pos radius)
  (let (selected)
    (chem:map-atoms nil
                    (lambda (a)
                      (let ((dist (geom:v3-magnitude (geom:sub (chem:get-position a) pos))))
                        (when (< dist radius)
                          (push a selected))))
                    matter)
    selected))

(defun atoms-outside-sphere (matter pos radius)
  (let (selected)
    (chem:map-atoms nil
                    (lambda (a)
                      (let ((dist (geom:v3-magnitude (geom:sub (chem:get-position a) pos))))
                        (when (> dist radius)
                          (push a selected))))
                    matter)
    selected))
