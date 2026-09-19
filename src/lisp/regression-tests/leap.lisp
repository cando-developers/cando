
(in-package #:clasp-tests)

(leap:setup-default-paths)
#+(or)
(test-type load-frcmod
           (leap:load-amber-params "frcmod.ff99SB")
           chem:force-field)
#+(or)
(test-type load-frcmod-scee-scnb
           (leap:load-amber-params "frcmod.fb15")
           chem:force-field)

(test-type load-frcmod-pho
           (leap:load-amber-params (namestring (translate-logical-pathname "sys:extensions;cando;src;lisp;regression-tests;data;frcmod.pho")))
           chem:force-field)

;;; Numbering gaps are segment boundaries, not molecule boundaries.
(defun pdb-gap--residue (number &key code (chain :A) (context :main) hetatm)
  (make-instance 'leap.pdb::pdb-residue
                 :name :GLY :residue-name :GLY :original-residue-name :GLY
                 :res-seq number :i-code code :chain-id chain
                 :context context :hetatmp hetatm))

(test-true pdb-gap-numbering
  (and (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 179) (pdb-gap--residue 188))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 179) (pdb-gap--residue 180)))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 179)
                                         (pdb-gap--residue 188 :chain :B)))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 179)
                                         (pdb-gap--residue 188 :hetatm t)))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 179) (pdb-gap--residue 1)))))

(test-true pdb-gap-insertion-codes
  (and (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 10)
                                         (pdb-gap--residue 10 :code :A)))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 10 :code :A)
                                         (pdb-gap--residue 10 :code :B)))
       (not (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 10 :code :B)
                                         (pdb-gap--residue 11)))
       (leap.pdb::pdb-numbering-gap-p (pdb-gap--residue 10 :code :A)
                                     (pdb-gap--residue 10 :code :C))))

(test-true pdb-gap-marks-and-molecule-grouping
  (let ((a (pdb-gap--residue 179))
        (b (pdb-gap--residue 188))
        (c (pdb-gap--residue 200)))
    (leap.pdb::mark-pdb-numbering-gap a b)
    (leap.pdb::mark-pdb-numbering-gap b c)
    (and (eq (leap.pdb::context a) :end-gap)
         (equal (leap.pdb::context b) '(:gap-begin :end-gap))
         (eq (leap.pdb::context c) :gap-begin)
         (null (leap.pdb::gap-before a))
         (null (leap.pdb::gap-after c))
         (eq (leap.pdb::gap-after a) (leap.pdb::gap-before b))
         (eq (leap.pdb::gap-after b) (leap.pdb::gap-before c))
         (= 1 (hash-table-count
                (leap.pdb::molecules-from-sequences (list (list a b c))))))))

(test-true pdb-gap-preserves-terminal-marks
  (let ((a (pdb-gap--residue 1 :context :head))
        (b (pdb-gap--residue 5)))
    (leap.pdb::mark-pdb-numbering-gap a b)
    (leap.pdb::add-pdb-context-mark b :tail)
    (and (equal (leap.pdb::context a) '(:head :end-gap))
         (equal (leap.pdb::context b) '(:gap-begin :tail)))))
