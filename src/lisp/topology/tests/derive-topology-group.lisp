(defpackage #:topology-derivation-tests
  (:use #:cl)
  (:export #:run-tests))

(in-package #:topology-derivation-tests)

(defun atom-graph (topology)
  (loop for atom across (topology:constitution-atoms (topology:constitution topology))
        collect (list (topology:atom-name atom) (topology:element atom)
                      (topology:index atom)
                      (loop for bond in (topology:bonds atom)
                            collect (list (topology:to-atom-index bond)
                                          (topology:order bond))))))

(defun expect-error (function)
  (handler-case (progn (funcall function) nil)
    (error () t)))

(defun stereo-data (topology)
  (mapcar (lambda (configuration)
            (list (topology:atom-name configuration)
                  (topology::configuration configuration)))
          (topology::stereoconfigurations (topology:stereoisomer topology))))

(defun run-tests ()
  "Exercise derivation in a private group registry with synthetic chemistry."
  (let ((topology:*topology-groups* (make-hash-table :test #'eq)))
    (topology:define-topology (dgt-chiral dgt-side)
        ((CM C :plugs (:-side) :stereochemistry-type :both)
         (- HM) (- (F :F)) - (C1 C :plugs (:+tail)) (- H11) (- H12))
      :dihedrals ((((:-side :+tail) :phi-1) :phi)
                  (:CM :chi0) (:C1 :chi1))
      :cluster-dihedrals (:chi1 :phi :chi0)
      :rotamer-limits (:grid-search ((C1 (-180.0 0.0))))
      :rings (:start (:+tail (:CM :C1)))
      :restraints ((:dihedral-restraints (CM C1 H11 H12 0.0 1.0)))
      :code original-code)
    (topology:define-topology (dgt-simple dgt-side)
        ((CM C :plugs (:-side) :adjust topology:internal-planar-adjustment)
         (- H1) (- H2) - H3)
      :dihedrals ((:CM :chi0)))
    (let* ((source-names (copy-list (gethash 'dgt-side topology:*topology-groups*)))
           (sources (mapcar #'chem:find-topology source-names))
           (source-groups
             (loop for key being the hash-keys of topology:*topology-groups*
                     using (hash-value members)
                   collect (cons key (copy-list members))))
           (source-data
             (mapcar (lambda (source)
                       (copy-tree
                        (getf (topology:property-list source) :%topology-definition)))
                     sources))
           (name-function
             (lambda (name)
               (intern (concatenate 'string "E" (symbol-name name))
                       (symbol-package name))))
           (copies
             (topology:derive-topology-group
              'dgt-eside :from 'dgt-side :name-transform name-function
              :plug-name-transform
              (lambda (name)
                (case name (:-side :-eside) (:+tail :+etail) (otherwise name))))))
      (assert (= 3 (length sources) (length copies)
                 (length (gethash 'dgt-eside topology:*topology-groups*))))
      (dolist (group source-groups)
        (assert (equal (cdr group) (gethash (car group) topology:*topology-groups*))))
      (loop for source in sources
            for definition in source-data
            for copy = (chem:find-topology (funcall name-function (topology:name source)))
            do (assert (eq source (chem:find-topology (topology:name source))))
               (assert (equal definition
                              (getf (topology:property-list source) :%topology-definition)))
               (assert (equal (atom-graph source) (atom-graph copy)))
               (assert (not (eq (topology:constitution source) (topology:constitution copy))))
               (assert (not (eq (topology:joint-template source) (topology:joint-template copy))))
               (assert (eq :-side (topology:name (topology:find-in-plug source))))
               (assert (eq :-eside (topology:name (topology:find-in-plug copy))))
               (assert (eq :-eside
                           (topology:name (topology::in-plug (topology:joint-template copy)))))
               (assert (equal (stereo-data source) (stereo-data copy))))
      (let* ((source (first (gethash 'dgt-chiral topology:*topology-groups*)))
             (copy (chem:find-topology (funcall name-function source)))
             (properties (topology:residue-properties (topology:constitution copy)))
             (phi (find :phi (getf properties :dihedrals) :key #'topology:name)))
        (assert (equal '(:-eside :+etail) (topology::plug-path phi)))
        (assert (equal '(:chi1 :phi :chi0) (getf properties :cluster-dihedrals)))
        (assert (eq 'original-code (getf properties :code)))
        (assert (topology:has-plug-named copy :+etail))
        (assert (not (topology:has-plug-named copy :+tail)))
        (assert (= 1 (length (topology::restraints copy))))
        (assert (= 1 (length (topology:rotamer-limits copy))))
        (assert (equal '(:start (:+etail (:CM :C1)))
                       (getf (getf (getf (topology:property-list copy)
                                        :%topology-definition)
                                  :options)
                             :rings))))
      ;; Refresh existing variants, remove only selected annotations, and keep
      ;; the explicit cluster order rather than treating it as a set.
      (let ((refreshed
              (topology:derive-topology-group
               'dgt-eside :from 'dgt-side :name-prefix "E"
               :rename-plugs '((:-side . :-eside) (:+tail . :+etail))
               :omit-dihedrals '(:phi))))
        (assert (= 3 (length refreshed)
                   (length (gethash 'dgt-eside topology:*topology-groups*))))
        (let* ((name (first (gethash 'dgt-chiral topology:*topology-groups*)))
               (copy (chem:find-topology (funcall name-function name)))
               (properties (topology:residue-properties (topology:constitution copy))))
          (assert (not (find :phi (getf properties :dihedrals) :key #'topology:name)))
          (assert (equal '(:chi1 :chi0) (getf properties :cluster-dihedrals)))))
      ;; A subset group must not publish the other stereoisomer.
      (let ((selected (first (gethash 'dgt-chiral topology:*topology-groups*))))
        (setf (gethash 'dgt-subset topology:*topology-groups*) (list selected))
        (assert (= 1 (length (topology:derive-topology-group
                              'dgt-derived-subset :from 'dgt-subset
                              :name-prefix "SUB-" :rename-plugs '((:-side . :-eside))))))
        (assert (= 1 (length (gethash 'dgt-derived-subset topology:*topology-groups*)))))
      (let* ((selected (first (gethash 'dgt-chiral topology:*topology-groups*)))
             (copy (first
                    (topology:derive-topology-group
                     'dgt-suffix :from 'dgt-subset
                     :name-transform
                     (lambda (name)
                       (intern (concatenate 'string (symbol-name name) "-SUFFIX")
                               (symbol-package name)))
                     :rename-plugs '((:-side . :-eside))))))
        (assert (string= (symbol-name (topology:name copy))
                         (concatenate 'string (symbol-name selected) "-SUFFIX")))
        (assert (eq (topology:name copy)
                    (topology:name (topology:stereoisomer copy)))))
      ;; Failure after another family was prepared must publish nothing.
      (topology:define-topology cdgt-simple ((C C) (- H1) (- H2) (- H3) - H4))
      (assert (expect-error
               (lambda ()
                 (topology:derive-topology-group
                  'dgt-collision :from 'dgt-side :name-prefix "C"))))
      (assert (not (gethash 'dgt-collision topology:*topology-groups*)))
      (dolist (name (gethash 'dgt-chiral topology:*topology-groups*))
        (assert (not (chem:find-topology
                      (intern (concatenate 'string "C" (symbol-name name))
                              (symbol-package name))
                      nil))))
      (assert (expect-error
               (lambda ()
                 (topology:derive-topology-group
                  'dgt-bad :from 'dgt-side :name-prefix "BAD-"
                  :name-transform #'identity))))
      (assert (expect-error
               (lambda ()
                 (topology:derive-topology-group
                  'dgt-bad :from 'dgt-side :name-transform (constantly nil)))))
      (assert (expect-error
               (lambda ()
                 (topology:derive-topology-group
                  'dgt-bad :from 'dgt-side :name-prefix "BAD-"
                  :plug-name-transform (lambda (name)
                                         (if (eq name :-side) :+eside name))))))
      (assert (not (gethash 'dgt-bad topology:*topology-groups*))))
    ;; Explicit omission of all dihedrals must survive further derivation.
    (let* ((copy (first
                  (topology:derive-topology-group
                   'dgt-no-dihedrals :from 'dgt-simple
                   :name-prefix "NONE-" :omit-dihedrals '(:chi0))))
           (copy-of-copy (first
                          (topology:derive-topology-group
                           'dgt-no-dihedrals-again :from 'dgt-no-dihedrals
                           :name-prefix "AGAIN-"))))
      (assert (null (getf (topology:residue-properties (topology:constitution copy))
                         :dihedrals)))
      (assert (null (getf (topology:residue-properties (topology:constitution copy-of-copy))
                         :dihedrals))))
    (topology:define-topology (dgt-ion dgt-ion-source)
        ((N N :plugs (:-side) :charge 1) (- H1) (- H2) - H3)
      :residue-properties (:residue-net-charge 1))
    (let* ((copy (first (topology:derive-topology-group
                        'dgt-ion-target :from 'dgt-ion-source :name-prefix "E"
                        :rename-plugs '((:-side . :-eside)))))
           (constitution (topology:constitution copy)))
      (assert (= 1 (getf (topology:residue-properties constitution) :residue-net-charge)))
      (assert (= 1 (getf (topology:properties
                          (topology:constitution-atom-named constitution :N))
                         :charge))))
    ;; Indexed two-bond incoming plugs keep both attachment atoms.
    (topology:define-topology (dgt-dual dgt-dual-source)
        ((C C :plugs (:-side.0)) (= O)
         - CA (- HA)
         - CB (- HB1) (- HB2)
         - CG (- HG1) (- HG2)
         - CD (- HD1) (- HD2)
         - (N N :plugs (:-side.1)) - CA))
    (let* ((copy (first (topology:derive-topology-group
                        'dgt-dual-target :from 'dgt-dual-source
                        :name-prefix "E" :rename-plugs '((:-side . :-eside)))))
           (plug (topology:find-in-plug copy)))
      (assert (eq :-eside (topology:name plug)))
      (assert (equal '(:C :N)
                     (map 'list #'topology:atom-name (topology:plug-bonds plug))))))
  (format t "~&Topology derivation regression checks passed.~%")
  t)
