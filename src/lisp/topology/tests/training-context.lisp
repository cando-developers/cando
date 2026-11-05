(in-package #:cl-user)

(defun test-training-context-parser ()
  (test-legacy-training-context-parser)
  (let* ((expression
           '(training-context :name "branched"
             :context (side :-side (body :as current)
                       (:+dkp (body :as following))
                       :-dkp (body :as preceding))
             :training-caps ((preceding :-dkp pros) (preceding :+side bnz))
             :shape-kind (sidechain :paired-with "pro->BODY"
                          :shape-key-sources ((current :first-ring) (preceding :second-ring)))))
         (context (topology::parse-training-context expression))
         (edges (topology::training-context-connections context))
         (current (topology::training-context-connection-target (first edges))))
    (assert (= 4 (length (topology::training-context-nodes context))))
    (assert (eq :sidechain (topology::training-context-shape-kind context)))
    (assert (equal "pro->BODY" (topology::training-context-paired-with context)))
    (assert (= 3 (length edges)))
    (assert (eq current (topology::training-context-connection-source (second edges))))
    (assert (eq current (topology::training-context-connection-source (third edges))))
    (assert (eq current (topology::training-context-selector-node
                        (first (topology::training-context-selectors context)))))
    (assert (= 2 (length (topology::training-context-caps context))))
    ;; Preserve shared node identity, not just node contents, on round-trip.
    (let* ((*readtable* cando.serialize::*cando-readtable*)
           (copy (cando.serialize:from-string (cando.serialize:as-string context)))
           (nodes (topology::training-context-nodes copy))
           (edges (topology::training-context-connections copy)))
      (assert (typep copy 'topology::training-context))
      (assert (equal expression (topology::training-context-expression copy)))
      (assert (equal "branched" (topology::training-context-name copy)))
      (assert (eq :sidechain (topology::training-context-shape-kind copy)))
      (assert (equal "pro->BODY" (topology::training-context-paired-with copy)))
      (assert (eq (first nodes) (topology::training-context-focus copy)))
      (assert (eq (second nodes) (topology::training-context-connection-source (third edges))))
      (assert (eq (fourth nodes) (topology::training-context-cap-node
                                  (first (topology::training-context-caps copy)))))
      (assert (eq (second nodes) (topology::training-context-selector-node
                                 (first (topology::training-context-selectors copy)))))))
  (dolist (path '( (side) ((body :as root))
                  (side :-side body (:+dkp body) :-dkp pro)))
    (topology::parse-training-context `(training-context :name "valid" :context ,path)))
  (dolist (options
            '((:name "x" :context (side :-side))
              (:name "x" :context (side (body :as wrong)))
              (:name "x" :context (side :-side (body :label x)))
              (:name "x" :context ((body :as x) :-dkp (body :as x)))
              (:name "x" :context (side) :training-caps ((missing :+side bnz)))
              (:name "x" :context ((body :as x) :+side side)
               :training-caps ((x :+side bnz)))
              (:name "x" :context ((body :as x))
               :shape-kind (sidechain :shape-key-sources ((x :ring :ring))))
              (:name "x" :context ((body :as x))
               :shape-kind (sidechain :shape-key-sources ((missing :ring))))
              (:name "x" :context (side) :unknown t)
              (:name "x" :name "y" :context (side))
              (:name "x" :context (side . bad))
              (:name "x" :context (side) :training-caps ((x :+dkp)))
              (:name "x" :context (side) :shape-kind (sidechain :paired-with 42))
              (:name "x" :context (side) :training-caps)
              (:context (side))))
    (assert (handler-case
                (progn (topology::parse-training-context (cons 'training-context options)) nil)
              (topology::training-context-parse-error () t))))
  (dolist (kind '((backbone) (:backbone) (backbone :paired-with-sidechains nil)
                  (backbone :paired-with-sidechains t)))
    (let ((context (topology::parse-training-context
                    `(training-context :name "backbone" :context (body) :shape-kind ,kind))))
      (assert (eq :backbone (topology::training-context-shape-kind context)))
      (assert (eq (getf (rest kind) :paired-with-sidechains)
                  (topology::training-context-paired-with-sidechains context)))
      (assert (null (topology::training-context-selectors context)))))
  (dolist (kind '((sidechain) (:sidechain)))
    (assert (eq :sidechain
                (topology::training-context-shape-kind
                 (topology::parse-training-context
                  `(training-context :name "side" :context (side) :shape-kind ,kind))))))
  (dolist (kind '(nil backbone (unknown) (backbone :paired-with "x")
                  (backbone :shape-key-sources nil)
                  (backbone :paired-with-sidechains 1)
                  (sidechain :paired-with-sidechains t)
                  (sidechain :paired-with "") (sidechain :paired-with :absent)
                  (sidechain :paired-with)
                  (sidechain :paired-with "a" :paired-with "b")
                  (sidechain :unknown t) (sidechain . bad)
                  (sidechain :shape-key-sources 42)))
    (assert (handler-case
                (progn (topology::parse-training-context
                        `(training-context :name "bad" :context (body) :shape-kind ,kind)) nil)
              (topology::training-context-parse-error () t))))
  (dolist (option '(:paired-with :shape-key-sources))
    (assert (handler-case
                (progn (topology::parse-training-context
                        `(training-context :name "bad" :context (body) ,option nil)) nil)
              (topology::training-context-parse-error () t))))
  ;; Foldamer-specific options must not be accepted by the default hook.
  (assert (handler-case
              (progn (topology::parse-training-context
                      '(training-context :name "unknown option" :context (body)
                        :shape-kind (backbone :domain-option t))) nil)
            (topology::training-context-parse-error () t)))
  (format t "Training-context parser tests passed.~%")
  t)

(defun test-legacy-training-context-parser ()
  (let* ((context (topology::parse-legacy-training-context
                   '(:description "source on cap"
                     :pattern (side :-side naa :+amide (:cap caa) :+side (:cap ala))
                     :shape-kind (:sidechain :shape-key-sources
                                  (:path (:-side :phi) (:+amide (:psi-1 :psi)))))))
         (selector (second (topology::training-context-selectors context))))
    (assert (= 2 (length (topology::training-context-nodes context))))
    (assert (eq (topology::training-context-selector-node selector)
                (first (topology::training-context-caps context)))))
  (let* ((expression
           '(pattern :description "legacy sidechain"
             :pattern (side :-side body (:-dkp pro) :+dkp (:cap cgly))
             :shape-kind (:sidechain :paired-with "PRO"
                          :shape-key-sources (:path (:-side :first-ring)
                                                   (:-dkp (:ring :second-ring))))))
         (before (copy-tree expression))
         (context (topology::parse-legacy-training-context expression))
         (nodes (topology::training-context-nodes context))
         (selectors (topology::training-context-selectors context)))
    (assert (equal before expression))
    (assert (= 3 (length nodes)))
    (assert (null (topology::training-context-node-label (first nodes))))
    (assert (eq 'topology::foo (topology::training-context-node-label (second nodes))))
    (assert (eq 'topology::bar (topology::training-context-node-label (third nodes))))
    (assert (eq (third nodes) (topology::training-context-selector-node (second selectors))))
    (assert (eq :ring (topology::training-context-selector-cache-field (second selectors))))
    (assert (eq :second-ring (topology::training-context-selector-output-name (second selectors)))))
  (let* ((context (topology::parse-legacy-training-context
                   '(:description "cap chain" :pattern (body :+amide (:cap caa) :+side (:cap ala))
                     :shape-kind (:backbone :backbone-shape-key-selector protein-backbone))))
         (caps (topology::training-context-caps context)))
    (assert (= 1 (length (topology::training-context-nodes context))))
    (assert (eq (first caps) (topology::training-context-cap-node (second caps))))
    (assert (eq 'protein-backbone
                (topology::shape-kind-backbone-shape-key-selector
                 (topology::training-context-shape-kind-properties context)))))
  (dolist (expression
            '((:description "bad" :pattern (side :-side))
              (:description "bad" :pattern (side :+a (:ring-cap x)))
              (:description "bad" :pattern (side) :shape-kind
               (:sidechain :shape-key-sources (:path (:-side :ring))))
              (:description "bad" :pattern (side :-side body) :shape-kind
               (:sidechain :shape-key-sources (:path (:-side :ring :ring))))))
    (assert (handler-case (progn (topology::parse-legacy-training-context expression) nil)
              (topology::training-context-parse-error () t))))
  t)
