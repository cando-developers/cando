(progn
(ql:quickload :myc)

(load (compile-file "~/common-lisp/tlt/design/src/design/blueprint.lisp"))
;; AFTER blueprint.lisp - blueprint-energies.lisp is in the DESIGN package and calls TOPOLOGY::
;; accessors that blueprint.lisp defines.  Neither file is in design.asd.
(load (compile-file "~/common-lisp/tlt/design/src/design/blueprint-energies.lisp"))
;; AFTER both - blueprint-debug.lisp is in DESIGN and calls TOPOLOGY:: accessors from
;; blueprint.lisp (MRKINDEX-INFO, BLUEPRINT-SLOT-RESIDUE-SET, BLUEPRINT-ATOM-I3).
(load (compile-file "~/common-lisp/tlt/design/src/design/blueprint-debug.lisp"))
)
;;; WHICH PACKAGE THIS FILE READS IN
;;;
;;; With :TESTS on it reads in CLASP-TESTS, so TEST-TRUE (framework.lisp:117) resolves and the
;;; sections below become assertions the regression runner can score.  Without it, PATH-DESIGN, as
;;; before.
;;;
;;; The consequence is that EVERY external symbol must be package-qualified - a bare PATH-DESIGN
;;; symbol resolves under :TESTS-off and is an unbound-variable/undefined-function under :TESTS-on,
;;; which is the worst kind of difference because only one of the two configurations shows it.
;;; The file's OWN symbols (*BLUEPRINT*, *PERSONA*, the BM- helpers) stay bare and simply intern
;;; into whichever package is current; they are consistent within a load either way.
#+tests(in-package #:clasp-tests)
#-tests(in-package #:path-design)

(progn

  ;; One call site per check, whether or not tests are on.  With :TESTS the form's value is
  ;; asserted; without it the form just runs.  Either way it runs ONCE - these checks take minutes
  ;; and print the diagnostics that make a failure readable, so running them twice, or running them
  ;; only under :TESTS, would both be wrong.
  (defmacro bench-check (name form)
    "Run FORM for its output; with :TESTS on, also assert it returned true."
    #+tests `(test-true ,name ,form)
    #-tests `(progn ,form))

  (defmacro with-tee ((path) &body body)
    "Run BODY with all standard output going to both the REPL and PATH."
    (let ((f (gensym)) (b (gensym)))
      `(with-open-file (,f ,path :direction :output
                                 :if-exists :supersede :if-does-not-exist :create)
         (let* ((,b (make-broadcast-stream *standard-output* ,f))
                (*standard-output* ,b)
                (*trace-output*    ,b)
                (*error-output*    ,b))
           (unwind-protect (progn ,@body)
             (finish-output ,f))))))

  (defmacro defpose (name form)
    "Evaluate FORM (which returns a pose), bind the special var NAME to it, and print
   wall-time + total-energy + the pose's rep-weight.  Re-evaluating re-runs FORM from its
   (fixed) input pose, so you can run a mover repeatedly and compare.  (Energy is omitted
   when the pose's total-energy is NIL — e.g. right after rescan-mover or freeze-and-rescan
   clears it.  The w= is always printed: it is the weight the pose is SCORED at, and
   total-energy is only comparable between poses that share it.)"
    `(progn
       (defvar ,name)
       (time
        (let ((t0 (get-internal-real-time)))
          (setf ,name ,form)
          (format t "defpose ~&~32a ~7,2fs   w=~,3f~@[   E=~,3f~]~%"
                  ',name
                  (/ (float (- (get-internal-real-time) t0) 1d0) internal-time-units-per-second)
                  (design:rep-weight ,name)
                  (design:total-energy ,name))
          (finish-output t)))
       ',name))
  )

;;; ---- optional: build *pose-start* / *protocol* from a starting-pose file + a task ----
;;; (path-design::build-initial-pose-and-protocol needs the same TASK the refine task uses; skip this
;;;  if you already have *pose-start* / *protocol* bound another way.)
(progn
  (defvar *pose-start*)
  (defvar *protocol*)
  ;; Declared here so BM-BUILD-GRAPH compiles cleanly; the values live in the
  ;; knob block below, which must be evaluated before the first BM-BUILD-GRAPH call.
  (defvar *bm-backbone-count*)
  (defvar *bm-backbone-seed*)
  (defun bm-setup (starting-poses-task refine-task output-file)
    "Run the upstream STARTING-POSES-TASK to produce the starting-pose file REFINE-TASK
     consumes, load it, and build *pose-start* / *protocol*."
    (ensure-directories-exist (task:file-pathname output-file))
    (format t "~&bm-setup: running starting-poses-task ...~%") (finish-output t)
    (task:execute starting-poses-task) ; produce the starting-pose file(s)
    (let* ((starting-pose-file (first (task:inputs refine-task)))
           (starting-pose (cando.serialize:load-cando
                           (task:file-pathname starting-pose-file)))) ; load-cando wants a pathname
      (multiple-value-bind (pose protocol)
          (path-design::build-initial-pose-and-protocol refine-task starting-pose output-file)
        (defvar *pose-start*)
        (setf *pose-start* pose
              *protocol*    protocol)
        (format t "~&*pose-start* ready.  total-energy=~,3f~%"
                (design:total-energy *pose-start*))
        (values))))

;;; ---- build the graph (don't run it) and pull the task pair for one combo ----
  (defun bm-build-graph (&key (backbone-count *bm-backbone-count*)
                           (backbone-seed  *bm-backbone-seed*))
    "Build the path-design graph in memory.  Requires *spiro-space* + the other params
     already set (your bench-run.lisp setfs).  A single EXPAND-GRAPH materializes the
     starting-poses-tasks, all refine tasks, and the sort task.

     BACKBONE-COUNT/BACKBONE-SEED limit the search to a stable random subset of the
     backbones: select-backbone-indices (graph.lisp:36) runs a Fisher-Yates over a
     random-state derived from the seed, so every node given the same seed picks the
     same subset.  BACKBONE-COUNT NIL searches all backbones."
    (unless path-design:*spiro-space*
      (error "path-design:*spiro-space* is not set"))
    (when (and backbone-count path-design:*limit-backbone-indices*)
      (error "path-design:*limit-backbone-indices* is ~s, which overrides the seeded ~
              subset (graph.lisp:245).  (setf path-design:*limit-backbone-indices* nil) ~
              to use :backbone-count ~d instead."
             path-design:*limit-backbone-indices* backbone-count))
    (let* ((target (make-instance 'myc:myc-target
                                  :spiro-space path-design:*spiro-space*
                                  :limit-backbone-indices
                                  (when backbone-count
                                    (list :random-backbone-count backbone-count
                                          :random-backbone-seed  backbone-seed))))
           (graph  (make-instance 'path-design::path-design-graph :target target)))
      (task:expand-graph graph)
      graph))

  (defun bm-backbone-indices (graph)
    "Sorted backbone indices actually present in GRAPH — i.e. what the seeded subset
     picked.  task:tasks is a fill-pointer vector, so MAP rather than LOOP..IN."
    (let (result)
      (map nil (lambda (tk)
                 (when (typep tk 'path-design::starting-poses-task)
                   (pushnew (path-design::backbone-index tk) result)))
           (task:tasks graph))
      (sort result #'<)))

  (defun bm-find-tasks (graph &key backbone (path 1) (offset-pm 450) (angle-deg 240))
    "Return (values starting-poses-task refine-task output-file) for one combo.
     output-file is the refine task's OWN output (refined-path-pose.cando), the same
     thing the refine task's execute uses (execute.lisp:732).

     BACKBONE defaults to the first backbone in the graph rather than 0, because a
     seeded random subset almost never contains index 0."
    (let* ((available (bm-backbone-indices graph))
           (backbone  (or backbone (first available))))
      (unless backbone
        (error "graph has no starting-poses-tasks at all"))
      (unless (member backbone available)
        (error "backbone ~d is not in this graph.  Available (~d): ~s"
               backbone (length available) available))
      (format t "~&bm-find-tasks: backbone=~d of ~d in graph ~s~%"
              backbone (length available) available)
      (let ((refine (find-if (lambda (tk)
                               (and (typep tk 'path-design::path-design-refine-mover-task)
                                    (eql (path-design::backbone-index tk) backbone)
                                    (eql (path-design::path-index tk) path)
                                    (let ((in (first (task:inputs tk))))
                                      (and (eql (path-design::offset-pm in) offset-pm)
                                           (eql (path-design::angle-deg in) angle-deg)))))
                             (task:tasks graph)))
            (sp     (find-if (lambda (tk)
                               (and (typep tk 'path-design::starting-poses-task)
                                    (eql (path-design::backbone-index tk) backbone)
                                    (eql (path-design::path-index tk) path)))
                             (task:tasks graph))))
        (unless refine (error "no refine task for backbone=~d path=~d offset=~d angle=~d"
                              backbone path offset-pm angle-deg))
        (unless sp (error "no starting-poses-task for backbone=~d path=~d" backbone path))
        (values sp refine (first (task:outputs refine))))))
  )


;;; =====================================================================
;;; Setup per-pass knobs (defaults mirror the small bench; tweak between passes)
;;; =====================================================================

(progn
  (defvar *protocol* nil "Shared protocol — bound once at setup, not changed by any mover.")
  (defparameter *bm-rep-weight*            0.1d0)
  (defparameter *bm-rigid-body-iterations* 2000)
  (defparameter *bm-design-iterations*     20000)
  (defparameter *bm-design-warm-up*        10000)
  (defparameter *bm-verbose*               nil)
  (defparameter *bm-backbone-count*        40
    "How many backbones to search.  NIL means all — far too many for *spiro8-space*.")
  (defparameter *bm-backbone-seed*         1234
    "Seed for the backbone subset.  Same seed => same subset on every node.")
  (setf path-design:*spiro-space* path-design:*spiro7-space*)
  ;; The explicit-index list takes priority over the seeded subset (graph.lisp:245),
  ;; and bench-run.lisp sets it — clear it so *bm-backbone-count* is what applies.
  (setf path-design:*limit-backbone-indices* nil))



(unwind-protect
     (progn
       (dribble "~/tmp/dump.log")

;;; =====================================================================
;;; PHASE 1 — rigid-body / hremc   (one rep-weight; repeat block per weight)
;;; =====================================================================

       (setf *bm-rep-weight* 0.1d0) ; set the rep-weight for this pass

       ;; 0) Get the starting point
       (ext:with-flame-profile (:path "~/public_html/a000-setup-pose.svg")
         (defparameter *graph* (bm-build-graph))
         ;; :backbone omitted — defaults to the first index the seeded subset actually picked.
         ;; (bm-backbone-indices *graph*) lists them if you want a different one.
         (multiple-value-call #'bm-setup
           (bm-find-tasks *graph* :path 1 :offset-pm 450 :angle-deg 240))
         (format t "*pose-start* and *protocol* are set~%"))

       ;; a1) compute-energies: build cached-energies for this rep-weight
       (ext:with-flame-profile (:path "~/public_html/a001-pose-p1-ce0.svg")
         (format t "Generating pose-p1-ce0~%")
         (defpose *pose-p1-ce0*
             (design:apply-mover (make-instance 'design:compute-energies-mover)
                                 *pose-start* *protocol*)))

       ;; The rotamers-database argument is gone: compute-blueprint-layout now derives it
       ;; per oligomer-shape from that shape's own foldamer, so the spiro loci stop being
       ;; looked up in the amber-protein database and stop coming back with one slot each.
       (format t "############################################ blueprint being built ~%")
       (defparameter *blueprint*
         (topology::make-blueprint (design:assembler path-design::*pose-p1-ce0*)))
       (format t "############################################ blueprint build is done~%")

       (let ((loci        (topology::loci *blueprint*))
             (total-slots (topology::total-slots *blueprint*)))
         ;; Pass 2 has not run, so these must still be unbound.  If any of them answers, the
         ;; blueprint quietly inherited the source assembler's materialization and every
         ;; atresidue-keyed lookup would miss without complaint.
         (loop for slot in '(topology::aggregate topology::energy-function
                             topology::ataggregate topology::joint-tree topology::adjustments)
               do (format t "~&  ~24a ~:[unbound (correct - pass 2)~;BOUND - unexpected!~]~%"
                          slot (slot-boundp *blueprint* slot)))
         (let ((sidechain (count :sidechain loci :key #'topology::kind))
               (backbone (count :backbone loci :key #'topology::kind))
               (fixed (count :fixed loci :key #'topology::kind)))
           (format t "~&loci = ~d   (~d sidechain, ~d backbone, ~d fixed)   total-slots = ~d~%"
                   (length loci) sidechain backbone fixed total-slots)
           ;; 16 bytes per pair entry - the interleaved fa_rest/fa_rep pair.
           (format t "~&pair table = ~:d entries = ~,2f GB~%~%"
                   (floor (* total-slots (1+ total-slots)) 2)
                   (/ (* 16d0 (floor (* total-slots (1+ total-slots)) 2)) 1d9)))
         ;; Show the first few of EACH tier - a run where the fixed tier is empty and one
         ;; where it swallowed everything both look fine if you only print the head.
         (flet ((show (label test limit)
                  (format t "~&~a:~%" label)
                  (loop with shown = 0
                        for l across loci
                        when (funcall test l)
                          do (when (>= shown limit) (loop-finish))
                             (incf shown)
                             (let ((*print-pretty* nil))
                               (format t "  locus ~3d  base ~8d  extent ~7d   ~a~%"
                                       (topology::locus l) (topology::base l)
                                       (topology::extent l)
                                       ;; NIL matcher means nothing matched the SPACE, which
                                       ;; silently yields one slot - it must not be mistaken
                                       ;; for a monomer that genuinely has no rotamers.
                                       (loop for m across (topology::monomers l)
                                             collect (list (topology::monomer-name m)
                                                           (topology::rotamer-slot-count m)
                                                           (if (topology::matcher m)
                                                               :matched
                                                               :NO-MATCHER))))))))
           (show "sidechain loci (own mrkindexes)"
                 (lambda (l) (eq (topology::kind l) :sidechain)) 15)
           (show "backbone loci (extent 0 - no mrkindex)"
                 (lambda (l) (eq (topology::kind l) :backbone)) 10)
           (show "fixed loci" (lambda (l) (eq (topology::kind l) :fixed)) 5))
         (let ((unmatched (loop for l across loci
                                append (loop for m across (topology::monomers l)
                                             unless (topology::matcher m)
                                               collect (topology::monomer-name m)))))
           (format t "~&~%monomers with NO matcher: ~d~@[  ~s~]~%"
                   (length unmatched) (remove-duplicates unmatched)))
         ;; Split by tier.  Tallying every monomer summed to 526 against a total-slots of 470 -
         ;; the difference being the five backbone monomers, whose rotamer-slot-count bounds a
         ;; rotamer INDEX rather than counting mrkindexes.  Two numbers that look like they
         ;; should agree and cannot.  The sidechain tally MUST equal total-slots, so it says so.
         (flet ((tally (test label expected)
                  (let ((counts (make-hash-table :test #'equal))
                        (sum 0))
                    (loop for l across loci
                          when (funcall test l)
                            do (loop for m across (topology::monomers l)
                                     for n = (topology::rotamer-slot-count m)
                                     do (incf (gethash (topology::monomer-name m) counts 0) n)
                                        (incf sum n)))
                    (format t "~&~%~a:~%" label)
                    (loop for k being the hash-keys of counts using (hash-value v)
                          do (format t "   ~20a ~8d~%" k v))
                    (format t "   ~20a ~8d~@[  ~a~]~%" "TOTAL" sum
                            (when expected
                              (if (= sum expected)
                                  (format nil "= total-slots")
                                  (format nil "*** MISMATCH - total-slots is ~d ***"
                                          expected)))))))
           (tally (lambda (l) (eq (topology::kind l) :sidechain))
                  "mrkindex slots by monomer-name" total-slots)
           (tally (lambda (l) (eq (topology::kind l) :backbone))
                  "backbone rotamer capacity by monomer-name (NOT mrkindex slots)" nil)))

       ;; ---- persona: load rotamers into the slots, and check the bound ----
       ;; No materialization needed - slot-to-rotamer holds rotamer objects from the database,
       ;; so this runs while the blueprint's assembler-base slots are still unbound.  OVERFLOWS
       ;; must be empty: a non-empty list means rotamer-slot-count was computed too small and
       ;; rotamers are being silently dropped.
       (format t "~&~%!!!! PERSONA: LOAD ROTAMERS INTO THE SLOTS, AND CHECK THE BOUND~%")
       (defparameter *persona* (topology::make-persona *blueprint* :fill nil))
       (multiple-value-bind (loaded overflows)
           (topology::fill-persona *persona* :verbose t)
         (declare (ignorable loaded))
         (format t "~&~s~%" *persona*)
         (format t "~&contexts resolved: ~d of ~d monomers~%"
                 (count-if-not #'null (topology::monomer-contexts *persona*))
                 (topology::monomer-count *blueprint*))
         (format t "~&shape-keys set:    ~d loci~%"
                 (count-if-not #'null (topology::shape-keys *persona*)))
         (format t "~&~%active/capacity per monomer:~%")
         (topology::report-persona *persona*)
         (when overflows
           (format t "~&*** THE BOUND IS WRONG - ~d monomer~:p overflowed ***~%"
                   (length overflows))))

       ;; ---- drive the selection through every combination the space allows ----
       ;; One fill proves the bound covers one backbone state; it claims to cover all of them.
       ;; This mutates monomer-indexes and restores them.
       (format t "~&~%!!!! DRIVE THE SELECTION THROUGH EVERY COMBINATION THE SPACE ALLOWS~%")
       (format t "~&~%")
       (defparameter *bound-failures*
         (topology::check-bound-over-selections *persona*))

       ;; ---- pass 2 stage 1: base structure, sidechain loci declined ----
       ;;
       ;; :ENSURE-TRAINED T.  Train the SMIRNOFF cache on real trainers FIRST, so the blueprint
       ;; hits the cache and real SMIRNOFF never touches the fanned-out graph.
       ;;
       ;; NIL was the other option and is not practical: it makes the blueprint miss, which runs
       ;; real SMIRNOFF over a graph where CA carries 58-79 substituents - including RING
       ;; PERCEPTION over that hub, which is both meaningless and enormous.  It would also harvest
       ;; the result, poisoning the cache for every later build in the session.
       ;;
       ;; The cost of skipping it: the four reference energies below were measured under exactly
       ;; that unpoisoned-cache condition, so they are no longer an apples-to-apples comparison.
       ;; Expect them to MOVE, and treat that as a correction rather than a regression - the
       ;; two-pass parameterization split is therefore NOT independently verified by this run.
       (format t "~&~%!!!! PASS 2 STAGE 1: BASE STRUCTURE, SIDECHAIN LOCI DECLINED~%")
       (defparameter *ensure-trained* t)
       (topology::materialize-blueprint *blueprint*
                                        :ensure-trained *ensure-trained*
                                        :verbose t)
       (multiple-value-bind (residues atresidues)
           (topology::blueprint-materialized-counts *blueprint*)
         ;; One residue per NON-sidechain locus - build-monomer-p declines the sidechain ones -
         ;; plus one per SLOT from the fan-out.  A sidechain locus contributes nothing directly;
         ;; its residues ARE the slots.
         (let* ((loci (topology::loci *blueprint*))
                (sidechain-loci (count :sidechain loci :key #'topology::kind))
                (base-residues (- (length loci) sidechain-loci))
                (slot-residues (topology::total-slots *blueprint*))
                (expected (+ base-residues slot-residues)))
           (format t "~&~%materialize-blueprint:~%")
           (format t "   residues   ~4d~%" residues)
           (format t "   atresidues ~4d~@[  ~a~]~%" atresidues
                   (unless (= residues atresidues)
                     "*** MISMATCH - every residue should have an atresidue ***"))
           ;; Which wrong number it is says which stage failed, so name them rather than
           ;; printing one generic marker.
           (format t "   expected   ~4d  (~d base + ~d slots)~@[  ~a~]~%"
                   expected base-residues slot-residues
                   (cond ((= residues expected) nil)
                         ((= residues (length loci))
                          "*** BUILD-MONOMER-P DID NOT FIRE ***")
                         ((= residues base-residues)
                          "*** FAN-OUT DID NOT RUN ***")
                         (t "*** UNEXPECTED COUNT ***")))))

       ;; every slot's root joint should have the backbone attachment joint as its parent,
       ;; and each fanned backbone joint should now carry 3 + N children
       (loop for bp-locus across (topology::loci *blueprint*)
             when (topology::owns-slots-p bp-locus)
               do (let ((expected (topology::blueprint-slot-parent-joint *blueprint* bp-locus))
                        (n 0) (wrong 0))
                    (loop for bp-monomer across (topology::monomers bp-locus)
                          do (loop for slot below (topology::rotamer-slot-count bp-monomer)
                                   for scan = (aref (topology::rotamer-scans *blueprint*)
                                                    (topology::blueprint-mrkindex bp-locus bp-monomer slot))
                                   for root = (aref (topology::joints (topology::atresidue scan)) 0)
                                   do (incf n)
                                      (unless (eq (kin:parent root) expected) (incf wrong))))
                    (format t "~&locus ~3d  ~3d slots  ~3d wrong parent  backbone joint children ~d~%"
                            (topology::locus bp-locus) n wrong
                            (length (kin:joint-children expected)))))

       ;; ---- residue/atresidue INDEX PAIRING ----
       ;;
       ;; WALK-ATOMS-JOINTS (assembler.lisp:922-933) pairs residue i with atresidue i by INDEX and
       ;; writes joint positions into atoms on that basis.  Equal COUNTS do not imply equal ORDER,
       ;; and the blueprint appends 470 residues in stage 2a and 470 atresidues in stage 2b as two
       ;; separate passes.  A mismatch writes every coordinate into the wrong atom - which looks
       ;; like a molecule that is scattered or missing rather than like an error.
       (format t "~&~%!!!! RESIDUE/ATRESIDUE INDEX PAIRING~%")
       (loop for molecule-index below (chem:content-size (topology::aggregate *blueprint*))
             for molecule = (chem:content-at (topology::aggregate *blueprint*) molecule-index)
             for atmol = (aref (topology::atmolecules (topology::ataggregate *blueprint*))
                               molecule-index)
             do (let ((bad (loop for i below (min (chem:content-size molecule)
                                                  (length (topology::atresidues atmol)))
                                 unless (eq (chem:content-at molecule i)
                                            (topology::residue
                                             (aref (topology::atresidues atmol) i)))
                                   collect i)))
                  (format t "~&molecule ~d  residues ~4d  atresidues ~4d  mispaired ~4d~@[  ~a~]~%"
                          molecule-index (chem:content-size molecule)
                          (length (topology::atresidues atmol)) (length bad)
                          (when bad
                            (format nil "*** first at index ~d ***" (first bad))))))

       ;; ---- does BLUEPRINT-SLOT-PARENT-JOINT pick the atom we think it does? ----
       ;;
       ;; It finds the attachment atom's position in the topology's CONSTITUTION-ATOMS and then
       ;; takes that index out of (JOINTS BACKBONE-ATRESIDUE).  That is only right if joints are
       ;; indexed parallel to constitution-atoms.  If they are indexed by residue atom order
       ;; instead, the joint belongs to a neighbouring atom - and hanging CB off the carbonyl
       ;; carbon instead of CA puts it one bond length from C, which is what a misplaced fan
       ;; looks like.  The existing "wrong parent" check below cannot see this: it compares
       ;; against the same possibly-wrong answer.
       (format t "~&~%!!!! DOES BLUEPRINT-SLOT-PARENT-JOINT PICK THE ATOM WE THINK IT DOES?~%")
       (loop for bp-locus across (topology::loci *blueprint*)
             when (topology::owns-slots-p bp-locus)
               do (multiple-value-bind (parent-joint index)
                      (topology::blueprint-slot-parent-joint *blueprint* bp-locus)
                    (declare (ignorable parent-joint))
                    (multiple-value-bind (sidechain-name backbone-name)
                        (topology::slot-attachment-atom-names
                         (aref (topology::monomers bp-locus) 0) bp-locus)
                      (declare (ignorable sidechain-name))
                      (let* ((monomer (topology::locus-backbone-attachment bp-locus))
                             (position (gethash monomer (topology::monomer-positions *blueprint*)))
                             (molecule (chem:content-at (topology::aggregate *blueprint*)
                                                        (topology::molecule-index position)))
                             (residue (chem:content-at molecule
                                                       (topology::residue-index position)))
                             (at-index (when (< index (chem:content-size residue))
                                         (chem:get-name (chem:content-at residue index)))))
                        (format t "~&locus ~3d  attach ~a  joint-index ~3d  residue atom there ~a~@[  ~a~]~%"
                                (topology::locus bp-locus) backbone-name index at-index
                                (unless (eq at-index backbone-name)
                                  "*** JOINT INDEX IS NOT THE ATTACHMENT ATOM ***"))))))

       ;; ---- SLOT ROOT JOINT ANCESTRY, the blueprint's own ANCHOR TRACE ----
       ;;
       ;; A BondedJoint places its atom from THREE ancestors: distance to parent, angle to
       ;; grandparent, dihedral to great-grandparent.  The working single-scan path prints
       ;;    ANCHOR TRACE mol1 res62  joint0=CB  p=CA  gp=N  ggp=C
       ;; where ggp is the PREVIOUS residue's carbonyl carbon (its i3x is lower than N's).  A
       ;; rotamer's internals are defined against exactly that chain.
       ;;
       ;; The blueprint hangs each slot off a receptor CA joint, and the receptor's base residues
       ;; are RESIDUE-SHAPEs, which get XYZJOINTS - joints that hold positions outright and have
       ;; no need to be chained N->CA.  If they are not chained, gp and ggp here are not N and
       ;; prev-C, the dihedral means something else, and CB lands somewhere different at every
       ;; residue.  That is exactly the reported symptom.
       ;;
       ;; JOINT ids are (molecule residue atom) - see FIND-ATOM-FOR-JOINT, assembler.lisp:966.
       (format t "~&~%!!!! SLOT ROOT JOINT ANCESTRY, THE BLUEPRINT'S OWN ANCHOR TRACE~%")
       (flet ((joint-atom-name (joint)
                (when joint
                  (destructuring-bind (m r a) (kin:id joint)
                    (let ((mol (chem:content-at (topology::aggregate *blueprint*) m)))
                      (chem:get-name (chem:content-at (chem:content-at mol r) a)))))))
         (loop for bp-locus across (topology::loci *blueprint*)
               when (topology::owns-slots-p bp-locus)
                 do (let* ((bp-monomer (aref (topology::monomers bp-locus) 0))
                           (scan (aref (topology::rotamer-scans *blueprint*)
                                       (topology::blueprint-mrkindex bp-locus bp-monomer 0)))
                           (root (aref (topology::joints (topology::atresidue scan)) 0))
                           (p (and root (kin:parent root)))
                           (gp (and p (kin:parent p)))
                           (ggp (and gp (kin:parent gp))))
                      (format t "~&BLUEPRINT TRACE locus ~3d  joint0=~a  p=~a  gp=~a  ggp=~a~@[  ~a~]~%"
                              (topology::locus bp-locus)
                              (joint-atom-name root) (joint-atom-name p)
                              (joint-atom-name gp) (joint-atom-name ggp)
                              (unless (and p gp ggp)
                                "*** ANCESTRY INCOMPLETE - dihedral reference undefined ***")))))

       ;; ---- SLOT INTERNAL ADJUSTMENTS ARE REGISTERED AND RESOLVABLE ----
       ;;
       ;; Ancestry being right is not enough.  A rotamer does not store CB's dihedral; AMBER-
       ;; PROTEIN's BETA-CARBON adjustment derives it from the sibling C joint as c-dih - 120 deg
       ;; (aminoacids.lisp:84-94).  With no adjustment registered - which is the state when the
       ;; :ADJUST BETA-CARBON clauses are commented out - CB keeps the raw rotamer dihedral and
       ;; lands on the carbonyl carbon.  That was the fan-misplacement bug.
       ;;
       ;; Two ways this silently reverts to broken, so check both:
       ;;   1. no adjustment registered on the slot atresidues at all
       ;;   2. registered but INITIALIZE-ADJUSTMENT could not resolve the sibling joints, because
       ;;      the parent CA is an XYZJOINT whose children were never linked.  BETA-CARBON's
       ;;      INITIALIZE-ADJUSTMENT does (kin:joint/child-with-name jca :C) and stores the result
       ;;      without checking it; a NIL there does not fail until INTERNAL-ADJUST walks it.
       ;;
       ;; The :C probe is BETA-CARBON-SPECIFIC and must be gated on the adjustment class.  Applying
       ;; it to every adjustment blows up on the ligand's spiro loci, whose parent joint has
       ;; children (:CAC :CM :CM ...) and correctly has no :C at all.
       ;;
       ;; Use KIN:JOINT/JOINT-CHILDREN, not KIN:JOINT/CHILD-WITH-NAME: the latter SIGNALS on a miss
       ;; (joint.cc:397) rather than returning NIL, which turns a diagnostic into a crash.
       ;; Collect the slot atresidues FLAT before testing them.  The nested locus/monomer/slot
       ;; walk plus a COND plus a DOLIST puts the innermost form eleven parens deep, which is how
       ;; the first version of this check ended up closing its own LET.
       ;; FIND-SYMBOL signals a PACKAGE-ERROR on a missing package, so gate it on FIND-PACKAGE.
       (format t "~&~%!!!! SLOT INTERNAL ADJUSTMENTS ARE REGISTERED AND RESOLVABLE~%")
       (let* ((bc-class (let ((p (find-package "AMBER-PROTEIN")))
                          (and p (find-symbol "BETA-CARBON" p))))
              (internal-adjustments (topology::internal-adjustments
                                     (topology::adjustments *blueprint*)))
              (slot-atresidues
                (loop for bp-locus across (topology::loci *blueprint*)
                      when (topology::owns-slots-p bp-locus)
                        nconc (loop for bp-monomer across (topology::monomers bp-locus)
                                    nconc (loop for slot below (topology::rotamer-slot-count bp-monomer)
                                                for scan = (aref (topology::rotamer-scans *blueprint*)
                                                                 (topology::blueprint-mrkindex
                                                                  bp-locus bp-monomer slot))
                                                when (and scan (topology::atresidue scan))
                                                  collect (topology::atresidue scan)))))
              (with-adj 0) (without-adj 0) (beta-carbon 0) (unresolvable 0))
         (flet ((child-named-p (joint name)
                  ;; KIN:PARENT SIGNALS when the parent is unbound (joint.h:162) - test with
                  ;; KIN:JOINT/PARENT-BOUND-P first.  KIN:JOINT/CHILD-WITH-NAME signals on a miss
                  ;; (joint.cc:397), so walk KIN:JOINT/JOINT-CHILDREN instead.
                  (and joint
                       (find name (kin:joint/joint-children joint) :key #'kin:joint/name)))
                (parent-of (adjustment)
                  (let ((j (topology::joint adjustment)))
                    (and j (kin:joint/parent-bound-p j) (kin:parent j)))))
           (dolist (atres slot-atresidues)
             (let ((adjs (gethash atres internal-adjustments)))
               (if (null adjs)
                   (incf without-adj)
                   (incf with-adj))
               (dolist (a adjs)
                 ;; The :C probe is BETA-CARBON-specific - gate it on the class.
                 (when (and bc-class (typep a bc-class))
                   (incf beta-carbon)
                   ;; A diagnostic must never abort the run it is diagnosing.  Anything unexpected
                   ;; here counts as unresolvable and gets reported, not signalled.
                   (unless (handler-case (child-named-p (parent-of a) :C)
                             (error (e)
                               (format t "~&  adjustment probe failed on ~s: ~a~%" a e)
                               nil))
                     (incf unresolvable)))))))
         (format t "~&slot adjustments: ~d slots with an internal adjustment, ~d without~
                    ~@[, ~d BETA-CARBON~]~@[, ~d UNRESOLVABLE (no sibling :C under the parent ~
                    joint)~]~%"
                 with-adj without-adj (and (plusp beta-carbon) beta-carbon)
                 (and (plusp unresolvable) unresolvable))
         (when (zerop beta-carbon)
           (format t "~&*** NO BETA-CARBON ADJUSTMENTS ON ANY SLOT - CB will keep the raw rotamer ~
                        dihedral and land on the carbonyl carbon.  Check that :ADJUST BETA-CARBON ~
                        is enabled in aminoacids.lisp. ***~%")))

       ;; ---- LK SOLVATION TYPES MUST BE ASSIGNED ON EVERY ATOM ----
       ;;
       ;; BUILD-RESIDUE-FOR-MONOMER-NAME stamps :LK-SOLVATION-ATOM-TYPE from the stereoisomer atom
       ;; and writes NIL when that slot is unbound (oligomer.lisp:27-31).  The C++ then reads it
       ;; with getPropertyOrDefault(...,nil) (energyRosettaLKSolvation.cc:281) and hands the NIL to
       ;; FIND-LKSOLVATION-TYPE, which errors - but only at EVALUATION time, long after the
       ;; residue was built, and only for whichever atom happens to be reached first.
       ;;
       ;; Grouped by residue NAME because that is the discriminating axis: the blueprint
       ;; materializes EVERY allowed monomer name at every locus, while the piece-assembler only
       ;; ever builds the currently selected one.  A name that is allowed but never selected can
       ;; go through the whole existing pipeline without anyone needing its LK types.  If the
       ;; untyped residues are exactly such names, the gap is in the foldamer's products file
       ;; (ASSIGN-ALL-LK-SOLVATION-ATOM-TYPES, build-rotamers.lisp:13) and not in the blueprint.
       ;;
       ;; If instead EVERY residue is untyped, the assignment never ran for this foldamer at all.
       (format t "~&~%!!!! LK SOLVATION TYPES MUST BE ASSIGNED ON EVERY ATOM~%")
       (let ((by-residue (make-hash-table :test 'equal))
             (total 0) (untyped 0))
         (chem:do-molecules (mol (topology::aggregate *blueprint*))
           (chem:do-residues (res mol)
             (chem:do-atoms (a res)
               (incf total)
               (unless (chem:matter-get-property-or-default a :lk-solvation-atom-type nil)
                 (incf untyped)
                 (incf (gethash (string (chem:get-name res)) by-residue 0))))))
         (format t "~&lk-solvation types: ~d of ~d atoms UNTYPED~@[ across ~d residue name~:p~]~%"
                 untyped total (and (plusp untyped) (hash-table-count by-residue)))
         (when (plusp untyped)
           (let (rows)
             (maphash (lambda (name n) (push (cons name n) rows)) by-residue)
             (loop for (name . n) in (sort rows #'> :key #'cdr)
                   do (format t "~&    ~16a ~6d untyped atom~:p~%" name n)))
           (format t "~&*** Rosetta LK solvation cannot be evaluated for these - ~
                        FIND-LKSOLVATION-TYPE errors on a NIL type. ***~%")))

       ;; ---- fold: internals from the persona, then coordinates, then into the atoms ----
       ;;
       ;; Until this existed the aggregate was chemically complete and geometrically meaningless -
       ;; FILL-PERSONA loads rotamer INDICES, never internal coordinates, so every rotamer of a
       ;; monomer sat superimposed at its template position.
       ;;
       ;; EVERY slot gets internals - all 470.  About 38 of them have no loaded rotamer (the
       ;; reserved refinement slot plus the slack between a monomer's rotamer count and the
       ;; locus-wide bound) and are folded with a stand-in rotamer from the same monomer.  They are
       ;; never selected, so their geometry is arbitrary; what matters is that they FOLD, because
       ;; the fold walks the whole tree and a slot left at zero internals puts its first atom on
       ;; top of its parent and kills KIN:UPDATE-XYZ-COORDS with a zero-length b-c.
       ;;
       ;; A stand-in count far from ~38 means the persona is not loading what it should.
       (format t "~&~%!!!! FOLD: INTERNALS FROM THE PERSONA, THEN COORDINATES, THEN INTO THE ATOMS~%")
       (format t "~&~%")
       ;; The LIGAND-ORIENTATION must be the POSE's, not a fresh (MAKE-ORIENTATION).  An identity
       ;; orientation is not a no-op placeholder - the ligand folds in its own local frame and
       ;; lands at the origin, ~150 A from the receptor, which looks like a missing ligand in a
       ;; viewer even though the molecule is in the aggregate and in the mol2.  It is POSE-rate
       ;; state (a rigid-body ligand move changes it), so the blueprint cannot hold it.
       (multiple-value-bind (coords internals)
           (topology::build-blueprint-coordinates
            *blueprint* *persona* :verbose t
            :ligand-orientation (design:orientation path-design::*pose-p1-ce0*))
         ;; Escapes the binding because FILL-BLUEPRINT-SINGLE-ENERGIES scores THIS fold - the
         ;; geometry the energies describe has to be the one this persona produced.
         (defparameter *blueprint-coords* coords)
         ;; Kept for CHECK-INTERNALS-VECTORS-AGAINST-POSE: the vector the fold actually consumed,
         ;; AFTER ADJUST-ALL-INTERNALS.  Regenerating it later would re-run the fill and silently
         ;; compare a fresh vector rather than the one these coordinates came from.
         (defparameter *blueprint-internals* internals)
         (format t "~&coordinate vector length ~d for ~d atoms~%"
                 (length coords) (chem:number-of-atoms (topology::aggregate *blueprint*)))

         ;; ---- LIGAND MUST LAND ON THE RECEPTOR ----
         ;;
         ;; A wrong LIGAND-ORIENTATION does not fail: TRANSFORM-EXTERNALS-TO-GLOBAL-FRAME applies
         ;; whatever it is handed and reports success.  The only symptom is the ligand sitting at
         ;; the origin while the receptor is near (-11 28 155) - visible by eye as "the ligand is
         ;; missing" and by nothing else.  Centroids need no atom-to-atom mapping, which is why
         ;; DESIGN::MOLECULE-CENTROID exists; reuse it rather than rolling one here.
         ;;
         ;; Molecule 0 is the ligand and molecule 1 the receptor, the same convention
         ;; REPORT-RECEPTOR-PLACEMENT uses (protein-binder.lisp:4147-4148), so this line is
         ;; directly comparable to the placement[...] lines earlier in the log.
         (format t "~&~%!!!! LIGAND MUST LAND ON THE RECEPTOR~%")
         (multiple-value-bind (lx ly lz ln) (design::molecule-centroid *blueprint* coords 0)
           (multiple-value-bind (rx ry rz rn) (design::molecule-centroid *blueprint* coords 1)
             (let ((d (sqrt (+ (expt (- lx rx) 2) (expt (- ly ry) 2) (expt (- lz rz) 2)))))
               (format t "~&placement[blueprint]: ligand n=~d (~,4f ~,4f ~,4f)  ~
                          receptor n=~d (~,4f ~,4f ~,4f)  separation ~,1f A~@[~%~a~]~%"
                       ln lx ly lz rn rx ry rz d
                       (when (> d 40.0)
                         "*** LIGAND IS NOT ON THE RECEPTOR - it folded in its own local frame.  ~
                          Check that :LIGAND-ORIENTATION is the POSE's orientation and not a ~
                          fresh (MAKE-ORIENTATION). ***"))))))

       ;; ---- write it out and LOOK AT IT ----
       ;;
       ;; The fan-out geometry has only ever been checked NUMERICALLY - parent joint identity and
       ;; child counts.  That confirms the tree's shape, not that atoms land anywhere sensible.  A
       ;; wrong attachment atom, a backwards bond, or a rotamer built off the wrong parent frame is
       ;; obvious in a viewer; through energies it is just "fa_rep is huge", which is
       ;; indistinguishable from twenty other causes.
       ;;
       ;; Expect a HEDGEHOG: 76 sidechains radiating from one CA at locus 53.  That is what a
       ;; fan-out IS, not a bug.
       (format t "~&~%!!!! WRITE IT OUT AND LOOK AT IT~%")
       (defparameter *blueprint-mol2* "/tmp/blueprint-fanout.mol2")
       (let ((agg (topology::aggregate *blueprint*)))
         (chem:save-mol2 (topology::aggregate *blueprint*) *blueprint-mol2* t))
       (format t "~&wrote ~a~%" *blueprint-mol2*)

       ;; ---- stage 3: backbone group + one bonded group per slot, by DETACH/REATTACH ----
       ;;
       ;; Breaks all 470 fan-out bonds, generates the backbone against a clean molecule, then
       ;; reattaches ONE slot at a time to generate its group, and restores every bond at the end.
       ;;
       ;; This exists because the filter-based version (BUILD-BLUEPRINT-SLOT-ENERGY-GROUPS, still
       ;; present) was spending its time enumerating terms it would then reject: with 79
       ;; substituents on a CA there are ~3,000 angles and tens of thousands of dihedrals through
       ;; that one hub, re-enumerated on all 470 passes.  Detaching removes them from the graph
       ;; instead of filtering them out of it, so each pass sees ordinary valence.
       ;;
       ;; BONDED ONLY - generateNonbondEnergyFunctionTables does not take a group yet, so these
       ;; carry no fa_rep.  Every rotamer of a monomer has nearly identical bonded energy, so this
       ;; verifies plumbing and is NOT yet usable for choosing rotamers.
       ;;
       ;; What to read: per-locus term counts (zero means the factory matched nothing), the
       ;; backbone group's term count, and the elapsed time against the filter version.  The two
       ;; drivers should produce the SAME per-slot term counts - if they disagree, one of the two
       ;; partitions is wrong and the detach version is the one to believe, since its exclusions
       ;; are structural rather than predicated.
       (format t "~&~%!!!! STAGE 3: BACKBONE GROUP + ONE BONDED GROUP PER SLOT, BY DETACH/REATTACH~%")
       (format t "~&~%")
       (multiple-value-bind (built backbone-group backbone-inter-group pair-group)
           (ext:with-flame-profile (:path "~/public_html/build-blueprint-energy-groups-by-detach.svg" :rate 97)
             (topology::build-blueprint-energy-groups-by-detach *blueprint* :verbose t))
         (format t "~&slot groups built ~d, backbone intra ~d terms, backbone inter ~d terms~%"
                 built
                 (chem:number-of-terms backbone-group)
                 (chem:number-of-terms backbone-inter-group))
         (defparameter *backbone-group* backbone-group)
         (defparameter *backbone-inter-group* backbone-inter-group)
         ;; The pair group is generated INSIDE the detach scan, with the fan-out bonds broken -
         ;; GENERATE-INTO-GROUP parameterizes every molecule before any keep-predicate runs, and a
         ;; fully attached spiro NG has 16 bonds, which SMIRNOFF rejects outright.
         (defparameter *pair-group* pair-group))

       ;; ---- the ENERGIES object, addressed by the blueprint's frozen mrkindexes ----
       ;;
       ;; Only the ADDRESSING is built here; every energy is left at zero.  Filling them needs the
       ;; per-slot component groups above AND the nonbond terms, which do not exist until cando is
       ;; rebuilt - a group holding only stretch/angle/dihedral cannot rank rotamers.
       ;;
       ;; ENERGY-SCALE is NIL deliberately: it is only read when energies are evaluated, and
       ;; passing a fake one would let a fill step run and produce numbers that look real.
       ;;
       ;; CHECK-BLUEPRINT-ENERGIES-ADDRESSING validates INTERNAL CONSISTENCY - that the parallel
       ;; vectors agree with each other and that SAMPLER's range subtraction reproduces each
       ;; locus's [base, base+extent).  It CANNOT catch a systematically wrong MONOMER-INDEX,
       ;; because a wrong-but-consistent index passes every one of its tests.  That is what the
       ;; mrkey dump below is for: it prints the monomer NAME each mrkey resolves to, the way
       ;; ENERGIES-LMKEY-INDEX-TO-MONOMER-NAME will read it, so the names can be eyeballed against
       ;; the locus table printed earlier in this log.
       (format t "~&~%!!!! THE ENERGIES OBJECT, ADDRESSED BY THE BLUEPRINT'S FROZEN MRKINDEXES~%")
       (format t "~&~%")
       (defparameter *blueprint-energies*
         (design::make-energies-for-blueprint *blueprint* :energy-scale nil))
       ;; ZEROP, not the value itself: this one returns a PROBLEM COUNT, so a true return means
       ;; failure.  Asserting it directly would invert the test and pass only when broken.
       (bench-check blueprint-energies-addressing-consistent
         (zerop (design::check-blueprint-energies-addressing
                 *blueprint* *blueprint-energies*)))
       (let ((mrkey-vector (design::mrkey-vector *blueprint-energies*))
             (monomer-vector (design::monomer-vector *blueprint-energies*)))
         (format t "~&mrkey spot-check (first mrkindex of each monomer-locus):~%")
         (loop for mlocus below (length monomer-vector)
               for mlmm = (design::monomer-locus-max-mrkindex *blueprint-energies*)
               for start = (if (zerop mlocus) 0 (1+ (aref mlmm (1- mlocus))))
               for mrkey = (aref mrkey-vector start)
               for design-monomer = (aref monomer-vector mlocus)
               for names = (topology:monomers design-monomer)
               for mindex = (design::mrkey-monomer-index mrkey)
               ;; Resolve the name only when the index is in range - reading it unconditionally
               ;; would signal on exactly the case this line exists to report.
               for in-range = (< mindex (length names))
               do (format t "~&  monomer-locus ~2d  mrkindex ~3d..~3d  monomer-index ~d -> ~a~%"
                          mlocus start (aref mlmm mlocus) mindex
                          (if in-range
                              (elt names mindex)
                              (format nil "*** OUT OF RANGE, only ~d name~:p: ~s ***"
                                      (length names) names)))))

       ;; ---- fill the single-scan energies from the per-slot groups ----
       ;;
       ;; Expected BEFORE the cando rebuild: nonzero fa_rest (the bonded terms are real) and
       ;; fa_rep exactly 0.0 everywhere, with "470 evaluated with NO nonbond component".  That is
       ;; the honest reading of a bonded-only group, not a fill bug.  AFTER the rebuild the
       ;; no-nonbond count should be 0 and fa_rep should vary between rotamers - that variation IS
       ;; the rotamer ranking signal, so a run where fa_rep stays 0 has learned nothing.
       ;;
       ;; The intermolecular channel has no groups at all yet: the detach scan sets only
       ;; INTRAMOLECULAR-ENERGY-COMPONENTS (blueprint.lisp:2221), so expect 470 there too.
       (format t "~&~%!!!! FILL THE SINGLE-SCAN ENERGIES FROM THE PER-SLOT GROUPS~%")
       (format t "~&~%")
       (design::fill-blueprint-single-energies *blueprint* *blueprint-energies*
                                               *blueprint-coords* :verbose t)
       (design::report-blueprint-single-energies *blueprint* *blueprint-energies*)

       ;; ---- backbone fill ----
       ;;
       ;; Without this both backbone channels read 0.0 and the total silently omits ~5900 terms.
       ;; FILL-BLUEPRINT-SINGLE-ENERGIES fills only the single-scan vectors.
       (format t "~&~%!!!! BACKBONE FILL~%")
       (design::fill-blueprint-backbone-energies *blueprint* *blueprint-energies*
                                                 *blueprint-coords*
                                                 *backbone-group* *backbone-inter-group*
                                                 :verbose t)

       ;; ---- pair-scan ----
       ;;
       ;; ONE group, re-pointed at each pair by SET-MATTERS - not one group per pair.  Generated
       ;; back in the detach scan while the bonds were broken; see that call site.
       ;;
       ;; What to read, in order of what it decides:
       ;;
       ;;   REJECTED vs EVALUATED - the locus bounding-sphere hit rate, and the whole question of
       ;;     whether pair-scan is affordable.  Nobody has measured it on this system; the
       ;;     experiment built to measure it (pairList.h:212-256) was left disabled.  A low
       ;;     rejection rate means one sphere per locus is too loose a bound and something tighter
       ;;     is needed before this is usable.
       ;;
       ;;   SKIPPED - same-locus pairs, which are mutually exclusive rotamers and must never be
       ;;     evaluated.  Reported so the exclusion is visible rather than implicit.
       ;;
       ;;   NON-ZERO - a fill that runs to completion and writes 110,685 zeros is indistinguishable
       ;;     from a correct one by any count, exactly as for the single scans.
       (format t "~&~%!!!! PAIR-SCAN~%")
       (format t "~&~%")
       (design::fill-blueprint-pair-energies *blueprint* *blueprint-energies*
                                             *blueprint-coords* *pair-group* :verbose t)
       (let* ((m (design::pair-scan-energy-lower-triangular-matrix *blueprint-energies*))
              (flags (design::intermolecular-p-pair-scan-energy-lower-triangular-matrix
                      *blueprint-energies*))
              (nonzero 0) (inter 0) (worst 0.0d0))
         (loop for k below (length flags)
               for rest = (aref m (* 2 k))
               for rep = (aref m (1+ (* 2 k)))
               do (unless (and (zerop rest) (zerop rep)) (incf nonzero))
                  (when (= 1 (aref flags k)) (incf inter))
                  (when (> rep worst) (setf worst rep)))
         (format t "~&pair energies: ~d non-zero of ~d slots, ~d flagged intermolecular, ~
                    max fa_rep ~,4f~@[~%~a~]~%"
                 nonzero (length flags) inter worst
                 (when (zerop nonzero)
                   "*** EVERY pair energy is 0.0 - either the sphere test rejected everything or ~
                    the group evaluated to nothing. ***")))

       ;; ---- total energy from the energies object ----
       ;;
       ;; Scored by MCSTATE-ENERGY, the same function the MC movers use - not a reimplementation.
       ;; The six-bucket breakdown is what makes a mismatch actionable: wrong by exactly the
       ;; backbone term points somewhere very different from wrong by the pair term.
       (format t "~&~%!!!! TOTAL ENERGY FROM THE ENERGIES OBJECT~%")
       (format t "~&~%")
       ;; The mcstate must describe THE POSE'S rotamers, not slot 0 everywhere - otherwise the
       ;; total describes a different structure and comparing it to the pose is meaningless.
       ;; The POSE is passed so the mcstate and the pose total below are provably the same
       ;; structure.  The selections come from the assembler's oligomer-shapes, which the blueprint
       ;; SHARES rather than copies - so without the pose there was nothing connecting the mcstate
       ;; to the pose CHECK-TOTAL-ENERGY scores, and a mismatch would have read as a missing term.
       (defparameter *blueprint-mcstate*
         (design::blueprint-mcstate-for-pose *blueprint* path-design::*pose-p1-ce0*))
       ;; ---- can this pose be compared at all? ----
       ;;
       ;; FIRST, because it decides whether anything below means what it appears to.
       (bench-check blueprint-pose-shapes-are-rotamers
         (design::check-pose-shapes-are-rotamers *blueprint* path-design::*pose-p1-ce0*))

       ;; ---- STAGE 1: same MOLECULE? ----
       ;;
       ;; BEFORE the energy comparison, deliberately.  Bonded terms are pure topology - a
       ;; stretch/angle/dihedral/1-4 exists because of the bond graph, not because of where the
       ;; atoms sit - so these two sets must match EXACTLY however far apart the two structures are
       ;; folded.  If they differ, the energy disagreement below is explained and needs no further
       ;; interpretation; if they match, the gap is geometry and the coordinate comparison is next.
       ;;
       ;; *BACKBONE-INTER-GROUP* is deliberately not passed: it holds 0 terms and its factory
       ;; refuses every bonded class, so it cannot contribute here.
       (format t "~&~%!!!! STAGE 1: SAME MOLECULE?~%")
       (bench-check blueprint-bonded-terms-match-pose
         (design::check-bonded-terms-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate* *backbone-group*))

       ;; ---- STAGE 2: same PLACE? ----
       ;;
       ;; Stage 1 covers only the BONDED classes; the nonbond components - which carry essentially
       ;; all of fa_rep - are never compared against the pose, because a pair is a term there only
       ;; if it falls within RPAIRLIST and two folds legitimately differ.  This is what
       ;; disambiguates: coordinates matching means the nonbond sets must match too, so a residual
       ;; fa_rep difference would be an enumeration bug; coordinates differing means geometry, and
       ;; the nonbond difference is downstream of it.
       ;; THE POSE COORDINATES EVERY CHECK BELOW USES.
       ;;
       ;; Built once here and threaded through, rather than each check defaulting to its own
       ;; refold - one fold, so every section is reading the same structure and a difference
       ;; between two sections cannot be an artefact of two separate folds.
       ;;
       ;; NOT (DESIGN:COORDS POSE): that vector was folded before
       ;; MAKE-ASSEMBLER-FOR-COMPLEX-SCANS swapped fresh ROTAMER-SHAPEs onto the moveable receptor
       ;; monomers (protein-binder.lisp:5186), and nothing re-folds afterwards, so it holds the
       ;; receptor's loaded sidechains while its own labels name rotamers the blueprint builds.
       (defparameter *pose-refold-coords*
         (design::refold-pose-coordinates path-design::*pose-p1-ce0*
                                          :stream *standard-output*))

       (bench-check blueprint-coordinates-match-pose
         (design::check-coordinates-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate* *blueprint-coords*
          :pose-coords *pose-refold-coords*))

       ;; ---- STAGE 2b: was the right ROTAMER selected? ----
       ;;
       ;; Stage 1 structurally cannot answer this: every rotamer of a monomer has identical bonded
       ;; terms over identically-named atoms, so selecting the wrong one passes it perfectly and is
       ;; wrong about every atom past CB.  Every rotamer is already materialized in the blueprint's
       ;; coordinates, so this scans them and names the slot that actually fits.
       (bench-check blueprint-rotamer-selection-matches-pose
         (design::check-rotamer-selection-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate* *blueprint-coords*
          :pose-coords *pose-refold-coords*))

       ;; ---- slot -> database rotamer index ----
       ;;
       ;; Every rotamer index in the system indexes ONE vector: a monomer-context's ROTAMER-VECTOR.
       ;; A slot number is a position in the ACTIVE subset, a monomer-shape's ROTAMER-INDEX is an
       ;; index into that full vector, and every confusion here has come from treating one as the
       ;; other.  This writes the mapping out, so the pose's index is a LOOKUP.
       (bench-check blueprint-slot-rotamer-index-audit
         (design::report-slot-rotamer-indexes
          *blueprint* *persona* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- same rotamer, different geometry: compare the JOINTS ----
       ;;
       ;; Internals are applied positionally, three per joint, and a no-op joint still advances the
       ;; index - so identical INTERNALS-VALUES give identical geometry only if the two atresidues
       ;; present the same joints in the same order.
       (bench-check blueprint-joint-order-matches-pose
         (design::check-joint-order-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- and what the first joint hangs OFF ----
       ;;
       ;; The joints inside the residue match; the first one's reference frame comes from its parent
       ;; chain, which runs out into the backbone where that check never looked.
       (bench-check blueprint-root-joint-ancestry-matches-pose
         (design::check-root-joint-ancestry-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- is the agreed-on rotamer index an index into the SAME vector? ----
       ;;
       ;; The slot/database audit says both sides hold rotamer 5605, and the geometry says they hold
       ;; different conformations.  Both can be true: a rotamer index only names a conformation
       ;; relative to a MONOMER-CONTEXT, and the two sides derive their context by different routes -
       ;; the pose from a table built with the assembler, the blueprint recomputed per persona.
       (bench-check blueprint-monomer-context-matches-pose
         (design::check-monomer-context-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- the angles themselves, in degrees ----
       ;;
       ;; Everything above compares a NAME for a conformation.  This compares the conformation: the
       ;; dihedral each side hands to each joint.  It is the last thing upstream of the coordinates,
       ;; so it decides between a LOOKUP fault (different angles) and a FOLD fault (same angles,
       ;; different coordinates) with no third possibility.
       (bench-check blueprint-internal-dihedrals-match-pose
         (design::check-internal-dihedrals-against-pose
          *blueprint* *persona* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- and which internals vector the coordinates actually came from ----
       ;;
       ;; The catalogue values agree, so the rewrite is in the internals VECTOR - the only stage
       ;; between the rotamer and the coordinates.  Each side's vector is compared against the
       ;; catalogue (did adjustments move it?) and against the angle measured from that side's own
       ;; built structure (did the fold honour it?).  The measurement is the arbiter.
       (bench-check blueprint-internals-vectors-match-pose
         (design::check-internals-vectors-against-pose
          *blueprint* *persona* path-design::*pose-p1-ce0* *blueprint-mcstate*
          *blueprint-coords* *blueprint-internals*
          :pose-coords *pose-refold-coords*))

       ;; ---- and WHY the root adjustment fires on one side only ----
       ;;
       ;; The angles show the pose adjusts the spiro root dihedral and the blueprint does not.
       ;; This reports the adjustment objects themselves - registered, initialized, or dead -
       ;; because "never registered" and "registered then disabled by the fan-out" need different
       ;; fixes and the angles cannot tell them apart.
       (bench-check blueprint-adjustments-match-pose
         (design::check-slot-adjustments-against-pose
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- and whether the pose's coordinates are its own pipeline's output ----
       ;;
       ;; CB's phi comes from the backbone via BETA-CARBON, so it cannot depend on the rotamer -
       ;; yet the pose's CB sits 3.84 deg off its own internals while the blueprint's matches.
       ;; Re-folding the pose from its own internals says whether the stored coordinates came from
       ;; this pipeline at all, or whether the internals I reconstruct for it are wrong.
       (bench-check blueprint-pose-refold-is-faithful
         (design::check-pose-refold
          *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate*))

       ;; ---- write both structures out and LOOK at them ----
       ;;
       ;; The blueprint file is TRIMMED to what the mcstate describes - backbone plus one rotamer
       ;; per locus - because its full aggregate holds every rotamer of every monomer at once, which
       ;; is what blueprint-fanout.mol2 already shows and is unreadable for a comparison.  The
       ;; backbones superimpose, so loading both in a viewer puts the misplaced sidechains straight
       ;; in front of you.
       (design::write-blueprint-and-pose-mol2
        *blueprint* path-design::*pose-p1-ce0* *blueprint-mcstate* *blueprint-coords*
        :pose-coords *pose-refold-coords*)

       (design::blueprint-total-energy *blueprint-energies* *blueprint-mcstate*)

       ;; ---- THE comparison: blueprint vs POSE ----
       ;;
       ;; The pose is the ground truth.  complex-scans is NOT a reference - it never matched a pose
       ;; either, and chasing agreement with it would be chasing its flaws.
       ;;
       ;; A disagreement here is real and points at the blueprint.  Read the bucket breakdown above
       ;; to place it: wrong by roughly the backbone term, the single term or the pair term each
       ;; implicate a different half of the machinery.
       ;; Score the pose DIRECTLY - its TOTAL-ENERGY slot is NIL, because
       ;; COMPUTE-ENERGIES-MOVER computes cached energies rather than a scored total.
       ;; MULTIPLE-VALUE-CALL so the pose's fa_rest and fa_rep arrive too, not just the combined
       ;; total.  With them CHECK-TOTAL-ENERGY prints the comparison split into its two halves,
       ;; which is what distinguishes a term-set difference from a geometry difference; without
       ;; them it silently falls back to comparing totals only.
       ;; SCORED ON A RE-FOLD, not on (DESIGN:COORDS POSE).  The pose's stored coordinates were
       ;; folded before MAKE-ASSEMBLER-FOR-COMPLEX-SCANS swapped in fresh ROTAMER-SHAPEs for the
       ;; moveable receptor monomers (protein-binder.lisp:5186), and nothing re-folds afterwards -
       ;; so those sidechains sit at their loaded conformation while their labels name a rotamer
       ;; the blueprint faithfully builds.  Comparing against them compares two structures.
       ;;
       ;; POSE-REFOLD above measured that gap: refold-vs-stored is 0.00 at the three ligand loci
       ;; and 132-180 deg at all nine receptor loci.  This scores what the labels describe, so a
       ;; remaining difference is the blueprint's rather than the pipeline's.
       ;;
       ;; *POSE-REFOLD-COORDS* is the one built before the coordinate checks - the SAME fold those
       ;; sections read, so the energy and the geometry describe one structure.
       (format t "~&~%!!!! THE COMPARISON: BLUEPRINT VS POSE~%")
       (bench-check blueprint-total-energy-matches-pose
         (multiple-value-call #'design::check-total-energy
           *blueprint* *blueprint-energies* *blueprint-mcstate*
           (design::pose-total-energy-direct path-design::*pose-p1-ce0*
                                             :coords *pose-refold-coords*)))

       ;; ---- term-by-term: DUPLICATE detection within the blueprint ----
       ;;
       ;; Every group that contributes to this mcstate, keyed by (class . sorted i3 list).  An
       ;; interaction appearing twice is counted twice in the total, and that is the one failure
       ;; mode no energy number reveals - it just makes the total wrong by a plausible amount.
       ;; The i3<j3 dedup rule and the backbone/slot partition are what prevent it.
       ;;
       ;; NOTE this compares the blueprint AGAINST ITSELF.  Comparing against a pose needs a
       ;; different key: the pose's assembler has its own atom table with ~1000 atoms where the
       ;; blueprint's has 6762, so i3 identifies different atoms in the two and the sets cannot be
       ;; diffed directly.  A name-based key - (molecule-index residue-index atom-name) mapped
       ;; through each side's atom table - is what that needs, plus the slot-residue-to-pose-
       ;; residue correspondence for whichever rotamer is installed.
       (format t "~&~%!!!! TERM-BY-TERM: DUPLICATE DETECTION WITHIN THE BLUEPRINT~%")
       (let ((table (make-hash-table :test 'equal)))
         (design::dump-group-terms *backbone-group* table :label :backbone-intra)
         (design::dump-group-terms *backbone-inter-group* table :label :backbone-inter)
         (loop for mrk across *blueprint-mcstate*
               for scan = (aref (topology::rotamer-scans *blueprint*) mrk)
               do (design::dump-group-terms
                   (topology::intramolecular-energy-components scan) table
                   :label (list :slot-intra mrk))
                  (design::dump-group-terms
                   (topology::intermolecular-energy-components scan) table
                   :label (list :slot-inter mrk)))
         ;; Pair terms for the selected rotamers only - 13 loci gives 78 pairs, not 110,215.
         (let ((vecs (design::blueprint-slot-atom-vectors *blueprint*)))
           (loop for a below (length *blueprint-mcstate*)
                 do (loop for b from (1+ a) below (length *blueprint-mcstate*)
                          for i = (aref *blueprint-mcstate* a)
                          for j = (aref *blueprint-mcstate* b)
                          ;; EVALUATE before dumping.  SET-MATTERS only invalidates the pair
                          ;; list; the terms are rebuilt by the next evaluation.  Dumping straight
                          ;; after set-matters reads whatever the PREVIOUS pair left behind, so
                          ;; every pair reports the same stale term list and each of those terms
                          ;; looks duplicated 78 times.
                          do (topology::set-pair-scan-matters
                              *pair-group* (aref vecs i) (aref vecs j))
                             (design::evaluate-group-alone
                              *blueprint* *pair-group* *blueprint-coords*
                              (chem:make-nvector 2 :initial-element 0.0d0)
                              (design::energy-scale *blueprint-energies*))
                             (design::dump-group-terms *pair-group* table
                                                       :label (list :pair i j)))))
         ;; DUPLICATION IS THE SAME KEY UNDER *DIFFERENT LABELS*, not a key appearing more than
         ;; once.  A label is the bucket, and the failure this exists to catch is one interaction
         ;; landing in two buckets - backbone AND a slot, or two slots - because that is what gets
         ;; counted twice in the total.
         ;;
         ;; Repetition under ONE label is legitimate and common: a proper dihedral is emitted once
         ;; per Fourier component (smirnoff-cache.lisp, the FOURIER-TERMS loop), so a four-component
         ;; torsion appears four times over the same four atoms, in one bucket.  Counting those as
         ;; duplicates reported 589 defects that were not defects - and did so only once this dump
         ;; began covering dihedrals at all, since the old allow-list restricted it to
         ;; ENERGY-ROSETTA-NONBOND where each pair genuinely appears once.
         (let ((dupes 0)
               (multiplicity 0))
           (maphash (lambda (k v)
                      (let ((buckets (remove-duplicates v :test #'equal)))
                        (cond ((> (length buckets) 1)
                               (incf dupes)
                               (when (<= dupes 10)
                                 (format t "~&  DUPLICATE ~s in ~d buckets: ~s~%"
                                         k (length buckets) buckets)))
                              ((> (length v) 1)
                               (incf multiplicity)))))
                    table)
           (format t "~&blueprint term dump: ~d distinct interactions, ~d in MULTIPLE BUCKETS, ~
                      ~d with multiple terms in one bucket (Fourier components - expected)~
                      ~@[~%~a~]~%"
                   (hash-table-count table) dupes multiplicity
                   (when (plusp dupes)
                     "*** An interaction counted in two buckets is counted twice in the total. ***"))))

       ;; ---- copy-remap equivalence ----
       ;;
       ;; The detach scan parameterizes the FIRST rotamer of each monomer and renames its bonded
       ;; terms into every other rotamer.  That is exact only if rotamers of a monomer really do hold
       ;; the same atoms in the same order; a wrong remap yields a full set of plausible terms with
       ;; some indexes on the wrong atoms and no error.
       ;;
       ;; Sample slots that were REMAPPED, never a monomer's first - the reference agrees with itself
       ;; and proves nothing.  TOPOLOGY::BLUEPRINT-LOCUS-OF-MRKINDEX returns the slot number, which
       ;; is what says whether an index is a reference.
       (format t "~&~%!!!! COPY-REMAP EQUIVALENCE~%")
       (format t "~&~%")
       (let ((sample (loop for mrk in '(1 2 95 235 236)
                           for slot = (nth-value 2 (topology::blueprint-locus-of-mrkindex
                                                    *blueprint* mrk))
                           when (and slot (plusp slot)) collect mrk)))
         (if (null sample)
             (format t "~&remap check: no sampled mrkindex is a non-reference rotamer - SKIPPED~%")
             (unless (topology::check-remapped-bonded-terms *blueprint* sample)
               (format t "~&*** Copy-remapped bonded terms do not match a real parameterization. ~
                          Every rotamer that is not a monomer's first is suspect. ***~%"))))

       ;; ---- drift regression ----
       ;;
       ;; The pair lists are cached across evaluations and re-enumerated only when the coordinates
       ;; have moved far enough.  Getting that threshold wrong is SILENT: the energy comes back a
       ;; plausible number computed from a stale list.  This jitters every atom by a sweep of
       ;; distances that straddles the thresholds, and checks the cached energy against a forced
       ;; rebuild at the same coordinates.
       ;;
       ;; Run BEFORE the query-set check, which is destructive - it leaves groups on the full-table
       ;; path, and the point here is to exercise the query-set enumeration's caching.
       (format t "~&~%!!!! DRIFT REGRESSION~%")
       (format t "~&~%")
       (let* ((scan (aref (topology::rotamer-scans *blueprint*) 234))
              (group (and scan (topology::intramolecular-energy-components scan))))
         (if (null group)
             (format t "~&drift check: mrkindex 234 has no intramolecular group - SKIPPED~%")
             (unless (design::check-drift-consistency
                      *blueprint* group *blueprint-coords*
                      :energy-scale (design::energy-scale *blueprint-energies*))
               (format t "~&*** Pair lists are being reused when they should be rebuilt.  Every ~
                          energy above that was computed after any coordinate change is suspect. ***~%"))))

       ;; ---- THE correctness check ----
       ;;
       ;; Everything above is "it runs and the numbers look plausible", which is exactly what a
       ;; subtly wrong enumeration also produces.  This is the one test that distinguishes them:
       ;; the query set changes only HOW candidate pairs are found, never WHICH survive, so the two
       ;; enumerations must agree exactly.
       ;;
       ;; Run LAST and on a few slots of different sizes - it leaves each group on the slow
       ;; full-table path, which is correct but expensive, and it is destructive to the query set.
       ;; mrk 234 is ARG (17 atoms, many intra-slot pairs), 94 is LEU, 0 is a spiro METHYL.
       (format t "~&~%!!!! THE CORRECTNESS CHECK: DOES THE FAST ENUMERATION FIND THE SAME PAIRS?~%")
       (format t "~&  Every energy above was computed with a QUERY SET: the pair enumeration only ~
                    considers atoms in the~%~
                  ~&  slot's index range against a target mask, instead of every atom in the ~
                    table.  That is a search~%~
                  ~&  restriction, NOT a modelling one - the keep-interaction predicate still ~
                    decides which candidate pairs~%~
                  ~&  survive, and it expresses the same restriction the query set does.  So both ~
                    enumerations must produce~%~
                  ~&  the SAME energy, and any difference is a bug in the fast path rather than a ~
                    choice.~%~
                  ~&~%~
                  ~&  Each slot below is evaluated twice: once as the run evaluated it, then again ~
                    with the query set~%~
                  ~&  cleared so candidates come from the whole table.  The second pass is slow by ~
                    construction - that is~%~
                  ~&  the reference it is being checked against.~%~
                  ~&~%~
                  ~&  What a disagreement would mean, by shape:~%~
                  ~&    slow path ~~2x the fast one    the i3 < j3 dedup rule is inverted, so pairs ~
                    are counted twice~%~
                  ~&    slow path MORE negative       the intra-slot brute-force loop is missing ~
                    pairs the full table finds~%~
                  ~&    slow path finds extra terms   the target mask is too narrow and the fast ~
                    path never looked there~%")
       (format t "~&~%")
       (let ((agree 0) (checked 0))
         (dolist (mrk '(0 94 234))
           (incf checked)
           (when (design::check-query-set-against-full-table
                  *blueprint* *blueprint-energies* *blueprint-coords* mrk)
             (incf agree)))
         (format t "~&query-set check: ~d of ~d slots agree~@[~%~a~]~%"
                 agree checked
                 (unless (= agree checked)
                   "*** The fast enumeration does not reproduce the full-table one.  Every energy ~
                    above is suspect.  Check the i3<j3 dedup rule, the intra-slot brute-force ~
                    loop, and the target mask. ***")))

       (format t "Done bench-blueprint.lisp~%")
       )
  (dribble))
