
(in-package :geom)

(progn
  ;; This sets things up so that the clasp-cleavir compiler will inline calls to
  ;; vx et al., to call cc_OVector3_vx et al. directly. The import of this is
  ;; that the result is an unboxed double. As such, if the result is only used
  ;; in ways that don't require boxing (e.g. being put in a double-specialized
  ;; vector, or used for arithmetic) the time and memory required for boxing
  ;; are skipped, improving efficiency.
  ;; Right now the interfaces to do this cleanly have not yet been developed
  ;; so we're kind of messing with compiler internals. BEWARE. FIXME.

  (if (eq (geom:vecreal-type) 'double-float)
      (progn
        (cmp:define-primitive "cc_OVector3_vx" :double-float '(:object))
        (cmp:define-primitive "cc_OVector3_vy" :double-float '(:object))
        (cmp:define-primitive "cc_OVector3_vz" :double-float '(:object))

        (clasp-cleavir::defvprimop-intrinsic %vx ((:double-float) :object) "cc_OVector3_vx")
        (clasp-cleavir::defvprimop-intrinsic %vy ((:double-float) :object) "cc_OVector3_vy")
        (clasp-cleavir::defvprimop-intrinsic %vz ((:double-float) :object) "cc_OVector3_vz")
        )
      (progn
        (cmp:define-primitive "cc_OVector3_vx" :single-float '(:object))
        (cmp:define-primitive "cc_OVector3_vy" :single-float '(:object))
        (cmp:define-primitive "cc_OVector3_vz" :single-float '(:object))

        (clasp-cleavir::defvprimop-intrinsic %vx ((:single-float) :object) "cc_OVector3_vx")
        (clasp-cleavir::defvprimop-intrinsic %vy ((:single-float) :object) "cc_OVector3_vy")
        (clasp-cleavir::defvprimop-intrinsic %vz ((:single-float) :object) "cc_OVector3_vz")
        ))

  (cc-bir-to-bmir::deftransform vx %vx t)
  (cc-bir-to-bmir::deftransform vy %vy t)
  (cc-bir-to-bmir::deftransform vz %vz t)

  ;; These are in place so that (a) type derivation doesn't insert stupid
  ;; checks, and (b) the bir-to-bmir transforms actually fire. (b) is a KLUDGE
  ;; due to how identities work in cleavir/setup.lisp.
  (clasp-cleavir::define-deriver geom:vx (o)
    (declare (ignore o))
    (let ((sys clasp-cleavir:*clasp-system*))
      (cleavir-ctype:single-value
       (cleavir-ctype:range (geom:vecreal-type) '* '* sys) sys)))
  (clasp-cleavir::define-deriver geom:vy (o)
    (declare (ignore o))
    (let ((sys clasp-cleavir:*clasp-system*))
      (cleavir-ctype:single-value
       (cleavir-ctype:range (geom:vecreal-type) '* '* sys) sys)))
  (clasp-cleavir::define-deriver geom:vz (o)
    (declare (ignore o))
    (let ((sys clasp-cleavir:*clasp-system*))
      (cleavir-ctype:single-value
       (cleavir-ctype:range (geom:vecreal-type) '* '* sys) sys)))
  )
