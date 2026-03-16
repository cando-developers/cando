;;;; -------------------------------
;;;; energy-kernels.lisp
;;;; -------------------------------

(in-package :mathkernel-user)


(defparameter *pipeline*
  (stmt-ir:make-optimization-pipeline
   :name :kernel-full
   :optimizations
   (list

    ;; linear canonicalization to unify linear forms before factoring + CSE
    (stmt-ir:make-optimization
     :name :linear-canonicalization
     :function 'stmt-ir:linear-canonicalization-optimization)

    ;; factor sums
    (stmt-ir:make-optimization
     :name :factor-sums
     :function 'stmt-ir:factor-sums-optimization
     :keyword-args (list :min-uses 2
                         :min-factors 1
                         :min-size 4))
    ;; CSE over full block
    (stmt-ir:make-optimization
     :name :cse-full
     :function 'stmt-ir:cse-block-multi-optimization
     :keyword-args (list :max-passes 50
                         :min-uses 2
                         :min-size 1))

    ;; inverse exponent optimization
    (stmt-ir:make-optimization
     :name :inverse-expt
     :function 'stmt-ir:inverse-expt-optimization
     :keyword-args (list :min-uses 2))

    ;; temp+param factoring
    (stmt-ir:make-optimization
     :name :factor-temp-param
     :function 'stmt-ir:factor-temp-param-products-optimization
     :keyword-args (list :min-uses    2
                         :min-factors 2
                         :max-factors 3 ))
    ;; copy-propagation
    (stmt-ir:make-optimization
     :name :copy-propagate
     :function 'stmt-ir:copy-propagate-optimization)

    ;; expt aliasing
    (stmt-ir:make-optimization
     :name :alias-assigned-exprs
     :function 'stmt-ir:alias-assigned-exprs-optimization)

    ;; sign normalization
    (stmt-ir:make-optimization
     :name :normalize-signs
     :function 'stmt-ir:normalize-signs-optimization))))

(defparameter kernels nil)


#+(or)
(progn
  (warn "!!!!!!!!!!!!!!!! Not building any kernels !!!!!!!!!!!!!!!!!")
  (test-nonbond-dd))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       ;; ================================================================
       ;; Flat-bottom anchor restraint
       ;;
       ;; Restrains atom to lie within a sphere of radius R0 centered at (xa,ya,za).
       ;; Energy is zero inside the sphere (the "flat bottom") and harmonic outside:
       ;;
       ;;   E = 0                         if r <= R0
       ;;   E = ka * (r - R0)^2           if r >  R0
       ;;
       ;; where r = distance from atom to center.
       ;;
       ;; Gradient:  dE/dq_i = 2*ka*(r-R0) * (q_i - q_i0) / r
       ;; Hessian:   full (outer-product + curvature terms)
       ;; ================================================================

(build-kernel-group (kernels "flat_bottom_anchor_restraint" (:energy :gradient :hessian))
  (:pipeline *pipeline*)
  (:coordinate-inputs 1)

  (:term-struct
   (:name "flat_bottom_anchor_restraint_term")
   (:var term)
   (:inputs ((|double| ka)
             (|double| R0)
             (|double| R02)
             (|double| xa)
             (|double| ya)
             (|double| za)))
   (:struct ((|double| ka)
             (|double| R0)
             (|double| R02)
             (|double| xa)
             (|double| ya)
             (|double| za)))
   (:derived ()))

  (:body
   (stmt-block :!base
     (=. dx "x1 - xa")
     (=. dy "y1 - ya")
     (=. dz "z1 - za")
     (=. r2 "dx*dx + dy*dy + dz*dz")
     (=. rr "sqrt(r2)")

     (stmt-cond
      ("rr<=R0"
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian))
       (=. energy "0.0"))
      (t
       (=. dr "rr - R0")
       (=. dE_dr "2.0*ka*dr" :modes (:gradient :hessian))
       (=. d2E_dr2 "2.0*ka" :modes (:hessian))
       (=. energy "ka * dr * dr")))))

  (:derivatives
   (:mode :manual
    :intermediates (rr)

    :intermediate->coord
    ((rr ((x1 "dx/rr")
          (y1 "dy/rr")
          (z1 "dz/rr"))))

    :intermediate->coord2
    ((rr (((x1 x1) "(rr*rr - dx*dx)/(rr^3)")
          ((x1 y1) "-(dx*dy)/(rr^3)")
          ((x1 z1) "-(dx*dz)/(rr^3)")
          ((y1 y1) "(rr*rr - dy*dy)/(rr^3)")
          ((y1 z1) "-(dy*dz)/(rr^3)")
          ((z1 z1) "(rr*rr - dz*dz)/(rr^3)"))))

    :energy->intermediate
    (:gradient ((rr "dE_dr"))
     :hessian  (((rr rr) "d2E_dr2")))

    :hessian-modes ((rr :full))
    :geometry-check :warn)))



(build-kernel-group (kernels "rosetta_elec_cutoff" (:gradient :energy :hessian))
  (:pipeline *pipeline*)
  (:coordinate-inputs 2)

  (:params-struct
   (:name "rosetta_elec_parameters")
   (:var "params")
   (:struct ((|double| elec_weight 1.0)
             (|double| eps_core 6.0)
             (|double| eps_solvent 80.0)
             (|double| rmin 1.45)
             (|double| rlow 1.85)
             (|double| rhi 4.5)
             (|double| rcut 5.5)
             (|double| rpairlist 7.5)
             )))

  (:term-struct
   (:name "rosetta_elec_term")
   (:var term)
   (:inputs ((|double| kqq)
             ))
   (:struct ((|double| kqq)
             (|double| E_rmin)
             (|double| AA_low)
             (|double| BB_low)
             (|double| CC_low)
             (|double| DD_low)
             (|double| AA_high)
             (|double| BB_high)
             (|double| CC_high)
             (|double| DD_high)))
   (:derived
    (
     ;; Eelec at rmin
     (rx_rmin "rmin*0.25")
     (ex_rmin "exp(-rx_rmin)")
     (gx_rmin "(1.0 + rx_rmin + 0.5*rx_rmin*rx_rmin) * ex_rmin")
     (ep_rmin "eps_solvent + (eps_core - eps_solvent)*gx_rmin")
     (E_rmin "kqq/ep_rmin")

     ;; Eelec + dE/dr at rlow
     (rx_rlow "rlow*0.25")
     (ex_rlow "exp(-rx_rlow)")
     (gx_rlow "(1.0 + rx_rlow + 0.5*rx_rlow*rx_rlow) * ex_rlow")
     (ep_rlow "eps_solvent + (eps_core - eps_solvent)*gx_rlow")
     (E_rlow "kqq/ep_rlow")
     (gx_p_rlow "(-0.5*rx_rlow*rx_rlow) * ex_rlow")
     (ep_p_rlow "(eps_core - eps_solvent) * gx_p_rlow * 0.25")
     (dE_rlow "-kqq * ep_p_rlow / (ep_rlow*ep_rlow)")

     ;; Eelec + dE/dr at rhi
     (rx_rhi "rhi*0.25")
     (ex_rhi "exp(-rx_rhi)")
     (gx_rhi "(1.0 + rx_rhi + 0.5*rx_rhi*rx_rhi) * ex_rhi")
     (ep_rhi "eps_solvent + (eps_core - eps_solvent)*gx_rhi")
     (E_rhi "kqq/ep_rhi")
     (gx_p_rhi "(-0.5*rx_rhi*rx_rhi) * ex_rhi")
     (ep_p_rhi "(eps_core - eps_solvent) * gx_p_rhi * 0.25")
     (dE_rhi "-kqq * ep_p_rhi / (ep_rhi*ep_rhi)")

     (dr_low  "rlow - rmin")
     (dr_high "rcut - rhi")

     ;; Low Hermite (shifted)
     (aa_low_shifted "(2*E_rmin - 2*E_rlow + (0.0 + dE_rlow)*dr_low) / (dr_low*dr_low*dr_low)")
     (bb_low_shifted "(-3*E_rmin + 3*E_rlow - (0.0 + dE_rlow)*dr_low) / (dr_low*dr_low)")
     (cc_low_shifted "0.0")
     (dd_low_shifted "E_rmin")

     ;; Low Hermite (unshifted)
     (AA_low "aa_low_shifted")
     (BB_low "bb_low_shifted - 3*aa_low_shifted*rmin")
     (CC_low "cc_low_shifted - 2*bb_low_shifted*rmin + 3*aa_low_shifted*rmin*rmin")
     (DD_low "dd_low_shifted - cc_low_shifted*rmin + bb_low_shifted*rmin*rmin - aa_low_shifted*rmin*rmin*rmin")

     ;; High Hermite (shifted)
     (aa_high_shifted "(2*E_rhi - 0.0 + (dE_rhi + 0.0)*dr_high) / (dr_high*dr_high*dr_high)")
     (bb_high_shifted "(-3*E_rhi + 0.0 - (2*dE_rhi + 0.0)*dr_high) / (dr_high*dr_high)")
     (cc_high_shifted "dE_rhi")
     (dd_high_shifted "E_rhi")

     ;; High Hermite (unshifted)
     (AA_high "aa_high_shifted")
     (BB_high "bb_high_shifted - 3*aa_high_shifted*rhi")
     (CC_high "cc_high_shifted - 2*bb_high_shifted*rhi + 3*aa_high_shifted*rhi*rhi")
     (DD_high "dd_high_shifted - cc_high_shifted*rhi + bb_high_shifted*rhi*rhi - aa_high_shifted*rhi*rhi*rhi"))))

  #+(or)(:layout ((1 . i3x1) (2 . i3x2))
           ((X . 0) (Y . 1) (Z . 2)))

  #+(or)(:coord-vars (x1 y1 z1
                   x2 y2 z2))
  #+(or)(:coord-load
         (coords-from-position
          ((x1 y1 z1 i3x1)
           (x2 y2 z2 i3x2))))

  (:body
   (stmt-block :!base
     (=. dx "x1 - x2")
     (=. dy "y1 - y2")
     (=. dz "z1 - z2")
     (=. r2 "dx*dx + dy*dy + dz*dz")
     (=. rr "sqrt(r2)")

     (stmt-cond
      ("rr<=rmin"
       (=. energy "elec_weight*E_rmin")
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian)))

      ("rr<=rlow"
       (=. energy "elec_weight*(AA_low*rr*rr*rr + BB_low*rr*rr + CC_low*rr + DD_low)")
       (=. dE_dr   "elec_weight*(3*AA_low*rr*rr + 2*BB_low*rr + CC_low)" :modes
           (:gradient :hessian))
       (=. d2E_dr2 "elec_weight*(6*AA_low*rr + 2*BB_low)" :modes (:hessian)))

      ("rr<=rhi"
       ;; runtime Eelec(rr)
       (=. rx "rr*0.25")
       (=. ex "exp(-rx)")
       (=. gx "(1.0 + rx + 0.5*rx*rx) * ex")
       (=. gx_p "(-0.5*rx*rx) * ex")
       (=. gx_pp "(0.5*rx*rx - rx) * ex")
       (=. ep "eps_solvent + (eps_core - eps_solvent)*gx")
       (=. ep_p "(eps_core - eps_solvent) * gx_p * 0.25")
       (=. ep_pp "(eps_core - eps_solvent) * gx_pp * 0.0625")
       (=. ee "kqq/ep")
       (=. energy "elec_weight*ee")
       (=. dE_dr   "elec_weight*(-kqq*ep_p/(ep*ep))" :modes (:gradient :hessian))
       (=. d2E_dr2 "elec_weight*(-kqq*(ep_pp/(ep*ep) - 2*ep_p*ep_p/(ep*ep*ep)))" :modes
           (:hessian)))

      ("rr<=rcut"
       (=. energy "elec_weight*(AA_high*rr*rr*rr + BB_high*rr*rr + CC_high*rr + DD_high)")
       (=. dE_dr   "elec_weight*(3*AA_high*rr*rr + 2*BB_high*rr + CC_high)" :modes
           (:gradient :hessian))
       (=. d2E_dr2 "elec_weight*(6*AA_high*rr + 2*BB_high)" :modes (:hessian)))

      (t
       (=. energy "0.0")
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian))))))

  (:derivatives
   (:mode :manual
    :intermediates (rr)

    :intermediate->coord
    ((rr ((x1 "dx / rr")
          (y1 "dy / rr")
          (z1 "dz / rr")
          (x2 "-dx / rr")
          (y2 "-dy / rr")
          (z2 "-dz / rr"))))

    :intermediate->coord2
    ((rr (
          ((x1 x1) "(rr*rr - dx*dx)/(rr^3)")
          ((x1 y1) "(-dx*dy)/(rr^3)")
          ((x1 z1) "(-dx*dz)/(rr^3)")
          ((x1 x2) "(dx*dx - rr*rr)/(rr^3)")
          ((x1 y2) "dx*dy/(rr^3)")
          ((x1 z2) "dx*dz/(rr^3)")

          ((y1 y1) "(rr*rr - dy*dy)/(rr^3)")
          ((y1 z1) "(-dy*dz)/(rr^3)")
          ((y1 x2) "dx*dy/(rr^3)")
          ((y1 y2) "(dy*dy - rr*rr)/(rr^3)")
          ((y1 z2) "dy*dz/(rr^3)")

          ((z1 z1) "(rr*rr - dz*dz)/(rr^3)")
          ((z1 x2) "dx*dz/(rr^3)")
          ((z1 y2) "dy*dz/(rr^3)")
          ((z1 z2) "(dz*dz - rr*rr)/(rr^3)")

          ((x2 x2) "(rr*rr - dx*dx)/(rr^3)")
          ((x2 y2) "(-dx*dy)/(rr^3)")
          ((x2 z2) "(-dx*dz)/(rr^3)")

          ((y2 y2) "(rr*rr - dy*dy)/(rr^3)")
          ((y2 z2) "(-dy*dz)/(rr^3)")

          ((z2 z2) "(rr*rr - dz*dz)/(rr^3)"))))

    :energy->intermediate
    (:gradient ((rr "dE_dr"))
     :hessian  (((rr rr) "d2E_dr2")))

    :hessian-modes ((rr :full))
    :geometry-check :warn)))


(build-kernel-group (kernels "rosetta_nonbond_dd_cutoff" (:gradient :energy :hessian))
  (:pipeline *pipeline*)
  (:coordinate-inputs 2)

  #+(or)(:layout ((1 . I3X1) (2 . I3X2))
           ((X . 0) (Y . 1) (Z . 2)))

  #+(or)(:coord-vars (x1 y1 z1
                   x2 y2 z2))

  #+(or)(:coord-load
         (coords-from-position
          ((x1 y1 z1 i3x1)
           (x2 y2 z2 i3x2))))

  (:params-struct
   (:name "rosetta_nonbond_parameters")
   (:var "params")
   (:struct ((|double| rep_weight 1.0)
             (|double| rswitch 4.5)
             (|double| rcut 6.0)
             (|double| rpairlist 9.0)
             )))

  (:term-struct
   (:name "rosetta_nonbond_term")
   (:var term)
   (:inputs ((|double| A)
             (|double| C)))
   (:struct ((|double| sigma)
             (|double| epsilon)
             (|double| mij)
             (|double| bij)
             (|double| aa)
             (|double| bb)
             (|double| cc)
             (|double| dd)))
   (:derived
    (
     ;; derive Rosetta sigma/epsilon from AMBER A,C
     (sigma6  "2.0*A/C")
     (sigma   "sigma6^(1.0/6.0)")
     (epsilon "C*C/(4.0*A)")

     ;; repulsive linearization at r = 0.6*sigma
     (rrep "0.6*sigma")
     (sr_rep "sigma/rrep")
     (sr2_rep "sr_rep*sr_rep")
     (sr6_rep "sr2_rep*sr2_rep*sr2_rep")
     (sr12_rep "sr6_rep*sr6_rep")
     (rep_rrep "epsilon*(sr12_rep - 2*sr6_rep + 1.0)")
     (drep_rrep "epsilon*((-12*sr12_rep + 12*sr6_rep)/rrep)")
     (mij "drep_rrep")
     (bij "rep_rrep - mij*rrep")

     ;; attractive cubic from rswitch to rcut (params.rswitch / params.rcut)
     (sr_sw "sigma/rswitch")
     (sr2_sw "sr_sw*sr_sw")
     (sr6_sw "sr2_sw*sr2_sw*sr2_sw")
     (sr12_sw "sr6_sw*sr6_sw")
     (atr_sw "epsilon*(sr12_sw - 2*sr6_sw)")
     (d_atr_sw "epsilon*((-12*sr12_sw + 12*sr6_sw)/rswitch)")
     (dr "rcut - rswitch")

     ;; shifted cubic coefficients
     (aa_shift "(2*atr_sw + d_atr_sw*dr)/(dr*dr*dr)")
     (bb_shift "(-3*atr_sw - 2*d_atr_sw*dr)/(dr*dr)")
     (cc_shift "d_atr_sw")
     (dd_shift "atr_sw")

     ;; unshifted cubic coefficients
     (aa "aa_shift")
     (bb "bb_shift - 3*aa_shift*rswitch")
     (cc "cc_shift - 2*bb_shift*rswitch + 3*aa_shift*rswitch*rswitch")
     (dd "dd_shift - cc_shift*rswitch + bb_shift*rswitch*rswitch - aa_shift*rswitch*rswitch*rswitch")
     )))

  (:body
   (stmt-block :!base
     (=. dx "x1 - x2")
     (=. dy "y1 - y2")
     (=. dz "z1 - z2")
     (=. r2 "dx*dx + dy*dy + dz*dz")
     (=. rr "sqrt(r2)")

     (stmt-cond
      ("rr<=0.6*sigma"
       (=. fa_rep "mij*rr + bij")
       (=. fa_atr "-epsilon")
       (=. dE_dr "rep_weight*mij" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian))
       (=. energy "rep_weight*fa_rep + fa_atr"))

      ("rr<=sigma"
       (=. sr   "sigma/rr")
       (=. sr2  "sr*sr")
       (=. sr6  "sr2*sr2*sr2")
       (=. sr12 "sr6*sr6")
       (=. fa_rep_epsilon "epsilon*(sr12 - 2*sr6)")
       (=. fa_rep "fa_rep_epsilon + epsilon")
       (=. fa_atr "-epsilon")
       (=. dE_dr   "rep_weight*epsilon*((-12*sr12 + 12*sr6)/rr)" :modes (:gradient :hessian))
       (=. d2E_dr2 "rep_weight*epsilon*((156*sr12 - 42*sr6)/(rr*rr))" :modes (:hessian))
       (=. energy "rep_weight*fa_rep + fa_atr"))

      ("rr<=rswitch"
       (=. sr   "sigma/rr")
       (=. sr2  "sr*sr")
       (=. sr6  "sr2*sr2*sr2")
       (=. sr12 "sr6*sr6")
       (=. fa_atr "epsilon*(sr12 - 2*sr6)")
       (=. dE_dr   "epsilon*((-12*sr12 + 12*sr6)/rr)" :modes (:gradient :hessian))
       (=. d2E_dr2 "epsilon*((156*sr12 - 42*sr6)/(rr*rr))" :modes (:hessian))
       (=. energy "fa_atr"))

      ("rr<=rcut"
       (=. fa_atr "aa*rr*rr*rr + bb*rr*rr + cc*rr + dd")
       (=. dE_dr   "(3*aa*rr*rr + 2*bb*rr + cc)" :modes (:gradient :hessian))
       (=. d2E_dr2 "(6*aa*rr + 2*bb)" :modes (:hessian))
       (=. energy "fa_atr"))

      (t
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian))
       (=. energy "0.0")))))

  (:derivatives
   (:mode :manual
    :intermediates (rr)

    :intermediate->coord
    ((rr ((x1 "dx / rr")
          (y1 "dy / rr")
          (z1 "dz / rr")
          (x2 "-dx / rr")
          (y2 "-dy / rr")
          (z2 "-dz / rr"))))

    :intermediate->coord2
    ((rr (
          ((x1 x1) "(rr*rr - dx*dx)/(rr^3)")
          ((x1 y1) "(-dx*dy)/(rr^3)")
          ((x1 z1) "(-dx*dz)/(rr^3)")
          ((x1 x2) "(dx*dx - rr*rr)/(rr^3)")
          ((x1 y2) "dx*dy/(rr^3)")
          ((x1 z2) "dx*dz/(rr^3)")

          ((y1 y1) "(rr*rr - dy*dy)/(rr^3)")
          ((y1 z1) "(-dy*dz)/(rr^3)")
          ((y1 x2) "dx*dy/(rr^3)")
          ((y1 y2) "(dy*dy - rr*rr)/(rr^3)")
          ((y1 z2) "dy*dz/(rr^3)")

          ((z1 z1) "(rr*rr - dz*dz)/(rr^3)")
          ((z1 x2) "dx*dz/(rr^3)")
          ((z1 y2) "dy*dz/(rr^3)")
          ((z1 z2) "(dz*dz - rr*rr)/(rr^3)")

          ((x2 x2) "(rr*rr - dx*dx)/(rr^3)")
          ((x2 y2) "(-dx*dy)/(rr^3)")
          ((x2 z2) "(-dx*dz)/(rr^3)")

          ((y2 y2) "(rr*rr - dy*dy)/(rr^3)")
          ((y2 z2) "(-dy*dz)/(rr^3)")

          ((z2 z2) "(rr*rr - dz*dz)/(rr^3)"))))

    :energy->intermediate
    (:gradient ((rr "dE_dr"))
     :hessian  (((rr rr) "d2E_dr2")))

    :hessian-modes ((rr :full))
    :geometry-check :warn)))

(build-kernel-group (kernels "rosetta_lk_solvation" (:gradient :energy :hessian))
  (:documentation "
Term parameters deltaGfreei, Vj, lambdai, sigmai - they are provided in the term-struct
 Pi is the physical constant PI (3.14159265358979...)
 c0, c1 are constants that will be defined in the params-struct
 f_solv_low_poly and f_solv_high_poly are cubic polynomials with form \"aa*x^3+bb*x^2+cc*x+dd\"
; Here is the text from the paper: The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design:

Solvation. Native-like protein conformations minimize the
exposure of hydrophobic side chains to the surrounding polar
solvent. Unfortunately, explicitly modeling all of the interactions between solvent and protein atoms is computationally
expensive. Instead, Rosetta represents the solvent as bulk water
based upon the Lazaridis−Karplus (LK) implicit Gaussian
exclusion model.36 Rosetta’s solvation model has two
components: an isotropic solvation energy, called fa_sol,
which assumes that bulk water is uniformly distributed around
the atoms (Figure 2A), and an anisotropic solvation energy,
called lk_ball_wtd, which accounts for specific waters near
polar atoms that form the solvation shell (Figure 2B).
The isotropic (LK) model36 is based on the function fdesolv
that describes the energy required to desolvate (remove contacting water)
an atom i when it is approached by a
neighboring atom j. In Rosetta, we exclude the LK deltaGref term
because we implement our own reference energy (discussed later).
The energy of the atom-pair interaction varies with the separation distance, dij,
the experimentally determined vaporto-water transfer free energy, deltaGi
free, the sum of the atomic radii, sigmaij,
the correlation length, lambdai , and the atomic volume of
the desolvating atom, Vj (eq 11(fdesolv)):

At short distances, fa_rep prevents atoms from overlapping;
however, many protocols briefly downweight or disable the
fa_rep term. To avoid scenarios where fdesolv encourages
atom-pair overlap in the absence of fa_rep, we smoothly
increase the value of the function to a constant at short
distances where the van der Waals spheres overlap (dij = sigmaij).
At large distances, the function asymptotically approaches zero;
therefore, we truncate the function at 6.0 Å for speed. We also
transition between the constants at short and long distances
using the distance-dependent cubic polynomials f_solv_low_poly and f_solv_high_poly
with constants c0 = 0.3 Å and c1 = 0.2 Å that define a
window for smoothing. The overall desolvation function is
given by eq 12(gdesolv(i,j)):
The total isotropic solvation energy, Efa_sol, is computed as a
sum of the energies for atom j desolvating atom i and vice versa,
scaled by the previously defined connectivity weight (eq 13(Efa_sol)):
")




  (:pipeline *pipeline*)
  (:coordinate-inputs 2)

  (:params-struct
   (:name "rosetta_lk_solvation_parameters")
   (:var "params")
   (:struct ((|double| c0 0.3)
             (|double| c1 0.2)
             (|double| r_solv_low 4.5)
             (|double| r_solv_high 6.0)
             (|double| rpairlist 8.0)
             )))

  (:term-struct
   (:name "rosetta_lk_solvation_term")
   (:var term)
   (:inputs ((|double| lk_dgfree_i)
             (|double| lk_lambda_i)
             (|double| lj_radius_i)
             (|double| lk_volume_i)
             (|double| lk_dgfree_j)
             (|double| lk_lambda_j)
             (|double| lj_radius_j)
             (|double| lk_volume_j)))
   (:struct ((|double| lj_radius_i)
             (|double| lj_radius_j)
             (|double| sigmaij)
             (|double| r0_low)
             (|double| r1_low)
             (|double| r2_high)
             (|double| r3_high)
             (|double| dr_low)
             (|double| dr_high)
             (|double| inv2pi32)
             (|double| lambda2_i)
             (|double| lambda2_j)
             (|double| pref_i_j)
             (|double| pref_j_i)

             (|double| f_contact_i)
             (|double| f_contact_j)
             (|double| f_r1_i)
             (|double| f_r1_j)
             (|double| df_r1_i)
             (|double| df_r1_j)
             (|double| f_r2_i)
             (|double| f_r2_j)
             (|double| df_r2_i)
             (|double| df_r2_j)

             (|double| AA_low_i)
             (|double| BB_low_i)
             (|double| CC_low_i)
             (|double| DD_low_i)
             (|double| AA_low_j)
             (|double| BB_low_j)
             (|double| CC_low_j)
             (|double| DD_low_j)

             (|double| AA_high_i)
             (|double| BB_high_i)
             (|double| CC_high_i)
             (|double| DD_high_i)
             (|double| AA_high_j)
             (|double| BB_high_j)
             (|double| CC_high_j)
             (|double| DD_high_j)))
   (:derived
    (
     (sigmaij "lj_radius_i + lj_radius_j")
     (r0_low "sigmaij - c0")
     (r1_low "sigmaij + c1")
     (r2_high "r_solv_low")
     (r3_high "r_solv_high")
     (dr_low "r1_low - r0_low")
     (dr_high "r3_high - r2_high")

     (inv2pi32 "1.0/(2.0*Pi*sqrt(Pi))")
     (lambda2_i "lk_lambda_i*lk_lambda_i")
     (lambda2_j "lk_lambda_j*lk_lambda_j")
     (pref_i_j "-inv2pi32 * lk_dgfree_i * (1.0/lk_lambda_i) * lk_volume_j")
     (pref_j_i "-inv2pi32 * lk_dgfree_j * (1.0/lk_lambda_j) * lk_volume_i")

     ;; fdesolv at contact (r = sigmaij)
     (delta_contact_i "sigmaij - lj_radius_i")
     (exp_contact_i "exp(-(delta_contact_i*delta_contact_i)/lambda2_i)")
     (f_contact_i "pref_i_j * exp_contact_i / (sigmaij*sigmaij)")
     (delta_contact_j "sigmaij - lj_radius_j")
     (exp_contact_j "exp(-(delta_contact_j*delta_contact_j)/lambda2_j)")
     (f_contact_j "pref_j_i * exp_contact_j / (sigmaij*sigmaij)")

     ;; fdesolv at r1_low
     (delta_r1_i "r1_low - lj_radius_i")
     (exp_r1_i "exp(-(delta_r1_i*delta_r1_i)/lambda2_i)")
     (f_r1_i "pref_i_j * exp_r1_i / (r1_low*r1_low)")
     (df_r1_i "f_r1_i * (-2.0*(delta_r1_i/lambda2_i + 1.0/r1_low))")
     (delta_r1_j "r1_low - lj_radius_j")
     (exp_r1_j "exp(-(delta_r1_j*delta_r1_j)/lambda2_j)")
     (f_r1_j "pref_j_i * exp_r1_j / (r1_low*r1_low)")
     (df_r1_j "f_r1_j * (-2.0*(delta_r1_j/lambda2_j + 1.0/r1_low))")

     ;; fdesolv at r2_high
     (delta_r2_i "r2_high - lj_radius_i")
     (exp_r2_i "exp(-(delta_r2_i*delta_r2_i)/lambda2_i)")
     (f_r2_i "pref_i_j * exp_r2_i / (r2_high*r2_high)")
     (df_r2_i "f_r2_i * (-2.0*(delta_r2_i/lambda2_i + 1.0/r2_high))")
     (delta_r2_j "r2_high - lj_radius_j")
     (exp_r2_j "exp(-(delta_r2_j*delta_r2_j)/lambda2_j)")
     (f_r2_j "pref_j_i * exp_r2_j / (r2_high*r2_high)")
     (df_r2_j "f_r2_j * (-2.0*(delta_r2_j/lambda2_j + 1.0/r2_high))")

     ;; Low Hermite (shifted)
     (aa_low_i_shift "(2*f_contact_i - 2*f_r1_i + (0.0 + df_r1_i)*dr_low) / (dr_low*dr_low*dr_low)")
     (bb_low_i_shift "(-3*f_contact_i + 3*f_r1_i - (0.0 + df_r1_i)*dr_low) / (dr_low*dr_low)")
     (cc_low_i_shift "0.0")
     (dd_low_i_shift "f_contact_i")

     (aa_low_j_shift "(2*f_contact_j - 2*f_r1_j + (0.0 + df_r1_j)*dr_low) / (dr_low*dr_low*dr_low)")
     (bb_low_j_shift "(-3*f_contact_j + 3*f_r1_j - (0.0 + df_r1_j)*dr_low) / (dr_low*dr_low)")
     (cc_low_j_shift "0.0")
     (dd_low_j_shift "f_contact_j")

     ;; Low Hermite (unshifted)
     (AA_low_i "aa_low_i_shift")
     (BB_low_i "bb_low_i_shift - 3*aa_low_i_shift*r0_low")
     (CC_low_i "cc_low_i_shift - 2*bb_low_i_shift*r0_low + 3*aa_low_i_shift*r0_low*r0_low")
     (DD_low_i "dd_low_i_shift - cc_low_i_shift*r0_low + bb_low_i_shift*r0_low*r0_low - aa_low_i_shift*r0_low*r0_low*r0_low")

     (AA_low_j "aa_low_j_shift")
     (BB_low_j "bb_low_j_shift - 3*aa_low_j_shift*r0_low")
     (CC_low_j "cc_low_j_shift - 2*bb_low_j_shift*r0_low + 3*aa_low_j_shift*r0_low*r0_low")
     (DD_low_j "dd_low_j_shift - cc_low_j_shift*r0_low + bb_low_j_shift*r0_low*r0_low - aa_low_j_shift*r0_low*r0_low*r0_low")

     ;; High Hermite (shifted)
     (aa_high_i_shift "(2*f_r2_i - 0.0 + (df_r2_i + 0.0)*dr_high) / (dr_high*dr_high*dr_high)")
     (bb_high_i_shift "(-3*f_r2_i + 0.0 - (2*df_r2_i + 0.0)*dr_high) / (dr_high*dr_high)")
     (cc_high_i_shift "df_r2_i")
     (dd_high_i_shift "f_r2_i")

     (aa_high_j_shift "(2*f_r2_j - 0.0 + (df_r2_j + 0.0)*dr_high) / (dr_high*dr_high*dr_high)")
     (bb_high_j_shift "(-3*f_r2_j + 0.0 - (2*df_r2_j + 0.0)*dr_high) / (dr_high*dr_high)")
     (cc_high_j_shift "df_r2_j")
     (dd_high_j_shift "f_r2_j")

     ;; High Hermite (unshifted)
     (AA_high_i "aa_high_i_shift")
     (BB_high_i "bb_high_i_shift - 3*aa_high_i_shift*r2_high")
     (CC_high_i "cc_high_i_shift - 2*bb_high_i_shift*r2_high + 3*aa_high_i_shift*r2_high*r2_high")
     (DD_high_i "dd_high_i_shift - cc_high_i_shift*r2_high + bb_high_i_shift*r2_high*r2_high - aa_high_i_shift*r2_high*r2_high*r2_high")

     (AA_high_j "aa_high_j_shift")
     (BB_high_j "bb_high_j_shift - 3*aa_high_j_shift*r2_high")
     (CC_high_j "cc_high_j_shift - 2*bb_high_j_shift*r2_high + 3*aa_high_j_shift*r2_high*r2_high")
     (DD_high_j "dd_high_j_shift - cc_high_j_shift*r2_high + bb_high_j_shift*r2_high*r2_high - aa_high_j_shift*r2_high*r2_high*r2_high")
     )))

  (:body
   (stmt-block :!base
     (=. dx "x1 - x2")
     (=. dy "y1 - y2")
     (=. dz "z1 - z2")
     (=. r2 "dx*dx + dy*dy + dz*dz")
     (=. rr "sqrt(r2)")

     (stmt-cond
      ("rr<=r0_low"
       (=. energy "f_contact_i + f_contact_j")
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian)))

      ("rr<=r1_low"
       (=. energy "AA_low_i*rr*rr*rr + BB_low_i*rr*rr + CC_low_i*rr + DD_low_i + AA_low_j*rr*rr*rr + BB_low_j*rr*rr + CC_low_j*rr + DD_low_j")
       (=. dE_dr   "3*AA_low_i*rr*rr + 2*BB_low_i*rr + CC_low_i + 3*AA_low_j*rr*rr + 2*BB_low_j*rr + CC_low_j" :modes (:gradient :hessian))
       (=. d2E_dr2 "6*AA_low_i*rr + 2*BB_low_i + 6*AA_low_j*rr + 2*BB_low_j" :modes (:hessian)))

      ("rr<=r2_high"
       (=. inv_rr "1.0/rr")
       (=. inv_rr2 "inv_rr*inv_rr")

       (=. di "rr - lj_radius_i")
       (=. exp_i "exp(-(di*di)/lambda2_i)")
       (=. fi "pref_i_j * exp_i * inv_rr2")
       (=. hi "-2.0*(di/lambda2_i + inv_rr)")
       (=. dfi "fi*hi" :modes (:gradient :hessian))
       (=. d2fi "fi*(hi*hi + (-2.0/lambda2_i + 2.0*inv_rr2))" :modes (:hessian))

       (=. dj "rr - lj_radius_j")
       (=. exp_j "exp(-(dj*dj)/lambda2_j)")
       (=. fj "pref_j_i * exp_j * inv_rr2")
       (=. hj "-2.0*(dj/lambda2_j + inv_rr)")
       (=. dfj "fj*hj" :modes (:gradient :hessian))
       (=. d2fj "fj*(hj*hj + (-2.0/lambda2_j + 2.0*inv_rr2))" :modes (:hessian))

       (=. energy "fi + fj")
       (=. dE_dr   "dfi + dfj" :modes (:gradient :hessian))
       (=. d2E_dr2 "d2fi + d2fj" :modes (:hessian)))

      ("rr<=r3_high"
       (=. energy "AA_high_i*rr*rr*rr + BB_high_i*rr*rr + CC_high_i*rr + DD_high_i + AA_high_j*rr*rr*rr + BB_high_j*rr*rr + CC_high_j*rr + DD_high_j")
       (=. dE_dr   "3*AA_high_i*rr*rr + 2*BB_high_i*rr + CC_high_i + 3*AA_high_j*rr*rr + 2*BB_high_j*rr + CC_high_j" :modes (:gradient :hessian))
       (=. d2E_dr2 "6*AA_high_i*rr + 2*BB_high_i + 6*AA_high_j*rr + 2*BB_high_j" :modes (:hessian)))

      (t
       (=. energy "0.0")
       (=. dE_dr "0.0" :modes (:gradient :hessian))
       (=. d2E_dr2 "0.0" :modes (:hessian))))))

  (:derivatives
   (:mode :manual
    :intermediates (rr)

    :intermediate->coord
    ((rr ((x1 "dx / rr")
          (y1 "dy / rr")
          (z1 "dz / rr")
          (x2 "-dx / rr")
          (y2 "-dy / rr")
          (z2 "-dz / rr"))))

    :intermediate->coord2
    ((rr (
          ((x1 x1) "(rr*rr - dx*dx)/(rr^3)")
          ((x1 y1) "(-dx*dy)/(rr^3)")
          ((x1 z1) "(-dx*dz)/(rr^3)")
          ((x1 x2) "(dx*dx - rr*rr)/(rr^3)")
          ((x1 y2) "dx*dy/(rr^3)")
          ((x1 z2) "dx*dz/(rr^3)")

          ((y1 y1) "(rr*rr - dy*dy)/(rr^3)")
          ((y1 z1) "(-dy*dz)/(rr^3)")
          ((y1 x2) "dx*dy/(rr^3)")
          ((y1 y2) "(dy*dy - rr*rr)/(rr^3)")
          ((y1 z2) "dy*dz/(rr^3)")

          ((z1 z1) "(rr*rr - dz*dz)/(rr^3)")
          ((z1 x2) "dx*dz/(rr^3)")
          ((z1 y2) "dy*dz/(rr^3)")
          ((z1 z2) "(dz*dz - rr*rr)/(rr^3)")

          ((x2 x2) "(rr*rr - dx*dx)/(rr^3)")
          ((x2 y2) "(-dx*dy)/(rr^3)")
          ((x2 z2) "(-dx*dz)/(rr^3)")

          ((y2 y2) "(rr*rr - dy*dy)/(rr^3)")
          ((y2 z2) "(-dy*dz)/(rr^3)")

          ((z2 z2) "(rr*rr - dz*dz)/(rr^3)"))))

    :energy->intermediate
    (:gradient ((rr "dE_dr"))
     :hessian  (((rr rr) "d2E_dr2")))

    :hessian-modes ((rr :full))
    :geometry-check :warn)))

(build-kernel-group (kernels "amber_nonbond14" (:gradient :energy :hessian))
  (:pipeline *pipeline*)
  (:coordinate-inputs 2)
  (:term-struct
   (:name "AmberNonbond14Term")
   (:var term)
   (:inputs ((|double| A)
             (|double| C)
             (|double| kqq)
             ))
   (:struct ((|double| A)
             (|double| C)
             (|double| kqq)))
   (:derived ()))

  (:body
   (stmt-block :!base
     (=. dx "x1 - x2")
     (=. dy "y1 - y2")
     (=. dz "z1 - z2")
     (=. r2 "dx*dx + dy*dy + dz*dz")
     (=. r  "sqrt(r2)")
     (=. invr  "1.0 / r")
     (=. invr2 "invr*invr")
     (=. invr6 "invr2*invr2*invr2")

     ;; E(r) = A/r^12 - C/r^6 + kqq/r
     (=. e_lj   "A*invr6*invr6 - C*invr6")
     (=. e_coul "kqq*invr")
     (=. energy "e_lj + e_coul")))

  (:derivatives
   (:mode :manual
    :intermediates (r)

    ;; radial geometry as in stretch
    :intermediate->coord
    ((r ((x1 "dx / r")
         (y1 "dy / r")
         (z1 "dz / r")
         (x2 "-dx / r")
         (y2 "-dy / r")
         (z2 "-dz / r"))))

    :intermediate->coord2
    ((r (
         ((x1 x1) "(r*r - dx*dx)/(r^3)")
         ((x1 y1) "(-dx*dy)/(r^3)")
         ((x1 z1) "(-dx*dz)/(r^3)")
         ((x1 x2) "(dx*dx - r*r)/(r^3)")
         ((x1 y2) "dx*dy/(r^3)")
         ((x1 z2) "dx*dz/(r^3)")

         ((y1 y1) "(r*r - dy*dy)/(r^3)")
         ((y1 z1) "(-dy*dz)/(r^3)")
         ((y1 x2) "dx*dy/(r^3)")
         ((y1 y2) "(dy*dy - r*r)/(r^3)")
         ((y1 z2) "dy*dz/(r^3)")

         ((z1 z1) "(r*r - dz*dz)/(r^3)")
         ((z1 x2) "dx*dz/(r^3)")
         ((z1 y2) "dy*dz/(r^3)")
         ((z1 z2) "(dz*dz - r*r)/(r^3)")

         ((x2 x2) "(r*r - dx*dx)/(r^3)")
         ((x2 y2) "(-dx*dy)/(r^3)")
         ((x2 z2) "(-dx*dz)/(r^3)")

         ((y2 y2) "(r*r - dy*dy)/(r^3)")
         ((y2 z2) "(-dy*dz)/(r^3)")

         ((z2 z2) "(r*r - dz*dz)/(r^3)"))))

    ;; E(r) = A/r^12 - C/r^6 + kqq/r
    ;; dE/dr   = -12*A/r^13 + 6*C/r^7 - kqq/r^2
    ;; d²E/dr² = 156*A/r^14 - 42*C/r^8 + 2*kqq/r^3
    :energy->intermediate
    (:gradient ((r "-12*A/(r^13) + 6*C/(r^7) - kqq/(r^2)"))
     :hessian  (((r r) "156*A/(r^14) - 42*C/(r^8) + 2*kqq/(r^3)")))

    :hessian-modes ((r :full))
    :geometry-check :warn)))


(write-all kernels)

#+(or)
(progn


;;; Run all the tests
     #+(or)(stmt-ir.tests:run-all)

     (build-kernel-group (kernels "stretch" (:hessian :gradient :energy))
       (:pipeline *pipeline*)
       (:coordinate-inputs 2)

       (:term-struct
        (:name "stretch_term")
        (:var term)
        (:inputs ((|double| kb)
                  (|double| r0)))
        (:struct ((|double| kb)
                  (|double| r0))))

       (:body
        (stmt-block :!base
          (=. dx "x2 - x1")
          (=. dy "y2 - y1")
          (=. dz "z2 - z1")
          (=. r2 "dx*dx + dy*dy + dz*dz")
          (=. r  "sqrt(r2)")
          (=. dr "r - r0")
          ;; E(r) = kb * (r - r0)^2
          (=. energy "kb*dr*dr")))
       )


     (build-kernel-group (kernels "angle" (:energy :gradient :hessian))
       (:pipeline *pipeline*)
       (:coordinate-inputs 3)
       (:term-struct
        (:name "angle_term")
        (:var term)
        (:inputs ((|double| kt)
                  (|double| t0)))
        (:struct ((|double| kt)
                  (|double| t0))))

       (:body
        (stmt-block :!base
          ;; bond vectors about atom 2
          (=. vx1 "x1 - x2")
          (=. vy1 "y1 - y2")
          (=. vz1 "z1 - z2")

          (=. vx2 "x3 - x2")
          (=. vy2 "y3 - y2")
          (=. vz2 "z3 - z2")

          ;; dot products and norms
          (=. dot   "vx1*vx2 + vy1*vy2 + vz1*vz2")
          (=. n1_sq "vx1*vx1 + vy1*vy1 + vz1*vz1")
          (=. n2_sq "vx2*vx2 + vy2*vy2 + vz2*vz2")
          (=. n1    "sqrt(n1_sq)")
          (=. n2    "sqrt(n2_sq)")

          (=. cos_theta "dot / (n1*n2)")
          (raw-c "if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();")
          (=. sin_theta "sqrt(1 - cos_theta^2)")
          (=. theta     "acos(cos_theta)")
          (=. dtheta    "theta - t0")

          ;; E(theta) = kt * (theta - t0)^2  (kt includes the 1/2 if desired)
          (=! energy "kt*dtheta*dtheta")))
       )


     (defparameter *dihedral-pipeline*
       (stmt-ir:make-optimization-pipeline
        :name :dihedral-opt
        :optimizations
        (list
         (stmt-ir:make-optimization
          :name :linear-canonicalization
          :function 'stmt-ir:linear-canonicalization-optimization)
         (stmt-ir:make-optimization
          :name :inverse-expt
          :function 'stmt-ir:inverse-expt-optimization
          :keyword-args (list :min-uses 3))
         (stmt-ir:make-optimization
          :name :factor-sums
          :function 'stmt-ir:factor-sums-optimization
          :keyword-args (list :min-uses 3 :min-factors 1 :min-size 6))
         (stmt-ir:make-optimization
          :name :cse-full
          :function 'stmt-ir:cse-block-multi-optimization
          :keyword-args (list :max-passes 8 :min-uses 3 :min-size 4))
         (stmt-ir:make-optimization
          :name :copy-propagate
          :function 'stmt-ir:copy-propagate-optimization)
         (stmt-ir:make-optimization
          :name :alias-assigned-exprs
          :function 'stmt-ir:alias-assigned-exprs-optimization)
         (stmt-ir:make-optimization
          :name :normalize-signs
          :function 'stmt-ir:normalize-signs-optimization))))

     (build-kernel-group (kernels "dihedral_fast" (:energy :gradient :hessian))
       (:pipeline *dihedral-pipeline*)
       (:coordinate-inputs 4)

       (:term-struct
        (:name "dihedral_term")
        (:var term)
        (:inputs ((|double| V)
                  (|double| n)
                  (|double| sinPhase)
                  (|double| cosPhase)
                  ))
        (:struct ((|double| V)
                  (|double| n)
                  (|double| sinPhase)
                  (|double| cosPhase)
                  ))
        )
       (:body
        (stmt-block :!base
          ;; bond vectors
          (=. v1x "x2 - x1")
          (=. v1y "y2 - y1")
         (=. v1z "z2 - z1")

          (=. v2x "x3 - x2")
          (=. v2y "y3 - y2")
          (=. v2z "z3 - z2")

          (=. v3x "x4 - x3")
          (=. v3y "y4 - y3")
          (=. v3z "z4 - z3")

          ;; plane normals
          (=. c1x "v1y*v2z - v1z*v2y")
          (=. c1y "v1z*v2x - v1x*v2z")
          (=. c1z "v1x*v2y - v1y*v2x")

          (=. c2x "v2y*v3z - v2z*v3y")
          (=. c2y "v2z*v3x - v2x*v3z")
          (=. c2z "v2x*v3y - v2y*v3x")

          ;; norms and dot products
          (=. v2_sq "v2x*v2x + v2y*v2y + v2z*v2z")
          (=. v2_len "sqrt(v2_sq)")

          ;; torsion via atan2
          (=. t1 "v2_len*(v1x*c2x + v1y*c2y + v1z*c2z)")
          (=. t2 "c1x*c2x + c1y*c2y + c1z*c2z")
          (=. phi "atan2(t1, t2)")

          ;; Amber-style torsion:
          ;; E(phi) = V * (1 + cos(n*phi - phase))
          ;; Use sinPhase/cosPhase inputs instead of the raw phase.
          (=. nphi      "n*phi")
          (=. sin_nphi  "sin(nphi)")
          (=. cos_nphi  "cos(nphi)")
          (=. cos_angle "cos_nphi*cosPhase + sin_nphi*sinPhase") ;; cos(nphi - phase)
          (=. sin_angle "sin_nphi*cosPhase - cos_nphi*sinPhase") ;; sin(nphi - phase)
          (=! energy "V*(1.0 + cos_angle)")))

       (:derivatives
        (:mode :manual
         :intermediates (phi)

         ;; no manual geometry; AD provides dphi/dcoord

         ;; E(phi) = V*(1 + cos_angle), where
         ;; cos_angle = cos(nphi - phase) = cos_nphi*cosPhase + sin_nphi*sinPhase
         ;; sin_angle = sin(nphi - phase) = sin_nphi*cosPhase - cos_nphi*sinPhase
         ;; dE/dphi  = -V*n*sin_angle
         ;; d²E/dphi² = -V*n*n*cos_angle
         :energy->intermediate
         (:gradient ((phi "-V*n*sin_angle"))
          :hessian  (((phi phi) "-V*n*n*cos_angle")))

         :hessian-modes ((phi :full)) ;; :outer-product-only))
         :geometry-check :warn)))


       (build-kernel-group (kernels "nonbond_dd_cutoff" (:gradient :energy :hessian))
         (:pipeline *pipeline*)

         (:template "template <typename HESSIAN>")

         (:params ((|double| A) ;; LJ A coefficient  (A / r^12)
                   (|double| B) ;; LJ B coefficient  (B / r^6)
                   (|double| qq) ;; Coulomb prefactor
                   (|double| invdd) ;; 1.0/epsilon(r) = 1.0/dd*r
                   (|double| rswitch) ;; switching start
                   (|double| rswitch2) ;; switching start squared
                   (|double| rcut)     ;; cutoff
                   (|double| rcut2)    ;; cutoff^2
                   (|double| inv_range) ;; 1.0/(rcut - rswitch)
                   (|size_t| i3x1)
                   (|size_t| i3x2)
                   (|double*| position)
                   (|double*| energy_accumulate)
                   (|double*| force)
                   (HESSIAN hessian)
                   (|double*| dvec)
                   (|double*| hdvec)))

         (:layout ((1 . I3X1) (2 . I3X2))
                  ((X . 0) (Y . 1) (Z . 2)))

         (:coord-vars (x1 y1 z1
                          x2 y2 z2))

         (:coord-load
          (coords-from-position
           ((x1 y1 z1 i3x1)
            (x2 y2 z2 i3x2))))

         (:body
          (stmt-block :!base
            (=. dx "x1 - x2")
            (=. dy "y1 - y2")
            (=. dz "z1 - z2")
            (=. r2 "dx*dx + dy*dy + dz*dz")
            (stmt-ir:make-if-stmt
             (expr-ir:parse-expr "r2 < rcut2")
             (stmt-block :!cut
               (=. invr  "r2^-0.5")     ; one pow/sqrt
               (=. r     "r2*invr")     ; reuse, no extra sqrt
               (=. invr2 "invr*invr")
               (=. invr3 "invr*invr2")
               (=. invr4 "invr2*invr2")
               (=. invr6 "invr2*invr2*invr2")
               (=. e_lj   "A*invr6*invr6 - B*invr6")
               (=. e_coul "qq*invr2*invdd")
               (=. e_base "e_lj + e_coul")
               ;; base radial derivatives
               (=. dE_base_dr   "(-12.0*A*invr6*invr6*invr) + (6.0*B*invr6*invr) + (-2.0*qq*invr3*invdd)"
                   :modes (:gradient :hessian))
               (=. d2E_base_dr2 "(156.0*A*invr6*invr6*invr2) + (-42.0*B*invr6*invr2) + (6.0*qq*invr2*invr2*invdd)"
                   :modes (:hessian))
               (stmt-ir:make-if-stmt
                (expr-ir:parse-expr "r2 < rswitch2")
                (stmt-block :!switch
                  (=! energy "e_base")
                  (=. dE_dr   "dE_base_dr" :modes (:gradient :hessian))
                  (=. d2E_dr2 "d2E_base_dr2" :modes (:hessian)))
                (stmt-block :!skin
                  (=. drs "r - rswitch")
                  (=. t1 "drs*inv_range")
                  (=. t2 "t1*t1")
                  (=. t3 "t2*t1")
                  (=. t4 "t2*t2")
                  (=. t5 "t3*t2")
                  (=. s "1.0 - 10.0*t3 + 15.0*t4 - 6.0*t5")
                  (=. ds_dt "(-30.0*t2 + 60.0*t3 - 30.0*t4)" :modes (:gradient :hessian))
                  (=. d2s_dt2 "(-60.0*t1 + 180.0*t2 - 120.0*t3)" :modes (:hessian))
                  (=. ds_dr "ds_dt*inv_range" :modes (:gradient :hessian))
                  (=. d2s_dr2 "d2s_dt2*inv_range*inv_range" :modes (:hessian))
                  (=! energy "s*e_base")
                  (=. dE_dr   "(s*dE_base_dr) + (ds_dr*e_base)" :modes (:gradient :hessian))
                  (=. d2E_dr2 "(s*d2E_base_dr2) + (2.0*ds_dr*dE_base_dr) + (d2s_dr2*e_base)" :modes (:hessian)))))))

          )

         (:derivatives
          (:mode :manual
           :intermediates (r)

           ;; radial geometry
           :intermediate->coord
           ((r ((x1 "dx / r")
                (y1 "dy / r")
                (z1 "dz / r")
                (x2 "-dx / r")
                (y2 "-dy / r")
                (z2 "-dz / r"))))

           :intermediate->coord2
           ((r (
                ((x1 x1) "(r*r - dx*dx)/(r^3)")
                ((x1 y1) "(-dx*dy)/(r^3)")
                ((x1 z1) "(-dx*dz)/(r^3)")
                ((x1 x2) "(dx*dx - r*r)/(r^3)")
                ((x1 y2) "dx*dy/(r^3)")
                ((x1 z2) "dx*dz/(r^3)")

                ((y1 y1) "(r*r - dy*dy)/(r^3)")
                ((y1 z1) "(-dy*dz)/(r^3)")
                ((y1 x2) "dx*dy/(r^3)")
                ((y1 y2) "(dy*dy - r*r)/(r^3)")
                ((y1 z2) "dy*dz/(r^3)")

                ((z1 z1) "(r*r - dz*dz)/(r^3)")
                ((z1 x2) "dx*dz/(r^3)")
                ((z1 y2) "dy*dz/(r^3)")
                ((z1 z2) "(dz*dz - r*r)/(r^3)")

                ((x2 x2) "(r*r - dx*dx)/(r^3)")
                ((x2 y2) "(-dx*dy)/(r^3)")
                ((x2 z2) "(-dx*dz)/(r^3)")

                ((y2 y2) "(r*r - dy*dy)/(r^3)")
                ((y2 z2) "(-dy*dz)/(r^3)")

                ((z2 z2) "(r*r - dz*dz)/(r^3)"))))

           ;; radial link energy ↔ r via branchwise dE_dr, d2E_dr2
           :energy->intermediate
           (:gradient ((r "dE_dr"))
            :hessian  (((r r) "d2E_dr2")))

           :hessian-modes ((r :full))
           :geometry-check :warn)))




       (build-kernel-group (kernels "chiral_restraint" (:energy :gradient :hessian))
         ;; Optimization pipeline
         (:pipeline *pipeline*)

         (:template "template <typename HESSIAN>")

         ;; Parameters: (ctype name)
         ;;   K   – force constant
         ;;   CO  – target chiral value
         (:params ((|double| K)
                   (|double| CO)
                   (|size_t| i3x1)
                   (|size_t| i3x2)
                   (|size_t| i3x3)
                   (|size_t| i3x4)
                   (|double*| position)
                   (|double*| energy_accumulate)
                   (|double*| force)
                   (HESSIAN hessian)
                   (|double*| dvec)
                   (|double*| hdvec)
                   ))

         ;; Layout: 4 atoms with 3D coords
         (:layout ((1 . I3X1) (2 . I3X2) (3 . I3X3) (4 . I3X4))
                  ((#\X . 0) (#\Y . 1) (#\Z . 2)))

         ;; Coordinate vars in expr-var package
         (:coord-vars (x1 y1 z1
                          x2 y2 z2
                          x3 y3 z3
                          x4 y4 z4))

         ;; How to load coords from position[]
         (:coord-load
          (coords-from-position
           ((x1 y1 z1 I3X1)
            (x2 y2 z2 I3X2)
            (x3 y3 z3 I3X3)
            (x4 y4 z4 I3X4))))

         ;; Manual derivative spec: we introduce one intermediate Q
         ;; and give dE/dQ. We deliberately DROP d²E/dQ² so that this
         ;; channel contributes only the outer‑product term to the Hessian.
         (:derivatives
          (:mode :hybrid
           :intermediates (Q)

           ;; No manual geometry; let AD fill:
           ;;   du/dq = dQ/dq, d²u/dqidqj = d²Q/dqidqj
           :intermediate->coord ()
           :intermediate->coord2 ()

           ;; Energy derivatives w.r.t Q:
           ;;   E = K*(CO + Q)^3
           ;;   dE/dQ   = 3*K*(CO + Q)^2
           ;;   We omit d²E/dQ² so that only the
           ;;   outer-product term (dE/dQ * dQ/dqi dQ/dqj) is kept.
           :energy->intermediate
           (:gradient ((Q "3*K*(CO + Q)^2"))
            :hessian  ())

           ;; Only outer‑product Hessian contribution from this channel
           :hessian-modes ((Q :outer-product-only))

           ;; Check geometry (AD vs auto-filled) and warn on mismatch
           :geometry-check :warn))

         ;; Kernel body
         (:body
          (stmt-block :!base
            ;; r13 = r1 - r3
            (=. dx13 "x1 - x3")
            (=. dy13 "y1 - y3")
            (=. dz13 "z1 - z3")

            ;; r23 = r2 - r3
            (=. dx23 "x2 - x3")
            (=. dy23 "y2 - y3")
            (=. dz23 "z2 - z3")

            ;; r43 = r4 - r3
            (=. dx43 "x4 - x3")
            (=. dy43 "y4 - y3")
            (=. dz43 "z4 - z3")

            ;; norms
            (=. r13_2 "dx13*dx13 + dy13*dy13 + dz13*dz13")
            (=. r23_2 "dx23*dx23 + dy23*dy23 + dz23*dz23")
            (=. r43_2 "dx43*dx43 + dy43*dy43 + dz43*dz43")

            (=. r13 "sqrt(r13_2)")
            (=. r23 "sqrt(r23_2)")
            (=. r43 "sqrt(r43_2)")

            ;; cross = r13 × r23
            (=. cx "dy13*dz23 - dz13*dy23")
            (=. cy "dz13*dx23 - dx13*dz23")
            (=. cz "dx13*dy23 - dy13*dx23")

            ;; triple product V = (r13 × r23) · r43
            (=. V "cx*dx43 + cy*dy43 + cz*dz43")

            ;; normalized chiral measure Q = V / (|r13| |r23| |r43|)
            (=. denom "r13*r23*r43")
            (=. Q "V/denom")

            ;; Energy: E = K * (CO + Q)^3
            (=! ENERGY "K * (CO + Q)^3")

            ;; Accumulate E/G/H at this point
            )))


       (build-kernel-group (kernels "anchor" (:energy :gradient :hessian))
         (:pipeline *pipeline*)
         (:layout ((1 . I3X1))
                  ((X . 0) (Y . 1) (Z . 2)))
         (:coord-vars (x1 y1 z1))
         (:coord-load
          (coords-from-position
           ((x1 y1 z1 i3x1))))
         (:body
          (stmt-block :!base
            (=. dx "x1 - xa")
            (=. dy "y1 - ya")
            (=. dz "z1 - za")
            (=. r2 "dx*dx + dy*dy + dz*dz")
            (=! energy "ka * r2")
            ))
         (:template "template <typename HESSIAN>")
         (:params ((|double| ka)
                   (|double| xa)
                   (|double| ya)
                   (|double| za)
                   (|size_t| i3x1)
                   (|double*| position)
                   (|double*| energy_accumulate)
                   (|double*| force)
                   (HESSIAN hessian)
                   (|double*| dvec)
                   (|double*| hdvec)
                   ))
         (:derivatives
          (:mode :manual
           :intermediates (dx dy dz r2)
           :intermediate->coord
           ((dx ((x1 "1")))
            (dy ((y1 "1")))
            (dz ((z1 "1")))
            (r2 ((x1 "2*dx")
                 (y1 "2*dy")
                 (z1 "2*dz"))))
           :intermediate->coord2
           ((dx (((x1 x1) "0")))
            (dy (((y1 y1) "0")))
            (dz (((z1 z1) "0")))
            (r2 (((x1 x1) "2")
                 ((y1 y1) "2")
                 ((z1 z1) "2"))))
           :energy->intermediate
           (:gradient
            ((r2 "ka"))
            :hessian
            (((r2 r2) "0")))
           :hessian-modes
           ((dx :none)
            (dy :none)
            (dz :none)
            (r2 :outer-product-only))
           :geometry-check :warn)))





       (build-kernel-group (kernels "dihedral_restraint_fast" (:energy :gradient :hessian))
         (:pipeline *dihedral-pipeline*)

         ;; Parameters:
         ;;   kdh  – force constant
         ;;   phi0 – target dihedral angle (radians)

         (:template "template <typename HESSIAN>")
         (:params ((|double| kdh)
                   (|double| phi0)
                   (|size_t| i3x1)
                   (|size_t| i3x2)
                   (|size_t| i3x3)
                   (|size_t| i3x4)
                   (|double*| position)
                   (|double*| energy_accumulate)
                   (|double*| force)
                   (HESSIAN hessian)
                   (|double*| dvec)
                   (|double*| hdvec)
                   ))

         ;; 4-atom dihedral 1–2–3–4
         (:layout ((1 . I3X1) (2 . I3X2) (3 . I3X3) (4 . I3X4))
                  ((#\X . 0) (#\Y . 1) (#\Z . 2)))

         (:coord-vars (x1 y1 z1
                          x2 y2 z2
                          x3 y3 z3
                          x4 y4 z4))

         (:coord-load
          (coords-from-position
           ((x1 y1 z1 I3X1)
            (x2 y2 z2 I3X2)
            (x3 y3 z3 I3X3)
            (x4 y4 z4 I3X4))))

         ;; Manual derivative spec:
         ;;   intermediates: phi, deltaPhi, deltaPhiModFn
         ;;   E = kdh * deltaPhiModFn^2
         ;;   dE/d(deltaPhiModFn)   = 2*kdh*deltaPhiModFn
         ;;   d²E/d(deltaPhiModFn)² = 2*kdh
         ;; Geometry of phi/deltaPhi/deltaPhiModFn is left to AD+auto-fill.
         ;; For truncated Newton, we keep only the dominant outer‑product
         ;; Hessian terms from deltaPhiModFn and drop curvature terms from
         ;; phi and deltaPhi.
         (:derivatives
          (:mode :hybrid
           :intermediates (phi deltaPhi deltaPhiModFn)

           :intermediate->coord ()
           :intermediate->coord2 ()

           :energy->intermediate
           (:gradient ((deltaPhiModFn "2*kdh*deltaPhiModFn"))
            :hessian  (((deltaPhiModFn deltaPhiModFn) "2*kdh")))

           ;; Only deltaPhiModFn contributes to the Hessian; phi and deltaPhi
           ;; are treated as having no second‑derivative contribution.
           :hessian-modes ((phi :none)
                           (deltaPhi :none)
                           (deltaPhiModFn :outer-product-only))

           :geometry-check :warn))

         (:body
          (stmt-block :!base
            ;; Bond vectors
            (=. dx12 "x2 - x1")
            (=. dy12 "y2 - y1")
            (=. dz12 "z2 - z1")

            (=. dx23 "x3 - x2")
            (=. dy23 "y3 - y2")
            (=. dz23 "z3 - z2")

            (=. dx34 "x4 - x3")
            (=. dy34 "y4 - y3")
            (=. dz34 "z4 - z3")

            ;; Planes normals: n1 = (r12 × r23), n2 = (r23 × r34)
            (=. n1x "dy12*dz23 - dz12*dy23")
            (=. n1y "dz12*dx23 - dx12*dz23")
            (=. n1z "dx12*dy23 - dy12*dx23")

            (=. n2x "dy23*dz34 - dz23*dy34")
            (=. n2y "dz23*dx34 - dx23*dz34")
            (=. n2z "dx23*dy34 - dy23*dx34")

            ;; Norms of normals
            (=. n1_2 "n1x*n1x + n1y*n1y + n1z*n1z")
            (=. n2_2 "n2x*n2x + n2y*n2y + n2z*n2z")

            ;; Cos(phi) = (n1 · n2) / (|n1|n2|)

            (=. inv_n1n2 "(n1_2*n2_2)^-0.5")
            (=. dot_n1n2 "n1x*n2x + n1y*n2y + n1z*n2z")
            (=. cosphi "dot_n1n2 * inv_n1n2")

            ;; phi in [0,pi]
            (=. phi "acos(cosphi)")

            ;; raw delta = phi - phi0
            (=. deltaPhi "phi - phi0")

            ;; For now, deltaPhiModFn == deltaPhi; wrapping is done externally
            (=. deltaPhiModFn "deltaPhi")

            ;; Harmonic restraint energy
            (=! ENERGY "kdh*deltaPhiModFn*deltaPhiModFn")

            )))

       (write-all kernels)

       #+(or)
       (mathkernel:emit-c-tests (loop for kernel in kernels
                                      when (member (mathkernel::kernel-group kernel) '("stretch" "angle" "dihedral") :test 'string=)
                                        collect kernel)
                                "/home/meister/Development/cando/extensions/cando/include/cando/chem/energyKernels/"
                                "~/tmp/tests/position.c"
                                "~/tmp/tests/calls.c"
                                "~/tmp/tests/all-tests.c")
       )
