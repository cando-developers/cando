
;;Term parameters deltaGfreei, Vj, lambdai, sigmai - they are provided in the term-struct
;; Pi is the physical constant PI (3.14159265358979...)
;; c0, c1 are constants that will be defined in the params-struct
;; f_solv_low_poly and f_solv_high_poly are cubic polynomials with form "aa*x^3+bb*x^2+cc*x+dd"
;;; Here is the text from the paper "The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design":

"""Solvation. Native-like protein conformations minimize the
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
"""

sigmai = "lj_radius_i"

sigmaij = "lj_radius_i + lj_radius_j"

(eq 11) fdesolv(i,j,dij) = "-(1/(2*pi^(3/2))*lk_dgfree_i*(1/lk_lambda_i)*lk_volume_j*(1/(dij^2)))*exp(-((dij-lj_radius_i)/lk_lambda_i)^2)"

(eq 12) gdesolv(i,j) =
(stmt-cond
 ("dij<=(sigmaij - c0)"
  "fdesolv(i,j,sigmaij)")
 ("dij<=(sigmaij + c1)"
  "f_solv_low_poly(i,j,dij)")
 ("dij<=4.5"
  "fdesolv(i,j,dij)")
 ("dij<=6.0"
  "f_solv_high_poly(i,j,dij)")
 (t 0.0))

(eq 13) Efa_sol = "gdesolv(i,j) + gdesolv(j,i)"

