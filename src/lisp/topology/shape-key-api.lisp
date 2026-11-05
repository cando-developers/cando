(in-package #:topology)

;;; Foldamer systems specialize these operations on their own classes.
;;; Topology owns the contract so geometry code need not depend on foldamer.

(defgeneric bin-shape-component (foldamer kind parameters)
  (:documentation
   "Convert physical parameters into a list of integer bin indices.
KIND is :DIHEDRAL, :RING5, or :RING6. PARAMETERS is respectively
(angle), (q phi), or (q theta phi). All angles are in radians;
amplitudes use the molecule's coordinate units. Methods validate the
component kind and values and apply the foldamer's binning policy."))

(defgeneric unbin-shape-component (foldamer kind bins)
  (:documentation
   "Decode canonical integer bins into representative physical parameters.
Return (angle), (q phi), or (q theta phi) for :DIHEDRAL, :RING5,
or :RING6. All angles are in radians. Use the foldamer policy that
produced BINS. Binning the result must reproduce BINS; the original
unrounded measurements cannot be recovered."))
