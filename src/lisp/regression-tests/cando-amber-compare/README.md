# Cando / exported Amber topology comparison

Load the comparison code (this loads `leap` and `spiros/spiros`, but does not run
the comparison):

```lisp
(load "~/Development/cando/extensions/cando/src/lisp/regression-tests/cando-amber-compare/compare.lisp")
(cando-amber-compare:compare-ligand)
(cando-amber-compare:compare-receptor)
(cando-amber-compare:compare-complex)
```

This rebuilds an Amber-component energy function from `data/ligand.cando`,
preserving its `:smirnoff-spiros` molecule force-field assignment. It separately
reads `data/ligand.prmtop` using Cando's existing topology reader. Both use the
same first-frame coordinates saved in the aggregate, after checking atom counts,
residue boundaries, atom names, elements and coordinate offsets.

The receptor and complex wrappers use `receptor.cando`/`receptor.prmtop` and
`complex.cando`/`complex.prmtop`, respectively. All three share the same comparison
implementation, keyword options, and five return values. Each molecule retains
its saved force-field assignment, including the mixed-force-field complex.
Use `:data-directory` to point at a directory containing freshly generated files.

The report compares BOND, ANGLE, DIHED, VDWAALS, EEL, VDW14, EEL14, OTHER and
TOTAL in kcal/mol. Restraints are disabled; dielectric and nonbond scales are
one. GB/SA are not calculated. All rows print before any mismatch signals.
Default tolerance is `1d-4 + 1d-6 * max(abs(left), abs(right))` per row.

For inspection without signaling on differences:

```lisp
(multiple-value-bind (rows rebuilt from-prmtop raw-rebuilt raw-prmtop)
    (cando-amber-compare:compare-ligand :error-on-mismatch nil)
  (values rows rebuilt from-prmtop raw-rebuilt raw-prmtop))
```

This is a comparison of two **Cando evaluations**, not a Sander comparison.
The regenerated energy function also depends on the currently loaded force-field
data; the aggregate does not embed the force-field database. No fixture files
are modified. External Amber comparisons are not yet wired in.

The local `#+(or)` disabling mismatch signaling has been preserved; currently
`:error-on-mismatch` does not cause failures to signal. PASS/FAIL rows are still
calculated and returned for all three species.

## Dump dihedral terms

To write both components automatically during the comparison:

```lisp
(cando-amber-compare:compare-ligand :dihedral-dump-directory #P"/tmp/")
(cando-amber-compare:compare-receptor :dihedral-dump-directory #P"/tmp/")
(cando-amber-compare:compare-complex :dihedral-dump-directory #P"/tmp/")
```

This writes `ligand-rebuilt-dihedrals.lisp` and `ligand-prmtop-dihedrals.lisp`
in the requested directory, creating it if necessary and overwriting previous
dumps. Both are written before energy evaluation, so a subsequent energy
mismatch error does not prevent the dumps. Omit the option (or pass `NIL`)
to leave files untouched. Return values are unchanged.
Receptor and complex dumps use the corresponding species prefix, so their files
do not overwrite the ligand dumps.

`dump-dihedral-terms` accepts the dihedral component, not the whole energy
function. It writes one readable plist per Fourier term, including proper/improper
status, one-based atom indices, names, zero-based coordinate offsets, amplitude,
periodicity, phase in radians, and the actual cached cosine/sine used in evaluation.
The differences between the cached values and `cos(phase)`/`sin(phase)` are also
printed. Proper terms are canonicalized by reversing all four atom indices,
names, and offsets when the first index exceeds the last; improper atom order
is preserved. Rows are sorted lexicographically by the four atom indices, then
numerically by periodicity. Equal sort keys retain their original relative order.
Terms are not evaluated or deduplicated, and the component is not modified.
The return value is
the term count. Atom names can repeat; use atom indices to match terms.

After reloading `compare.lisp`, dump both sides (these explicit filenames are
overwritten if they already exist):

```lisp
(multiple-value-bind (rows rebuilt from-prmtop)
    (cando-amber-compare:compare-ligand :error-on-mismatch nil)
  (declare (ignore rows))
  (with-open-file (out "/tmp/ligand-rebuilt-dihedrals.lisp"
                       :direction :output :if-exists :supersede)
    (cando-amber-compare:dump-dihedral-terms
     (chem:get-dihedral-component rebuilt) :stream out))
  (with-open-file (out "/tmp/ligand-prmtop-dihedrals.lisp"
                       :direction :output :if-exists :supersede)
    (cando-amber-compare:dump-dihedral-terms
     (chem:get-dihedral-component from-prmtop) :stream out)))
```

Term IDs are omitted to simplify comparison. The dumped `:v` amplitude is
rounded to three decimal places; the energy component itself is unchanged.
