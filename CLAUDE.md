# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Scope of this file

This is the **cando** extension subtree. The enclosing Clasp tree has its own `CLAUDE.md` at `../../CLAUDE.md` covering the Clasp build system (koga + ninja), GC variants, the static analyzer / `.sif` descriptors, and how cando plugs in as an extension. Read that file for the big picture; this file covers what's specific to the cando subtree itself.

## Build convention

**Do not run `ninja`, `make`, or `./koga` yourself. The user drives builds and reports results.** This is a project-wide rule (also stated in the parent `CLAUDE.md`); the practical implication for work inside this subtree is:

- After editing `.cc`/`.h` files in `src/` or `include/`, just stop and tell the user what changed — they'll rebuild.
- After **adding or removing** a `.cc`/`.h` file, you must update the relevant `cscript.lisp` (see below); flag this in your response so the user knows to run `./koga --deep-clean`.

## cscript.lisp — the source manifest

cando's source files are not discovered automatically; they're listed in `cscript.lisp` files that the koga build driver reads. There is a tree of them:

- `./cscript.lisp` (this directory) — top-level: declares the `:cando-user` system and the install scripts.
- `src/cscript.lisp` — recurses into the source subdirectories (`adapt/`, `geom/`, `units/`, `chem/`, `kinematics/`, `lisp/`, `main/`, `data/`, `examples/`) and declares the required `expat` library.
- `src/chem/cscript.lisp`, `src/kinematics/cscript.lisp`, etc. — list every `.cc` file in that directory that gets compiled into libclasp.

When adding a new C++ class:
1. Header goes under `include/cando/<subsystem>/`.
2. Source goes under `src/<subsystem>/`.
3. Add the new `.cc` to the matching `src/<subsystem>/cscript.lisp` (look at the existing `(k:sources :libclasp #~"foo.cc" ...)` block).
4. Tell the user to run `./koga --deep-clean` — without that, the new file won't be picked up.

The `k:recurse`, `k:sources`, `k:systems`, and `k:library` helpers are koga primitives; the parent Clasp tree's koga script defines them.

## Code layout

- `include/cando/` — C++ public headers, organized by subsystem (`chem/`, `geom/`, `kinematics/`, `adapt/`, `units/`, `energy-functions/`, `main/`). Mirrors `src/`.
- `src/` — C++ sources:
  - `chem/` — the bulk of cando: atoms, bonds, residues, aggregates, forcefield databases (`ffStretchDb`, `ffNonbondDb`, `ffAngleDb`, `ffItorDb`, `ffPtorDb`, `ffVdwDb`, `ffTypesDb`), energy components (`energyAngle`, `energyDihedral`, `energyAtomTable`, etc.), the minimizer, SMARTS/chemInfo matching, AMBER and Rosetta forcefield glue, topology, conformation explorer.
  - `geom/`, `kinematics/`, `units/`, `adapt/`, `main/`, `data/`, `examples/`, `analysis/`, `common/`, `mathematica/`, `openmmPackage/`, `tests/` — supporting C++.
  - `lisp/` — Lisp side of the extension. Each subdirectory is typically an ASDF system. Important ones:
    - `cando/` (system `:cando`) — the main cando Lisp library: dynamics, energy, geometry, modelling, SMARTS, SMILES, charges, restraints, topology helpers.
    - `cando-user/` (system `:cando-user`) — top-level user image; `extension.lisp` wires up AMBER/Rosetta paths, lparallel threads, and the `cando-user` package as the default REPL package.
    - `leap/`, `tirun/`, `topology/`, `smarts/`, `smiles/`, `lomap/`, `rosetta/`, `psi4/`, `xtal/`, `sketch2d/`, `structure-editor/`, `cando-jupyter/`, `tirun-jupyter/`, `regression-tests/`, etc. — feature systems.
- `chem-lisp/` — a single auxiliary Lisp file (`print-read.lisp`).
- `docker/`, `documentation/`, `patch/`, `sandbox/` — packaging, docs, patches, and developer scratch.
- `cando-user-install`, `cando-apptainer-install`, `scando-apptainer-install`, `scando-zeus-install` — installer shell scripts (also listed in the top-level `cscript.lisp` under `:install-bin`).

## Runtime environment knobs

`src/lisp/cando-user/extension.lisp` reads these at startup:

- `AMBERHOME` / `AMBER_HOME` — root for AMBER tools; sets the `AMBER:` logical host. Falls back to `/opt/amber/`.
- `ROSETTA_HOME` — root for Rosetta; sets the `ROSETTA:` logical host. Falls back to `/opt/rosetta/`.
- `CANDO_THREADS` — integer; sizes the lparallel kernel. Defaults to `sys:num-logical-processors`.
- `:no-auto-lparallel` in `*features*` — skip lparallel kernel creation.

Profiling: `CLASP_FLAME_PROFILE` (parent Clasp tree) — set to `1` or `path=...:duration=...:rate=...`, then send `SIGUSR2` to dump a flame-graph SVG. See `../../src/lisp/kernel/lsp/flamegraph.lisp`.

## Style

`.clang-format` here is LLVM-based with `ColumnLimit: 132`, `SortIncludes: false`, and a handful of whitespace-sensitive Clasp macros (`CL_LAMBDA`, `CL_LISPIFY_NAME`, `CL_DECLARE`, `CL_NAME`) declared so they aren't reformatted. C++ style is LLVM; Lisp style follows whatever is already in the file you're editing.
