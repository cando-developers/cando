#!/bin/sh
#
# Regenerate the cando energy kernel C sources from the mathkernel DSL.
#
# This is NOT part of the cando build.  The generated .c/.h files under
# include/cando/chem/energyKernels/ are checked in and treated as source;
# ninja only copies them.  Run this by hand -- via `make kernels` from the
# top of the tree -- when you change a kernel definition in
# energy-kernels.lisp or change the mathkernel generator, then review the
# diff, commit it, and rebuild normally.
#
# It deliberately runs in SBCL rather than cando: the generator produces the
# sources cando is compiled from, so running it inside cando would be a
# bootstrap loop.
#
# Usage:
#   regenerate-kernels.sh            regenerate in place
#   regenerate-kernels.sh --check    regenerate to a temp dir and diff against
#                                    the committed files; exit 1 if they differ
#
# Environment:
#   MATHKERNEL_DIR   location of the mathkernel repo
#                    (default: $HOME/common-lisp/tlt/mathkernel)
#   SBCL             sbcl binary to use (default: sbcl)

set -eu

here=$(cd "$(dirname "$0")" && pwd)
target="$here/../../../include/cando/chem/energyKernels"
: "${MATHKERNEL_DIR:=$HOME/common-lisp/tlt/mathkernel}"
: "${SBCL:=sbcl}"

check=no
case "${1:-}" in
  --check) check=yes ;;
  "")      ;;
  *)       echo "usage: $(basename "$0") [--check]" >&2; exit 64 ;;
esac

if [ ! -d "$MATHKERNEL_DIR" ]; then
  echo "regenerate-kernels: mathkernel not found at $MATHKERNEL_DIR" >&2
  echo "  set MATHKERNEL_DIR to its location (see kernels/README.md)" >&2
  exit 2
fi

if ! command -v "$SBCL" >/dev/null 2>&1; then
  echo "regenerate-kernels: $SBCL not found on PATH" >&2
  exit 2
fi

if [ ! -d "$target" ]; then
  echo "regenerate-kernels: target directory not found: $target" >&2
  exit 2
fi
target=$(cd "$target" && pwd)

if [ "$check" = yes ]; then
  outdir=$(mktemp -d "${TMPDIR:-/tmp}/cando-kernels.XXXXXX")
  trap 'rm -rf "$outdir"' EXIT
else
  outdir="$target"
fi

# write-all reads KERNEL_PATH; it needs a trailing slash to be treated as a
# directory by merge-pathnames.
MATHKERNEL_DIR="$MATHKERNEL_DIR" \
KERNELS_DIR="$here/" \
KERNEL_PATH="$outdir/" \
  "$SBCL" --dynamic-space-size 8192 --noinform --non-interactive \
          --load "$here/regenerate-kernels.lisp"

if [ "$check" = no ]; then
  echo "regenerate-kernels: wrote to $target"
  echo "  review with: git -C \"$target\" diff ."
  exit 0
fi

# --check: compare only the files the generator produced.  The committed
# directory also holds hand-written and legacy kernels that this generator
# does not emit, so a whole-directory comparison would be misleading.
status=0
for f in "$outdir"/*; do
  base=$(basename "$f")
  if [ ! -f "$target/$base" ]; then
    echo "check: MISSING from tree: $base"
    status=1
  elif ! cmp -s "$f" "$target/$base"; then
    echo "check: OUT OF DATE: $base"
    status=1
  fi
done

if [ "$status" -eq 0 ]; then
  echo "check: generated kernels are up to date"
else
  echo "check: run 'make kernels' and commit the result" >&2
fi
exit "$status"