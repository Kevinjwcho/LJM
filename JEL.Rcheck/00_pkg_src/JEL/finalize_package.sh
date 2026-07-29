#!/usr/bin/env bash
# ======================================================================
# finalize_package.sh  — run from Package_working/JEL_2.1
#   1. regenerate man/ + NAMESPACE from roxygen (picks up the new @param)
#   2. R CMD check   3. install
# ======================================================================
set -e
PKG_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PKG_DIR/.."

# 1. regenerate documentation (updates man/*.Rd usage with update_c, h_init,
#    min_window_obs, h, y_vars, Bs) and NAMESPACE
Rscript -e 'roxygen2::roxygenise("JEL_2.1")'

# 2. build + check
R CMD build JEL_2.1
TARBALL="$(ls -t JEL_*.tar.gz | head -1)"
R CMD check "$TARBALL"            # add --as-cran for a stricter pass

# 3. install (so library(JEL) in the simulations picks up the c-update build)
R CMD INSTALL "$TARBALL"

echo "DONE. Review the check log:  ${TARBALL%.tar.gz}.Rcheck/00check.log"
echo "Run tests only:  Rscript -e 'devtools::test(\"JEL_2.1\")'"
