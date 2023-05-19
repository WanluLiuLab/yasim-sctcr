#!/usr/bin/env bash
set -ue
SHDIR="$(dirname "$(readlink -f "${0}")")"
cd "${SHDIR}/" || exit 1
rm -rf "${SHDIR}/git_dist"
mkdir -p "${SHDIR}/git_dist"

git clone --mirror "$(git remote get-url origin)" "${SHDIR}/git_dist/yasim-sctcr.git"

for dist in labw_utils yasim; do
    cd "${SHDIR}/deps/${dist}" || exit 1
    git clone --mirror "$(git remote get-url origin)" "${SHDIR}/git_dist/${dist}.git"
    cd "${SHDIR}/" || exit 1
done

tar cvf git_dist.tar git_dist
rm -rf "git_dist"
