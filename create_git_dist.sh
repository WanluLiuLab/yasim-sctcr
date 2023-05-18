#!/usr/bin/env bash
set -ue
SHDIR="$(dirname "$(readlink -f "${0}")")"
cd "${SHDIR}/" || exit 1
rm -rf "${SHDIR}/git_dist"
mkdir -p "${SHDIR}/git_dist"

git clone --mirror "$(git remote get-url origin)" "${SHDIR}/git_dist/yasim.git"

cd "${SHDIR}/deps/labw_utils" || exit 1
git clone --mirror "$(git remote get-url origin)" "${SHDIR}/git_dist/labw_utils.git"

cd "${SHDIR}/" || exit 1

tar cvf "${SHDIR}/git_dist.tar" "${SHDIR}/git_dist"
rm -rf "${SHDIR}/git_dist"
