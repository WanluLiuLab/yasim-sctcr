#!/usr/bin/env bash
set -ue

git ls-files |
    grep \.py$ |
    while read -r fn; do
        sed 's;\#.*;;' <"${fn}" |
            sed 's;^[[:blank:]]*;;' |
            sed 's;[[:blank:]]*$;;' |
            sed 's;$;'" # ${fn}"';'
    done |
    grep 'import ' |
    grep -v '^\#' |
    grep -v 'typing' |
    sort | uniq
