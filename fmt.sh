#!/usr/bin/env bash

git ls-files |
    grep -v 'deps' |
    grep -v '.idea/' |
    grep '\.sh$' |
    while read -r line; do
        if [ -e "${line}" ]; then
            {
                echo SHFMT "${line}"
                shfmt -i 4 --write "${line}"
            } &
        fi
    done

git ls-files |
    grep -v 'deps' |
    grep -v '.idea/' |
    grep -e '\.py$' -e '\.pyi$' |
    while read -r line; do
        if [ -e "${line}" ]; then
            {
                echo BLACK "${line}"
                black --line-length 120 --target-version py38 --quiet "${line}"
            } &
        fi
    done
wait

#FormatR has bugs.
#git ls-files |
#    grep -v 'deps' |
#    grep -v '.idea/' |
#    grep -e '\.R$' |
#    while read -r line; do
#        if [ -e "${line}" ]; then
#            {
#                echo FormatR "${line}"
#                echo "formatR::tidy_file(\"${line}\")" | Rscript /dev/stdin
#            } &
#        fi
#    done
#wait
