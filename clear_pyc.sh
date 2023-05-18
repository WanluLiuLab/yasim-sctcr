#!/usr/bin/env bash
find . | grep \.pyc$ | while read -r line; do rm -v "${line}"; done
find . | grep \.pytest_cache$ | while read -r line; do rm -vrf "${line}"; done
