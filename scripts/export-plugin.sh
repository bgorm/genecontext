#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 1 || -z $1 ]]; then
    echo 'Usage: bash scripts/export-plugin.sh NEW_DIRECTORY' >&2
    exit 1
fi
case "$1" in
    /*) destination=$1 ;;
    *) destination=$PWD/$1 ;;
esac
if [[ -e $destination || -L $destination ]]; then
    echo "Destination already exists: $destination" >&2
    exit 1
fi

cd "$(dirname "$0")/.."
# Publish only these files; local data, chat exports, notes, and build products
# stay in the development workspace.
files=(
    README.md LICENSE THIRD_PARTY.md Makefile .gitignore
    src/annotations.hpp src/common.hpp src/genes.hpp src/model.hpp
    src/options.hpp src/plugin.cpp src/store.hpp
    scripts/bootstrap.sh scripts/export-plugin.sh
    tests/helpers.py tests/test_genecontext.py tests/test_transcripts.py
    .github/workflows/build.yml
)
for file in "${files[@]}"; do
    [[ -f $file && ! -L $file ]] || { echo "Missing or symlinked source file: $file" >&2; exit 1; }
done

mkdir -p "$destination"
for file in "${files[@]}"; do
    mkdir -p "$destination/$(dirname "$file")"
    cp "$file" "$destination/$file"
    chmod 644 "$destination/$file"
done
printf 'Exported genecontext sources to %s\n' "$destination"
