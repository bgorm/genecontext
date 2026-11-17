#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    echo 'Usage: bash scripts/bootstrap.sh (set PLUGIN_HTSLIB to reuse matching headers)'
    exit 0
fi
if [[ $# != 0 ]]; then echo 'Unexpected bootstrap argument' >&2; exit 1; fi
headers=${PLUGIN_HTSLIB:-vendor/htslib-1.25}
if [[ ! -f $headers/htslib/hts.h && -z ${PLUGIN_HTSLIB:-} && -f ../locusfinder/vendor/htslib-1.25/htslib/hts.h ]]; then
    headers=../locusfinder/vendor/htslib-1.25
fi
if [[ ! -f $headers/htslib/hts.h ]]; then
    if [[ $headers != vendor/htslib-1.25 ]]; then
        echo "HTSlib headers not found in $headers" >&2; exit 1
    fi
    mkdir -p vendor
    archive=vendor/htslib-1.25.tar.bz2
    curl -fL --retry 2 https://github.com/samtools/htslib/releases/download/1.25/htslib-1.25.tar.bz2 -o "$archive"
    if command -v sha256sum >/dev/null 2>&1; then
        checksum=$(sha256sum "$archive")
    elif command -v shasum >/dev/null 2>&1; then
        checksum=$(shasum -a 256 "$archive")
    else
        echo 'A SHA-256 checker is required' >&2; exit 1
    fi
    if [[ ${checksum%% *} != 158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730 ]]; then
        echo "SHA-256 mismatch: $archive" >&2; exit 1
    fi
    tar -xjf "$archive" -C vendor
fi
make -j"${JOBS:-4}" plugin PLUGIN_HTSLIB="$headers"
