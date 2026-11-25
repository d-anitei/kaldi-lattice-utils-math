#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  printf 'Usage: %s <lattice.gz> <words.txt>\n' "$0" >&2
  exit 2
fi

if [[ ! -f $1 || ! -r $1 ]]; then
  printf 'Error: cannot read lattice file %s\n' "$1" >&2
  exit 1
fi

if [[ ! -f $2 || ! -r $2 ]]; then
  printf 'Error: cannot read vocabulary %s\n' "$2" >&2
  exit 1
fi

mkdir -p TMP-IDX
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
index_tmp_dir=$(mktemp -d TMP-IDX/index.XXXXXX)
index_text="$index_tmp_dir/idx.txt"
index_gzip="$index_tmp_dir/idx.gz"
trap 'rm -f -- "$index_text" "$index_gzip"; rmdir -- "$index_tmp_dir"' EXIT

gzip -dc -- "$1" |
  lattice-char-index-position-math --print-args=true --verbose=2 \
    --nbest=25 \
    --min-posterior=1e-10 \
    --symbol-table="$2" \
    --open-delimiters="{ {@" \
    --close-delimiters="} }@" \
    --operator-segment-prefixes="@bin @rel @punct" \
    ark:- ark,t:- 2>TMP-IDX/log |
  VOCAB=$2 gawk '
    BEGIN {
      vocab = ENVIRON["VOCAB"]
      while ((getline < vocab) > 0) {
        words[$2] = $1
      }
      close(vocab)
    }
    {
      delete scores
      delete rows
      for (i = 2; i + 4 <= NF; i += 6) {
        count = split($i, ids, "_")
        expression = ""
        for (j = 1; j <= count; j++) {
          expression = expression words[ids[j]] "¬"
        }
        start = $(i + 1)
        key = expression "_" start
        score = exp($(i + 4))
        if (!(key in scores) || scores[key] < score) {
          scores[key] = score
          rows[key] = $1 " " expression " " score " " $(i + 2) " " $(i + 3) " " start
        }
      }
      for (key in rows) {
        print rows[key]
      }
    }
  ' |
  sed 's/\.png / /g; s/  \+/ /g; s/¬¬\+/¬/g; s/ ¬/ /g; s/¬ / /g' > "$index_text"

python3 "$script_dir/filter_inconsistent_MEs.py" "$index_text"
gzip -9 -c -- "$index_text" > "$index_gzip"
mv -f -- "$index_gzip" TMP-IDX/idx.gz
