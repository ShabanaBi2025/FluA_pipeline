#!/bin/bash

REF_DIR="reference/h1n1"
INPUT="$REF_DIR/reference.fna"

# ✅ Check input exists
if [[ ! -f "$INPUT" ]]; then
    echo "❌ reference.fna not found in $REF_DIR"
    exit 1
fi

echo "📦 Splitting segments from $INPUT"

# Split segments directly into root directory
awk -v outdir="$REF_DIR" '
  /^>/ {
    if (out) close(out)
    if (match($0, /segment[[:space:]]+([0-9]+)/, seg)) {
      out = outdir "/segment_" seg[1] ".fasta"
    } else {
      out = outdir "/segment_unknown.fasta"
    }
  }
  { if (out) print >> out }
' "$INPUT"

echo "🧼 Removing original reference file: $INPUT"
rm "$INPUT"

echo "✅ Segments saved in $REF_DIR and reference.fna deleted"
