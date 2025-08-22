#!/bin/bash

# Path to segment FASTA files
REF_DIR="reference/h1n1"
cd "$REF_DIR" || {
    echo "❌ Could not enter $REF_DIR"
    exit 1
}

# Rename each segment to its standard gene name
for f in segment_*.fasta; do
  case "$f" in
    segment_1.fasta) mv "$f" PB2.fasta ;;
    segment_2.fasta) mv "$f" PB1.fasta ;;
    segment_3.fasta) mv "$f" PA.fasta ;;
    segment_4.fasta) mv "$f" HA.fasta ;;
    segment_5.fasta) mv "$f" NP.fasta ;;
    segment_6.fasta) mv "$f" NA.fasta ;;
    segment_7.fasta) mv "$f" M.fasta ;;
    segment_8.fasta) mv "$f" NS.fasta ;;
    *) echo "⚠️ Unexpected segment file: $f" ;;
  esac
done

echo "✅ Segment files renamed successfully in $REF_DIR"
