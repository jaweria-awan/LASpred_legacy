#!/bin/sh
# Compile and run the original 2011 code in Docker without editing it.
#
# Main.java hard-codes Windows paths (E:/...). This script copies the source
# into the container, rewrites only those three paths, then compiles and runs.
#
# Usage: ./run.sh <input.fasta> <output.txt>
set -eu

if [ $# -ne 2 ]; then
  echo "usage: $0 <input.fasta> <output.txt>" >&2
  exit 1
fi
input=$1
output=$2

[ -f "$input" ] || { echo "error: $input not found" >&2; exit 1; }
case "$input$output" in
  /*|*..*) echo "error: use paths inside this folder, e.g. example.fasta" >&2; exit 1 ;;
esac
# The original loops forever if the file does not start with a FASTA header.
head -c 1 "$input" | grep -q '>' || {
  echo "error: $input must start with a '>' header line" >&2; exit 1; }

docker run --rm -v "$PWD":/w eclipse-temurin:21-jdk sh -c "
  set -e
  mkdir -p /tmp/build/finalcode && cd /tmp/build
  sed -e 's#E:/e-coli15pos.txt#/w/$input#' \
      -e 's#E:/outpute-coli.txt#/w/$output#' \
      -e 's#E:/InputFileForHumanPosi.txt#/w/InputFileForHumanPosi.txt#' \
      /w/src/finalcode/Main.java > finalcode/Main.java
  javac -nowarn finalcode/Main.java
  java finalcode.Main
"
echo "wrote $output"
