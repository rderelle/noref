#!/usr/bin/env bash
# Run the release binary on the golden dataset (directory and table inputs) and
# compare all output files with the expected ones.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

BIN="$(pwd)/target/release/noref"
DATA="$(pwd)/tests/data"
OUTPUT_FILES=(
  "noref_raw_alignment.fas"
  "noref_filtered_alignment.fas"
  "noref_raw_missing.tsv"
  "noref_filtered_missing.tsv"
  "noref_raw_NJ.tree"
  "noref_filtered_NJ.tree"
)

run_case() {
  local name="$1"
  shift
  local work_dir="target/golden_${name}"
  rm -rf "$work_dir"
  mkdir -p "$work_dir"
  (cd "$work_dir" && "$BIN" "$@" --nj)
  for file_name in "${OUTPUT_FILES[@]}"; do
    diff -u "${DATA}/expected/${file_name}" "${work_dir}/${file_name}"
  done
  # Recombination plot (--plot), without its comment line (it names the version).
  local plot="noref_recombination.gff"
  if [[ -f "${work_dir}/${plot}" ]]; then
    diff -u <(sed 3d "${DATA}/expected/${plot}") <(sed 3d "${work_dir}/${plot}")
  fi
  echo "Golden test passed: ${name}."
}

# The expected outputs were made with -k 15 -f 0.85 (as in tests/golden.rs).
run_case "directory" -i "${DATA}/input" -k 15 -f 0.85
run_case "table" -I "${DATA}/input.tsv" -t 2 -k 15 -f 0.85
run_case "plot" -i "${DATA}/input" -p G7 -k 15 -f 0.85
