#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if (($# != 0)); then
    printf 'Usage: %s\n' "$0" >&2
    exit 2
fi

gff_file="$repo_root/fixtures/gencode_sorted.gff.gz"
gff_coordinate_index="$repo_root/fixtures/gencode_sorted.gff.gz.tbi"
bed_file="$repo_root/fixtures/sorting/gencode_v46.bed.gz"
bed_coordinate_index="$repo_root/fixtures/sorting/gencode_v46.bed.gz.tbi"

for fixture in "$gff_file" "$gff_coordinate_index" "$bed_file" "$bed_coordinate_index"; do
    if [[ ! -f "$fixture" ]]; then
        printf 'Required benchmark fixture not found: %s\n' "$fixture" >&2
        exit 1
    fi
done

target_dir="$(cargo metadata --no-deps --format-version 1 \
    | sed -n 's/.*"target_directory":"\([^"]*\)".*/\1/p')"
if [[ -z "$target_dir" ]]; then
    printf 'Could not determine the Cargo target directory.\n' >&2
    exit 1
fi
gai_bin="$target_dir/release/gai"

printf 'Building the release gai binary (outside benchmark timings)...\n' >&2
cargo build --release --bin gai >&2
if [[ ! -x "$gai_bin" ]]; then
    printf 'Expected release binary not found or not executable: %s\n' "$gai_bin" >&2
    exit 1
fi

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/gai-benchmark.XXXXXXXX")"
cleanup() {
    rm -rf -- "$tmp_dir"
}
trap cleanup EXIT

now_ns() {
    local timestamp
    timestamp="$(date +%s%N)"
    if [[ ! "$timestamp" =~ ^[0-9]+$ ]]; then
        printf 'date must support nanosecond timestamps (date +%%s%%N)\n' >&2
        return 1
    fi
    printf '%s' "$timestamp"
}

format_duration() {
    local nanoseconds="$1"
    printf '%d.%03ds' "$((nanoseconds / 1000000000))" "$(((nanoseconds % 1000000000) / 1000000))"
}

file_size_bytes() {
    local size
    size="$(wc -c < "$1")"
    size="${size//[[:space:]]/}"
    printf '%s' "$size"
}

format_megabytes() {
    local bytes="$1"
    local whole="$((bytes / 1000000))"
    local fraction="$(((bytes % 1000000 + 5000) / 10000))"
    if ((fraction == 100)); then
        whole="$((whole + 1))"
        fraction=0
    fi
    printf '%d.%02d' "$whole" "$fraction"
}

build_index() {
    local label="$1"
    local source="$2"
    local coordinate_index="$3"
    local output="$4"
    local reported_attributes="$5"
    shift 5
    local -a attributes=("$@")
    local start_ns end_ns duration raw_size index_size

    printf 'Building %s index with attributes %s...\n' "$label" "$reported_attributes" >&2
    start_ns="$(now_ns)"
    "$gai_bin" build-index "$source" \
        --coordinate-index "$coordinate_index" \
        --output "$output" \
        "${attributes[@]}" >&2
    end_ns="$(now_ns)"
    duration="$((end_ns - start_ns))"
    raw_size="$(file_size_bytes "$source")"
    index_size="$(file_size_bytes "$output")"

    printf '| %s | %s | %s | %s | %s |\n' \
        "$label" \
        "$reported_attributes" \
        "$(format_megabytes "$raw_size")" \
        "$(format_megabytes "$index_size")" \
        "$(format_duration "$duration")"
}

query_index() {
    local label="$1"
    local attributes="$2"
    local source="$3"
    local coordinate_index="$4"
    local gai_index="$5"
    local query_type="$6"
    local query_term="$7"
    local query_output="$tmp_dir/query-output"
    local start_ns end_ns duration record_count

    printf 'Querying %s (%s, %s: %s)...\n' \
        "$label" "$attributes" "$query_type" "$query_term" >&2
    start_ns="$(now_ns)"
    "$gai_bin" query-index "$source" "$query_term" \
        --coordinate-index "$coordinate_index" \
        --gai "$gai_index" \
        --match "$query_type" > "$query_output"
    end_ns="$(now_ns)"
    duration="$((end_ns - start_ns))"
    record_count="$(wc -l < "$query_output")"
    record_count="${record_count//[[:space:]]/}"

    printf '| %s | %s | %s | `%s` | %s | %s |\n' \
        "$label" "$attributes" "$query_type" "$query_term" \
        "$record_count" "$(format_duration "$duration")"
}

gff_attributes_1='gene_name'
gff_attributes_2='gene_name, transcript_name'
gff_attributes_3='ID, gene_name, transcript_name'
bed_attributes='name'

gff_index_1="$tmp_dir/gencode-gene-name.gai"
gff_index_2="$tmp_dir/gencode-gene-transcript-name.gai"
gff_index_3="$tmp_dir/gencode-id-gene-transcript-name.gai"
bed_index="$tmp_dir/gencode-bed-name.gai"

printf '## Index build results\n\n'
printf '| Index / file | Attributes indexed | Raw compressed file size (MB) | Index size (MB) | Build time |\n'
printf '| :--- | :--- | ---: | ---: | ---: |\n'
printf 'Timing windows: index-build rows exclude cargo compilation; query rows include the full gai query invocation and captured output write.\n' >&2

build_index 'gencode_sorted.gff.gz' "$gff_file" "$gff_coordinate_index" \
    "$gff_index_1" "$gff_attributes_1" \
    --attribute gene_name
build_index 'gencode_sorted.gff.gz' "$gff_file" "$gff_coordinate_index" \
    "$gff_index_2" "$gff_attributes_2" \
    --attribute gene_name --attribute transcript_name
build_index 'gencode_sorted.gff.gz' "$gff_file" "$gff_coordinate_index" \
    "$gff_index_3" "$gff_attributes_3" \
    --attribute ID --attribute gene_name --attribute transcript_name
build_index 'gencode_v46.bed.gz' "$bed_file" "$bed_coordinate_index" \
    "$bed_index" "$bed_attributes"

printf '\n## Query results\n\n'
printf '| Index queried | Attributes indexed | Query type | Query | Records matched | Query time |\n'
printf '| :--- | :--- | :--- | :--- | ---: | ---: |\n'

for index_path in "$gff_index_1" "$gff_index_2" "$gff_index_3"; do
    case "$index_path" in
        "$gff_index_1") attributes="$gff_attributes_1" ;;
        "$gff_index_2") attributes="$gff_attributes_2" ;;
        "$gff_index_3") attributes="$gff_attributes_3" ;;
    esac
    query_index 'gencode_sorted.gff.gz' "$attributes" "$gff_file" \
        "$gff_coordinate_index" "$index_path" exact 'brca1'
    query_index 'gencode_sorted.gff.gz' "$attributes" "$gff_file" \
        "$gff_coordinate_index" "$index_path" contains 'orf'
    query_index 'gencode_sorted.gff.gz' "$attributes" "$gff_file" \
        "$gff_coordinate_index" "$index_path" prefix 'brca'
    query_index 'gencode_sorted.gff.gz' "$attributes" "$gff_file" \
        "$gff_coordinate_index" "$index_path" regex 'c\d+orf'
done

bed_name='ENST00000607096.1'
query_index 'gencode_v46.bed.gz' "$bed_attributes" "$bed_file" \
    "$bed_coordinate_index" "$bed_index" exact "$bed_name"
query_index 'gencode_v46.bed.gz' "$bed_attributes" "$bed_file" \
    "$bed_coordinate_index" "$bed_index" contains '607096'
query_index 'gencode_v46.bed.gz' "$bed_attributes" "$bed_file" \
    "$bed_coordinate_index" "$bed_index" prefix 'ENST000006070'
query_index 'gencode_v46.bed.gz' "$bed_attributes" "$bed_file" \
    "$bed_coordinate_index" "$bed_index" regex '^ENST00000607096\.1$'
