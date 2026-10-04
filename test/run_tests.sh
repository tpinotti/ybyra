#!/usr/bin/env bash
# Run the ybyra tests on the example data. Requires the ybyra environment to be activated.
# Use `--update` to overwrite the expected results in `data/` with the results of this run.

set -u

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "${TEST_DIR}")"
OUT_DIR="${TEST_DIR}/out"
CORES="${CORES:-2}"
UPDATE=false
if [ "${1:-}" = "--update" ]; then
    UPDATE=true
fi

FAILED=0
pass() { echo "PASS    $1"; }
fail() { echo "FAIL    $1: $2"; FAILED=$((FAILED + 1)); }

# Set up the working dir for a test case: units, config, and links to the example data.
setup_case() {
    local dir="${OUT_DIR}/$1"
    mkdir -p "${dir}"
    cp "${TEST_DIR}/units.tsv" "${dir}/units.tsv"
    cp "$2" "${dir}/config.yaml"
    ln -s "${ROOT_DIR}/example/bams" "${dir}/bams"
    ln -s "${ROOT_DIR}/example/ref" "${dir}/ref"
}

# Run snakemake in the given working dir, with additional arguments, logging to that dir.
run_snakemake() {
    local dir="$1"
    shift
    snakemake --snakefile "${ROOT_DIR}/workflow/Snakefile" --directory "${dir}" "$@" \
        > "${dir}/snakemake.log" 2>&1
}

# If the step size exploration was run, check that its outputs exist, and that its placements at
# the configured step size are the ones in aggregate.yplace.
check_step_size() {
    local dir="$1"
    grep -q "^step_size_exploration: true" "${dir}/config.yaml" || return 0
    for f in step_size.tsv overview.pdf samples.pdf; do
        [ -s "${dir}/step_size/${f}" ] || return 1
    done
    diff <(awk -F'\t' 'NR > 1 { print $1, $2, $3 }' "${dir}/aggregate.yplace" | sort) \
         <(awk -F'\t' 'NR > 1 && $4 == "True" && $5 != "" { print $1, $5, $6 }' \
            "${dir}/step_size/step_size.tsv" | sort) > /dev/null
}

rm -rf "${OUT_DIR}"

# Full runs with each tree, compared to the expected results.
for name in hg37_ftdna hg37_isogg hg37_yfull; do
    dir="${OUT_DIR}/${name}"
    expected="${TEST_DIR}/data/expected_aggregate.${name}.yplace"
    setup_case "${name}" "${TEST_DIR}/configs/${name}.yaml"

    if ! run_snakemake "${dir}" --cores "${CORES}"; then
        fail "${name}" "snakemake failed, see ${dir}/snakemake.log"
    elif [ ! -s "${dir}/aggregate.pdf" ] || [ ! -s "${dir}/score_ties.pdf" ]; then
        fail "${name}" "plots missing or empty"
    elif ! check_step_size "${dir}"; then
        fail "${name}" "step size exploration missing, or inconsistent with aggregate.yplace"
    elif ${UPDATE}; then
        cp "${dir}/aggregate.yplace" "${expected}"
        echo "UPDATE  ${name}"
    elif ! diff "${expected}" "${dir}/aggregate.yplace"; then
        fail "${name}" "aggregate.yplace differs from expected"
    else
        pass "${name}"
    fi
done

# Invalid configs, which need to be rejected by the config validation.
for name in bad_typo bad_strbool bad_enum; do
    dir="${OUT_DIR}/${name}"
    setup_case "${name}" "${TEST_DIR}/configs/${name}.yaml"

    if run_snakemake "${dir}" --dry-run; then
        fail "${name}" "invalid config was accepted"
    elif ! grep -q "Error validating config file" "${dir}/snakemake.log"; then
        fail "${name}" "snakemake failed, but not on config validation, see ${dir}/snakemake.log"
    else
        pass "${name}"
    fi
done

if [ "${FAILED}" -gt 0 ]; then
    echo "${FAILED} test(s) failed"
    exit 1
fi
echo "All tests passed"
