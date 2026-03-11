#!/bin/bash
set -u -o pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${ROOT_DIR}/.." && pwd)"
CONFIG_PATH="${1:-${ROOT_DIR}/config.yml}"
JULIA_BIN="${JULIA_BIN:-julia}"

if [ ! -f "${CONFIG_PATH}" ]; then
    echo "Config file not found: ${CONFIG_PATH}"
    exit 1
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"

RESULTS_DIR="$(cd "${REPO_DIR}" && "${JULIA_BIN}" --project=. -e 'using YAML; cfg = YAML.load_file(ARGS[1]); println(abspath(joinpath(dirname(ARGS[1]), String(cfg["results_dir"]))))' "${CONFIG_PATH}")"
mkdir -p "${RESULTS_DIR}/raw_logs" "${RESULTS_DIR}/csv"

mapfile -t N_VALUES < <(cd "${REPO_DIR}" && "${JULIA_BIN}" --project=. -e 'using YAML; cfg = YAML.load_file(ARGS[1]); for n in cfg["n_values"]; println(Int(n)); end' "${CONFIG_PATH}")

BACKENDS=("vec_views" "vec_vectors" "mat_opt_cacheauto" "mat_opt_cacheoff" "mat_scan_cacheauto" "mat_scan_cacheoff")
BENCHMARK_KINDS=("oracle" "fw" "afw")

run_case() {
    local benchmark_kind="$1"
    local n="$2"
    local backend="$3"
    local case_id="${benchmark_kind}_n${n}_${backend}"
    local raw_log="${RESULTS_DIR}/raw_logs/${case_id}.log"

    echo "Running ${case_id}"
    if [ "${benchmark_kind}" = "oracle" ]; then
        if (cd "${REPO_DIR}" && /usr/bin/time -v "${JULIA_BIN}" --project=. "${ROOT_DIR}/run_oracle_bench.jl" "${CONFIG_PATH}" "${n}" "${backend}" > "${raw_log}" 2>&1); then
            status=0
        else
            status=$?
        fi
    else
        if (cd "${REPO_DIR}" && /usr/bin/time -v "${JULIA_BIN}" --project=. "${ROOT_DIR}/run_algorithm_probe.jl" "${CONFIG_PATH}" "${benchmark_kind}" "${n}" "${backend}" > "${raw_log}" 2>&1); then
            status=0
        else
            status=$?
        fi
    fi
    echo "Exit status: ${status}" >> "${raw_log}"
}

for n in "${N_VALUES[@]}"; do
    for backend in "${BACKENDS[@]}"; do
        for benchmark_kind in "${BENCHMARK_KINDS[@]}"; do
            run_case "${benchmark_kind}" "${n}" "${backend}"
        done
    done
done

(cd "${REPO_DIR}" && "${JULIA_BIN}" --project=. "${ROOT_DIR}/summarize_results.jl" "${CONFIG_PATH}")
