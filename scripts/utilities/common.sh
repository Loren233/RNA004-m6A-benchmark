#!/usr/bin/env bash

require_var() {
    local name="$1"
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: required variable is not set: ${name}" >&2
        exit 1
    fi
}

require_file() {
    local path="$1"
    [[ -s "$path" ]] || { echo "ERROR: missing input: $path" >&2; exit 1; }
}

require_dir() {
    local path="$1"
    [[ -d "$path" ]] || { echo "ERROR: missing directory: $path" >&2; exit 1; }
}
