#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALLOWLIST_FILE="${ROOT_DIR}/scripts/private_api_allowlist.txt"
CURRENT_FILE="$(mktemp)"
trap 'rm -f "${CURRENT_FILE}"' EXIT

cd "${ROOT_DIR}"

rg -n --no-heading "\._seq|\._data" biokit tests -S > "${CURRENT_FILE}" || true

if ! diff -u "${ALLOWLIST_FILE}" "${CURRENT_FILE}"; then
  echo
  echo "Private API guard failed."
  echo "If these changes are intentional, update ${ALLOWLIST_FILE}."
  exit 1
fi

echo "Private API guard passed."
