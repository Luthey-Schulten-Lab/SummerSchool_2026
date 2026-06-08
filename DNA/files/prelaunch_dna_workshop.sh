#!/bin/bash

# Copy the full-cell DNA workshop template into the user's bgvl directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_BASE="/projects/bgvl/${USER}"

# Template may live next to this script (git clone) or on the shared bgvl copy.
for TEMPLATE_DIR in \
    "${SCRIPT_DIR}/DNA_SummerSchool_2026" \
    "/projects/bgvl/SummerSchool_2026/DNA/files/DNA_SummerSchool_2026" \
    "/projects/bgvl/SummerSchool_2026/DNA/DNA_SummerSchool_2026"
do
    if [[ -d "${TEMPLATE_DIR}/scripts" ]]; then
        break
    fi
    TEMPLATE_DIR=""
done

if [[ -z "${TEMPLATE_DIR}" ]]; then
    echo "ERROR: workshop template not found." >&2
    echo "Expected DNA_SummerSchool_2026/scripts under:" >&2
    echo "  ${SCRIPT_DIR}/DNA_SummerSchool_2026" >&2
    echo "  /projects/bgvl/SummerSchool_2026/DNA/files/DNA_SummerSchool_2026" >&2
    exit 1
fi

for LAUNCH_SCRIPT in \
    "${SCRIPT_DIR}/launch_simulation.sh" \
    "/projects/bgvl/SummerSchool_2026/DNA/files/launch_simulation.sh"
do
    if [[ -f "${LAUNCH_SCRIPT}" ]]; then
        break
    fi
    LAUNCH_SCRIPT=""
done

if [[ -z "${LAUNCH_SCRIPT}" ]]; then
    echo "ERROR: launch_simulation.sh not found." >&2
    exit 1
fi

mkdir -p "${DEST_BASE}"
cp -p "${LAUNCH_SCRIPT}" "${DEST_BASE}/"
cp -rp "${TEMPLATE_DIR}" "${DEST_BASE}/"

echo "Workshop files copied to ${DEST_BASE}"
echo "Template source: ${TEMPLATE_DIR}"
