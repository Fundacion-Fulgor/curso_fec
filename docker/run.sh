#!/usr/bin/env bash
# Run a command inside the course container, with the repository mounted at /work.
# The image is built the first time. Works on Linux, WSL and macOS (bash 3.2).
#
#   docker/run.sh                      interactive shell
#   docker/run.sh pytest -q            any command, run from the current directory
#   docker/run.sh --rebuild            rebuild the image, then open a shell
#
# FEC_IMAGE sets the image name and FEC_JOBS the parallel build jobs (default: CPUs, at most 4).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
image="${FEC_IMAGE:-fec-env:1}"

if ! command -v docker > /dev/null 2>&1; then
    echo "docker was not found in PATH" >&2
    exit 1
fi

rebuild=0
if [ "${1:-}" = "--rebuild" ]; then
    rebuild=1
    shift
fi

if [ "$rebuild" -eq 1 ] || ! docker image inspect "$image" > /dev/null 2>&1; then
    cpus="$(getconf _NPROCESSORS_ONLN 2> /dev/null || echo 2)"
    jobs="${FEC_JOBS:-$(( cpus < 4 ? cpus : 4 ))}"
    docker build \
        --build-arg JOBS="$jobs" \
        --build-arg USER_UID="$(id -u)" \
        --build-arg USER_GID="$(id -g)" \
        --tag "$image" \
        "$here"
fi

# keep the current directory if it is inside the repository
workdir=/work
case "$PWD/" in
    "$repo"/*) workdir="/work${PWD#"$repo"}" ;;
esac

tty_flags=""
if [ -t 0 ] && [ -t 1 ]; then
    tty_flags="-it"
fi

# shellcheck disable=SC2086
exec docker run --rm $tty_flags \
    --volume "$repo":/work \
    --workdir "$workdir" \
    "$image" "$@"
