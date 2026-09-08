#!/bin/bash -l
#SBATCH --job-name=oifsbuild
#SBATCH --partition=compute
#SBATCH --account=ab0246
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
#
# Build an OpenIFS tree with esm_master on a COMPUTE node.
#
# WHY NOT ON THE LOGIN NODE.  The IFS build runs `make -j32`, and several of those
# targets each spawn `fcm make -j 256` -- FCM_PARALLEL comes from CMake's
# ProcessorCount(), which reports 256 on a Levante login node (2x EPYC 7763).  That is
# fine from a plain ssh shell, and has been for years.  It is NOT fine from inside a
# VSCode-remote / agent session: on 2026-08-21 the IDE tooling alone held ~1310 threads
# (devin 399, language_server 300, github-mcp 218, claude 206) against RLIMIT_NPROC of
# 2048 -- ulimit -u counts THREADS on Linux -- leaving ~700 for the build.  Every attempt
# died with "ifort: error #10103: can't fork process: Resource temporarily unavailable",
# and the first failure orphaned ~250 fcm workers that then broke the retries.
#
# A batch job gets its own allocation and does not share the login session's budget, so
# the as-released FCM_PARALLEL=256 is left untouched.
#
# USAGE:  sbatch build_oifs_compute.sh <model_dir> <esm_master_target>
#   e.g.  sbatch build_oifs_compute.sh /work/.../oifsamip-cy48        recomp-oifs-48r1
#         sbatch build_oifs_compute.sh /work/.../awiesm3-develop      recomp-oifs-48r1v5
set -eu
MODEL_DIR="${1:?need model dir}"
TARGET="${2:?need esm_master target, e.g. recomp-oifs-48r1}"

cd "$MODEL_DIR"
echo "host=$(hostname)  cpus=$(nproc --all)  ulimit-u=$(ulimit -u)"
echo "threads in this job at start: $(ps -u "$USER" -L --no-headers | wc -l)"
echo "building $TARGET in $MODEL_DIR"
esm_master "$TARGET"
echo "esm_master exit=$?"
