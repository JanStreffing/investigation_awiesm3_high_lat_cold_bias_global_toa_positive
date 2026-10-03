# Handoff: esm_tools experiment venv build went from ~3.5 min to ~18 min (albedo, 2026-10-01)

## The problem

Every new experiment that `esm_runscripts` submits on albedo builds its own Python venv
(`<exp>/.venv_esmtools`, the `use_venv: True` contained-run mode). Until this afternoon that
took about 3.5 minutes. At 20:05 on 2026-10-01 the build for `PICAL_crunveg_gm1500` took
**17 min 48 s**. Nothing in the runscript or in `~/esm_tools` changed between the fast builds
and the slow one. Find out why, and make the build fast again, or make it not happen per
experiment.

It matters because each 2-member test now costs 40 minutes of wall time before anything
reaches SLURM. A user waiting on `squeue` sees nothing for that time.

## Timings (from `/albedo/work/projects/p_awiesm3_cmip7/jstreffi/submit_*.log`, line `...finished H:MM:SS`)

| submitted | experiment | venv build |
|---|---|---|
| 2026-09-29 21:25 | PICAL_crunveg_ihf0 | 3:14 |
| 2026-09-30 21:09 | PICAL_crunveg_ob1200 | 3:12 |
| 2026-10-01 11:23-12:31 | TEST_dt* and TEST_xd* (10 runs, sequential) | 3:15-4:16 |
| 2026-10-01 13:50 | PICAL_crunveg_nx1800 | 3:19 |
| 2026-10-01 20:05 | **PICAL_crunveg_gm1500** | **17:48** |
| 2026-10-01 20:25 | PICAL_crunveg_gmramp | **21:02** (slow again) |

## What was observed during the slow build

- The log stops at `Building virtual env, please be patient (this takes about 3 minutes)...`
  for the whole time. The submit log is block-buffered, so set `PYTHONUNBUFFERED=1` to see output while it runs.
- For about 15 minutes (20:07 to about 20:22) the running step was
  `<venv>/bin/pip install -q git+https://github.com/JanStreffing/ocp-tool`.
  That is `_install_required_plugins` in `virtual_env_builder.py:216-230`. ocp-tool is a
  required plugin of the awiesm3 setup.
- That pip process was in state `D` (disk sleep) at about 110 MB RSS. Its only open file in
  `/tmp` was a wheel being unpacked (`/tmp/pip-unpack-*/jedi-0.20.0-py2.py3-none-any.whl`, later
  `anyio-4.15.1-...whl`). A wheel file did not grow over 12 s. The pip processes then
  turned over every few seconds, so it was working through many packages slowly, not hanging on one.
- The finished venv has 415 entries in site-packages; nx1800's from 13:50 has 418.

## Already ruled out

- **Local disk:** `/tmp` on the login node is `/dev/sdb1`, 1.7 TB, 1 % used. `dd` with `oflag=direct` wrote 466 MB/s.
- **PyPI bandwidth:** `curl` of a jedi wheel from files.pythonhosted.org returned 200 at 5.6 MB/s in 0.28 s.
- **pip version:** both venvs have pip 26.2.1. The system python has pip 24.0.
- **Package set:** `pip freeze` of the gm1500 and nx1800 venvs is identical except for one line:
  `esm-tools` was installed from `esm_tools@135b682c` at 13:50 and from `@a376cf18` at 20:05.
  The branch `feat/awiesm3-v3.4-co2` moved on GitHub in between; the local `~/esm_tools`
  checkout is at `322f0ed1f`.
- **Login node load:** load average 5.9-8.3 on 28 users. That is normal for albedo, but compare it with the next fast build.

## Leads, roughly in order of how cheap they are to test

1. **Reproduce and time each step.** In a scratch directory, create a venv with the system python
   (`/albedo/soft/sw/spack-sw/python/3.10.4-daos5pe/bin/python3.10 -m venv x`), then run the same
   steps as `venv_bootstrap` (`virtual_env_builder.py:233-268`) one by one under `time`, with
   `pip -v`: `pip install -U pip wheel`, the esm_tools install (`git+https://github.com/esm-tools/esm_tools@feat/awiesm3-v3.4-co2`, see
   `_install_tools_general:111-213`), then `pip install git+https://github.com/JanStreffing/ocp-tool`.
   `pip -v` shows whether the time goes to resolution (repeated `Collecting` or backtracking messages),
   download, build or install.
2. **Commit `a376cf18` on esm_tools `feat/awiesm3-v3.4-co2`.** Check whether it changed
   `setup.py`, `pyproject.toml` or the requirements. A loosened or new pin there could send pip into long resolution.
3. **ocp-tool's dependencies.** It installs from a git URL with no pins visible here. A new release
   of one of its dependencies that day could cause backtracking. jedi 0.20.0 and anyio 4.15.1 went
   past, so check their release dates.
4. **Wheel cache.** The builder passes `--find-links=$HOME/.cache/pip/wheels` and
   `--wheel-dir=$HOME/.cache/pip/wheels` (`virtual_env_builder.py:163,175,197,209`). Check the
   size and number of files there, and whether 116 leftover `/tmp/pip-*` directories from earlier builds point to a cache that is failing.
5. **Concurrency.** Earlier the same evening, two `esm_runscripts` submissions started at the same
   moment (gm1500 and gmramp). gmramp failed with
   `FileNotFoundError: Could not determine where configs's path is inside the esm-tools installation`
   (`esm_tools/__init__.py:208`, from `~/.local/bin/esm_runscripts`), and gm1500 stalled and was killed by
   a 15-minute `timeout`. Afterwards the host install was fine again (`esm_tools.get_config_filepath()` returns
   `~/esm_tools/configs/`). Find out whether parallel builds share a cache or the editable install unsafely.
   Until that is known, submit one experiment at a time.

## Possible fixes, once the cause is known

- Pin or cache whatever ocp-tool pulls, or install it from a prebuilt wheel.
- Build one template venv per esm_tools commit and copy or `--system-site-packages` it, instead of rebuilding per experiment.
  Copying a venv is not trivial: `bin/` shebangs and `activate` hold absolute paths.
- Install the plugin with `--no-deps` if its dependencies are already satisfied by esm_tools.

Any fix to `virtual_env_builder.py` belongs on the esm_tools branch the runs use
(`feat/awiesm3-v3.4-co2`). Do not push or open a PR without asking Jan.

## Do not touch

- The `.venv_esmtools` of running or queued experiments under
  `/albedo/work/projects/p_awiesm3_cmip7/jstreffi/runtime/awiesm3-v3.4/`. Later legs reuse them.
- `/albedo/pool` (shared; the auto-mode classifier refuses changes there).
- Do not submit model runs to test this. A venv build alone does not need SLURM.

## Environment notes

- Login nodes: RLIMIT_NPROC is 2048. Export `OPENBLAS_NUM_THREADS=1` before any numpy work.
- `ssh levante` works if a comparison with a levante build helps; there, `module` is not available in non-login shells.
