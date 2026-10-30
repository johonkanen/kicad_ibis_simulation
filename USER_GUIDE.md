# User guide: running sweeps without hangs

How to set up an IBIS sweep with `kicad_ibis` so it finishes in minutes, how to tell a slow run from a hung one, and
what to do about each. The numbers come from the LPDDR4 routed-board sweep in
https://github.com/johonkanen/agilex3_soc_control_card (`sim/lpddr4_si/routed_si.py`): 54 nets of real copper, 426
netlists, on a 32-core, 19 GB WSL2 machine. The method itself is in `CLAUDE.md`; the module list is in `README.md`.

## The short version

1. **One net per netlist** when the nets are uncoupled. Never put a whole byte or bus into one netlist.
2. **`.save` only the vectors you analyse** (driver die, receiver pin, receiver die).
3. **Delete each `wrdata` file** once it has been loaded and reduced to numbers.
4. **Time one netlist alone** before starting a sweep, then size `-j` from its peak RAM.
5. **Keep the machine awake**, and write results as they come in, not only at the end.

## Why a run is slow

ngspice prints nothing while a transient runs, and the `.log` that `spice.ngspice()` writes appears only at the
end. So a run that will take an hour and a run that will never finish look the same: one process at 100 % CPU and no
output. The causes we have hit, and their fixes:

| cause | symptom | fix | what it gave |
|---|---|---|---|
| **Many nets in one netlist** | the run takes far longer than the sum of its nets alone. ngspice has one timestep for the whole circuit, so every net is solved at the step the busiest edge of any net needs, and the matrix grows with every route graph | one netlist per net (`build(..., [net])`), run in parallel | a byte of 11 DQ/DM/DQS nets: not finished after 75 min. The same nets alone: 11–20 s each |
| **All nodes stored** | RAM grows to GBs per run; parallel runs get OOM-killed | `n.directive('.save v(a) v(b) ...')` with only the probed nodes | a byte netlist: GBs → 283 MB |
| **Big output files** | `work/` reaches 10+ GB; disk fills | `os.remove(dat)` after `spice.load()` | 130–290 MB per byte netlist, kept |
| **Very short lossless lines (T)** | "timestep too small", or a run that crawls | `fastnet.Netlist.tl` lumps lines < 15 ps into damped L/C (automatic) | CKE trunk: stopped at 20 of 102 ns → 5 s |
| **Undamped lumped L/C** (your own vias, stubs) | 100+ GHz tank rings on every edge, sub-ps steps all run long | R across each L (≈ Z0, or a few 100 Ω for a via) | as above |
| **Many IBIS drivers in one netlist** | ngspice > 9 GB | few drivers per netlist (one per net) | ~0.3 GB each at 150 bits |
| **KIBIS conversion** | the first run of a new sweep spends minutes in kicad-cli | `KibisCache` (automatic): each buffer, corner and rate is converted once | 748 conversions once; 0 s on later runs |

Rule of thumb: if a netlist takes much longer than (number of nets) × (one net alone), it is the shared timestep.
Split it.

## Before a long sweep

1. **Time one netlist of the worst kind** (the longest route, the min corner, the highest rate) with
   `/usr/bin/time -v` and note the wall time and "Maximum resident set size". Min corners run slower than typ.
2. **Choose `-j`** = min(free cores, free RAM / peak RSS, with 2–3 GB left over). Each worker is a separate
   `python3 -I -c ...` process that loads libngspice. The parent also holds the eye data of every row.
3. **Estimate the total**: netlists × time ÷ `-j`. If it comes out at hours, something above has not been applied.
4. **Stop the machine sleeping.** On WSL2, Windows sleep pauses the VM. The workers resume afterwards, but the run
   loses that time and every estimate is off. Set the Windows power plan to "never sleep" while plugged in, or on
   native Linux run the sweep under `systemd-inhibit`.
5. **Write results incrementally** (append each finished netlist's rows to a CSV) so a killed or stopped sweep can
   resume and skip what is done. A sweep that writes only at the end loses everything when it is killed.

## Is it hung or just slow?

- **Find the workers.** ngspice runs in-process through ctypes, so `pgrep ngspice` finds nothing. Look for the
  library instead: `pgrep -af libngspice`. The netlist name is in the command line (`source .../X.cir`).
- **Check CPU time against wall time**: `ps -o pid,etime,time,rss -p <pid>`. TIME close to ELAPSED and RSS still
  changing means it is computing. TIME far below ELAPSED means it was paused (sleep) or is waiting for memory.
- **Check which outputs finished**: `ls -lt work/*.dat`. If no new file has appeared for many times the single-net
  time, the netlists still running are too big. Do not wait them out: stop the sweep and split.
- **Find the slow net** by running each net of the slow netlist alone with a timeout. Remember to export
  `KICAD_APPRUN` in that shell, or every run fails in 0 s with "libngspice not found".

## Stopping a sweep cleanly

- `ProcessPoolExecutor` workers each start an ngspice child. Killing the parent leaves the children running and
  still using RAM. Kill the children too: `pgrep -f libngspice` → `kill <pids>`.
- Avoid `pkill -f '<pattern>'` from a shell whose own command line contains the pattern: it kills that shell
  first. Kill by PID instead.
- Remove the partial `.dat` / `.cir` files of the stopped points before you rerun.

## Reference sweep

`sim/lpddr4_si/routed_si.py` in the board repo follows all of the above. It extracts copper with `extract_routes.py`,
builds one netlist per net, saves only the probed vectors and deletes each `.dat` after reading it. It runs 12
workers.
