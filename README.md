# kicad_ibis_simulation

Headless IBIS signal-integrity simulation with KiCad 10 and ngspice: Python writes a KiCad schematic (IBIS buffers,
TLINEs, vias), `kicad-cli` converts it to SPICE (KiCad's IBIS converter, KIBIS), libngspice simulates it, and the
eye tools fold the result into eye diagrams with Rx-mask margins. The schematics remain ordinary KiCad files that
open and simulate in Eeschema.

`CLAUDE.md` is the full method (for Claude Code, and a good checklist for people too), including the KiCad 10 IBIS
parser problems found so far and their workarounds.

## Contents

| module | what |
|---|---|
| `kicad_ibis/sch.py` | `Sch`: write a simulation schematic (IBIS driver/receiver, TLINE, R/L/C, pi vias, directives) |
| `kicad_ibis/spice.py` | `export()`: kicad-cli netlist + KIBIS cache clean-up; `ngspice()`: run through libngspice; `load()` |
| `kicad_ibis/eye.py` | KIBIS PRBS7 bit labels, folding, eye height/width, mask margins, plots |
| `kicad_ibis/show_eye.py` | interactive eye viewer for any ngspice `wrdata` file |
| `kicad_ibis/p2p.py` | one point-to-point link in one command (schematic, run, eye) |
| `kicad_ibis/ibisfix.py` | make vendor IBIS files load in KiCad 10; read V-t tables for driver validation |
| `kicad_ibis/zsolve.py` | 2D field solver: trace width -> Z0, ps/mm; edge-coupled pairs -> Zodd/Zdiff/Zeven (stripline, microstrip with mask) |

## Set-up

```sh
kicad-10.0.7-x86_64.AppImage --appimage-extract          # once, anywhere
export KICAD_APPRUN=$PWD/squashfs-root/AppRun
pip install numpy scipy matplotlib
```

As a submodule of a board project: `git submodule add https://github.com/johonkanen/kicad_ibis_simulation.git
tools/kicad_ibis_simulation`, then run the modules with `PYTHONPATH=tools/kicad_ibis_simulation python3 -m ...` or
`sys.path.insert` in your scripts.

## Examples

```sh
# fix a vendor file for KiCad's parser (keep vendor files out of git)
python3 -m kicad_ibis.ibisfix vendor.ibs models/vendor_kicad.ibs

# trace width -> impedance
python3 -m kicad_ibis.zsolve strip --below 0.357 --above 0.1 --t 0.0152 --er-below 4.16 --er-above 4.6 -w 0.1 0.12
# ... and a differential pair (adds --gap): Zodd, Zdiff, Zeven, Zcm, odd/even ps/mm
python3 -m kicad_ibis.zsolve strip --below 0.357 --above 0.1 --t 0.0152 --er-below 4.16 --er-above 4.6 -w 0.1 --gap 0.1 0.2

# one link, LPDDR4 write at 1600 MT/s
python3 -m kicad_ibis.p2p --drv models/fpga.ibs:FPGA:T62:out_40 --rx models/dram.ibs:DRAM:C11:ODT_48 \
    --line 55:1.5:7.0 --line 40:17:7.0 --line 55:1.5:7.0 --via 0.5n:0.4p --rate 1600 --mask rect:140m:0.22

# look at an eye from any run
python3 -m kicad_ibis.show_eye work/x.dat --list
python3 -m kicad_ibis.show_eye work/x.dat --drv xu1.die0 --rx dq_rx xu2.die0 --rate 1600 --mask rect:140m:0.22
```

A full project using it (LPDDR4, FPGA to DRAM, sweeps over impedance, length, corners and data rate):
`sim/lpddr4_si/` in https://github.com/johonkanen/agilex3_soc_control_card.

## Limits

Lossless lines, no coupling in the simulation (no crosstalk; a differential pair is simulated as two uncoupled legs at its odd-mode impedance from `zsolve --gap`), no SSO/power-integrity, no AMI/jitter.
Good for termination, drive-strength, topology and trace-width decisions on short parallel buses; not a substitute
for a full SI tool on SerDes links.
