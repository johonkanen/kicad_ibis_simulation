# Delay-matching rules in KiCad 10 (time domain): instructions for Claude

Use this when a bus (DDR/LPDDR, a parallel interface) must be matched in time, so that pcbnew's tuning patterns and
DRC hold the targets. The goal is that **green means good**: a net inside its rule meets the timing checks of the
part's design guide. The IBIS simulation then checks the routed result. Getting there needs four things to agree:
KiCad's delay model, your own model of the board, the DRC rules and the tuning patterns.

Reference implementation (LPDDR4, Efinix Ti135 + one x32 die, T-branched CA/CK):
[johonkanen/efinix_soc_board](https://github.com/johonkanen/efinix_soc_board) `tools/ddr/`: `kicad_delays.sh`,
`ddr_timing.py`, `delay_rules.sh`.

## 1. How KiCad counts delay

- **Tuning profile** (Board Setup > Tuning Profiles, saved in the `.kicad_pro` under `tuning_profiles`): one
  propagation delay per signal layer (ps/mm) and one `via_prop_delay` (ps per mm of barrel). Take the layer values
  from a field solve of the real stack-up (`python3 -m kicad_ibis.zsolve`), not from a rule of thumb.
- **A via** counts the stack-up height between the outermost layers it joins, times `via_prop_delay`. On a 1.6 mm
  8-layer board at 7.2 ps/mm that is F.Cu->In2 1.8 ps, F.Cu->In5 9.4 ps and F.Cu->B.Cu 11.1 ps. A flat "x ps per via"
  is wrong by several ps, and wrong by a different amount per layer pair, which matters when a byte's DQS and DQ run
  on different layers.
- **Pad-to-die delay** (pad property) adds the package delay of a pad. Set it on the device pads from the vendor's
  package data. A part without public package data gets 0; say so.
- **DRC** (`length` constraint in ps) sums **all** of the net's copper: every branch, stub and loop.
- **The tuner** (tuning pattern) reads the **pad-to-pad path** through the pattern's track.
- On a point-to-point net with clean copper the two agree. On a **branched** net (a DDR CA/CK T or fly-by) DRC reads
  the path plus the other branches. A stray stub or a duplicated segment makes DRC read more than the tuner. That
  can hide a real stub: look for it, don't widen the rule.

## 2. Read KiCad's own numbers

Get the delay per net from KiCad itself with a probe rule on a copy of the board. Add a rule that every matched net
fails, run DRC headless, and parse the `actual` value of each `length_out_of_range` item:

```
(rule "probe_delay"
	(condition "A.hasNetclass('lpddr_dq') || A.hasNetclass('lpddr_ca') || A.hasNetclass('lpddr_diff')")
	(constraint length (min 1ps) (max 2ps)))
```

`kicad-cli pcb drc --severity-all -o out.rpt copy.kicad_pcb` then gives "actual N ps" per net (`kicad_delays.sh`).
The copy needs the project's `.kicad_pro` and `.kicad_dru` next to it under the same base name, or the tuning
profile is missing.

## 3. Your model must reproduce KiCad

Build your own per-net model from the board (pcbnew Python): a copper graph per net, with per-layer delays from the
same profile, vias by stack-up height (layer depths from the board's `(stackup ...)`), and pad-to-die on the pads.
Then take the shortest path from the driver pad to each receiver pad. Use it to compute targets and the path of
branched nets.

**Check it before using it.** Print KiCad readout minus model path for every net. On point-to-point nets it must be
~0 (within about 1 ps). A constant offset on all nets of one layer means the via model is wrong. A single net far off
means extra copper (a stub or loop) on that net: fix the board. On a branched net the difference is that net's other
branches.

## 4. Rules that make green mean good

- **Generate the rules, never hand-edit them.** Write one `length` rule per net or group into a marked block at the
  end of the `.kicad_dru` (`# ---- BEGIN generated ... ----` / `END`), and regenerate after every reroute. Later
  rules win in KiCad, so a hand exception placed before the block can be overridden by it; put net-specific
  exceptions for other constraint types (widths) before the block.
- **Targets follow the longest net of a group.** Meanders only add delay, so the group target is at least the
  longest member. Add the guide's ordering constraints, for example Efinix: every DQS shorter than CK, so the CA/CK
  target is at least the latest DQS + 5 ps.
- **Point-to-point nets:** `(min T+lo) (opt T+opt) (max T+hi)`, all in ps. The tuner sizes the meander to opt and
  shows the pattern green anywhere between min and max.
- **Branched nets (T or fly-by):** one rule must serve the tuner (path) and DRC (whole net). Put **min and opt in
  path terms** and **max = path target + tol + that net's other branches** (KiCad readout minus model path, measured
  on the routed board). The tuner, and pcbnew's "Update all tuning patterns", aim at opt, so the meander lands on the
  path target. A path below min shows in the tuner (pattern not green). DRC (whole net = path + branch) shows a net
  that is too long. **Good = the pattern is green after an update and DRC is green.**
  **Check the tuner against your model on every branched net.** On some nets the tuner's path takes in part of the
  other branch (seen: +11 ps, about one pad via, and +3 ps; it depends on how the T is drawn), so those patterns read
  "tuned" while the real path is short. Each pattern stores its reading in `last_tuning_length` (mm); map it to ps
  with the pattern's own `target_length_min/max` and `target_delay_min/max` pairs. Shift that net's min/opt by
  (tuner - path). Readings are as fresh as the last "Update all tuning patterns": regenerate the rules right after an
  update and save, then update again. Do not put the rule in whole-net
  terms and override the tuning patterns with the path target (`override_custom_rules`): when the GUI regenerates
  the patterns it tunes them to the rule and ignores the override (KiCad 10.0.7, seen 2026-10-31: CK and CKE grew by
  ~24 ps).
- **Skew rules** (`skew ... within_diff_pairs`) compare whole nets. On a branched pair (CK) the branches differ, so
  the rule can never pass. Leave that pair out of the skew rule and hold P and N with their own per-net windows:
  +-tol each gives P/N <= 2 tol, so choose tol below half the guide's limit (CK 0.9 ps -> +-0.4 ps).
- **Centre what should be centred.** Group members around a common target (CA vs CK): opt at the target, +-tol
  symmetric, with 2 tol inside the guide's limit. For a window that is one-sided in the guide (DQ - DQS 0..66 ps),
  aim near the low edge (opt DQS + 10) and keep the max tight (DQS + 15). Then every DQ in a byte sits within 10 ps
  of the others. A wide max lets meanders drift and turns into DQ-to-DQ skew while everything shows green.
- **Tolerances are the margin between the guide's limit and your model error.** Prove the model error small
  (section 3) before you spend the margin.

## 5. Setting up the tuner in pcbnew

What the user does in the GUI, and what to tell them. Observed with KiCad 10.0.7.

### Configuration that worked (LPDDR4, 8-layer 1.6 mm, 52 patterns)

| Where | Setting |
|---|---|
| Board Setup > Tuning Profiles | One profile per impedance class, time-domain tuning enabled: `LPDDR4 SE 40R` (single-ended) and `LPDDR4 DIFF 80R` (pairs, with the pair gap per layer). Per signal layer: its reference planes and its delay in ps/mm from the field solve (here F/B.Cu 6.11, In1/In6 6.97, inner striplines 7.00-7.01). Via delay 7.2 ps/mm of barrel. |
| Net classes (`net_settings.classes[]`) | `tuning_profile` set on every bus class (`lpddr_dq`, `lpddr_ca` -> SE profile; `lpddr_diff` -> DIFF profile). Without it the class has no time-domain profile (no ps targets). `clearance` = the crosstalk spacing (0.2 mm), with DRC exceptions below it (section 5, spacing). |
| Board Setup > tuning pattern defaults (`tuning_pattern_settings`) | Single track and skew: min amplitude 0.2 mm, max 1.0 mm, spacing 0.6 mm, rounded corners 80 %, two-sided. Diff pair: same, spacing 1.0 mm. |
| Each pattern | Time domain on; Override custom rules off (targets from the generated rules); rounded, two-sided. Where space was short the user lowered max amplitude (0.5-0.8 mm) and spacing (0.3-0.5 mm) on that pattern; long pair runs needed max amplitude up to 2.6 mm. |
| Pattern types | Single-track pattern on each single-ended net and on each leg's trunk where needed; diff-pair *length* pattern for a pair's delay; skew patterns only for small P/N corrections (they do not set the pair length). |
| Interactive router | Walk around (or Shove). |
| Rules file | Generated `length` rules per net/group (section 4); spacing rules whose BGA exceptions test both items (below). |

- **One tuning pattern per net to tune.** Single-ended nets: a single-track pattern. Differential pairs: a diff-pair
  *length* pattern for the pair's delay. A diff-pair *skew* pattern only matches P to N and never sets the pair's
  length (DQS stayed at its old length while its rule asked for +30 ps). An update also rebuilt one skew pattern
  wrongly: its meander was removed and P/N skew went from 1 ps to 9 ps. Prefer per-leg length rules (+-tol each) to
  skew patterns.
- **Leave "Override custom rules" off.** The patterns then take min/opt/max from the generated rules. An override is
  ignored by "Update all tuning patterns" (section 4).
- **The work loop:** regenerate the rules (script) -> the user runs "Update all tuning patterns" and saves -> check
  with your model and DRC (`--refill-zones`) -> regenerate the rules if a tuner offset or target moved -> update again.
  Rules generated from a board whose patterns are stale (copper changed after the last update) carry stale tuner
  readings.
- **A net with no pattern is not touched by an update.** List those nets with how much to add or remove (ps and mm
  on their layer). The user adds a pattern or edits the route. A route that is too long and hard to shorten: raise
  the reference instead (for a byte, the DQS target, so the longest DQ sits at opt), and lengthen the rest.
- **"too_short" after an update means the meander ran out of room.** Raise the pattern's max amplitude, stretch it
  over a longer piece of track, or move it.
- **Spacing to other nets: the router and the pattern updates hold custom clearance rules, but check every
  exception's condition with A and B swapped.** A BGA exception written as "A touches the courtyard" also covers a
  long segment that only starts under the BGA, along its whole length. The router evaluated the rule with that
  segment as A and used the 0.1 mm exception; DRC took the other order and flagged 0.2 mm. Every meander next to such
  a segment came out 0.12-0.19 mm from it, whatever the severity of the spacing rule. Require both items on the same
  courtyard:

  ```
  (rule "bga_clearance"
  	(condition "<bus netclasses on A> && ((A.intersectsCourtyard('U1') && B.intersectsCourtyard('U1')) || (A.intersectsCourtyard('U2') && B.intersectsCourtyard('U2')))")
  	(constraint clearance (min 0.1mm)))
  ```

  and the same "both" test in the exemption of the spacing rules. After that fix the meanders held both a 0.2 mm
  rule inside a group and a 0.3 mm rule between groups. A tuning pattern whose start lies inside a BGA courtyard
  still got the BGA value: move the pattern's start outside the courtyard. Putting the spacing into the netclass
  clearance as well (`net_settings.classes[].clearance` in the `.kicad_pro`, edited with pcbnew closed) is a
  backstop; it then needs DRC rules that allow less under the BGAs (as above), between P and N
  (`A.isCoupledDiffPair()`) and to copper outside the bus. Run DRC before and after such a change.
- **Router mode:** Walk around (or Shove), not "Highlight collisions", which places a meander into a clearance and
  only marks it.

## 6. Pitfalls

- A tuning pattern's colour is the status saved when it was last updated (`last_status`), not live. Run "Update
  all tuning patterns" and DRC after any change before calling the bus tuned.
- Commit (or copy) the board before running an update or DRC in the GUI. In one session (KiCad 10.0.7) every tuning
  pattern was turned into plain tracks with new IDs; the patterns could not be restored onto the new tracks.
- KiCad DRC reports at most 199 items of one type.
- Unfilled zones give pages of clearance errors on a copy: run `kicad-cli pcb drc --refill-zones` when comparing
  two boards.
- Tiny leftover track pieces (0.05 mm) next to a pattern or via can make the tuner read a different path. Look at
  the copper of any net whose readout and model disagree.
- Before writing the `.kicad_pcb`: the user may have it open. Work on a copy, check that the real file is unchanged
  since you copied it (`cmp`), and only then copy back.
- Report what is off as numbers per net (now -> target, ps), plus how much meander that is on its layer (ps / ps/mm).
