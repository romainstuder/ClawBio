# Stability Prediction Report

## Input

- **Structure**: `cftr.pdb` (chain A, 1480 residues)
- **Mutations**: 3
- **Methods requested**: rasp, thermompnn, foldx
- **Methods used**: rasp, thermompnn, foldx
- **Generated**: 2026-05-21 22:33 UTC

## Summary

ΔΔG values in kcal/mol. Positive = destabilizing.

| Mutation | rasp ΔΔG | thermompnn ΔΔG | foldx ΔΔG | Consensus | Direction | Confidence |
|---|---|---|---|---|---|---|
| A:F508A | +3.20 | +2.90 | +3.40 | +3.17 | **destabilizing** | high (3/3 methods) |
| A:R553Q | -0.80 | -0.40 | -1.10 | -0.77 | **neutral** | low (3/3 methods) |
| A:R1070W | -0.20 | -0.10 | -0.30 | -0.20 | **neutral** | high (3/3 methods) |

## Interpretation Guide

| ΔΔG (kcal/mol) | Meaning |
|---|---|
| > +1.0 | **Destabilizing** — protein likely misfolds or has reduced fold stability |
| -1.0 to +1.0 | **Neutral** — no major effect on folding stability |
| < -1.0 | **Stabilizing** — improves folding stability |

These thresholds approximate the energy of one hydrogen bond and follow common conventions in the directed-evolution literature.

## Method Agreement

- **Direction agreement**: 2/3 mutations
- **Magnitude agreement** (within 1.0 kcal/mol): 3/3
- **High-confidence calls**: 2/3

Disagreement between methods is informative, not a failure. Physics-based (FoldX) and ML-based (RaSP, ThermoMPNN) methods make different assumptions; where they diverge, the prediction warrants closer inspection.

## Flags

The following mutations have notes worth reviewing:

### A:R553Q
- Methods disagree on direction; review individual predictions

## Limitations

- Predictions assume a **rigid backbone**. Mutations affecting hinge motions, loop flexibility, or large conformational rearrangements may be under-predicted.
- All three methods predict **folding stability only**. Activity, binding affinity, and allosteric effects are not addressed.
- ML methods (RaSP, ThermoMPNN) reflect their training distributions; predictions on protein families heavily under-represented in PDB may be less reliable.
- FoldX is calibrated against directed-evolution datasets dominated by small, well-folded proteins.
- Deletions, insertions, and multi-residue changes are **not supported** in v1.
- Predictions are hypotheses, not experimental measurements.

## Methods

- **RaSP** — CNN trained on Rosetta ΔΔG across ~10M variants of 1,400 proteins. Fast, CPU-friendly, open source (Apache 2.0).
- **ThermoMPNN** — Graph transformer trained on the Megascale cDNA display dataset (~770K experimental ΔΔG measurements). State-of-the-art on the Megascale benchmark. Open source (MIT).
- **FoldX** — Empirical force field with terms for van der Waals, electrostatics, solvation, and backbone strain. Calibrated against directed-evolution data. Free for academic use; commercial license required otherwise.

## References

- Blaabjerg LM et al. (2023) eLife 12:e82593. doi:10.7554/eLife.82593
- Dieckhaus H et al. (2024) PNAS 121(6):e2314853121. doi:10.1073/pnas.2314853121
- Schymkowitz J et al. (2005) Nucleic Acids Res 33(W):W382-W388. doi:10.1093/nar/gki387
- Validation context: Studer RA, Christin P-A, Williams MA, Orengo CA (2014) Stability-activity tradeoffs constrain the adaptive evolution of RubisCO. PNAS 111(6):2223-2228. doi:10.1073/pnas.1310811111

## Reproducibility

This run produced the following replay artifacts:

- `reproducibility/commands.sh` — exact CLI invocation (missing)
- `reproducibility/environment.yml` — method versions and dependencies (missing)
- `reproducibility/output_checksum.txt` — SHA-256 of all outputs (missing)

**Replay**: `bash reproducibility/commands.sh`
**Verify**: `sha256sum -c reproducibility/output_checksum.txt`

## Disclaimer

This skill provides computational predictions for research planning only. It does not constitute clinical advice, and predictions should not be used as the sole basis for clinical, diagnostic, or therapeutic decisions. Experimental validation is required for any high-stakes use.

