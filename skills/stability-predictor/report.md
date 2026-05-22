# Stability Prediction Report

## Input

- **Structure**: `cli_demo_pdb.pdb` (chain A (of A), 1 residues)
- **Mutations**: 3
- **Methods requested**: mock
- **Methods used**: mock
- **Generated**: 2026-05-21 22:40 UTC

## Summary

ΔΔG values in kcal/mol. Positive = destabilizing.

| Mutation | mock ΔΔG | Consensus | Direction | Confidence |
|---|---|---|---|---|
| A:A328S | -3.00 | -3.00 | **stabilizing** | medium (1/1 methods) |
| A:L270I | +3.00 | +3.00 | **destabilizing** | medium (1/1 methods) |
| A:F508A | +3.00 | +3.00 | **destabilizing** | medium (1/1 methods) |

## Interpretation Guide

| ΔΔG (kcal/mol) | Meaning |
|---|---|
| > +1.0 | **Destabilizing** — protein likely misfolds or has reduced fold stability |
| -1.0 to +1.0 | **Neutral** — no major effect on folding stability |
| < -1.0 | **Stabilizing** — improves folding stability |

These thresholds approximate the energy of one hydrogen bond and follow common conventions in the directed-evolution literature.

## Method Agreement

- **Direction agreement**: 3/3 mutations
- **Magnitude agreement** (within 1.0 kcal/mol): 3/3
- **High-confidence calls**: 0/3

Disagreement between methods is informative, not a failure. Physics-based (FoldX) and ML-based (RaSP, ThermoMPNN) methods make different assumptions; where they diverge, the prediction warrants closer inspection.

## Flags

No mutations flagged for review.

## Limitations

- Predictions assume a **rigid backbone**. Mutations affecting hinge motions, loop flexibility, or large conformational rearrangements may be under-predicted.
- All three methods predict **folding stability only**. Activity, binding affinity, and allosteric effects are not addressed.
- ML methods (RaSP, ThermoMPNN) reflect their training distributions; predictions on protein families heavily under-represented in PDB may be less reliable.
- FoldX is calibrated against directed-evolution datasets dominated by small, well-folded proteins.
- Deletions, insertions, and multi-residue changes are **not supported** in v1.
- Predictions are hypotheses, not experimental measurements.

## Methods


## References

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
