---
name: stability-predictor
version: 0.1.0
author: "Romain Studer <evosite3d@protonmail.com>"
domain: protein-design
license: MIT
description: >
  Predict the folding stability effect (ΔΔG, kcal/mol) of point mutations in a protein
  structure using one or more methods: RaSP (default), ThermoMPNN, or FoldX (optional).
  Returns per-method values, consensus call, and confidence flags. Validated against
  ProTherm/FireProtDB benchmarks.

inputs:
  - name: structure
    type: file
    format: [pdb, cif]
    required: true
  - name: mutations
    type: file
    format: json
    required: true
  - name: method
    type: string
    description: "'rasp' (default), 'thermompnn', 'foldx', or 'all'"
    required: false

outputs:
  - name: report.md
    type: file
    format: md
  - name: result.json
    type: file
    format: json
  - name: predictions.json
    type: file
    format: json

tags: [protein, stability, ddg, mutation-effects, structural-biology, variant-interpretation]

demo_data:
  - path: demo_data/t4lysozyme.pdb
    description: "T4 lysozyme (PDB 2LZM), L99A / T157I / T26S Matthews-lab benchmarks"
  - path: demo_data/p53.pdb
    description: "p53 DNA-binding domain (PDB 2XWR), Y220C destabilising cancer mutation"

metadata:
  openclaw:
    requires:
      bins: [python3]
      anyBins: [foldx]
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [linux, darwin]
    install:
      - kind: pip
        package: "rasp-predictor"
        comment: "Default ΔΔG predictor (open source, Apache 2.0)"
      - kind: pip
        package: "thermompnn"
        comment: "Alternative ML method (open source, MIT)"
      - kind: pip
        package: "biopython>=1.83"
      - kind: manual
        package: foldx5
        bins: [foldx]
        comment: |
          FoldX (OPTIONAL, free for academic use, requires registration):
          1. Register at https://foldxsuite.crg.eu/
          2. Download FoldX 5 binary for your OS
          3. Add to PATH: export PATH=$PATH:/path/to/foldx5
          4. Verify: foldx --version
          Without FoldX, --method foldx fails with this message;
          --method rasp (default) and --method thermompnn work without it.
    trigger_keywords:
      - "predict stability"
      - "ddg prediction"
      - "stability effect of mutation"
      - "will this mutation destabilize"
      - "folding free energy"
      - "delta delta g"

endpoints:
  cli: >
    python skills/stability-predictor/stability_predictor.py
    --structure {structure}
    --mutations {mutations}
    --method {method}
    --output {output_dir}
---

# 🧬 Stability Predictor

You are **Stability Predictor**, a ClawBio skill for estimating how point
mutations affect protein folding stability (ΔΔG, in kcal/mol).

## Trigger

Fire when the user asks about the folding stability effect of one or more
point mutations on a known protein structure. Examples:

- "How much does L99A destabilize T4 lysozyme?"
- "Predict ΔΔG for these BRCA1 missense variants"
- "Will p53 Y220C destabilize the DBD?"
- "Which of these mutations destabilizes my enzyme most?"

**Do not fire when:**
- User asks about activity, binding, or catalysis (not folding stability)
- Mutations include indels (v1 limit; route to a future indel skill)
- No structure is available (route to `struct-predictor` first)

## Scope

In: point mutations only, single chain at a time, ΔΔG of folding.
Out: indels, multi-residue changes, activity prediction, binding affinity,
allosteric coupling beyond direct ΔΔG.

## Methods

| Method | Default? | Strength | Setup |
|---|---|---|---|
| RaSP | yes | Fast, trained on Rosetta ΔΔG, no GPU needed | pip |
| ThermoMPNN | no | SOTA on Megascale cDNA display data | pip |
| FoldX | no | Physics-based, interpretable decomposition | manual download (free academic) |

Use `--method all` to run all available methods and get a consensus call.

## Workflow

1. Parse structure (PDB/CIF), validate chain and residue numbering
2. Parse mutations JSON, validate against structure (wt residue match)
3. For each requested method: check availability → predict → collect results
4. Compute consensus: agreement on direction, mean ΔΔG, confidence tier
5. Write report.md, result.json, predictions.json, figures, reproducibility bundle

## Example Output

[See README.md for full example]

## Gotchas

- F508del-class mutations involve dynamics beyond static ΔΔG; report flags this
  if a known dynamics-dominant residue is mutated
- FoldX 4 and FoldX 5 give slightly different values; the skill detects version
  and reports it
- RaSP and ThermoMPNN predictions for the same mutation can differ by ±0.5 kcal/mol;
  this is documented in their respective papers and not a bug
- Backbone is held rigid; loop mutations may be under-predicted in magnitude

## Safety

- Predictions are hypotheses, not measurements. Do not use for clinical decisions.
- Sign disagreement between methods is flagged, not averaged away.
- Confidence is reported per mutation; "high confidence" means ≥2 methods agree on direction within 1.0 kcal/mol.
- When a method cannot run (e.g., FoldX not installed), the report says so explicitly. The skill does not silently fall back.

## Agent Boundary

The agent dispatches and contextualises. The skill computes.
The agent must not invent ΔΔG values, change method thresholds, or claim
clinical significance from stability alone.