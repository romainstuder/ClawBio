# 🧬 Stability Predictor

Predict the folding stability effect (ΔΔG, kcal/mol) of point mutations on a
protein structure using one or more of:

| Method     | Default | Source                                                                    | Setup                                    |
|------------|---------|---------------------------------------------------------------------------|------------------------------------------|
| RaSP       | ✅      | Blaabjerg et al. (2023) eLife. CNN trained on Rosetta ΔΔG (Apache 2.0)    | `pip install rasp-predictor`             |
| ThermoMPNN |         | Dieckhaus et al. (2024) PNAS. Graph transformer on Megascale data (MIT)   | `pip install thermompnn`                 |
| FoldX 5    |         | Schymkowitz et al. (2005) NAR. Empirical force field (free for academics) | manual download — see install_instructions |

Positive ΔΔG = destabilising. Negative = stabilising. The skill returns
per-method values, a consensus call, a confidence tier, and a reproducibility
bundle.

## 30-second quick start

```bash
# T4 lysozyme demo — L99A / T157I / T26S Matthews-lab benchmarks (default)
python skills/stability-predictor/stability_predictor.py --demo --output /tmp/sp_demo

# p53 demo — Y220C destabilising cancer mutation
python skills/stability-predictor/stability_predictor.py \
  --demo --demo-set p53 --output /tmp/sp_p53

# Real run — your own structure + mutations
python skills/stability-predictor/stability_predictor.py \
  --structure my_protein.pdb \
  --mutations my_mutations.json \
  --method rasp \
  --output /tmp/sp_out
```

Demo runs complete in seconds: predictions are served from a bundled
literature-grounded cache, so RaSP / ThermoMPNN / FoldX do not need to be
installed to see the pipeline end-to-end.

## Mutation input format

`mutations.json`:

```json
{
  "chain": "A",
  "mutations": [
    {"position": 99,  "wt": "L", "mt": "A"},
    {"position": 157, "wt": "T", "mt": "I"}
  ]
}
```

`chain` at the top level is the default; each mutation may override it.

## Output

```
output_dir/
├── report.md                   # Primary deliverable (markdown)
├── result.json                 # Top-level summary
├── predictions.json            # Per-mutation × per-method detail
├── figures/
│   ├── ddg_heatmap.png         # Mutations × methods, colour = ΔΔG
│   ├── per_mutation_bars.png   # Grouped bars per mutation
│   └── method_agreement.png    # Scatter of method-vs-method ΔΔG (if ≥2 methods)
└── reproducibility/
    ├── commands.sh             # Exact replay command
    ├── environment.yml         # Method versions
    ├── input_checksum.txt      # SHA-256 of structure + mutations
    └── output_checksum.txt     # SHA-256 of every output file
```

## Methods cheatsheet

```bash
# Default (RaSP)
--method rasp

# Single alternative
--method thermompnn
--method foldx

# Run everything available and consensus across methods
--method all
```

Methods that aren't installed are skipped with a warning; if none are
available the run exits with code 4 and prints install instructions.

## Demo data

| Demo set     | Structure | Mutations              | Source                                                |
|--------------|-----------|------------------------|-------------------------------------------------------|
| `t4lysozyme` | PDB 2LZM  | L99A, T157I, T26S      | Matthews-lab benchmark (Eriksson 1992; Pjura 1993)    |
| `p53`        | PDB 2XWR  | Y220C                  | Bullock et al. 1997; Joerger & Fersht 2007            |

Both demos ship the full crystal structures (downloaded from the RCSB
PDB) plus per-method prediction caches under `demo_data/`. The cache
values are literature-grounded illustrations — they reproduce the
published direction and magnitude (L99A ≈ +5 kcal/mol, Y220C ≈ +3
kcal/mol) but are not new measurements.

## Validation

- **Tier 1** (every commit): the bundled `tests/` suite runs the CLI on
  both demo sets and verifies all expected outputs.
- **Tier 2** (every commit): consensus logic, IO parsing, method
  availability are unit-tested.
- **Tier 3** (`pytest -m slow`): correlation against ProTherm /
  FireProtDB benchmarks. Targets per upstream publications: RaSP
  Spearman ≥ 0.55, ThermoMPNN ≥ 0.60, FoldX ≥ 0.50. (Benchmark fetcher
  not implemented in the v1 scaffold — see `tests/test_validation.py`.)

## Citations

- Blaabjerg LM et al. (2023) "Rapid protein stability prediction using
  deep learning representations." eLife 12:e82593.
  doi:[10.7554/eLife.82593](https://doi.org/10.7554/eLife.82593)
- Dieckhaus H et al. (2024) "Transfer learning to leverage larger
  datasets for improved prediction of protein stability changes." PNAS
  121(6):e2314853121.
  doi:[10.1073/pnas.2314853121](https://doi.org/10.1073/pnas.2314853121)
- Schymkowitz J et al. (2005) "The FoldX web server: an online force
  field." Nucleic Acids Res 33(W):W382-W388.
  doi:[10.1093/nar/gki387](https://doi.org/10.1093/nar/gki387)
- Eriksson AE, Baase WA, Zhang X-J, Heinz DW, Blaber M, Baldwin EP,
  Matthews BW (1992) "Response of a protein structure to cavity-creating
  mutations and its relation to the hydrophobic effect." Science
  255(5041):178-183.
  doi:[10.1126/science.1553543](https://doi.org/10.1126/science.1553543)
- Joerger AC, Fersht AR (2007) "Structure-function-rescue: the diverse
  nature of common p53 cancer mutants." Oncogene 26(15):2226-2242.
  doi:[10.1038/sj.onc.1210291](https://doi.org/10.1038/sj.onc.1210291)
- Studer RA, Christin P-A, Williams MA, Orengo CA (2014)
  "Stability-activity tradeoffs constrain the adaptive evolution of
  RubisCO." PNAS 111(6):2223-2228.
  doi:[10.1073/pnas.1310811111](https://doi.org/10.1073/pnas.1310811111)

*ClawBio is a research and educational tool. It is not a medical device
and does not provide clinical diagnoses.*
