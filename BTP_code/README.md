# Computational Evaluation and Prioritisation of De Novo Protein Binders Against α-Cobratoxin

This repository contains computational workflows for evaluating inherited de novo protein-binder designs against α-cobratoxin from *Naja kaouthia*. The analyses combine predicted complex confidence, receptor-footprint overlap, sequence-based specificity, cross-reactivity predictions, and interface solvent-accessible surface area.

> **Scope:** This project performs computational evaluation and candidate prioritisation. Predicted structures and confidence scores do not establish experimental binding, toxin neutralisation, or therapeutic efficacy.

## Project background

α-Cobratoxin is a long-chain three-finger neurotoxin that binds postsynaptic nicotinic acetylcholine receptors (nAChRs) at the neuromuscular junction. This interferes with acetylcholine-mediated muscle activation and can cause paralysis.

The project evaluates whether inherited designed binders are computationally compatible with the toxin and whether their predicted interfaces overlap the receptor-associated toxin surface.

The inherited design library was generated before this evaluation work. The analysis in this repository should not be interpreted as having generated those original binder sequences.

## Analysis overview

The major workflow is:

1. **Receptor-footprint reconstruction**
   - Uses PDB structure `1YI5`, containing α-cobratoxin bound to the acetylcholine-binding protein (AChBP), a receptor-like structural model.
   - Treats chains A–E as receptor-like subunits and chains F–J as toxin copies.
   - Counts heavy-atom contacts at a 4.0 Å cutoff.
   - Identifies toxin positions that repeatedly contact the receptor-like structure.

2. **Sequence conservation and hotspot analysis**
   - Uses multiple sequence alignment (MSA) to compare long-chain toxin sequences.
   - Compares long-chain and short-chain toxin sequence groups to estimate position-wise conservation differences.
   - Conservation differences are used as an exploratory measure of potential long-chain selectivity.

3. **Binder–toxin contact and receptor-footprint overlap**
   - Reads predicted binder–toxin complex PDB files.
   - Uses chain A as the binder and chain B as the toxin in the predicted complexes.
   - Detects heavy-atom contacts at a 4.0 Å cutoff.
   - Pairwise-aligns the complex toxin sequence against chain F of `1YI5` to map residue positions to the reference structure.
   - Compares mapped binder-contact positions with the receptor-associated footprint.

4. **Contact-weighted specificity**
   - Combines binder contact counts with the difference between long-chain and short-chain conservation.
   - Calculates a contact-weighted score for each candidate.
   - This score is a computational prioritisation metric, not an experimentally measured specificity or binding probability.

5. **Boltz-2 cross-reactivity evaluation**
   - The completed analysis comprised 16 binders evaluated against a 19-member long-chain toxin panel (304 predictions), plus 16 inherited-target control predictions.
   - Total: 320 predicted complexes, processed in 10 batches of 32.
   - The results were summarised using interface confidence metrics, including ipTM.
   - ipTM describes confidence in the predicted relative arrangement of interacting chains; it is not a direct measurement of binding affinity.

6. **Interface SASA/BSA**
   - Calculates solvent-accessible surface area (SASA) for isolated binder and toxin chains and for the complete complex.
   - Uses the Shrake–Rupley method with a 1.4 Å probe radius and 200 sphere points.
   - Calculates total buried surface area (BSA) as:

     `BSA_total = SASA_binder + SASA_toxin - SASA_complex`

   - Reports BSA per partner as `BSA_total / 2`.
   - The recorded BSA analysis was performed on three representative candidates, not the entire shortlist.

## Repository contents

The consolidated script is:

- `btp_major_analyses.py` — contains the principal analysis functions and command-line entry points.

The script includes the following sections:

| Section | Purpose |
|---|---|
| Common helpers | PDB parsing, residue conversion, heavy-atom selection |
| Receptor footprint | Reconstruct toxin–AChBP contact positions in 1YI5 |
| Clustal Omega wrapper | Run MSA on a supplied FASTA file |
| Pairwise alignment and overlap | Map binder-contact residues to the 1YI5 reference and calculate footprint overlap |
| Long/short conservation | Calculate position-wise conservation and specificity differences from aligned FASTA files |
| Contact-weighted specificity | Combine contact counts and conservation differences |
| SASA/BSA | Calculate buried surface area for representative complexes |
| Boltz-2 manifest | Create a 320-job prediction manifest |
| Boltz-2 result processing | Create cross-reactivity matrices and candidate/toxin summaries |

## Biological and structural references

- **1YI5:** α-Cobratoxin bound to acetylcholine-binding protein. The AChBP complex is used as a receptor-like structural reference; it is not the complete human neuromuscular nAChR.
- **1CTX:** Free α-cobratoxin structure, used as an additional structural reference in the project.
- **Target toxin:** α-Cobratoxin, UniProt accession `P01391`.

## Recorded computational environment

The Boltz-2 predictions were run on the `megatron` server using an NVIDIA L40S GPU.

| Component | Recorded information |
|---|---|
| GPU | NVIDIA L40S |
| GPU memory | Approximately 46 GB |
| NVIDIA driver | 580.178.04 |
| CUDA version reported by `nvidia-smi` | 13.0 |
| Python in Boltz environment | 3.10 |
| Boltz environment | `/data/guest_anirban/boltz_env` |
| MSA server | `https://api.colabfold.com` |
| MSA pairing strategy | `greedy` |
| Maximum MSA sequences | 16,384 |

The exact Boltz-2 package version, Biopython version, PyTorch version, and total wall-clock time should be recorded from the active environment and run logs. The CUDA version shown by `nvidia-smi` is the driver-supported CUDA version and does not necessarily equal the CUDA runtime against which PyTorch was built.

## Input data and expected layout

The script was developed around the following server-side layout. Update the paths near the top of `btp_major_analyses.py` if your directory structure differs.

```text
/data/guest_anirban/
├── structures/
│   ├── 1YI5.pdb
│   └── 1CTX.pdb
└── protein_design/
    └── long_chain_design/
        ├── AlphaFold_designs_121_A28_33/
        ├── AlphaFold_designs_121_A30_34/
        ├── AlphaFold_designs_125_A36_39/
        ├── AlphaFold_designs_125_A53_56/
        ├── long_chain_panel_19.fasta
        ├── long_chain_MSA.fa
        └── short_chain_mature_sequences.fasta
```

Expected predicted-complex convention:

- Chain A: binder
- Chain B: toxin
- File suffix: `_af2pred.pdb`

The script's default paths are project-specific. Confirm all paths and chain identifiers before execution.

## Requirements

Core Python packages used by the script:

- Python 3.10 or compatible
- Biopython
- NumPy
- pandas
- SciPy

Additional requirements depend on the selected workflow:

- Clustal Omega executable (`clustalo`) for MSA generation.
- Boltz-2 installed in its intended environment for prediction generation. The consolidated script does not recreate the original Boltz-2 YAML files or prediction command because the exact installed-version YAML schema and original invocation were not preserved in the recovered code.

Install common Python dependencies in an appropriate environment:

```bash
python -m pip install biopython numpy pandas scipy
```

Check installed versions:

```bash
python - <<'PY'
import sys
import Bio
import numpy
import pandas
import scipy

print("Python:", sys.version)
print("Biopython:", Bio.__version__)
print("NumPy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("SciPy:", scipy.__version__)
PY
```

For the Boltz-2 environment, also check:

```bash
python - <<'PY'
import torch
import importlib.metadata as md

print("PyTorch:", torch.__version__)
print("PyTorch CUDA runtime:", torch.version.cuda)

for package in ("boltz", "biopython"):
    try:
        print(package, md.version(package))
    except md.PackageNotFoundError:
        print(package, "not found")
PY

nvidia-smi
```

## Running the analysis script

Run commands from the directory containing `btp_major_analyses.py`, after editing the configured paths if needed.

### 1. Reconstruct receptor footprint

```bash
python btp_major_analyses.py footprint
```

Default output:

- `receptor_contacts_1YI5.csv`

### 2. Run Clustal Omega

```bash
python btp_major_analyses.py clustalo input_sequences.fasta aligned_sequences.fasta
```

This calls the installed `clustalo` executable with FASTA output.

### 3. Calculate binder–receptor-footprint overlap

```bash
python btp_major_analyses.py overlap
```

Default outputs:

- `epitope_overlap_summary.csv`
- `binder_residue_contacts.csv`

**Important before running:** Verify the consensus footprint residue set in the script against the final `consensus_receptor_footprint.csv`. The script contains a provisional hard-coded set based on the recovered notes, and the notes indicate that one of the 19 consensus positions needs confirmation from the original CSV.

Also verify the PairwiseAligner scoring parameters against the original analysis script before treating them as exact historical settings. The current values in the consolidated script are documented defaults for this reconstruction, not independently verified original-run provenance.

### 4. Calculate contact-weighted specificity

```bash
python btp_major_analyses.py specificity
```

Expected inputs:

- `binder_residue_contacts.csv`
- `long_vs_short_specificity_profile_final.csv`

Default outputs:

- `candidate_specificity_summary.csv`
- `candidate_footprint_specificity_details.csv`

The long-chain and short-chain alignments must already be mapped to the same reference residue coordinate system before calculating the conservation profile. Do not compare unrelated alignment-column numbers directly.

### 5. Calculate representative candidate BSA

```bash
python btp_major_analyses.py sasa
```

Default output:

- `representative_bsa_results.csv`

The script uses three representative PDB paths defined in the function. Check that these files exist before running.

## Boltz-2 prediction workflow and result processing

The recorded completed prediction campaign consisted of:

- 16 binder sequences
- 19 long-chain toxin sequences
- 304 binder–panel combinations
- 16 inherited-target controls
- 320 total predictions
- 10 batches of 32

The consolidated script contains two supporting functions:

1. `create_boltz_manifest(...)` — creates a CSV manifest describing the binder–toxin combinations and batch assignments.
2. `summarize_boltz_results(...)` — reads a consolidated prediction-results CSV and generates the ipTM matrix and summary tables.

The manifest function deliberately does **not** generate Boltz YAML files. The exact YAML input schema and original successful CLI invocation were not recoverable from the available project notes. Use the schema documented for the installed Boltz-2 version rather than assuming a generic YAML format.

The result-processing function expects a CSV with at least these columns:

```text
binder_id,toxin_id,iptm
```

It optionally accepts:

```text
iplddt,ipde,job_type
```

The `job_type` column should identify panel predictions and inherited-target controls, for example:

```text
panel_cross_reactivity
inherited_target_control
```

Default outputs:

- `cross_reactivity_iptm_matrix.csv`
- `candidate_crossreactivity_summary.csv`
- `toxin_crossreactivity_summary.csv`
- `inherited_reference_scores.csv`

If your actual Boltz output uses different column names, rename or map them before calling the function.

## Interpretation and limitations

### Structural confidence metrics

- pLDDT estimates local structural confidence.
- PAE estimates uncertainty in relative residue positioning.
- ipTM estimates confidence in the relative arrangement of interacting chains.
- These are model-confidence measures, not experimental binding measurements.

### Contact analysis

A 4.0 Å heavy-atom distance cutoff identifies close structural contacts. Atom-pair counts are not interaction energies and do not independently establish chemical complementarity or binding strength.

### Sequence specificity

The long-chain versus short-chain conservation difference is an exploratory sequence-based discriminator. It depends on the selected sequence panel, sequence quality, alignment register, and sample size. It is not an experimentally measured specificity value.

### BSA

Buried surface area estimates how much solvent-accessible surface becomes buried when the two chains form the predicted complex. It does not independently establish affinity or interface stability.

### Cross-reactivity predictions

Boltz-2 ipTM values across a toxin panel describe predicted complex-interface confidence. They should not be interpreted as binding probabilities or direct affinity measurements. MSA depth and model limitations can affect prediction confidence.

### Experimental validation still required

The computational workflow prioritises candidates for subsequent validation. Useful future steps include:

- Expression and purification of selected binders.
- Direct binding measurements such as SPR, BLI, or ITC.
- Receptor-competition assays to test whether a binder interferes with toxin–receptor interaction.
- Functional assays to test whether toxin activity is reduced.
- MD simulations to investigate the time-dependent behaviour of selected predicted complexes.
- Rosetta interface-energy calculations, such as mutation-based `flex ddG`, to explore possible interface optimisation.

A predicted binder–toxin complex is not evidence of toxin neutralisation. Functional neutralisation requires appropriate experimental testing.

## Reproducibility checklist

Before reporting or rerunning results, record:

- [ ] Exact input FASTA and PDB filenames.
- [ ] Exact chain IDs used for each structure.
- [ ] Exact Clustal Omega version and command.
- [ ] Exact Biopython version and PairwiseAligner scoring parameters.
- [ ] Final consensus receptor-footprint residue list and contact cutoff.
- [ ] Exact Boltz-2 package version and prediction configuration.
- [ ] PyTorch version and `torch.version.cuda`.
- [ ] NVIDIA driver and GPU model.
- [ ] MSA settings and server.
- [ ] Prediction batch logs and total wall-clock time.
- [ ] Output CSVs and prediction files.
- [ ] Any missing, failed, or excluded predictions.

## Citation notes

Relevant methodological references include:

- Bourne et al. (2005), structural basis of α-cobratoxin binding to acetylcholine-binding protein, *EMBO Journal*. DOI: `10.1038/sj.emboj.7600620`.
- Jumper et al. (2021), highly accurate protein structure prediction with AlphaFold, *Nature*. DOI: `10.1038/s41586-021-03819-2`.
- Watson et al. (2023), de novo design of protein structure and function with RFdiffusion, *Nature*. DOI: `10.1038/s41586-023-06415-8`.
- Dauparas et al. (2022), robust deep learning-based protein sequence design using ProteinMPNN, *Science*. DOI: `10.1126/science.add2187`.
- Shrake and Rupley (1973), solvent accessibility calculations, *Journal of Molecular Biology*. DOI: `10.1016/0022-2836(73)90011-9`.
- Passaro et al. (2025), Boltz-2 preprint. DOI: `10.1101/2025.06.14.659707`.

---

**Project status:** Computational screening, interface analysis, sequence-based specificity analysis, cross-reactivity prediction, and representative BSA calculations were completed as described in the project records. Experimental binding and functional neutralisation remain future validation stages.
