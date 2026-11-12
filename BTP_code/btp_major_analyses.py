#!/usr/bin/env python3
"""BTP major computational analyses: receptor footprint, sequence mapping,
contact-weighted specificity, SASA/BSA, and Boltz-2 result organisation.

The core structural-analysis settings below reproduce the recorded workflow.
The exact original Boltz-2 YAML schema/CLI invocation was not recoverable, so
this file creates its job manifest and post-processes exported scores rather
than inventing a supposedly exact YAML format or run command.

Recorded environment: Python 3.10; Boltz environment
/data/guest_anirban/boltz_env; NVIDIA L40S (~46 GB); ColabFold MSA server
https://api.colabfold.com; greedy pairing; max MSA sequences 16384.
"""
from pathlib import Path
from collections import Counter
import argparse
import subprocess
import numpy as np
import pandas as pd
from Bio import Align
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

# -----------------------------------------------------------------------------
# SECTION 0 — PROJECT PATHS AND ANALYSIS SETTINGS
# -----------------------------------------------------------------------------
PROJECT = Path('/data/guest_anirban/protein_design/long_chain_design')
REFERENCE_PDB = Path('/data/guest_anirban/structures/1YI5.pdb')
CAMPAIGN_DIRS = [
    PROJECT / 'AlphaFold_designs_121_A28_33',
    PROJECT / 'AlphaFold_designs_121_A30_34',
    PROJECT / 'AlphaFold_designs_125_A36_39',
    PROJECT / 'AlphaFold_designs_125_A53_56',
]
CUTOFF = 4.0  # Angstrom; heavy-atom contact threshold
RECEPTOR_CHAINS = list('ABCDE')
REFERENCE_TOXIN_CHAINS = list('FGHIJ')
BINDER_CHAIN, COMPLEX_TOXIN_CHAIN, REFERENCE_TOXIN_CHAIN = 'A', 'B', 'F'


def _load(path, sid='structure'):
    return PDBParser(QUIET=True).get_structure(sid, str(path))


def _heavy(atom):
    element = (getattr(atom, 'element', '') or '').strip().upper()
    return element != 'H' if element else not atom.get_name().upper().startswith('H')


def _aa_residues(chain):
    for residue in chain:
        if is_aa(residue, standard=False):
            atoms = [a for a in residue.get_atoms() if _heavy(a)]
            if atoms:
                yield residue, atoms


def _letter(residue):
    try:
        return seq1(residue.resname, custom_map={'UNK': 'X'})
    except Exception:
        return 'X'


def _sequence_and_residues(chain):
    residues = [r for r in chain if is_aa(r, standard=False) and _letter(r) != 'X']
    return ''.join(_letter(r) for r in residues), residues


# -----------------------------------------------------------------------------
# SECTION 1 — RECEPTOR-FOOTPRINT RECONSTRUCTION IN 1YI5
# Chains A-E are AChBP; chains F-J are the five toxin copies. Count heavy-atom
# pairs within 4.0 A using cKDTree. Output receptor_contacts_1YI5.csv.
# -----------------------------------------------------------------------------
def receptor_footprint(pdb=REFERENCE_PDB, output='receptor_contacts_1YI5.csv', cutoff=CUTOFF):
    model = _load(pdb, '1YI5')[0]
    receptor_xyz, receptor_meta = [], []
    for cid in RECEPTOR_CHAINS:
        if cid not in model:
            continue
        for residue, atoms in _aa_residues(model[cid]):
            for atom in atoms:
                receptor_xyz.append(atom.coord)
                receptor_meta.append((cid, residue))
    if not receptor_xyz:
        raise RuntimeError('No receptor heavy atoms found in chains A-E')
    tree = cKDTree(np.asarray(receptor_xyz, dtype=float))
    rows = []
    for tcid in REFERENCE_TOXIN_CHAINS:
        if tcid not in model:
            continue
        for tres, tatoms in _aa_residues(model[tcid]):
            n_pairs, r_hits, r_chains = 0, set(), set()
            for atom in tatoms:
                hits = tree.query_ball_point(atom.coord, cutoff)
                n_pairs += len(hits)
                for idx in hits:
                    rcid, rres = receptor_meta[idx]
                    r_hits.add((rcid, rres.id[1], rres.id[2], rres.resname))
                    r_chains.add(rcid)
            if n_pairs:
                details = ';'.join(f'{c}:{name}{n}{ic.strip()}' for c,n,ic,name in sorted(r_hits))
                rows.append({
                    'toxin_chain': tcid, 'toxin_residue_number': tres.id[1],
                    'insertion_code': tres.id[2].strip(), 'toxin_residue': _letter(tres),
                    'atom_pair_contacts': n_pairs, 'unique_receptor_residues': len(r_hits),
                    'receptor_chains': ','.join(sorted(r_chains)),
                    'receptor_residue_details': details,
                })
    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)
    print(f'Saved {len(df)} rows: {output}')
    return df


# -----------------------------------------------------------------------------
# SECTION 2 — CLUSTAL OMEGA MSA
# Used for multiple-sequence alignment of the toxin panels before conservation
# analysis. This calls the installed clustalo executable with FASTA output.
# -----------------------------------------------------------------------------
def run_clustalo(input_fasta, output_alignment, executable='clustalo'):
    cmd = [executable, '-i', str(input_fasta), '-o', str(output_alignment),
           '--outfmt=fasta', '--force']
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    print('Alignment saved:', output_alignment)


# -----------------------------------------------------------------------------
# SECTION 3 — PAIRWISE ALIGNMENT AND BINDER-FOOTPRINT OVERLAP
# Recorded PairwiseAligner settings: global; match +2; mismatch -1;
# gap-open -8; gap-extension -0.5. Each predicted complex: chain A binder,
# chain B toxin. Align complex chain B to 1YI5 chain F, then map binder-contact
# residues to reference residue numbering. Heavy-atom contact cutoff: 4 A.
# Outputs: epitope_overlap_summary.csv and binder_residue_contacts.csv.
# -----------------------------------------------------------------------------
def _aligner():
    a = Align.PairwiseAligner()
    a.mode = 'global'
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -8.0
    a.extend_gap_score = -0.5
    return a


def _query_to_reference_map(query, reference):
    # Biopython alignment target is reference; query is inherited complex toxin.
    aln = _aligner().align(reference, query)[0]
    mapping = {}
    ref_blocks, query_blocks = aln.aligned
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        for off in range(min(r1-r0, q1-q0)):
            mapping[q0+off] = r0+off
    return mapping, aln


def _toxin_contact_counts(binder_chain, toxin_chain, cutoff=CUTOFF):
    binder_xyz = [a.coord for _, atoms in _aa_residues(binder_chain) for a in atoms]
    if not binder_xyz:
        return {}
    tree = cKDTree(np.asarray(binder_xyz, dtype=float))
    counts = {}
    for residue, atoms in _aa_residues(toxin_chain):
        count = sum(len(tree.query_ball_point(atom.coord, cutoff)) for atom in atoms)
        if count:
            counts[(residue.id[1], residue.id[2].strip(), residue.resname)] = count
    return counts


def binder_footprint_overlap(campaign_dirs=CAMPAIGN_DIRS, reference_pdb=REFERENCE_PDB,
                             footprint_positions=None, summary_out='epitope_overlap_summary.csv',
                             contacts_out='binder_residue_contacts.csv', cutoff=CUTOFF):
    # IMPORTANT: supply the complete consensus footprint position set from
    # consensus_receptor_footprint.csv. The 19-position set from the recorded
    # analysis must be checked against that file; do not silently omit a site.
    if footprint_positions is None:
        raise ValueError('Pass the verified full consensus footprint residue-number set.')
    ref_chain = _load(reference_pdb, 'reference')[0][REFERENCE_TOXIN_CHAIN]
    ref_seq, ref_residues = _sequence_and_residues(ref_chain)
    ref_numbers = [r.id[1] for r in ref_residues]
    summary, contact_rows = [], []
    for folder in campaign_dirs:
        if not folder.exists():
            print('Skipping missing campaign folder:', folder)
            continue
        for pdb in sorted(folder.glob('*_af2pred.pdb')):
            model = _load(pdb, pdb.stem)[0]
            if BINDER_CHAIN not in model or COMPLEX_TOXIN_CHAIN not in model:
                print('Skipping missing chain A/B:', pdb)
                continue
            binder, toxin = model[BINDER_CHAIN], model[COMPLEX_TOXIN_CHAIN]
            query_seq, query_residues = _sequence_and_residues(toxin)
            qmap, alignment = _query_to_reference_map(query_seq, ref_seq)
            q_lookup = {(r.id[1], r.id[2].strip(), r.resname): i
                        for i, r in enumerate(query_residues)}
            counts = _toxin_contact_counts(binder, toxin, cutoff)
            mapped, overlap = set(), set()
            for key, ncontacts in counts.items():
                qi = q_lookup.get(key)
                if qi is None or qi not in qmap:
                    continue
                ri = qmap[qi]
                if ri >= len(ref_numbers):
                    continue
                refnum = ref_numbers[ri]
                mapped.add(refnum)
                is_fp = refnum in footprint_positions
                if is_fp:
                    overlap.add(refnum)
                contact_rows.append({
                    'candidate': pdb.stem, 'campaign': folder.name,
                    'complex_toxin_residue_number': key[0],
                    'complex_toxin_residue': _letter(query_residues[qi]),
                    'reference_residue_number': refnum,
                    'atom_pair_contacts': ncontacts,
                    'is_receptor_footprint': is_fp,
                })
            fp_contacts = sum(x['atom_pair_contacts'] for x in contact_rows
                              if x['candidate'] == pdb.stem and x['campaign'] == folder.name
                              and x['is_receptor_footprint'])
            summary.append({
                'candidate': pdb.stem, 'campaign': folder.name,
                'mapped_contact_residues': len(mapped),
                'footprint_residues_contacted': len(overlap),
                'footprint_overlap_fraction': len(overlap)/len(footprint_positions),
                'footprint_atom_pair_contacts': fp_contacts,
                'contacted_footprint_positions': ','.join(map(str, sorted(overlap))),
                'alignment_score': alignment.score,
            })
    sdf, cdf = pd.DataFrame(summary), pd.DataFrame(contact_rows)
    sdf.to_csv(summary_out, index=False)
    cdf.to_csv(contacts_out, index=False)
    print('Saved:', summary_out, contacts_out)
    return sdf, cdf


# -----------------------------------------------------------------------------
# SECTION 4 — LONG-CHAIN VS SHORT-CHAIN CONSERVATION PROFILE
# Compute conservation by alignment column for two MSA FASTA files already
# placed in the same reference register. Delta = long conservation - short
# conservation. Inputs used in project: six long-chain and five short-chain
# sequences. Output: long_vs_short_specificity_profile_final.csv.
# -----------------------------------------------------------------------------
def read_aligned_fasta(path):
    from Bio import SeqIO
    records = list(SeqIO.parse(str(path), 'fasta'))
    if not records:
        raise ValueError(f'No FASTA records in {path}')
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise ValueError(f'FASTA is not aligned: {path}')
    return {r.id: str(r.seq).upper() for r in records}


def _column_conservation(aligned):
    seqs = list(aligned.values())
    rows = []
    for col in range(len(seqs[0])):
        residues = [s[col] for s in seqs if s[col] not in '-.']
        if not residues:
            rows.append(('-', 0.0))
        else:
            aa, n = Counter(residues).most_common(1)[0]
            rows.append((aa, n/len(seqs)))
    return rows


def long_short_conservation(long_msa, short_msa, output='long_vs_short_specificity_profile_final.csv',
                            reference_positions=None):
    long_data, short_data = read_aligned_fasta(long_msa), read_aligned_fasta(short_msa)
    lc, sc = _column_conservation(long_data), _column_conservation(short_data)
    if len(lc) != len(sc):
        raise ValueError('Both alignments must be mapped into the same reference coordinate system.')
    rows = []
    for i, ((la, lv), (sa, sv)) in enumerate(zip(lc, sc)):
        pos = reference_positions[i] if reference_positions is not None else i+1
        rows.append({'reference_position': pos, 'long_consensus': la,
                     'long_conservation': lv, 'short_consensus': sa,
                     'short_conservation': sv, 'specificity_delta': lv-sv})
    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)
    print('Saved:', output)
    return df


# -----------------------------------------------------------------------------
# SECTION 5 — CONTACT-WEIGHTED SPECIFICITY
# Formula: S_b = sum_i(n_bi * Delta_i) / sum_i(n_bi)
# n_bi = binder atom-pair contact count at footprint position i;
# Delta_i = long-chain conservation minus short-chain conservation.
# Inputs: binder_residue_contacts.csv and long_vs_short_specificity_profile_final.csv
# Outputs: candidate_specificity_summary.csv and candidate_footprint_specificity_details.csv
# -----------------------------------------------------------------------------
def contact_weighted_specificity(contacts_csv='binder_residue_contacts.csv',
                                 profile_csv='long_vs_short_specificity_profile_final.csv',
                                 summary_out='candidate_specificity_summary.csv',
                                 details_out='candidate_footprint_specificity_details.csv'):
    contacts, profile = pd.read_csv(contacts_csv), pd.read_csv(profile_csv)
    needed_c = {'candidate', 'reference_residue_number', 'atom_pair_contacts', 'is_receptor_footprint'}
    needed_p = {'reference_position', 'specificity_delta'}
    if not needed_c.issubset(contacts.columns):
        raise ValueError(f'Contacts file missing columns: {needed_c-set(contacts.columns)}')
    if not needed_p.issubset(profile.columns):
        raise ValueError(f'Profile file missing columns: {needed_p-set(profile.columns)}')
    contacts = contacts[contacts.is_receptor_footprint.astype(bool)].copy()
    merged = contacts.merge(profile, left_on='reference_residue_number',
                            right_on='reference_position', how='left').dropna(subset=['specificity_delta'])
    merged['weighted_contact'] = merged.atom_pair_contacts * merged.specificity_delta
    merged.to_csv(details_out, index=False)
    rows = []
    for candidate, group in merged.groupby('candidate'):
        total = group.atom_pair_contacts.sum()
        rows.append({
            'candidate': candidate,
            'footprint_positions_contacted': group.reference_residue_number.nunique(),
            'positive_specificity_positions': group.loc[group.specificity_delta > 0, 'reference_residue_number'].nunique(),
            'strong_specificity_positions_delta_ge_0_5': group.loc[group.specificity_delta >= 0.5, 'reference_residue_number'].nunique(),
            'total_footprint_atom_pair_contacts': total,
            'mean_specificity_delta': group.specificity_delta.mean(),
            'contact_weighted_specificity': group.weighted_contact.sum()/total if total else np.nan,
        })
    result = pd.DataFrame(rows).sort_values('contact_weighted_specificity', ascending=False)
    result.to_csv(summary_out, index=False)
    print('Saved:', summary_out, details_out)
    return result, merged


# -----------------------------------------------------------------------------
# SECTION 6 — SASA / BSA (BIOPYTHON SHRAKE-RUPLEY)
# Recorded: chain A binder, chain B toxin; probe radius 1.4 A; 200 points.
# BSA_total = SASA_A + SASA_B - SASA_complex; BSA_per_partner = BSA_total/2.
# Calculated for three representative candidates, not all 16 candidates.
# -----------------------------------------------------------------------------
def calculate_bsa(pdb_path, probe_radius=1.4, n_points=200, binder_chain='A', toxin_chain='B'):
    parser = PDBParser(QUIET=True)
    complex_struct = parser.get_structure('complex', str(pdb_path))
    ShrakeRupley(probe_radius=probe_radius, n_points=n_points).compute(complex_struct, level='S')
    complex_sasa = float(complex_struct[0].sasa)
    isolated = parser.get_structure('isolated', str(pdb_path))[0]
    if binder_chain not in isolated or toxin_chain not in isolated:
        raise ValueError(f'{pdb_path}: expected chains {binder_chain} and {toxin_chain}')
    binder, toxin = isolated[binder_chain].copy(), isolated[toxin_chain].copy()
    ShrakeRupley(probe_radius=probe_radius, n_points=n_points).compute(binder, level='S')
    ShrakeRupley(probe_radius=probe_radius, n_points=n_points).compute(toxin, level='S')
    sasa_a, sasa_b = float(binder.sasa), float(toxin.sasa)
    bsa = sasa_a + sasa_b - complex_sasa
    return {'pdb': str(pdb_path), 'binder_sasa_A2': sasa_a, 'toxin_sasa_A2': sasa_b,
            'complex_sasa_A2': complex_sasa, 'bsa_total_A2': bsa,
            'bsa_per_partner_A2': bsa/2, 'probe_radius_A': probe_radius, 'n_points': n_points}


def representative_bsa(output='representative_bsa_results.csv'):
    pdbs = [
        PROJECT/'AlphaFold_designs_121_A28_33/rfdesigns_len_121_99_dldesign_2_af2pred.pdb',
        PROJECT/'AlphaFold_designs_121_A28_33/rfdesigns_len_121_47_dldesign_3_af2pred.pdb',
        PROJECT/'AlphaFold_designs_125_A53_56/rfdesigns_len_118_3_dldesign_0_af2pred.pdb',
    ]
    rows = [calculate_bsa(p) for p in pdbs if p.exists()]
    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)
    print('Saved:', output)
    return df


# -----------------------------------------------------------------------------
# SECTION 7 — BOLTZ-2 CROSS-REACTIVITY MANIFEST
# Recorded design: 16 binders x 19 long-chain toxins = 304; plus 16 controls
# against inherited target = 320 total; 10 batches of 32.
# The precise YAML schema and exact original Boltz CLI command were not
# available in the recovered history. This creates the auditable manifest;
# generate YAMLs using the schema documented by the installed Boltz release.
# -----------------------------------------------------------------------------
def create_boltz_manifest(binder_fasta, toxin_fasta, inherited_target_fasta,
                          manifest='boltz_batch_manifest.csv', batch_size=32):
    from Bio import SeqIO
    binders = list(SeqIO.parse(str(binder_fasta), 'fasta'))
    toxins = list(SeqIO.parse(str(toxin_fasta), 'fasta'))
    controls = list(SeqIO.parse(str(inherited_target_fasta), 'fasta'))
    if len(binders) != 16 or len(toxins) != 19 or len(controls) != 1:
        raise ValueError(f'Expected 16 binders, 19 toxins, and 1 inherited target; got {len(binders)}, {len(toxins)}, {len(controls)}')
    rows, job = [], 0
    for binder in binders:
        for toxin in toxins:
            job += 1
            rows.append({'job_id': f'job_{job:03d}', 'job_type': 'panel_cross_reactivity',
                         'binder_id': binder.id, 'toxin_id': toxin.id,
                         'binder_sequence': str(binder.seq), 'toxin_sequence': str(toxin.seq),
                         'batch': (job-1)//batch_size+1})
    for binder in binders:
        job += 1
        rows.append({'job_id': f'job_{job:03d}', 'job_type': 'inherited_target_control',
                     'binder_id': binder.id, 'toxin_id': controls[0].id,
                     'binder_sequence': str(binder.seq), 'toxin_sequence': str(controls[0].seq),
                     'batch': (job-1)//batch_size+1})
    df = pd.DataFrame(rows)
    df.to_csv(manifest, index=False)
    print(f'Saved {len(df)} jobs in {df.batch.nunique()} batches: {manifest}')
    return df


# -----------------------------------------------------------------------------
# SECTION 8 — BOLTZ-2 RESULT MATRIX AND SUMMARY
# Supply a consolidated CSV with binder_id, toxin_id, iptm and optionally
# job_type, iplddt, ipde. Output matrix and candidate/toxin summary tables.
# -----------------------------------------------------------------------------
def summarize_boltz_results(results_csv, matrix_out='cross_reactivity_iptm_matrix.csv',
                            candidate_out='candidate_crossreactivity_summary.csv',
                            toxin_out='toxin_crossreactivity_summary.csv',
                            controls_out='inherited_reference_scores.csv'):
    df = pd.read_csv(results_csv)
    required = {'binder_id', 'toxin_id', 'iptm'}
    if not required.issubset(df.columns):
        raise ValueError(f'Results CSV must contain {required}')
    if 'job_type' in df.columns:
        panel = df[df.job_type != 'inherited_target_control'].copy()
        controls = df[df.job_type == 'inherited_target_control'].copy()
    else:
        panel, controls = df.copy(), pd.DataFrame()
    matrix = panel.pivot_table(index='binder_id', columns='toxin_id', values='iptm', aggfunc='first')
    matrix.to_csv(matrix_out)
    cand = panel.groupby('binder_id').iptm.agg(mean_panel_iptm='mean', median_panel_iptm='median',
                                               min_panel_iptm='min', std_panel_iptm='std',
                                               n_predictions='count').reset_index().sort_values('mean_panel_iptm', ascending=False)
    tox = panel.groupby('toxin_id').iptm.agg(mean_candidate_iptm='mean', median_candidate_iptm='median',
                                             n_predictions='count').reset_index().sort_values('mean_candidate_iptm', ascending=False)
    cand.to_csv(candidate_out, index=False)
    tox.to_csv(toxin_out, index=False)
    if not controls.empty:
        controls[['binder_id', 'toxin_id', 'iptm']].to_csv(controls_out, index=False)
    print('Saved Boltz matrix and summary CSVs.')
    return matrix, cand, tox


# -----------------------------------------------------------------------------
# SECTION 9 — COMMAND LINE
# Examples:
#   python btp_major_analyses.py footprint
#   python btp_major_analyses.py overlap --footprint-positions 6,7,8,...
#   python btp_major_analyses.py specificity
#   python btp_major_analyses.py sasa
#   python btp_major_analyses.py clustalo input.fasta aligned.fasta
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='BTP major computational analyses')
    parser.add_argument('task', choices=['footprint', 'overlap', 'specificity', 'sasa', 'clustalo'])
    parser.add_argument('input', nargs='?')
    parser.add_argument('output', nargs='?')
    parser.add_argument('--footprint-positions', default=None,
                        help='Comma-separated verified 1YI5 consensus footprint residue numbers')
    args = parser.parse_args()
    if args.task == 'footprint':
        receptor_footprint()
    elif args.task == 'overlap':
        if not args.footprint_positions:
            parser.error('overlap requires --footprint-positions from the complete consensus footprint CSV')
        positions = {int(x) for x in args.footprint_positions.split(',')}
        binder_footprint_overlap(footprint_positions=positions)
    elif args.task == 'specificity':
        contact_weighted_specificity()
    elif args.task == 'sasa':
        representative_bsa()
    elif args.task == 'clustalo':
        if not args.input or not args.output:
            parser.error('clustalo requires input and output FASTA paths')
        run_clustalo(args.input, args.output)

if __name__ == '__main__':
    main()
