#!/usr/bin/env python3
"""Summarize DNA alignment site classes across a directory.

Gaps, missing characters, and IUPAC ambiguities are excluded from state counts.
Coordinates are 1-based in the original alignment.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path
from Bio import AlignIO

FORMATS = {'.fa':'fasta','.fas':'fasta','.fasta':'fasta','.fna':'fasta',
           '.phy':'phylip-relaxed','.phylip':'phylip-relaxed',
           '.nex':'nexus','.nexus':'nexus'}
BASES = frozenset('ACGT')
SUMMARY_GROUPS = [
    ('ORIGINAL MATRIX — composition and missing data', [
        'sequences', 'sites', 'sites_with_gaps', 'sites_with_ambiguities',
        'complete_acgt_sites', 'sites_with_missing', 'total_gap_cells',
        'total_ambiguous_cells', 'gap_cells_%', 'ambiguous_cells_%',
        'missing_data_%',
    ]),
    ('FILTERING — configuration and retained sites', [
        'deletion_mode', 'coverage_threshold_%', 'sites_retained',
        'sites_excluded', 'retained_sites_%',
    ]),
    ('RETAINED SITES — nucleotide variability and informativeness', [
        'constant', 'variable', 'variable_non_pi', 'singleton_sites',
        'parsimony_informative', 'all_missing', 'constant_%',
        'variable_%', 'parsimony_informative_%',
    ]),
]

DETAIL = ['file','position','classification','A','C','G','T','gaps',
          'ambiguities','observed_acgt','missing_total','distinct_states',
          'singleton_states','coverage_%','retained']


def analyze(path, detail_writer, deletion, coverage):
    aln = AlignIO.read(str(path), FORMATS[path.suffix.lower()])
    n, length = len(aln), aln.get_alignment_length()
    if n == 0 or length == 0:
        raise ValueError('Empty alignment')
    if any(len(record.seq) != length for record in aln):
        raise ValueError('Sequences have inconsistent lengths')
    seqs = [str(record.seq).upper().replace('U','T') for record in aln]
    counts = Counter()
    for pos, column in enumerate(zip(*seqs), 1):
        frequencies = Counter(c for c in column if c in BASES)
        gap = column.count('-')
        ambiguous = sum(c not in BASES and c != '-' for c in column)
        observed = sum(frequencies.values())
        distinct = len(frequencies)
        singletons = sum(v == 1 for v in frequencies.values())
        retained = (observed == n if deletion == 'complete' else
                    100 * observed >= coverage * n if deletion == 'partial' else True)
        counts['sites_retained'] += int(retained)
        if not observed:
            classification = 'all_missing'
        elif distinct == 1:
            classification = 'constant'
        elif sum(v >= 2 for v in frequencies.values()) >= 2:
            classification = 'parsimony_informative'
        else:
            classification = 'variable_non_pi'
        if retained:
            counts[classification] += 1
        counts['singleton_sites'] += int(retained and distinct > 1 and singletons > 0)
        counts['sites_with_gaps'] += int(gap > 0)
        counts['sites_with_ambiguities'] += int(ambiguous > 0)
        counts['complete_acgt_sites'] += int(observed == n)
        counts['sites_with_missing'] += int(observed < n)
        counts['total_gap_cells'] += gap
        counts['total_ambiguous_cells'] += ambiguous
        detail_writer.writerow({
            'file':path.name,'position':pos,'classification':classification,
            **{b:frequencies[b] for b in 'ACGT'},'gaps':gap,
            'ambiguities':ambiguous,'observed_acgt':observed,
            'missing_total':n-observed,'distinct_states':distinct,
            'singleton_states':singletons,'coverage_%':round(100*observed/n,2),
            'retained':int(retained)
        })
    variable = counts['parsimony_informative'] + counts['variable_non_pi']
    return {'file':path.name,'sequences':n,'sites':length,
            'constant':counts['constant'],'variable':variable,
            'variable_non_pi':counts['variable_non_pi'],
            'singleton_sites':counts['singleton_sites'],
            'parsimony_informative':counts['parsimony_informative'],
            'all_missing':counts['all_missing'],
            'constant_%':round(100*counts['constant']/counts['sites_retained'],2) if counts['sites_retained'] else 0,
            'variable_%':round(100*variable/counts['sites_retained'],2) if counts['sites_retained'] else 0,
            'parsimony_informative_%':round(100*counts['parsimony_informative']/counts['sites_retained'],2) if counts['sites_retained'] else 0,
            **{k:counts[k] for k in ('sites_with_gaps','sites_with_ambiguities',
                'complete_acgt_sites','sites_with_missing','total_gap_cells',
                'total_ambiguous_cells')},
            'gap_cells_%':round(100*counts['total_gap_cells']/(n*length),2),
            'ambiguous_cells_%':round(100*counts['total_ambiguous_cells']/(n*length),2),
            'missing_data_%':round(100*(counts['total_gap_cells']+counts['total_ambiguous_cells'])/(n*length),2),
            'deletion_mode':deletion, 'coverage_threshold_%':coverage if deletion == 'partial' else (100 if deletion == 'complete' else 0),
            'sites_retained':counts['sites_retained'],
            'sites_excluded':length-counts['sites_retained'],
            'retained_sites_%':round(100*counts['sites_retained']/length,2)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path, help='Folder containing alignments')
    p.add_argument('-o','--output',type=Path,default=None,help='Output folder (default: DIRECTORY/site_stats)')
    p.add_argument('-r','--recursive',action='store_true',help='Search subfolders recursively')
    p.add_argument('--deletion', choices=['none','complete','partial'], default='none',
                   help='Site filtering policy (default: none)')
    p.add_argument('--coverage', type=float, default=95,
                   help='Minimum percent A/C/G/T coverage for partial deletion (default: 95)')
    args = p.parse_args()
    if not 0 <= args.coverage <= 100:
        p.error('--coverage must be between 0 and 100')
    if args.deletion != 'partial' and args.coverage != 95:
        p.error('--coverage only applies with --deletion partial')
    if not args.directory.is_dir():
        p.error(f'Not a directory: {args.directory}')
    out = args.output or args.directory / 'site_stats'
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in (args.directory.rglob('*') if args.recursive else args.directory.iterdir())
                   if p.is_file() and p.suffix.lower() in FORMATS and out.resolve() not in p.resolve().parents)
    if not files:
        p.error('No supported alignment files found')
    results = []
    with (out/'sites.tsv').open('w', newline='', encoding='utf-8') as d, \
         (out/'errors.tsv').open('w', newline='', encoding='utf-8') as e:
        dw = csv.DictWriter(d, fieldnames=DETAIL, delimiter='\t')
        dw.writeheader()
        ew = csv.writer(e, delimiter='\t')
        ew.writerow(['file', 'error'])
        for file in files:
            try:
                row = analyze(file, dw, args.deletion, args.coverage)
                results.append(row)
                print(f"OK {file.name}: {row['sequences']} taxa, {row['sites']} sites, {row['parsimony_informative']} PI")
            except Exception as ex:
                ew.writerow([str(file), str(ex)])
                print(f'ERROR {file}: {ex}')

    # Transposed summary with clearly labeled groups; all values unchanged.
    with (out/'summary.tsv').open('w', newline='', encoding='utf-8') as s:
        sw = csv.writer(s, delimiter='\t')
        sw.writerow([''] + [row['file'] for row in results])
        for index, (section, fields) in enumerate(SUMMARY_GROUPS):
            if index:
                sw.writerow([])
            sw.writerow([section] + [''] * len(results))
            for statistic in fields:
                sw.writerow([statistic] + [row[statistic] for row in results])
    documentation = """ALIGNMENT SITE STATISTICS — METHODOLOGY AND DATA DICTIONARY

INPUT AND PROCESSING
- Processes each supported alignment in the input folder (optionally recursively).
- Supported formats: FASTA (.fa, .fas, .fasta, .fna), relaxed PHYLIP (.phy, .phylip), NEXUS (.nex, .nexus).
- Every alignment must contain sequences of equal length. Positions are numbered from 1 in the original alignment.
- A, C, G, T are the only unambiguous nucleotide states. U is normalized to T.
- '-' is a gap and is treated as missing data, NOT as a fifth state.
- '?' and all IUPAC ambiguity codes (including N, R, Y, S, W, K, M, B, D, H, V) are treated as missing data.
- Any other character besides A/C/G/T/- is also counted as ambiguous/missing; check the source alignment if unexpected symbols occur.
- Missing/ambiguous cells do not contribute to state frequencies or parsimony-informativeness.
- Site inclusion is controlled by --deletion: none (all sites), complete (100% A/C/G/T), partial (minimum --coverage percent A/C/G/T).
- Gaps and ambiguous characters never contribute nucleotide states, even in retained columns.
- sites.tsv preserves original coordinates and marks retained sites with 1 or 0.
- Raw missingness metrics are always calculated on the complete original matrix, regardless of filtering.
- Consequently, a site with only one observed A and all other cells missing is classified as constant, but has low coverage.

SITE CLASSIFICATION (mutually exclusive)
- all_missing: no unambiguous A/C/G/T in the column.
- constant: at least one A/C/G/T and exactly one distinct observed nucleotide.
- parsimony_informative: at least two distinct observed nucleotides, with at least two occurrences of EACH of at least two states (e.g. A,A,G,G,-).
- variable_non_pi: at least two distinct observed nucleotides but not parsimony-informative (e.g. A,A,A,G,-).
- variable = parsimony_informative + variable_non_pi.
- singleton_sites is NOT mutually exclusive with parsimony_informative: it counts variable columns containing any observed nucleotide state occurring exactly once (e.g. A,A,G,G,T).

SUMMARY.TSV: ONE MATRIX PER COLUMN, ORGANIZED IN THREE LABELED SECTIONS
- Section headings and blank separator rows are visual organization only; they are not statistics.
- ORIGINAL MATRIX contains raw alignment dimensions, missingness counts and cell percentages, independent of filtering.
- FILTERING contains the deletion mode, coverage threshold and retained/excluded site counts.
- RETAINED SITES contains variability and informativeness calculated only over retained columns.
- The first cell of the header row is empty; subsequent header cells are input filenames.

DATA DICTIONARY FOR SUMMARY.TSV
- sequences: number of sequence records (rows) in the alignment.
- sites: number of aligned positions (columns).
- constant: number of columns classified as constant among observed nucleotides.
- variable: number of columns with two or more distinct observed nucleotides.
- variable_non_pi: variable columns that are not parsimony-informative.
- singleton_sites: variable columns with at least one nucleotide state occurring once.
- parsimony_informative: columns with at least two observed states each occurring at least twice.
- all_missing: columns with no observed A/C/G/T.
- constant_%: 100 * constant / sites_retained (0 if none retained).
- variable_%: 100 * variable / sites_retained (0 if none retained).
- parsimony_informative_%: 100 * parsimony_informative / sites_retained (0 if none retained).
- sites_with_gaps: columns containing one or more '-' characters.
- sites_with_ambiguities: columns containing one or more ambiguous/missing non-gap characters.
- complete_acgt_sites: columns where EVERY sequence has an A/C/G/T character.
- sites_with_missing: columns with at least one gap or ambiguous/missing character; this is the UNION of sites_with_gaps and sites_with_ambiguities, not their sum.
- total_gap_cells: total number of '-' characters in the entire alignment.
- total_ambiguous_cells: total number of ambiguous/missing non-gap characters in the entire alignment.
- gap_cells_%: 100 * total_gap_cells / (sequences * sites).
- ambiguous_cells_%: 100 * total_ambiguous_cells / (sequences * sites).
- missing_data_%: 100 * (total_gap_cells + total_ambiguous_cells) / (sequences * sites).
- deletion_mode: none, complete or partial.
- coverage_threshold_%: minimum coverage used for filtering (0 for none, 100 for complete).
- sites_retained: number of columns satisfying the filter.
- sites_excluded: sites - sites_retained.
- retained_sites_%: 100 * sites_retained / sites.

PERCENTAGE DENOMINATORS
- Constant/variable/PI percentages use RETAINED columns; retained_sites_% uses all original columns.
- Cell percentages use ALL cells (number of sequences multiplied by alignment length).
- Values are rounded to two decimal places; totals may differ by 0.01 due to rounding.
- gap_cells_% + ambiguous_cells_% = missing_data_% before rounding.

SITES.TSV: ONE RECORD PER ORIGINAL ALIGNMENT POSITION
- file: input filename; position: original 1-based column number.
- classification: one of constant, variable_non_pi, parsimony_informative, all_missing.
- A, C, G, T: observed counts of each unambiguous nucleotide.
- gaps: number of '-' cells; ambiguities: number of non-gap missing/ambiguous cells.
- observed_acgt: sum of A+C+G+T; missing_total: gaps + ambiguities.
- distinct_states: number of distinct observed A/C/G/T states.
- singleton_states: number of observed nucleotide states that occur exactly once.
- coverage_%: 100 * observed_acgt / sequences.
- retained: 1 if included in classification totals, 0 otherwise.
- Classification of excluded sites is still shown for inspection, but excluded from summary counts.

ERRORS.TSV
- Lists input files that could not be read or analyzed. A failed file is omitted from summary.tsv.

INTERPRETATION CAUTIONS
- A 'constant' site may contain extensive missing data; it is not necessarily identical in all taxa.
- The number of parsimony-informative sites can change under different missing-data policies.
- IUPAC ambiguity codes are not resolved probabilistically or as sets of possible nucleotides.
- Gaps are not coded as indel characters, and adjacent gaps are not modeled as one evolutionary event.
- These are descriptive alignment statistics, not a phylogenetic signal or homoplasy test.
"""
    threshold = args.coverage if args.deletion == 'partial' else (100 if args.deletion == 'complete' else 0)
    if args.deletion == 'none':
        policy_sentence = (
            "All alignment columns were included irrespective of site coverage. "
            "Only unambiguous A, C, G and T observations contributed to site-state counts."
        )
        policy_label = "available-data analysis (no site deletion)"
    elif args.deletion == 'complete':
        policy_sentence = (
            "Complete deletion was applied: only columns with unambiguous A, C, G or T "
            "in every sequence (100% site coverage) were included in site statistics."
        )
        policy_label = "complete deletion (100% site coverage)"
    else:
        policy_sentence = (
            f"Partial deletion was applied with a minimum site coverage of {threshold:g}%: "
            f"only columns with at least {threshold:g}% unambiguous A, C, G or T "
            "observations were included in site statistics."
        )
        policy_label = f"partial deletion ({threshold:g}% minimum site coverage)"

    documentation += f"""
RUN SETTINGS
- deletion_mode: {args.deletion}
- coverage_threshold_%: {threshold:g}
- coverage parameter provided: {args.coverage:g}% (used only when deletion_mode=partial)
- selected policy: {policy_label}
- command-line arguments: --deletion {args.deletion}{f' --coverage {threshold:g}' if args.deletion == 'partial' else ''}

HOW TO CITE AND REPORT THIS ANALYSIS

A. READY-TO-ADAPT METHODS TEXT (ENGLISH; SETTINGS USED IN THIS RUN)
"Alignment statistics were calculated using the alignment_site_stats_citation.py
Python script (software version or archived release identifier: [INSERT]).
{policy_sentence} Gaps ('-'), missing characters ('?') and IUPAC ambiguity
codes (e.g. N, R, Y) were treated as missing data and did not contribute
nucleotide states. A site was classified as parsimony-informative when at
least two distinct unambiguous nucleotide states each occurred in at least
two sequences. Variable sites lacking this condition were classified as
variable but not parsimony-informative. Original-matrix gap and missing-data
percentages were calculated before site filtering; site-class percentages
were calculated relative to the number of retained sites. Alignments were
parsed with Biopython [CITE THE BIOPYTHON PAPER OR SOFTWARE RELEASE]."

B. SHORT METHODS TEXT (ENGLISH; SETTINGS USED IN THIS RUN)
"Constant, variable and parsimony-informative sites were counted using
alignment_site_stats_citation.py [CITE SCRIPT RELEASE] with {policy_label}.
Gaps and ambiguous nucleotides were treated as missing; missing-data
percentages refer to the original alignments."

C. SOFTWARE CITATION: WHAT TO REFERENCE
- Cite the actual archived version of this script if it has a persistent
  release (e.g. a Zenodo DOI); provide the exact version, author(s), year,
  title, repository/release URL and DOI if available.
- If no archived release exists, cite the script as custom, unpublished
  analysis code and provide a stable repository URL and commit hash or an
  archived copy as supplementary material. Do NOT invent a DOI, author,
  version number or peer-reviewed publication.
- Suggested reference template (fill in verified metadata):
  [AUTHOR(S)]. [YEAR]. alignment_site_stats_citation.py: DNA alignment site
  statistics with configurable missing-data filtering (version [VERSION],
  commit [COMMIT]). [REPOSITORY/ARCHIVE]. [DOI OR URL].
- Biopython is a dependency for alignment parsing, not the author of the
  custom site-classification logic. Cite Biopython separately if required
  by the journal, using verified bibliographic metadata.
- The MEGA-like terms 'complete deletion' and 'partial deletion' describe
  filtering policies; do not claim that MEGA was run or that results are
  identical to MEGA unless separately validated against that program.

D. REPORTING AND REPRODUCIBILITY CHECKLIST
- State the exact script version/commit, Python and Biopython versions,
  alignment filenames and whether subdirectories were searched.
- State the filtering mode explicitly; for partial deletion, report the
  numerical coverage threshold (e.g. --deletion partial --coverage 95).
- State that A/C/G/T alone count as observed nucleotides; U is mapped to T;
  gaps, '?' and all IUPAC ambiguities are missing.
- Report original alignment length and number of sequences, retained and
  excluded columns, and the denominator for each percentage.
- For full reproducibility, preserve input alignments, summary.tsv,
  sites.tsv, methodology.txt and the executed script/release.
- In none mode, a column with no A/C/G/T is retained and classified as
  all_missing; in complete/partial modes it is excluded (except that a
  partial threshold of 0% retains all columns).
- Note that singleton_sites can overlap parsimony_informative.
- These counts summarize site patterns; they do not establish phylogenetic
  resolution, homoplasy, or evolutionary rates.
"""

    (out/'methodology.txt').write_text(documentation, encoding='utf-8')
    print(f'Output: {out}')

if __name__ == '__main__':
    main()
