# Affine Gapped Sequence Alignment Tool

A Python implementation of the **Smith-Waterman algorithm** with **Affine Gap Penalties**. This tool performs local nucleotide sequence alignment between two FASTA files, allowing for more realistic biological modeling by penalizing gap openings and extensions differently.

### Features
* **Local Alignment:** Finds the most similar region between two sequences.
* **Affine Gaps:** Supports user-defined scoring for matches, mismatches, gap openings, and gap extensions.
* **FASTA Support:** Directly reads and parses `.fna` or `.fasta` files.

---

### Setup & Usage

#### Prerequisites
* Python 3.x

#### Execution
Run the script from the terminal using the following positional arguments:

```bash
python local_affine_alignment.py <file1> <file2> <match> <mismatch> <gap_open> <gap_extend>

Example Commands
To test with the sequences in this repository, use the following examples (assuming typical scoring parameters like match=2, mismatch=-1, gap_open=-10, gap_extend=-1):

1. Human vs. Mouse (1st Pair):
python local_affine_alignment.py human.fna mouse.fna 2 -1 -10 -1

2. Test Sequence 1 vs. Test Sequence 2 (2nd Pair):
python local_affine_alignment.py seq1.fna seq2.fna 2 -1 -10 -1

Implementation DetailsScoring Logic & FormulaThe tool calculates the cost of a gap of length $L$ using the affine model:$$W = o + (L \cdot e)$$o: Gap Open Penaltye: Gap Extend PenaltyThe alignment score at each cell $(i, j)$ is derived from the maximum of three states:$M(i, j)$: Residues are aligned (Match/Mismatch).$I_x(i, j)$: Sequence 1 is aligned to a gap (Deletion).$I_y(i, j)$: Sequence 2 is aligned to a gap (Insertion).Troubleshooting & FAQNegative vs. Positive Inputs: Standard biological scoring uses negative values for penalties. Ensure you input -10 and -1 for penalties rather than 10 and 1 to avoid rewarding gaps.Local vs. Global Alignment: This tool is strictly for Local Alignment. If comparing with external tools like EBI LALIGN, ensure they are set to "Local" mode; otherwise, scores will differ.Memory Usage: The implementation uses $O(n \cdot m)$ space. It is optimized for gene-sized sequences (like the provided test files) rather than full-genome scales.TestingAlignment results can be verified against the EBI LALIGN Tool using the same input parameters.
