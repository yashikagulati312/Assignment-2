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

Troubleshooting & FAQ:
*Negative vs. Positive Inputs: Standard biological scoring uses negative values for penalties. Ensure you input -10 and -1 for penalties rather than 10 and 1 to avoid rewarding gaps.
*Local vs. Global Alignment: This tool is strictly for Local Alignment. If comparing with external tools like EBI LALIGN, ensure they are set to "Local" mode; otherwise, scores will differ.
*Memory Usage: The implementation uses O(nm) space. It is optimized for gene-sized sequences (like the provided test files) rather than full-genome scales.
