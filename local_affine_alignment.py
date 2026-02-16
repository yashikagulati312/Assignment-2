import sys

NEG_INF = -10**9


def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()


def local_affine_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):

    n = len(seq1)
    m = len(seq2)

    # DP matrices
    M = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # Traceback matrices
    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Ix = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Iy = [[None] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0, "M")

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            # M matrix (local alignment: include 0)
            candidates = [
                (0, None),
                (M[i - 1][j - 1] + score, "M"),
                (Ix[i - 1][j - 1] + score, "Ix"),
                (Iy[i - 1][j - 1] + score, "Iy")
            ]

            M[i][j], prev_matrix = max(candidates, key=lambda x: x[0])

            if prev_matrix:
                trace_M[i][j] = (prev_matrix, i - 1, j - 1)

            # Ix matrix (gap in seq2)
            open_gap = M[i - 1][j] + gap_open
            extend_gap = Ix[i - 1][j] + gap_extend
            Ix[i][j] = max(open_gap, extend_gap)

            if Ix[i][j] == open_gap:
                trace_Ix[i][j] = ("M", i - 1, j)
            else:
                trace_Ix[i][j] = ("Ix", i - 1, j)

            # Iy matrix (gap in seq1)
            open_gap = M[i][j - 1] + gap_open
            extend_gap = Iy[i][j - 1] + gap_extend
            Iy[i][j] = max(open_gap, extend_gap)

            if Iy[i][j] == open_gap:
                trace_Iy[i][j] = ("M", i, j - 1)
            else:
                trace_Iy[i][j] = ("Iy", i, j - 1)

            # Track global maximum
            for matrix_name, matrix_val in [("M", M[i][j]),
                                            ("Ix", Ix[i][j]),
                                            ("Iy", Iy[i][j])]:
                if matrix_val > max_score:
                    max_score = matrix_val
                    max_pos = (i, j, matrix_name)

    # Traceback
    align1 = []
    align2 = []
    i, j, matrix = max_pos

    while True:

        if matrix == "M":
            if M[i][j] == 0:
                break
            prev = trace_M[i][j]
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])

        elif matrix == "Ix":
            prev = trace_Ix[i][j]
            align1.append(seq1[i - 1])
            align2.append("-")

        else:  # Iy
            prev = trace_Iy[i][j]
            align1.append("-")
            align2.append(seq2[j - 1])

        matrix, i, j = prev

    align1.reverse()
    align2.reverse()

    return max_score, "".join(align1), "".join(align2)


if __name__ == "__main__":

    if len(sys.argv) != 7:
        print("Usage:")
        print("python local_affine_alignment.py seq1.fasta seq2.fasta match mismatch gap_open gap_extend")
        sys.exit(1)

    seq1 = read_fasta(sys.argv[1])
    seq2 = read_fasta(sys.argv[2])

    match = int(sys.argv[3])
    mismatch = int(sys.argv[4])
    gap_open = int(sys.argv[5])
    gap_extend = int(sys.argv[6])

    score, alignment1, alignment2 = local_affine_alignment(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    print("\n===== BEST LOCAL ALIGNMENT (Affine Gap) =====\n")
    print(alignment1)
    print(alignment2)
    print("\nAlignment Score:", score)
