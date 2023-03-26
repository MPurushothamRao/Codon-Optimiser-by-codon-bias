def generate_dot_matrix(seq1, seq2, window_size=3):
    matrix = []
    for i in range(len(seq1) - window_size + 1):
        row = []
        for j in range(len(seq2) - window_size + 1):
            if seq1[i:i + window_size] == seq2[j:j + window_size]:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)

    return matrix


def gc(seq):
    i = 0
    for x in seq:
        if x == "G" or x == "C":
            i += 1
    cc = i / len(seq)
    return cc
