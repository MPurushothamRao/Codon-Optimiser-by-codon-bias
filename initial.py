import aa2codon
import refiner
import update_codon_freq
import cai2


def initiator(protein, reference):
    codon_freq = {
        'A': {'GCT': 0.27, 'GCC': 0.17, 'GCA': 0.24, 'GCG': 0.32},
        'C': {'TGT': 0.35, 'TGC': 0.65},
        'D': {'GAT': 0.45, 'GAC': 0.55},
        'E': {'GAA': 0.78, 'GAG': 0.22},
        'F': {'TTT': 0.29, 'TTC': 0.71},
        'G': {'GGT': 0.51, 'GGC': 0.44, 'GGA': 0.02, 'GGG': 0.03},
        'H': {'CAT': 0.28, 'CAC': 0.72},
        'I': {'ATT': 0.32, 'ATC': 0.67, 'ATA': 0.01},
        'K': {'AAA': 0.80, 'AAG': 0.20},
        'L': {'TTA': 0.03, 'TTG': 0.04, 'CTT': 0.04, 'CTC': 0.07, 'CTA': 0.01, 'CTG': 0.79},
        'M': {'ATG': 1.0},
        'N': {'AAT': 0.16, 'AAC': 0.84},
        'P': {'CCT': 0.09, 'CCC': 0.01, 'CCA': 0.14, 'CCG': 0.76},
        'Q': {'CAA': 0.17, 'CAG': 0.83},
        'R': {'CGT': 0.64, 'CGC': 0.35, 'CGA': 0.005, 'CGG': 0.002, 'AGA': 0.002, 'AGG': 0.001},
        'S': {'TCT': 0.31, 'TCC': 0.27, 'TCA': 0.04, 'TCG': 0.07, 'AGT': 0.05, 'AGC': 0.26},
        'T': {'ACT': 0.28, 'ACC': 0.57, 'ACA': 0.04, 'ACG': 0.11},
        'V': {'GTT': 0.26, 'GTC': 0.22, 'GTA': 0.15, 'GTG': 0.37},
        'W': {'TGG': 1.0},
        'Y': {'TAT': 0.35, 'TAC': 0.65},
        '*': {'TAA': 0.87, 'TAG': 0.01, 'TGA': 0.12}
    }

    optimised = []
    sett = {}
    for x in protein:
        if x in sett:
            sett[x] += 1
        else:
            sett[x] = 1
        optimised.append(aa2codon.aa2codon(x, codon_freq))
        reduction_factor = (protein.count(x) - sett[x]) / protein.count(x)
        update_codon_freq.update_codon_freq(x, codon_freq, reduction_factor)
    predicted = "".join(optimised)
    cai = (cai2.CAI(predicted, reference=reference))
    print(predicted)
    print(cai)
    predicted= refiner.refine(predicted, cai, reference)
    cai = (cai2.CAI(predicted, reference=reference))
    print(predicted)
    print(cai)
    return predicted, cai
