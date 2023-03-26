import numpy as np


def aa2codon(aa, codon_freq):
    aa_codons = codon_freq[aa]
    keys = list(codon_freq[aa].keys())
    values = list(codon_freq[aa].values())
    sorted_value_index = np.argsort(values)
    sorted_dict = {keys[i]: values[i] for i in sorted_value_index}
    codon_freq.update({aa:sorted_dict})
    return max(aa_codons, key=aa_codons.get)