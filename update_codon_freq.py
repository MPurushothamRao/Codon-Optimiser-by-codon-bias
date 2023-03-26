def update_codon_freq(x, codon_freq, reduction_factor):
    max_codon = max(codon_freq[x], key=codon_freq[x].get)
    max_value = codon_freq[x][max_codon]
    for codon in codon_freq[x]:
        if codon == max_codon:
            if x != "M" and x != "W":
                codon_freq[x][codon] = max_value* reduction_factor
    total_freq = sum(codon_freq[x].values())
    codon_freq = {x: {codon: (freq / float(total_freq))} for codon, freq in codon_freq[x].items() if total_freq != 0}
    return codon_freq

