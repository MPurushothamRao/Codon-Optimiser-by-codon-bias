import codon_freq
import random
import cai2
import codon2aa


def refine(predicted, cai, reference):
    temp = cai
    codons = [predicted[i:i + 3] for i in range(0, len(predicted), 3)]
    copied_codons = codons.copy()
    i=0
    for codon in codons:
        amino = codon2aa.codon2aa(codon)
        for x in codon_freq.codon_freq[amino].values():
            if codon_freq.codon_freq[amino][codon] != max(codon_freq.codon_freq[amino].values()):
                shuff = random.choice([k for k, v in codon_freq.codon_freq[amino].items()])
                copied_codons[i] = shuff
                cai = (cai2.CAI("".join(copied_codons), reference=reference))
                if temp >= cai:
                    copied_codons[i] = codon
        i+=1
    return "".join(copied_codons)
