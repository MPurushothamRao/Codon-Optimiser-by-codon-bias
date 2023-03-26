import argparse
from Bio import SeqIO
import json

parser = argparse.ArgumentParser(description='to create codon bias dictionary from reference nucleotide itself')
parser.add_argument('-ref', '--reference', type=str, help='reference file in nucleotide fasta')
parser.add_argument('-out', '--output', type=str, help='output file in dict.py')
args = parser.parse_args()
out = args.output
ref = args.reference

# Create a dictionary to store the codon frequencies
codon_prop = {}
codon_table = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

codon_freq = {aa: {} for aa in set(codon_table.values())}
# Loop over the sequence, counting the occurrence of each codon
reference = []
with open(ref, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        nucleotide = record.seq
        nucleotide = nucleotide.replace("\n", "")
        nucleotide = str(nucleotide)
        reference.append(nucleotide)
for seq in reference:
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in codon_prop:
            codon_prop[codon] += 1
        else:
            codon_prop[codon] = 1
# Calculate the total number of codons
total_codons = sum(codon_prop.values())

# Normalize the counts to obtain the frequencies
for codon in codon_prop:
    codon_prop[codon] /= total_codons
    codon_freq[codon_table[codon]][codon] = float("{:.3f}".format(codon_prop[codon]))

for aa in codon_freq:
    total_aa_freq = sum(codon_freq[aa].values())
    for codon in codon_freq[aa]:
        codon_freq[aa][codon] = float("{:.3f}".format(codon_freq[aa][codon] / total_aa_freq))
# Print the resulting dictionary and saving it in a file

with open(out, 'w') as output_file:
    output_file.write(json.dumps(codon_freq))

print(codon_freq)
