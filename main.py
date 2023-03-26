import argparse
import codon2aa
import initial
import generate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

parser = argparse.ArgumentParser(description='to optimise codon either from protein or nucleotide itself')
parser.add_argument('-in', '--inputnucl', type=str, help='input file in nucleotide fasta')
parser.add_argument('-ip', '--inputprot', type=str, help='input file in protein fasta')
parser.add_argument('-on', '--outputnucl', type=str, help='output file in nucleotide fasta')
parser.add_argument('-op', '--outputprot', type=str, help='output file in protein fasta')
parser.add_argument('-ref', '--reference', type=str, help='reference file in nucleotide fasta')
parser.add_argument('-repeat', '--repeat', type=bool, help='to access repeat', default=False)

args = parser.parse_args()

nucl = args.inputnucl
prot = args.inputprot
out = args.outputnucl
ref = args.reference
out_protein = args.outputprot
to_generate = args.repeat
proteins = []

if nucl is not None:
    with open(nucl) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            nucleotide = record.seq
            nucleotide = nucleotide.replace("\n", "")
            if len(nucleotide) % 3 != 0:
                print("Sequence should multiple's of 3")
                raise EOFError
            residue = [nucleotide[i:i + 3] for i in range(0, len(nucleotide), 3)]
            aa = []
            for x in residue:
                if x != '':
                    aa.append(codon2aa.codon2aa(x))
            protein = "".join(aa)
            proteins.append(protein)
else:
    with open(prot) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein = record.seq
            protein = protein.replace("\n", "")
            proteins.append(protein)
reference = []
with open(ref, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        nucleotide = record.seq
        nucleotide = nucleotide.replace("\n", "")
        reference.append(nucleotide)
for index in range(len(proteins)):
    protein = proteins[index]
    predicted, cai = initial.initiator(protein, reference)
    if to_generate:
        matrix = generate.generate_dot_matrix(predicted, predicted)
        df = pd.DataFrame(matrix)
        for x in df.columns:
            repe += df[x].value_counts()[1]
        repe = repe / (len(df.columns) * len(df.index))
        print(f"repeat {repe}")
    gc = generate.gc(predicted)
    with open(out, 'a') as output:
        record = SeqRecord(Seq(predicted), id=str(cai), description=str(gc))
        SeqIO.write(record, output, "fasta")
    if out_protein is not None:
        with open(out_protein, "a") as ort:
            record = SeqRecord(Seq(predicted), id=str(cai), description=str(len(predicted)))
            SeqIO.write(record, ort, "fasta")
    print("Random is not so random, CO-Incidence is purely hope of Incidence")
    index += 1
