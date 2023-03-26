# Codon-Optimiser-by-codon-bias
It is simple codon optimiser by using codon bias and codon usage

Use requiremnts.txt file for required packages
To create codon bias dict:

frequencer.py [-h] [-ref REFERENCE] [-out OUTPUT]

to create codon bias dictionary from reference nucleotide itself

optional arguments:
  -h, --help            show this help message and exit
 
  -ref REFERENCE, --reference REFERENCE reference file in nucleotide fasta
  
  -out OUTPUT, --output OUTPUT output file in dict.py
  
  Output of this should be used in case for different expression system in initiator.py and whenever codon_freq is used.
  
  
Usage:
main.py [-h] [-in INPUTNUCL] [-ip INPUTPROT] [-on OUTPUTNUCL] [-op OUTPUTPROT] [-ref REFERENCE] [-repeat REPEAT]

to optimise codon either from protein or nucleotide itself

optional arguments:
  -h, --help            show this help message and exit
  
  -in INPUTNUCL, --inputnucl INPUTNUCL
                        input file in nucleotide fasta
                        
  -ip INPUTPROT, --inputprot INPUTPROT
                        input file in protein fasta
                        
  -on OUTPUTNUCL, --outputnucl OUTPUTNUCL
                        output file in nucleotide fasta
                        
  -op OUTPUTPROT, --outputprot OUTPUTPROT
                        output file in protein fasta
                        
  -ref REFERENCE, --reference REFERENCE
                        reference file in nucleotide fasta
                        
  -repeat REPEAT, --repeat REPEAT
                        to access repeat
                        
example: python main.py -in references.fasta -on ref_70.fasta -ref references.fasta -repeat True

Reference file is from http://genomes.urv.cat/HEG-DB/ which is highly expressed gene repository
