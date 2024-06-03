'''
NAME
    Reading Frames

VERSION
    1

AUTHOR
    Carlos García González
        

DESCRIPTION
    Este script de python recibe como argumento el nombre de un archivo con el cual trabajar, una vez que sea así
    dibidira la secuencia por codones y por cada marco de lectura analizara esos codones, guardando la informacion
    de cada marco de lectura en un archivo FASTA, por lo tanto se tendran 6 archivos FASTA como salida.
        

CATEGORY
        

USAGE

    python template_program -i INPUT_FILE
    python template_program -i seq.nt.fa   
    

ARGUMENTS

    -i --INPUT_FILE

METHOD
    biopythonclass/archivos_FASTA

SEE ALSO


        
'''


# ===========================================================================
# =                            imports
# ===========================================================================

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# ===========================================================================
# =                            functions
# ===========================================================================

def get_codons(sequence, frame):
    """Regresa los codones por cada marco de lectura."""
    codons = []
    start = frame - 1
    for i in range(start, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

def write_codons_to_fasta(codons, frame, output_file):
    """Escribe los codones en un archivo FASTA."""
    with open(output_file, 'w') as file:
        file.write(f">Frame{frame}\n")
        file.write(' '.join(codons) + '\n')

def process_frames(sequence, output_prefix):
    """Procesa la secuencia y regresa los codones en un archivo ppor marco de lectura"""
    # Marcos de lectura del 1 al 3
    for frame in range(1, 4):
        codons = get_codons(sequence, frame)
        output_file = f"{output_prefix}_Frame{frame}.fa"
        write_codons_to_fasta(codons, frame, output_file)

    # Marcos de lectura del 4 al 6
    reverse_complement = sequence.reverse_complement()
    for frame in range(1, 4):
        codons = get_codons(reverse_complement, frame)
        output_file = f"{output_prefix}_Frame{frame+3}.fa"
        write_codons_to_fasta(codons, frame+3, output_file)



# ===========================================================================
# =                            main
# ===========================================================================


def main():
    parser = argparse.ArgumentParser(description="Procesa un archivo FASTA y como salida genera archivos FASTA usando los 6 marcos de lectura.")
    
    parser.add_argument("-i", "--input_file",
                        help = "El nombre del archivo FASTA a evaluar",
                        required = True)
    
    args = parser.parse_args()

    # Lectura de una secuencia de un archivo FASTA
    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_id = record.id
        sequence = record.seq
        process_frames(sequence, seq_id)

if __name__ == "__main__":
    main()
