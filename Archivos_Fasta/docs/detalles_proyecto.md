# Reading Frames

Fecha: 02/06/2024

**Participantes**:
García González Carlos

## Descripción del Problema
El propósito de este script es leer una secuencia de nucleótidos desde un archivo FASTA y generar archivos de salida que contengan los codones correspondientes a los seis posibles marcos de lectura. Un marco de lectura es una manera específica de dividir la secuencia de nucleótidos en tripletes consecutivos, que luego se traducen en aminoácidos para formar proteínas. Hay tres marcos de lectura posibles en la dirección 5'->3' y tres marcos en la dirección 3'->5' (utilizando la cadena complementaria inversa).


## Especificación de Requisitos

Requisitos Funcionales
Lectura de archivos FASTA:

El script debe ser capaz de leer archivos de secuencias en formato FASTA.
Procesamiento de secuencias:

Para cada secuencia en el archivo FASTA, el script debe extraer y procesar los codones en los seis posibles marcos de lectura (tres en la dirección 5'->3' y tres en la dirección 3'->5').
Generación de archivos de salida:

El script debe crear archivos de salida en formato FASTA para cada marco de lectura de cada secuencia, con nombres específicos que indiquen el marco de lectura.
Manejo de argumentos de línea de comandos:

El script debe aceptar como argumento de línea de comandos la ruta al archivo FASTA de entrada.
Formato de salida:

Los archivos de salida deben estar en formato FASTA, con los codones separados por espacios y el encabezado indicando el marco de lectura.
Requisitos No Funcionales
Eficiencia:

El script debe ser eficiente en términos de tiempo de ejecución, capaz de manejar archivos FASTA de gran tamaño con múltiples secuencias.
Compatibilidad:

El script debe ser compatible con Python 3 y usar las bibliotecas Biopython y argparse.
Portabilidad:

El script debe ser ejecutable en diferentes sistemas operativos donde Python esté instalado, incluyendo Linux, macOS y Windows.
Mantenibilidad:

El código debe ser claro y bien documentado, con funciones específicas para tareas individuales para facilitar la comprensión y mantenimiento.
Robustez:

El script debe manejar adecuadamente posibles errores, como la falta de archivo de entrada o archivos con formato incorrecto.
Uso de memoria:

El script debe ser eficiente en el uso de memoria, especialmente al manejar archivos de secuencias grandes.


## Análisis y Diseño
Para resolver este problema, se utilizarán varias funciones incorporadas en Python. A continuación, se muestra un pseudocódigo simple para ilustrar la lógica básica del script:

'''
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def get_codons(sequence, frame):
    """Returns the codons for a given reading frame."""
    codons = []
    start = frame - 1
    for i in range(start, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

def write_codons_to_fasta(codons, frame, output_file):
    """Writes the codons to a FASTA formatted file."""
    with open(output_file, 'w') as file:
        file.write(f">Frame{frame}\n")
        file.write(' '.join(codons) + '\n')

def process_frames(sequence, output_prefix):
    """Processes the sequence for all 6 frames and writes them to separate FASTA files."""
    # Frames 1 to 3: Forward frames
    for frame in range(1, 4):
        codons = get_codons(sequence, frame)
        output_file = f"{output_prefix}_Frame{frame}.fa"
        write_codons_to_fasta(codons, frame, output_file)

    # Frames 4 to 6: Reverse complement frames
    reverse_complement = sequence.reverse_complement()
    for frame in range(1, 4):
        codons = get_codons(reverse_complement, frame)
        output_file = f"{output_prefix}_Frame{frame+3}.fa"
        write_codons_to_fasta(codons, frame+3, output_file)

def main():
    parser = argparse.ArgumentParser(description="Process a FASTA file and output codons for all 6 reading frames.")
    parser.add_argument("input_file", help="Input FASTA file")
    args = parser.parse_args()

    # Read sequences from FASTA file
    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_id = record.id
        sequence = record.seq
        process_frames(sequence, seq_id)

if __name__ == "__main__":
    main()

'''




