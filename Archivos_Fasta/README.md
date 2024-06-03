# Reading Frames

El propósito de este script es leer una secuencia de nucleótidos desde un archivo FASTA y generar archivos de salida que contengan los codones correspondientes a los seis posibles marcos de lectura. Un marco de lectura es una manera específica de dividir la secuencia de nucleótidos en tripletes consecutivos, que luego se traducen en aminoácidos para formar proteínas. Hay tres marcos de lectura posibles en la dirección 5'->3' y tres marcos en la dirección 3'->5' (utilizando la cadena complementaria inversa).


## Uso

El script se ejecuta en una terminal que pueda usar python y solo se le debe de dar la siguiente informacion: python /biopython-class/archivos_fasta/src/template_program

## Salida

Genera seis archivos diferentes los cuales continen la secuencia de codones la cual va acambiar dependiendo del
marco de lectura usado, son 6 archivos de salida, uno por cada marco de lectura existente.


## Pruebas
El script incluye un archivo conjunto de pruebas unitarias

## Metadatos y documentacion
Este README ofrece información de uso básico. Para obtener información más detallada sobre el diseño y la implementación del script.


## Codigo fuente
El código fuente está disponible en este repositorio. Se acoge con satisfacción cualquier contribución o sugerencia a través de solicitudes pull request.

## Terminos de uso
Este script está disponible bajo la licencia MIT license. Consulte el archivo LICENSE para obtener más detalles.

## Como citar
Si utiliza este script en su trabajo, por favor citelo.

## Contactenos
email: <carlosgg@lcg.unam.mx>