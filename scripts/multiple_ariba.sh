#!/bin/sh
# Este script se llama con un directorio con fastq files, una base de datos y un directorio de salida
# Dependencias: ariba
#  https://github.com/sanger-pathogens/ariba
#  ariba debe estar en PATH
# module load intelpython/3.6
# spades.py 3.13
#  http://cab.spbu.ru/software/spades/
# module load SPAdes/3.13.0
# cd-hit
#  module load cdhit
#  # compilado en ~/bin/cd-hit-est -> ~/software/cdhit/cd-hit-est
#  module load gcc/7.2.0
# nucmer (en mummer)
# module load mummer
# module load bowtie2/2.2.9

# Verificar que directorio existe
# $1 es el primer parametro del script


THREADS=$(grep -c ^processor /proc/cpuinfo)
# Tiempo inicial
T0=`date +%s`

if test -z "$1"; then
	echo "ERROR. Este script necesita 3 parametros:"
	echo " 1: Directorio con archivos fastq"
	echo " 2: Una base de datos: argannot, card, ncbi, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core, vfdb_full, virulencefinder "
	echo " 3: Directorio de salida"
	exit
fi

# Verificar que el primer argumento del script sea un directorio
SAMPLE_DIR="$1"
if test ! -d "$SAMPLE_DIR"; then
	echo "ERROR $1 no es un directorio"
	exit
fi

# Revisar que el segundo argumento sea un nombre de base de datos
if test -z "$2"; then
	echo "ERROR: Falta poner la base de datos como segundo parametro"
	exit
# Verificar si la base de datos es una de las de ariba
elif echo $2 | grep -v -q -E 'argannot|card|ncbi|megares|plasmidfinder|resfinder|srst2_argannot|vfdb_core|vfdb_full|virulencefinder'; then
	echo La base de datos $2 no es una de la lista: argannot, card, ncbi, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core, vfdb_full, virulencefinder
fi
DATABASE=$2


# Validar que el directorio tenga archivos *_1.fastq.gz o *_1.fastq
for sample in $SAMPLE_DIR/*_1.fastq* ; do
	test -f "$sample" && _directorio_valido=1 ; break
done
if test -z $_directorio_valido ; then
	echo "El directorio $SAMPLE_DIR no tiene archivos fastq o fastq.gz"
	exit
fi

if test -z "$3"; then
	echo "ERROR: Falta poner un directorio de salida como tercer parametro"
	exit
fi


# Crear dir de outputs
OUT_DIR=$3/ariba_out_"$DATABASE"
mkdir -p $OUT_DIR/
mkdir -p $OUT_DIR/samples

# ARIBA: Descargar la base de datos e info de version en un directorio de referencia
if test ! -d $OUT_DIR/downloaded.$DATABASE; then
	echo "[$((`date +%s`-T0))s] Descargando base de datos $DATABASE con ariba..."
	ariba getref $DATABASE $OUT_DIR/downloaded.$DATABASE || exit
fi

# ARIBA: Preparar la base de datos
if test ! -d $OUT_DIR/$DATABASE.prepareref; then
	echo "[$((`date +%s`-T0))s] Preparando referencia con ariba..."
	ariba prepareref -f $OUT_DIR/downloaded.$DATABASE.fa -m $OUT_DIR/downloaded.$DATABASE.tsv $OUT_DIR/$DATABASE.prepareref
fi


# Ciclo: por cada muestra
#   ARIBA: run
for sample in $SAMPLE_DIR/*_1.fastq*; do
	sample_name=$( basename $sample | sed 's/_1.fastq.gz//' | sed 's/_1.fastq//' )
	if test -d $OUT_DIR/samples/$sample_name/ ; then
		echo "[$((`date +%s`-T0))s] $sample_name ya esta procesada"
	else
		echo "[$((`date +%s`-T0))s] Procesando sample $sample_name con ariba ..."
		ariba run --threads $THREADS $OUT_DIR/$DATABASE.prepareref $sample $SAMPLE_DIR/"$sample_name"_2.fastq* $OUT_DIR/samples/$sample_name/
	fi
done

# ARIBA: Summary
ariba summary $OUT_DIR/ariba_summary $OUT_DIR/samples/*/report.tsv
echo "Ariba summary en $OUT_DIR/ariba_summary.csv"
echo "Tiempo total: $((`date +%s`-T0)) segundos"
