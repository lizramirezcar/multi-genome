#!/bin/sh
#Este script se llama con una ruta a directorio con fastq, la ruta con el directorio para output, el accesion del cromosoma de referencia y opcional el accesion de los plasmidos de referencia 
#Para psudomolecula de solo cromosoma solo incluir el primer parametro

# Dependencias: 
#bwa -1.7.17
#samtools 1.7
#bcftools 1.7

# module load snp-site/2.4.1. 
# module load bwa/0.7.17
# module load bcftools/1.9
# module load raxml/8.2.12-SSE3-gcc-pthreads
# module load samtools/1.9

# CPU count
THREADS=$(grep -c ^processor /proc/cpuinfo)
#THREADS=4
MPILEUP_DEPTH=500
# En kabre
# raxml/8.2.12
#RAXML=raxmlHPraxml/8.2.12-SSE3-gcc-mpi
# raxml/git
#RAXML=raxmlHPC-SSE3
RAXML=raxmlHPC-PTHREADS-SSE3

SNPSITES=snp-sites
#SNPSITES=~/software/snp-sites/src/snp-sites

function say()   { echo "$1" ; };
function abort() { echo "$1" >&2 ; exit $2 ; };

# Verificar que directorio existe
# $1 es el primer parametro del script
if test -z "$1"; then
	abort "ERROR: Falta poner el directorio con archivos fastq como primer parametro"
fi

# Verificar que el primer argumento del script sea un directorio
SAMPLE_DIR="$1"
if test ! -d "$SAMPLE_DIR"; then
	abort "ERROR $1 no es un directorio"
fi

# Verificar que directorio existe
# $2 es el primer parametro del script
if test -z "$2"; then
	abort "ERROR: Falta poner el directorio de salida como segundo parametro"
fi

# Verificar que el segundo argumento del script sea un directorio
OUT_DIR="$2"
if test -d "$OUT_DIR"; then
	say "OUT_DIR es $OUT_DIR"
else
	say "Creando OUT_DIR: $OUT_DIR"
	mkdir -p $OUT_DIR 
fi
   
# Revisar que el tercer argumento se haya escrito, con test -z revisa si la longitud del texto escrito en la variable 3 es zero
if test -z "$3"; then
	abort "ERROR: Falta poner el accesion ID de la referencia como tercer parametro"
fi

CROMOSOMA=$3
REFERENCIA=$CROMOSOMA

#OUT_DIR=$SAMPLE_DIR/mapping_$REFERENCIA

mkdir -p $OUT_DIR

say "Script ejecutado: $0 $* "

# Descargar el fasta de la referencia si no existe
if test ! -s $OUT_DIR/$CROMOSOMA.fasta || head -1 $OUT_DIR/$CROMOSOMA.fasta | grep -v -q "^>" ; then
	say "Downloading chromosome $CROMOSOMA..."
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$CROMOSOMA&rettype=fasta&retmode=text" > $OUT_DIR/$CROMOSOMA.fasta 
fi
# Verificar que la descarga sea un fasta
if test ! -s $OUT_DIR/$CROMOSOMA.fasta || head -1 $OUT_DIR/$CROMOSOMA.fasta | grep -v -q "^>" ; then
	abort "ERROR al descargar la referencia"
fi

# Revisar que el cuarto argumento, si existe, se haya escrito
if test "$4"; then
	ACCESORIOS=$4
	# Descargar el fasta de la referencia si no existe, con test -f revisa que el file exista
	if test ! -s $OUT_DIR/$ACCESORIOS.fasta || head -1 $OUT_DIR/$ACCESORIOS.fasta | grep -v -q "^>" ; then
		say "Downloading $ACCESORIOS..."
		curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ACCESORIOS&rettype=fasta&retmode=text" > $OUT_DIR/$ACCESORIOS.fasta 
		# Verificar que la descarga sea un fasta
		if test ! -s $OUT_DIR/$ACCESORIOS.fasta || head -1 $OUT_DIR/$ACCESORIOS.fasta | grep -v -q "^>" ; then
			abort "ERROR al descargar los genes accesorios"
		fi
		#Concatenar cromosoma y accesorios
		cat $OUT_DIR/$ACCESORIOS.fasta >> $OUT_DIR/$CROMOSOMA.fasta
	fi
fi

# Indexar la referencia para bwa
say "Indexing reference file $OUT_DIR/$REFERENCIA.fasta"
bwa index $OUT_DIR/$REFERENCIA.fasta

rm -f $OUT_DIR/consol_aln_$REFERENCIA.fasta
touch $OUT_DIR/consol_aln_$REFERENCIA.fasta

# El patron del prefix se debe ajustar a los output files del secuenciador

# sanger files
#miseq files 
#for sample_1 in $SAMPLE_DIR/*_R1_001.fastq.gz; do
for sample_1 in $SAMPLE_DIR/*_1.fastq.gz $SAMPLE_DIR/*_R1_001.fastq.gz ; do
	test -f $sample_1 || continue

	#Get prefix del achivo. El patron 

	# sanger files
	prefix1=$( echo $sample_1 | sed 's/_1.fastq.gz//' );
	prefix1=$( echo $sample_1 | sed 's/_R1_001.fastq.gz//' );
	prefix=$(basename $prefix1)
	sample_2=$( echo $sample_1 | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/' );
	# sanger files
	sample_2=$( echo $sample_2 | sed 's/_1.fastq.gz/_2.fastq.gz/' );

	prefix_dir=$OUT_DIR/$prefix
	mkdir -p $prefix_dir

	#Mapping usando BWA, sort e index
	echo "[$prefix] Mapping sample $prefix with bwa against $REFERENCIA and $ACCESORIOS using $THREADS threads..."
	test -f $prefix_dir/"$prefix".sam || bwa mem -t $THREADS $OUT_DIR/$REFERENCIA.fasta "$sample_1" "$sample_2" > $prefix_dir/"$prefix".sam;
	echo "[$prefix] Using samtools view to format alignment file"
	test -f $prefix_dir/"$prefix".bam || samtools view -bS -q 15 $prefix_dir/"$prefix".sam > $prefix_dir/"$prefix".bam;
	echo "[$prefix] Using samtools sort..."
	test -f $prefix_dir/"$prefix"_sort.bam || samtools sort $prefix_dir/"$prefix".bam -o $prefix_dir/"$prefix"_sort.bam;
	echo "[$prefix] Using samtools index..."
	samtools index $prefix_dir/"$prefix"_sort.bam;

	#SNP calling original
#	echo "Calling SNPs"
#	samtools mpileup -DSugBf $SAMPLE_DIR/$REFERENCIA.fa $OUT_DIR/"$prefix"_sort.bam > $OUT_DIR/"$prefix"_tmp.bcf;
#	bcftools view -bcg $OUT_DIR/"$prefix"_tmp.bcf > $OUT_DIR/"$prefix".bcf;
#	bcftools index $OUT_DIR/"$prefix".bcf;
#	Build pseudo genome
#	echo "Building pseudo molecule..."

#	vcfutils.pl version modificada por sanger ftp://ftp.sanger.ac.uk/pub/resources/coursesandconferences/Mst5_2017/vcfutils.pl
#	bcftools view -cg $OUT_DIR/"$prefix".bcf "$CROMOSOMA" | vcfutils-fa.pl vcf2fa -d 5 -o $prefix - > $OUT_DIR/"$prefix".fa 	

#	SNP calling  con bcftools y Building pseudo genome
	echo "[$prefix] Calling SNPs on $prefix..."
	# --threads para call
	test -f $prefix_dir/"$prefix"_tmp.vcf.gz || bcftools mpileup --threads $THREADS --max-depth $MPILEUP_DEPTH -f $OUT_DIR/$REFERENCIA.fasta $prefix_dir/"$prefix"_sort.bam | bcftools call -mv -Oz -o $prefix_dir/"$prefix"_tmp.vcf.gz

	echo "[$prefix] Using bcftools index..."
	bcftools index --threads $THREADS $prefix_dir/"$prefix"_tmp.vcf.gz 
	echo "[$prefix] Using bcftools norm..."
	bcftools norm -f $OUT_DIR/$REFERENCIA.fasta $prefix_dir/"$prefix"_tmp.vcf.gz -Ob -o $prefix_dir/"$prefix"_norm.bcf
	echo "[$prefix] Using bcftools filter..."
	bcftools filter -Oz -s LowQual -e '%QUAL>=20 || DP>5' $prefix_dir/"$prefix"_norm.bcf > $prefix_dir/"$prefix"_flt.vcf.gz
	#bgzip -c $OUT_DIR/"$prefix"_norm_flt.vcf > $OUT_DIR/"$prefix"_flt.vcf.gz
	echo "[$prefix] Using bcftools index..."
	bcftools index $prefix_dir/"$prefix"_flt.vcf.gz
	echo "[$prefix] Normalizing, filtering and indexing variants on $prefix..."
	cat $OUT_DIR/$REFERENCIA.fasta | bcftools consensus -i 'type="snp"' $prefix_dir/"$prefix"_flt.vcf.gz | sed "s/$REFERENCIA/$prefix/" > $prefix_dir/"$prefix"-cons.fasta

	#generar consolidado de pseudomoleculas
	cat $prefix_dir/"$prefix"-cons.fasta >> $OUT_DIR/consol_aln_$REFERENCIA.fasta 
done
# correr snp-sites a consol_aln_$REFERENCIA.fa
echo "Using $SNPSITES ..."
$SNPSITES -o $OUT_DIR/consol_aln_snps_$REFERENCIA.fasta $OUT_DIR/consol_aln_$REFERENCIA.fasta

#numero de lineas
#cat $OUT_DIR/consol_aln_$REFERENCIA.fasta | while read line
#	do
#	count=$(echo $line | wc -c)
#	echo $count > $OUT_DIR/consol_aln_$REFERENCIA-counts.txt 
#	done

#correr RAXML
# $RAXML -T $THREADS -s $OUT_DIR/consol_aln_snps_$REFERENCIA.fasta -w $(readlink -f $OUT_DIR/ ) -n consol_aln_snps_"$REFERENCIA"_raxml -m GTRGAMMA -f d
raxml_dir=$(readlink -f $OUT_DIR/raxml_$REFERENCIA )
mkdir -p $raxml_dir
rm -f $raxml_dir/RAxML_*
$RAXML -T $THREADS -s $OUT_DIR/consol_aln_snps_$REFERENCIA.fasta -p 1723 -w $raxml_dir -n consol_aln_snps_"$REFERENCIA"_raxml -m GTRGAMMA -f d 


# Abrir con figtree
# java -jar "/windows/Users/lisbeth/Downloads/FigTree v1.4.3/FigTree v1.4.3/lib/figtree.jar"
