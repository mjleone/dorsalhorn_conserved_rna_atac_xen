#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --mem=24G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=chromap
#SBATCH --error=logs/chromap_%A_%a_out.txt
#SBATCH --output=logs/chromap_%A_%a_out.txt
#SBATCH --array 1-5
#SBATCH --no-requeue

ulimit -n 16000; source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq
CODEDIR=$PROJDIR/code/raw_code/preprocess_mouse_snATACseq
DATADIR=$PROJDIR/data/raw_data; TMPDIR=/scratch/bnphan
FASTQDIR=$DATADIR/fastq/Mouse_DH_snATAC-seq
BEDDIR=$DATADIR/bed
cd $PROJDIR/code/raw_code/preprocess_mouse_snATACseq

# sample-specific files
DIRLABEL=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 10) {print $1}' ${DATADIR}/tables/SampleSheet.csv )
SAMPLE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 10) {print $2}' ${DATADIR}/tables/SampleSheet.csv )

# mm10 reference genome indexed for chromap
BARCODE_LIST=/home/bnphan/src/cellranger-atac-1.2.0/cellranger-atac-cs/1.2.0/lib/python/barcodes/737K-cratac-v1_revcomp.txt
GENOME_FASTA=/home/bnphan/resources/genomes/mm10/mm10.fa
CHROMAP_IDX=/home/bnphan/resources/genomes/chromap/mm10/index
ALIGN_BED=${SAMPLE}.aln.bed

if [ ! -f ${BEDDIR}/$ALIGN_BED ]; then
	cd $TMPDIR; 
	## copy FASTQ to scratch dir
	rsync -Paq ${FASTQDIR}/${DIRLABEL}/${SAMPLE}*_R?_001.fastq.gz .

	## get the fastq files
	FQ1=$(ls ${SAMPLE}*_R1_001.fastq.gz | tr '\n' ',')
	BC1=$(ls ${SAMPLE}*_R2_001.fastq.gz | tr '\n' ',')
	FQ2=$(ls ${SAMPLE}*_R3_001.fastq.gz | tr '\n' ',')

	## align w/ chromap
	~/src/chromap-0.1_x64-linux/chromap --preset atac \
	-x $CHROMAP_IDX -r ${GENOME_FASTA} -t 8 \
	--remove-pcr-duplicates-at-cell-level \
	--bc-error-threshold 2 \
	--bc-probability-threshold .9 \
	-1 ${FQ1} -2 ${FQ2} -b ${BC1} -o ${ALIGN_BED} \
	--barcode-whitelist ${BARCODE_LIST}

	## move trimmed fastq to de-multiplex dir
	bgzip -c ${ALIGN_BED} > ${ALIGN_BED}.gz
	tabix -p bed ${ALIGN_BED}.gz
	rsync -Paq --remove-source-files ${ALIGN_BED}.gz* ${BEDDIR}
	rm -f ${SAMPLE}*_R?_001.fastq.gz ${ALIGN_BED} # clean up any old trim files

else echo 'chromap alignment bed file exists.'
fi



