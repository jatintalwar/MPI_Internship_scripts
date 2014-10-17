## for STAR run :

for file in `ls -1 *_trimmed.fastq`
do
PAIR=$file
echo "$PAIR"
~/installed_applications/STAR_2.3.0e/STAR --genomeDir ~/RNAseq_analysis/Genomes_Celegans/Star_genome_C_Elegans --readFilesIn $PAIR --runThreadN 15 --outReadsUnmapped Fastx --sjdbGTFfile ~/RNAseq_analysis/Genomes_Celegans/Caenorhabditis_elegans.WBcel235.76.gtf --outFileNamePrefix $PAIR
done


### for sam to bam :

for file in `ls -1 *.sam`
do
PAIR=$file
echo "$PAIR"
samtools view -S -b -o $PAIR.bam $PAIR
done

#### for sorting the bam:

for file in `ls -1 *.bam`
do
PAIR=$file
echo "$PAIR"
samtools sort $PAIR $PAIR.sorted
done