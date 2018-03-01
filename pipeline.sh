#This is a simple RNA-Seq QC and analysis pipeline. Packages such as FastQC, are installed by conda.

#view dataset
head Read_1.fastq
head Read_2.fastq

#Assessment of data quality
fastqc --nogroup Read_1.fastq Read_2.fastq
firefox Read_*_fastqc/*.html &

#Pre-processing of data
fastq-mcf -o Read_1_q30l50.fastq -o Read_2_q30l50.fastq -q 30 -l 50 --qual-mean 30 adapters.fasta Read_1.fastq Read_2.fastq

#Reassessment of data quality
fastqc --nogroup Read_1_q30l50.fastq Read_2_q30l50.fastq firefox Read*q30l50*fastqc/*.html

#build up the reference using bowtie2
cd Reference
bowtie2-build mm10_chr19-1-20000000.fasta mm10_chr19-1-20000000


#alignment using tophat2
cd ..
tophat2 -o tophat2_with_filtered_data --no-mixed \
--rg-id Lane-1 --rg-sample sample1 --rg-center XYZ --rg-platform Illumina \
-g 1 \
-G Reference/mm10_chr19-1-20000000_Ensembl.gtf Reference/mm10_chr19-1-20000000 \
/home/training/Data/03_Quality_control_and_data_preprocessing/Read_1_q30l50.fastq \
/home/training/Data/03_Quality_control_and_data_preprocessing/Read_2_q30l50.fastq

#QA and QC of the sequence alignment
cat align_summary.txt

#view mapping quality
samtools view SRR769316.bam | less -S

#mark duplicates
java -jar /Users/zheliu/Desktop/ngsPackages/picard.jar MarkDuplicates I=SRR769316.bam O=SRR769316_dupMark.bam M=SRR769316_dupMarkMetrics.xls
open SRR769316_dupMarkMetrics.xls


#create folder test and put the results in it
mkdir test

#inspect the inner distance between two reads to reconfirm the library size selection process is correctly carried out.
inner_distance.py -i SRR769316.bam -r Reference/mm10_chr19-1-20000000_Emsembl.bed -o test/SRR76931
open test/SRR769316.inner_distance_plot.pdf

#inspect whether there is 3' or 5' end degradation, especially useful for single cell data
geneBody_coverage.py -i SRR769316.bam -r Reference/mm10_chr19-1-20000000_Emsembl.bed -o test/SRR769316
open test/SRR769316.geneBodyCoverage.curves.pdf

#inspect the splicing junction categories (known, partial novel, novel) 
junction_annotation.py -i SRR769316.bam -r Reference/mm10_chr19-1-20000000_Emsembl.bed -o test/SRR769316
open test/SRR769316.splice*.pdf

#inspect whether additional sequencing is needed to verify a splicing juntion
junction_saturation.py -i SRR769316.bam -r Reference/mm10_chr19-1-20000000_Emsembl.bed -o test/SRR769316 -v 5
open test/SRR769316.junctionSaturation_plot.pdf

#missing

#Visualising mapping reads
samtools index 
samtools index SRR769316_duplicates_marked.bam

igv.sh

#Load the reference:
#Genomes => Create .genome File
#Enter mm10-chr19-1-20000000 for unique Id and descriptive name
#Select the fasta file in Data => 06_Visualisation_of_mapped_reads => Reference => mm10-chr19-1-20000000.fasta Select the Gene file => Data => 06_#Visualisation_of_mapped_reads => Reference => mm10-chr19-1-20000000.gtf Save the .genome file

#Load both sample BAM files:
##File => Load from File
##Navigate to the file: Data => 06_Visualisation_of_mapped_reads => SRR769314_duplicates_marked.bam File => Load from File
##Navigate to the file: Data => 06_Visualisation_of_mapped_reads => SRR769316_duplicates_marked.bam

#Estimate gene read counts
#perform gene read estimate with parameter 'intersection-nonempty'
htseq-count -f bam -r pos -t exon -s no -m intersection-nonempty SRR769314_duplicates_marked.bam Reference/mm10_chr19-1-20000000_Emsembl.gtf > ./test/SRR769314_duplicates_marked.htseq.count
htseq-count -f bam -r pos -t exon -s no -m intersection-nonempty SRR769316_duplicates_marked.bam Reference/mm10_chr19-1-20000000_Emsembl.gtf > ./test/SRR769316_duplicates_marked.htseq.count

#perform gene read estimate with parameter 'union'
htseq-count -f bam -r pos -t exon -s no -m union SRR769314_duplicates_marked.bam Reference/mm10_chr19-1-20000000_Emsembl.gtf > ./test/SRR769314_duplicates_marked.htseq.union.count
htseq-count -f bam -r pos -t exon -s no -m union SRR769316_duplicates_marked.bam Reference/mm10_chr19-1-20000000_Emsembl.gtf > ./test/SRR769316_duplicates_marked.htseq.union.count


