# ngs-project
contains rare disease genomics project about Trio WES Identification of a De Novo Variant in a simulated Paediatric Rare Disease Patient Using Trio WES Analysis
# following block of commands was ran on the each chromosome (shown here by using the example of chromosome 22)
CHR_NUM=22
CHR_ACC=NC_000022.11
GENOME=/mnt/c/Users/hp/ngs_project/a/GCF_000001405.40_GRCh38.p14_genomic.fna

mkdir -p /mnt/d/rare_disease_genomics/chr${CHR_NUM} && cd /mnt/d/rare_disease_genomics/chr${CHR_NUM}

samtools faidx $GENOME $CHR_ACC > reference.fa
bwa index reference.fa
samtools faidx reference.fa

samtools view -b -f 2 "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam" 22 -o proband_subset.bam
samtools view -b -f 2 "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam" 22 -o father_subset.bam
samtools view -b -f 2 "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam" 22 -o mother_subset.bam

samtools sort -n proband_subset.bam -o proband_sorted.bam
samtools fastq -F 0x900 -1 proband_R1.fastq.gz -2 proband_R2.fastq.gz -s proband_s.fastq.gz proband_sorted.bam
samtools sort -n father_subset.bam -o father_sorted.bam
samtools fastq -F 0x900 -1 father_R1.fastq.gz -2 father_R2.fastq.gz -s father_s.fastq.gz father_sorted.bam
samtools sort -n mother_subset.bam -o mother_sorted.bam
samtools fastq -F 0x900 -1 mother_R1.fastq.gz -2 mother_R2.fastq.gz -s mother_s.fastq.gz mother_sorted.bam

java -jar /usr/share/java/trimmomatic-0.39.jar PE -phred33 proband_R1.fastq.gz proband_R2.fastq.gz trimmed_proband_R1.fastq.gz unpaired_proband_R1.fastq.gz trimmed_proband_R2.fastq.gz unpaired_proband_R2.fastq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
java -jar /usr/share/java/trimmomatic-0.39.jar PE -phred33 father_R1.fastq.gz father_R2.fastq.gz trimmed_father_R1.fastq.gz unpaired_father_R1.fastq.gz trimmed_father_R2.fastq.gz unpaired_father_R2.fastq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
java -jar /usr/share/java/trimmomatic-0.39.jar PE -phred33 mother_R1.fastq.gz mother_R2.fastq.gz trimmed_mother_R1.fastq.gz unpaired_mother_R1.fastq.gz trimmed_mother_R2.fastq.gz unpaired_mother_R2.fastq.gz ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

bwa mem reference.fa trimmed_proband_R1.fastq.gz trimmed_proband_R2.fastq.gz > proband.sam
bwa mem reference.fa trimmed_father_R1.fastq.gz trimmed_father_R2.fastq.gz > father.sam
bwa mem reference.fa trimmed_mother_R1.fastq.gz trimmed_mother_R2.fastq.gz > mother.sam

samtools sort proband.sam -o proband.bam && samtools index proband.bam
samtools sort father.sam -o father.bam && samtools index father.bam
samtools sort mother.sam -o mother.bam && samtools index mother.bam

gatk AddOrReplaceReadGroups -I proband.bam -O proband_rg.bam -RGID proband -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM proband
gatk AddOrReplaceReadGroups -I father.bam -O father_rg.bam -RGID father -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM father
gatk AddOrReplaceReadGroups -I mother.bam -O mother_rg.bam -RGID mother -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM mother
samtools index proband_rg.bam && samtools index father_rg.bam && samtools index mother_rg.bam

gatk CreateSequenceDictionary -R reference.fa
gatk HaplotypeCaller -R reference.fa -I proband_rg.bam -O proband.vcf
gatk HaplotypeCaller -R reference.fa -I father_rg.bam -O father.vcf
gatk HaplotypeCaller -R reference.fa -I mother_rg.bam -O mother.vcf

bgzip -k proband.vcf && bgzip -k father.vcf && bgzip -k mother.vcf
bcftools index proband.vcf.gz && bcftools index father.vcf.gz && bcftools index mother.vcf.gz

bcftools merge proband.vcf.gz father.vcf.gz mother.vcf.gz -o trio.vcf
bcftools view -i 'FORMAT/GT[0]="0/1" || FORMAT/GT[0]="1/1"' trio.vcf | bcftools view -i 'FORMAT/GT[1]="0/0" && FORMAT/GT[2]="0/0"' -o denovo_candidates.vcf
grep -v "^#" denovo_candidates.vcf | wc -l

# result

of all the chromosomes in which the commands were run
no denovo variant was found

# discussion
while the Giab data is known for showing denovo variants, the lack of variants found is the project is theorised because of not going through entire genome because of hardware constraints.
