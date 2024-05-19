# Ancient_Fungi
Work performed for ancient Botrytis fungi analysis of two strains

## Production of whole-genome sequencing data

Genomic DNA was extracted from the colonies using the ????? protocol previously described in ref. ???? or the ???? Kit . The genomes were sequenced either at the ????? Platform of Institut ????, the Department ???? using the Illumina sequencing technology. Paired-end reads of ???? bp were obtained. Reads have been deposited at the NCBI Sequence Read Archive under BioProject ID PRJN######.

Each set of paired-end reads was mapped against the ??? reference genome ??? haplotype A or haplotype B53 downloaded from the ???? database (version ?? ??Date) using the Bowtie2 Alignment tool, version ???. SAMtools version ?? and Picard tools version ??? (http://broadinstitute.github.io/picard) were then used to filter, sort and convert SAM files. 

SNPs were called using Genome Analysis Toolkit version 3.1–157,58,59, according to the GATK Best Practices. SNPs and indels were filtered using these following parameters: VariantFiltration, QD < 2.0, LowQD, ReadPosRankSum < −8.0, LowRankSum, FS > 60.0, HightFS, MQRankSum < −12.5, MQRankSum, MQ < 40.0, LowMQ, HaplotypeScore > 13.0, HaploScore. Coverages were also calculated using the Genome Analysis Toolkit.


<details>
<summary>Raw data</summary>

```


```
 </details>




<details>
<summary>Creating the variant calling files</summary>



```
#!/usr/bin/


#MAPPING:
bowtie2-build B05_REF_normalized.fasta B05_REF_normalized.fasta;
#java -jar /root/Software/picard/build/libs/picard.jar  NormalizeFasta -I B05_REF.fasta -O B05_REF_normalized.fasta
#samtools faidx B05_REF_normalized.fasta
#java -jar /root/Software/picard/build/libs/picard.jar CreateSequenceDictionary -R /root/Desktop/Dagan_Fungi_reads/20210317_Re-Run/B05_REF_normalized.fasta -O /root/Desktop/Dagan_Fungi_reads/20210317_Re-Run/B05_REF_normalized.dict
bowtie2 -x B05_REF_normalized.fasta -U G2.4_S2_R1_001.fastq.gz.2.fq.gz,SAD_2.4_S8_R1_001.fastq.gz -S B05normalized_bowtie_vs_2.4.sam --no-unal -p 20;
bowtie2 -x B05_REF_normalized.fasta -U G2.3_S3_R1_001.fastq.gz.2.fq.gz,SAD_2.3_S9_R1_001.fastq.gz -S B05normalized_bowtie_vs_2.3.sam --no-unal -p 20;
bowtie2 -x B05_REF_normalized.fasta -U GBO5_S4_R1_001.fastq.gz.2.fq.gz,SAD_B05_S10_R1_001.fastq.gz -S B05normalized_bowtie_vs_GB05.sam --no-unal -p 20;
bowtie2 -x B05_REF_normalized.fasta -U G7B3_S5_R1_001.fastq.gz.2.fq.gz,SAD_7B3_S11_R1_001.fastq.gz -S B05normalized_bowtie_vs_7B.sam --no-unal -p 20;

#FORMAT SAM > BAM:
for i in B05normalized*.sam; do samtools view -S -b $i > $i.bam; samtools sort -@ 20 $i.bam -f $i.sorted.bam; samtools index $i.bam.sorted.bam; done
echo '\n';
echo '\n';
echo '\n';
#HaplotypeCaller(vcf):
for i in B05normalized*.sorted.bam;
do echo $i;
java -jar /root/Software/picard/build/libs/picard.jar ValidateSamFile -I $i -MODE SUMMARY;
java -jar /root/Software/picard/build/libs/picard.jar AddOrReplaceReadGroups -I $i -O $i.out.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM 20;
java -jar /root/Software/picard/build/libs/picard.jar ValidateSamFile -I $i.out.bam -MODE SUMMARY;
done
echo '\n';
echo '\n';
echo HaplotypeCaller;
echo '\n';
echo '\n';
for i in B05normalized*.out.bam; do samtools sort -@ 20 $i -f $i.sorted2.bam; samtools index $i.sorted2.bam;

/root/Software/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /root/Desktop/Dagan_Fungi_reads/20210317_Re-Run/B05_REF_normalized.fasta  -I $i.sorted2.bam -O /root/Desktop/Dagan_Fungi_reads/VCF_OUTPUT/$i.g.vcf.gz;
done

```
 </details>


<details>
<summary>VCF Correction</summary>

```
gatk --java-options "-Xmx4g" VariantFiltration --reference T4_REF.fa --variant 054_T4.g.vcf.gz --filter-expression "QD < 2.0" --filter-name "SNP_QD" --filter-expression "FS > 60.0" --filter-name "SNP_FS" --filter-expression "SOR > 4.0" --filter-name "SNP_SOR" --filter-expression "MQ < 40.0" --filter-name "SNP_MQ" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_MQRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum" --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz

gatk --java-options "-Xmx4g" SelectVariants --reference T4_REF.fa --variant 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz -exclude-filtered --exclude-non-variants --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
Using GATK jar /usr/local/hurcs/gatk4/4.2.5.0/share/gatk4-4.2.5.0-0/gatk-package-4.2.5.0-local.jar





```
 </details>




<details>
<summary>Results location of VCF identified on chromosomes </summary>




```
Using Rplot https://www.bioinformatics.com.cn/plot_basic_SNP_density_by_CMplot_107_en

```
 </details>




<details>
<summary>Results T4 (ALOC0100000) vs. strains 903053 and 903054 </summary>


Command:

```
bcftools stats                      054_T4.g.vcf.gz 053_T4.g.vcf.gz > joined_T4_2.3_vs_7B.stats.txt

plot-vcfstats                      joined_T4_2.3_vs_7B.stats.txt -p outdir_T4
```

Figure Total counts for indels and SNPs:


More comparisons can be found in the folder 



 </details>

