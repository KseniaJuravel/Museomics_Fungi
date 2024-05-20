# Ancient_Fungi
Work performed for ancient Botrytis fungi analysis of two strains

## Production of whole-genome sequencing data

Genomic DNA was extracted from the colonies using the ????? protocol previously described in ref. ???? or the ???? Kit . The genomes were sequenced either at the ????? Platform of Institut ????, the Department ???? using the Illumina sequencing technology. Paired-end reads of ???? bp were obtained. Reads have been deposited at the NCBI Sequence Read Archive under BioProject ID PRJN######.

Each set of paired-end reads was mapped against the ??? reference genome ??? haplotype ???  downloaded from the ???? database (version ?? ??Date) using the Bowtie2 Alignment tool, version 2.4.1 (Built on Fri Feb 28 17:23:43 UTC 2020). SAMtools samtools 1.9, Using htslib 1.9 and Picard tools version 1.70 (http://broadinstitute.github.io/picard) to filter, sort and convert SAM files. 

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

SNPs were called using Genome Analysis Toolkit version 4.2.0.0 and 4.2.5.0 according to the GATK Best Practices. SNPs and indels were filtered using the following parameters: VariantFiltration, QD < 2.0, LowQD, ReadPosRankSum < −8.0, LowRankSum, FS > 60.0, HightFS, MQRankSum < −12.5, MQRankSum, MQ < 40.0, LowMQ, HaplotypeScore > 13.0, HaploScore. Coverages were also calculated using the Genome Analysis Toolkit.

The parameters are shown in other fungi recent research:

https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000979#cited

Li X, Muñoz JF, Gade L, Argimon S, Bougnoux M-E, Bowers JR, Chow NA, Cuesta I, Farrer RA, Maufrais C, et al. 2023. Comparing genomic variant identification protocols for Candida auris. Microbial Genomics 9:1–19.

https://www.nature.com/articles/s41467-018-04787-4#citeas

Ropars, J., Maufrais, C., Diogo, D. et al. Gene flow contributes to diversification of the major fungal pathogen Candida albicans. Nat Commun 9, 2253 (2018). https://doi.org/10.1038/s41467-018-04787-4


<details>
<summary>VCF Correction command:</summary>

```
gatk --java-options "-Xmx4g" VariantFiltration --reference T4_REF.fa --variant 054_T4.g.vcf.gz --filter-expression "QD < 2.0" --filter-name "SNP_QD" --filter-expression "FS > 60.0" --filter-name "SNP_FS" --filter-expression "SOR > 4.0" --filter-name "SNP_SOR" --filter-expression "MQ < 40.0" --filter-name "SNP_MQ" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_MQRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum" --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz

gatk --java-options "-Xmx4g" SelectVariants --reference T4_REF.fa --variant 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz -exclude-filtered --exclude-non-variants --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
Using GATK jar /usr/local/hurcs/gatk4/4.2.5.0/share/gatk4-4.2.5.0-0/gatk-package-4.2.5.0-local.jar

```
 </details>


To ensure that the correction impacts the results significantly, the location of the variants was plotted to their location and manually inspected for variant T4 (ALOC0100000) vs. strains 903053. 

<details>
<summary>Results location of VCF identified on chromosomes comparison for correction</summary>


```
Using Rplot https://www.bioinformatics.com.cn/plot_basic_SNP_density_by_CMplot_107_en

```

Plot of locations without correction:



Plot of locations with correction:

T4_vs_903054
![T4_vs_903054](https://github.com/KseniaJuravel/Ancient_Fungi/blob/main/VCF_output/T4_vs_903054/13f5f469b08a2ec5.png)


 </details>




<details>
<summary>Results T4 (ALOC0100000) vs. strains 903053 and 903054 </summary>

<details>
<summary>Tool #1 for VCF analysis - bcftools</summary>


Command:

```
bcftools stats                      054_T4.g.vcf.gz 053_T4.g.vcf.gz > joined_T4_2.3_vs_7B.stats.txt

plot-vcfstats                      joined_T4_2.3_vs_7B.stats.txt -p outdir_T4
```

Figure Total counts for indels and SNPs:


More comparisons can be found in the folder. 


 </details>

 
<details>
<summary>Tool #2 for VCF analysis - vt peek</summary>


```
vt/vt peek 053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz
peek v0.5

options:     input VCF file            053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz


stats: no. of samples                     :          1
       no. of chromosomes                 :        118

       ========== Micro variants ==========

       no. of SNP                         :      67489
           2 alleles                      :           67484 (2.97) [50470/17014]
           3 alleles                      :               5 (0.67) [4/6]

       no. of INDEL                       :       3328
           2 alleles                      :            3312 (0.90) [1568/1744]
           3 alleles                      :              16 (0.23) [6/26]

       no. of SNP/INDEL                   :          5
           3 alleles                      :               5 (1.50) [3/2] (inf) [5/0]

       no. of micro variants              :      70822

       ++++++ Other useful categories +++++

        no. of complex substitutions      :          5
           3 alleles                      :               5 (1.50) [3/2] (inf) [5/0]


       ========= General summary ==========

       no. of VCF records                        :      70822
```

```
vt/vt peek 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
peek v0.5

options:     input VCF file            054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz


stats: no. of samples                     :          1
       no. of chromosomes                 :        118

       ========== Micro variants ==========

       no. of SNP                         :      89430
           2 alleles                      :           89424 (2.92) [66608/22816]
           3 alleles                      :               6 (0.33) [3/9]

       no. of INDEL                       :       4749
           2 alleles                      :            4724 (0.86) [2190/2534]
           3 alleles                      :              25 (0.79) [22/28]

       no. of SNP/INDEL                   :          1
           3 alleles                      :               1 (0.00) [0/1] (inf) [1/0]

       no. of micro variants              :      94180

       ++++++ Other useful categories +++++

        no. of complex substitutions      :          1
           3 alleles                      :               1 (0.00) [0/1] (inf) [1/0]


       ========= General summary ==========

       no. of VCF records                        :      94180
```

The following results represent the T4 (ALOC0100000) genome vs. the local B05 sequenced (1. no correction of parameters for the VCF output, 2. with correction):


```
(Not corrected HapplotypeCaller)
vt/vt peek B05_T4.g.vcf.gz
peek v0.5

options:     input VCF file            B05_T4.g.vcf.gz


stats: no. of samples                     :          1
       no. of chromosomes                 :        118

       ========== Micro variants ==========

       no. of SNP                         :     147695
           2 alleles                      :          147647 (2.66) [107267/40380]
           3 alleles                      :              48 (0.48) [31/65]

       no. of INDEL                       :      15900
           2 alleles                      :           15766 (0.90) [7475/8291]
           3 alleles                      :             134 (1.76) [171/97]

       no. of SNP/INDEL                   :         25
           3 alleles                      :              25 (0.92) [12/13] (1.00) [16/16]

       no. of micro variants              :     163620

       ++++++ Other useful categories +++++

        no. of complex substitutions      :         25
           3 alleles                      :              25 (0.92) [12/13] (1.00) [16/16]


       ========= General summary ==========

       no. of VCF records                        :     163620
```


```
(Corrected HapplotypeCaller)
vt/vt peek B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
peek v0.5

options:     input VCF file            B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz


stats: no. of samples                     :          1
       no. of chromosomes                 :        117

       ========== Micro variants ==========

       no. of SNP                         :      83016
           2 alleles                      :           83010 (2.92) [61844/21166]
           3 alleles                      :               6 (0.50) [4/8]

       no. of INDEL                       :       4046
           2 alleles                      :            4027 (0.88) [1886/2141]
           3 alleles                      :              19 (0.90) [18/20]

       no. of micro variants              :      87062

       ++++++ Other useful categories +++++


       ========= General summary ==========

       no. of VCF records                        :      87062
```
 </details>


<details>
<summary>Results</summary>

```


```
 </details>
 </details>

 The results shown by the tools indicate


<details>
<summary>Results B05.10 genome (AAID02000000) vs. strains 903053 and 903054 </summary>

<details>
<summary>Tool #1 for VCF analysis - bcftools</summary>


Command:

```


```
 </details>


<details>
<summary>Tool #2 for VCF analysis - vt peek</summary>


Command:

```


```
 </details>

The results shown by the tools indicate

 
 </details>



 ## Phylogenetic analyses

Phylogenetic analysis of multiple internal fragments of house-keeping genes (n=10) of the two
museum genomes with contemporary Botrytis genus genomes using the Multilocus sequence
typing (MLST) (Plesken et al., 2021) combined to create a phylogeny. The MLST database
consisted of 10 gene fragments representing published _B. cinerea_ strains that were sampled
from 14 different plants from 12 countries over 12 years.
The obtained cladogram shows the two museum strains as sister species with very high support
(100 bootstrap value) and both were placed as sisters cladding to Bcin sp. D12_BH20_4 from 
raspberry and C12_S_E7_4 from strawberry strains with high support (94 bootstrap value). T4 GCA
000292645.1 reference sequence is placed in the clade, which brunches off before the museum
strains sampled without high support (belove 85 bootstrap value), and the B05.10 NCBI
reference, which cluster together with the in-house B05.10 brunches off in a more basal location
farther along the cladogram. These findings show the museum strains to be closer to T4 than
B05.10 based on the analyzed genes in this dataset.

Command line:

The resulting cladogram:

The results shown by the cladogram indicate

## Selected genes analysis

