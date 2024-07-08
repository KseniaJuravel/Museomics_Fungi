# Museomics_Fungi
Work performed for ancient Botrytis fungi analysis of two strains

## Production of whole-genome sequencing data

Genomic DNA was extracted from the colonies using the ????? protocol previously described in ref. ???? or the ???? Kit . The genomes were sequenced either at the ????? Platform of Institut ????, the Department ???? using the Illumina ??? sequencing technology. Paired-end reads of ???? bp were obtained. Reads have been deposited at the NCBI Sequence Read Archive under BioProject ID PRJN######.

## The analysis of the whole-genome sequencing data

Each set of paired-end reads was mapped against the ??? reference genome ??? haplotype ???  downloaded from the ???? database (version ?? ??Date) using the Bowtie2 Alignment tool, version 2.4.1 (Built on Fri Feb 28 17:23:43 UTC 2020). SAMtools samtools 1.9, Using htslib 1.9 and Picard tools version 1.70 (http://broadinstitute.github.io/picard) to filter, sort and convert SAM files. 

SNPs were called using Genome Analysis Toolkit version 4.2.0.0 and 4.2.5.0 according to the GATK Best Practices. SNPs and indels were filtered using the following parameters: VariantFiltration, QD < 2.0, LowQD, ReadPosRankSum < −8.0, LowRankSum, FS > 60.0, HightFS, MQRankSum < −12.5, MQRankSum, MQ < 40.0, LowMQ, HaplotypeScore > 13.0, HaploScore. Coverages were calculated using the Samtools mpileup toolkit. Several tools were used to ensure the accuracy of the genomic variants detected from the data.

The parameters are shown in other fungi recent research:

https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000979#cited

Li X, Muñoz JF, Gade L, Argimon S, Bougnoux M-E, Bowers JR, Chow NA, Cuesta I, Farrer RA, Maufrais C, et al. 2023. Comparing genomic variant identification protocols for Candida auris. Microbial Genomics 9:1–19.

https://www.nature.com/articles/s41467-018-04787-4#citeas

Ropars, J., Maufrais, C., Diogo, D. et al. Gene flow contributes to diversification of the major fungal pathogen Candida albicans. Nat Commun 9, 2253 (2018). https://doi.org/10.1038/s41467-018-04787-4


<details>
<summary>Alignment of data</summary>

Sequencing data obtained for each of the strains:

```
515M Feb  3  2021 GBOS_S4_R1_001.fastq.gz.2.fq.gz
108M Feb  3  2021 G7B3_S5_R1_001.fastq.gz.2.fq.gz
410M Feb  3  2021 G2307_S3_R1_001.fastq.gz.2.fq.gz
3.2G Feb  1  2021 SAD_B05_S10_R1_001.fastq.gz
663M Feb  1  2021 SAD_7B3_S11_R1_001.fastq.gz
2.8G Feb  1  2021 SAD_234_S9_R1_001.fastq.gz
774M Nov 23  2020 GBOS_S4_R1_001.fastq.gz

```

Files concatenated to represent the strains of interest as follows:

```
cat GBOS_S4_R1_001.fastq.gz GBOS_S4_R1_001.fastq.gz.2.fq.gz SAD_B05_S10_R1_001.fastq.gz > B05.fastq.gz
cat SAD_7B3_S11_R1_001.fastq.gz G7B3_S5_R1_001.fastq.gz.2.fq.gz > 903054.fastq.gz
cat SAD_234_S9_R1_001.fastq.gz G2307_S3_R1_001.fastq.gz.2.fq.gz > 903053.fastq.gz

```

The concatenated files were checked for quality using fastqc tool (FastQC v0.11.8).



Reference genome obtained from NCBI for [_Botrytis cinerea B05.10_ GCF_000143535.2 from Feb 5, 2015](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000143535.2/) and [_Botrytis cinerea T4_ GCA_000292645.1 from Aug 22, 2012](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000292645.1/).

Bowtie2 (bowtie2-build-s version 2.3.4.3 64-bit) was used for alignment with the following parameters:

```
bowtie2 -p 70 -x $1 --very-fast --preserve-tags --no-unal -q -U $2 -2 $3 -S $3.sam;
samtools view -@ 70 -h -S -b -o $3.bam $3.sam;
samtools sort --threads 70 $3.bam -O BAM -o $3.sort.bam;
samtools index -@ 70 $3.sort.bam;

```


</details>


<details>
<summary>Calling variants from alignment</summary>

```


```
</details>



<details>
<summary>Variants comparison for the different strains against the reference genomes</summary>

<details>
<summary>Comparison to _Botrytis cinerea T4_ GCA_000292645.1</summary>



<details>
<summary>Comparison to _Botrytis cinerea B05.10_ GCF_000143535.2</summary>



</details>
</details>
</details>
























<details>
<summary>Alignments coverage</summary>

<details>
<summary>Command</summary>
 
```
samtools mpileup B05_bowtie_vs_2.3.sam.bam.sorted.bam | awk '{ count++ ; SUM += $4 } END { print "Total: " SUM "\t" "Nucleotides: " count "\t" "Average_coverage: " SUM/count }'
[mpileup] 1 samples in 1 input files
 
```
 </details>

```
T4 (ALOC0100000) ref with 903053 reads alignment has 

Total: 3014160672       Nucleotides: 37443825   Average_coverage: 80.4982
```
```
T4 (ALOC0100000) ref with 903054 reads alignment has Total:

Total: 1057337721       Nucleotides: 37370932   Average_coverage: 28.2931
```
```
T4 (ALOC0100000) ref with B05.10 local reads alignment has Total:

Total: 4986300933       Nucleotides: 37481885   Average_coverage: 133.032
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
<summary>VCF Correction command</summary>

```
sbatch -A gila.kahila GATK_correction.sh B05_REF_normalized.fasta B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.vcf.gz

gatk --java-options "-Xmx4g" VariantFiltration --reference T4_REF.fa --variant 054_T4.g.vcf.gz --filter-expression "QD < 2.0" --filter-name "SNP_QD" --filter-expression "FS > 60.0" --filter-name "SNP_FS" --filter-expression "SOR > 4.0" --filter-name "SNP_SOR" --filter-expression "MQ < 40.0" --filter-name "SNP_MQ" --filter-expression "MQRankSum < -12.5" --filter-name "SNP_MQRankSum" --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum" --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz

gatk --java-options "-Xmx4g" SelectVariants --reference T4_REF.fa --variant 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz -exclude-filtered --exclude-non-variants --output 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
Using GATK jar /usr/local/hurcs/gatk4/4.2.5.0/share/gatk4-4.2.5.0-0/gatk-package-4.2.5.0-local.jar

```
 </details>

<details>
<summary>Results for VCF comparisons after correction with multiple tools</summary>

<details>
<summary>Results T4 (ALOC0100000) vs. strains 903053 and 903054 </summary>

<details>
<summary>Plot of locations of variants</summary>
 
<details>
<summary>Tool </summary>
```
Using Rplot https://www.bioinformatics.com.cn/plot_basic_SNP_density_by_CMplot_107_en
```
 </details>



<details>
<summary>T4_vs_903054</summary>

![T4_vs_903054](https://github.com/KseniaJuravel/Ancient_Fungi/blob/main/VCF_output/T4_vs_903054/13f5f469b08a2ec5.png)
 </details>

<details>
<summary>T4_vs_903053</summary>

![T4_vs_903053](https://github.com/KseniaJuravel/Ancient_Fungi/blob/main/VCF_output/T4_vs_903053/0b452b346e3a3e9c.png)
 </details>
<details>
<summary>T4_vs_B05.10 genome (AAID02000000)</summary>

  </details>

 </details>



<details>
<summary>Tool #1 for VCF analysis - bcftools</summary>


<details>
<summary>Command</summary>

```
bcftools stats                      054_T4.g.vcf.gz 053_T4.g.vcf.gz > joined_T4_2.3_vs_7B.stats.txt

plot-vcfstats                      joined_T4_2.3_vs_7B.stats.txt -p outdir_T4
```

 </details>


Figure Total counts for indels and SNPs:

![](https://github.com/KseniaJuravel/Ancient_Fungi/blob/main/Figures%26Data/corrected_outdir_T4/venn_bars.snps.png)

More comparisons can be found in the folder. 


 </details>

 
<details>
<summary>Tool #2 for VCF analysis - vt peek & multi-partition </summary>


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
<summary>vt multi-partition</summary>


```
(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$ ../../VCF_OUTPUT/vt/vt multi_partition 053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       70822 variants
    B:       94180 variants

               no  [ts/tv] [ins/del]
    A-       1889  [1.98]  [1.37]
    -B      25247  [2.72]  [0.90]
    AB      68933  [3.00]  [0.84]

    Unique variants     :      96069
    Overall concordance :      71.75% (#intersection/#union)

Time elapsed: 0.38s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$ ../../VCF_OUTPUT/vt/vt multi_partition 053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file c   B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       70822 variants
    B:       94180 variants
    C:       87062 variants

                no  [ts/tv] [ins/del]
    A--       1534  [1.89]  [1.44]
    -B-      16481  [2.64]  [0.87]
    AB-      36874  [2.95]  [0.91]
    --C      45882  [2.86]  [0.94]
    A-C        355  [2.33]  [0.78]
    -BC       8766  [2.85]  [1.01]
    ABC      32059  [3.04]  [0.75]

    Unique variants     :     141951
    Overall concordance :      22.58% (#intersection/#union)

Time elapsed: 0.58s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$ ../../VCF_OUTPUT/vt/vt multi_partition 053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   053_T4.g.vcf.gz.RGsorteer.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       70822 variants
    B:       87062 variants

               no  [ts/tv] [ins/del]
    A-      38408  [2.91]  [0.99]
    -B      54648  [2.86]  [0.95]
    AB      32414  [3.04]  [0.75]

    Unique variants     :     125470
    Overall concordance :      25.83% (#intersection/#union)

Time elapsed: 0.37s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$ ../../VCF_OUTPUT/vt/vt multi_partition 054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   054_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05_T4.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       94180 variants
    B:       87062 variants

               no  [ts/tv] [ins/del]
    A-      53355  [2.86]  [0.90]
    -B      46237  [2.85]  [0.94]
    AB      40825  [3.00]  [0.81]

    Unique variants     :     140417
    Overall concordance :      29.07% (#intersection/#union)

Time elapsed: 0.42s
```

</details>

 </details> 


<details>
<summary>Results B05.10 genome (AAID02000000) vs. strains 903053 and 903054 </summary>

<details>
<summary>Plot of locations of variants</summary>

<details>
<summary>B05_vs_B05.10 genome (AAID02000000)</summary>

![](https://github.com/KseniaJuravel/Ancient_Fungi/blob/main/VCF_output/B05NCBI_vs_B05Local/8293d5a8309d2e5b.png)

  </details>

<details>
<summary>903053_vs_B05.10 genome (AAID02000000)</summary>

  </details>


<details>
<summary>903054_vs_B05.10 genome (AAID02000000)</summary>

  </details>
  
 </details>


<details>
<summary>Tool #1 for VCF analysis - bcftools</summary>

```

```

 </details>


<details>
<summary>Tool #2 for VCF analysis - vt peek & multi-partition</summary>

```
vt peek B05 vs 903053

stats: no. of samples                     :          1
       no. of chromosomes                 :         17

       ========== Micro variants ==========

       no. of SNP                         :      72862
           2 alleles                      :           72859 (2.92) [54257/18602]
           3 alleles                      :               3 (1.00) [3/3]

       no. of INDEL                       :       3439
           2 alleles                      :            3419 (0.99) [1698/1721]
           3 alleles                      :              20 (0.67) [16/24]

       no. of SNP/INDEL                   :          5
           3 alleles                      :               5 (1.50) [3/2] (inf) [5/0]

       no. of micro variants              :      76306

       ++++++ Other useful categories +++++

        no. of complex substitutions      :          5
           3 alleles                      :               5 (1.50) [3/2] (inf) [5/0]


       ========= General summary ==========

       no. of VCF records                        :      76306

```

```
vt peek B05 vs 903054

stats: no. of samples                     :          1
       no. of chromosomes                 :         17

       ========== Micro variants ==========

       no. of SNP                         :      96520
           2 alleles                      :           96515 (2.89) [71673/24842]
           3 alleles                      :               5 (0.11) [1/9]

       no. of INDEL                       :       4948
           2 alleles                      :            4912 (0.94) [2386/2526]
           3 alleles                      :              36 (0.95) [35/37]

       no. of SNP/INDEL                   :          6
           3 alleles                      :               6 (0.20) [1/5] (1.00) [4/4]

       no. of micro variants              :     101474

       ++++++ Other useful categories +++++

        no. of complex substitutions      :          6
           3 alleles                      :               6 (0.20) [1/5] (1.00) [4/4]


       ========= General summary ==========

       no. of VCF records                        :     101474

```

```

vt peek B05 vs B05 in house
stats: no. of samples                     :          1
       no. of chromosomes                 :         18

       ========== Micro variants ==========

       no. of SNP                         :        147
           2 alleles                      :             146 (1.28) [82/64]
           3 alleles                      :               1 (0.00) [0/2]

       no. of INDEL                       :        305
           2 alleles                      :             290 (0.54) [102/188]
           3 alleles                      :              15 (0.30) [7/23]

       no. of SNP/INDEL                   :          3
           3 alleles                      :               3 (0.00) [0/3] (0.25) [1/4]

       no. of micro variants              :        455

       ++++++ Other useful categories +++++

        no. of complex substitutions      :          3
           3 alleles                      :               3 (0.00) [0/3] (0.25) [1/4]


       ========= General summary ==========

       no. of VCF records                        :        455

```

 </details>


<details>
<summary>vt multi-partition</summary>
 
```

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$  ../../VCF_OUTPUT/vt/vt multi_partition B05normalized_bowtie_vs_*
multi_partition v0.5

Options:     input VCF file a   B05normalized_bowtie_vs_2.3.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05normalized_bowtie_vs_7B.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file c   B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       76306 variants
    B:      101474 variants
    C:         455 variants

                no  [ts/tv] [ins/del]
    A--       1948  [2.13]  [1.21]
    -B-      27116  [2.74]  [0.93]
    AB-      74288  [2.94]  [0.96]
    --C        372  [0.76]  [0.49]
    A-C         13  [0.00]  [0.86]
    -BC         13  [1.00]  [0.17]
    ABC         57  [5.67]  [0.89]

    Unique variants     :     103807
    Overall concordance :       0.05% (#intersection/#union)

Time elapsed: 0.41s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$  ../../VCF_OUTPUT/vt/vt multi_partition B05normalized_bowtie_vs_2.3.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz B05normalized_bowtie_vs_7B.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   B05normalized_bowtie_vs_2.3.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05normalized_bowtie_vs_7B.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       76306 variants
    B:      101474 variants

               no  [ts/tv] [ins/del]
    A-       1961  [2.13]  [1.19]
    -B      27129  [2.74]  [0.92]
    AB      74345  [2.94]  [0.96]

    Unique variants     :     103435
    Overall concordance :      71.88% (#intersection/#union)

Time elapsed: 0.41s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$  ../../VCF_OUTPUT/vt/vt multi_partition B05normalized_bowtie_vs_2.3.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   B05normalized_bowtie_vs_2.3.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:       76306 variants
    B:         455 variants

               no  [ts/tv] [ins/del]
    A-      76236  [2.92]  [0.99]
    -B        385  [0.77]  [0.48]
    AB         70  [4.86]  [0.88]

    Unique variants     :      76691
    Overall concordance :       0.09% (#intersection/#union)

Time elapsed: 0.18s

(base) ksenia.juravel@glacier-12:/sci/labs/gila.kahila/ksenia.juravel/aDNA_Fungi/VCF_OUTPUT/PASS$  ../../VCF_OUTPUT/vt/vt multi_partition B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz B05normalized_bowtie_vs_7B.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
multi_partition v0.5

Options:     input VCF file a   B05normalized_bowtie_vs_GB05.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz
             input VCF file b   B05normalized_bowtie_vs_7B.sam.sorted.bam.out.bam.sorted2.bam.g.vcf.gz.RGsorted.HaplotypeCaller.all.snp.filtered.vcf.gz.pass.vcf.gz

    A:         455 variants
    B:      101474 variants

               no  [ts/tv] [ins/del]
    A-        385  [0.75]  [0.50]
    -B     101404  [2.88]  [0.95]
    AB         70  [4.11]  [0.60]

    Unique variants     :     101859
    Overall concordance :       0.07% (#intersection/#union)

Time elapsed: 0.24s

```

 </details>


 </details>
 </details>
 </details>




![miltu-peek](https://github.com/KseniaJuravel/Museomics_Fungi/blob/main/Figures%26Data/vt_results.png)

The results shown by the tools indicate

 



 ## Phylogenetic analyses

Phylogenetic analysis of multiple internal fragments of house-keeping genes (n=10) of the two
museum genomes with contemporary Botrytis genus genomes using the Multilocus sequence
typing (MLST) (Plesken et al., 2021) combined to create a phylogeny. The MLST database
consisted of 10 gene fragments representing published _B. cinerea_ strains that were sampled
from 14 different plants from 12 countries over 12 years.
The obtained cladogram for 9266 nucleotide positions in the alignment of 95 strains shows the two museum strains
as sister species with very high support (100 bootstrap value) 
and both were placed as sisters cladding to Bcin sp. D12_BH20_4 from 
raspberry and C12_S_E7_4 from strawberry strains with high support (94 bootstrap value). T4 GCA
000292645.1 reference sequence is placed in the clade, which brunches off before the museum
strains sampled without high support (belove 75 bootstrap value), and the B05.10 NCBI
reference, which cluster together with the in-house B05.10 brunches off in a more basal location
farther along the cladogram. These findings show the museum strains to be closer to T4 than
B05.10 based on the analyzed genes in this dataset.

The command line:

```
for i in T4_MLST_10_Genes_subalign_subalign.fa; do iqtree -s $i -pre $i -m MFP -bb 1000 -alrt 5000 -nt AUTO; done
```

The alignment:

![Aln](https://github.com/KseniaJuravel/Museomics_Fungi/blob/main/Phylogeny/Screenshot%202024-06-04%20102029.png)

The resulting cladogram:


![Tree](https://github.com/KseniaJuravel/Museomics_Fungi/blob/main/Phylogeny/Phylogeny.png)

The results shown by the cladogram indicate museum strains to be closer to T4 than
B05.10 based on the analyzed genes in this dataset. Arrows indicate the datasets used for the analysis in this study.

## Selected genes analysis

<br>




</br>

The command line:

```


```

