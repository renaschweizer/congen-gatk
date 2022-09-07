# ConGen 2022: Reference-based genotyping with GATK
## September 7, 2022

This session is meant to guide the user through alignment of sequence reads to a reference followed by reference-based genotyping with GATK. Why align to a reference? Shafer et al. (2016) compared _de novo_ and reference-based pipelines for RADSeq data and found that using reference-based approaches with a closely-related genome led to summary statistics closest to the expected values. 

Why GATK? Many folks choose to follow the GATK “Best Practices Pipeline”, for many reasons: 1) GATK is freely available and widely used, 2) GATK is under continuous improvement by a team of bioinformaticians at MIT’s Broad Institute, and 3) GATK is well-documented and has a moderated Q&A forum.

Caveats of GATK? Like many software, there's a learning curve. GATK is not designed with non-model organisms in mind, so not all tools function as intended, and you have to make sure you're using it appropriately for your specific study. 

<img width="1200" alt="wholepipeline" src="https://user-images.githubusercontent.com/10552484/132532035-3c7981d9-0728-41b5-a112-5bddcd8de6e4.png">


# Goals

By the end of this session, you should know how to do the following: 
- adequately trim and filter raw sequence reads prior to alignment
- prepare a reference genome for the pipeline (e.g. by indexing)
- align sequence reads to a reference genome using BWA
- call genotypes using the GATK pipeline
- decide on appropriate hard filtering of genotypes

# A note on "pipelines"

As we discussed during the whiteboard session, and as you'll hear throughout ConGen, there is no single pipeline that works for everyone. I put together this particular set of commands after thinking carefully about my particular data set and downstream analysis plans, and after reading carefully through other resources. If someone has shared their scripts with you, make sure those scripts work for your data! There are tons of available sources, such as the GATK Best practice workflows (https://software.broadinstitute.org/gatk/best-practices/) and blogs (example pipeline for non-model organism: https://evodify.com/gatk-in-non-model-organism/). 

## Getting organized

Starting from our **home** directory, let's make a new folder for our genotyping pipeline then change to it. Substitute in your own user ID number, e.g. userN.data

```{bash}
cd ~
mkdir user10.data/genotyping_pipeline
cd user10.data/genotyping_pipeline
```

When I start a project, I like to make a series of directories to organize my data at various stages. This may not apply to us today, but I might also make directories to hold all of my scripts for the project, and other directories to hold bed files or analytical results. 

```{bash}
mkdir reference
mkdir raw_fastq
mkdir pass_qc_fastq
mkdir bam
mkdir bqsr
mkdir raw_genotypes_HC
mkdir called_genotypes_HC
mkdir filtered_genotypes_HC
mkdir quality_metrics
```
Today, we're going to be mapping a subset of the deer mouse exome (the same reads we used for last week's Intro to Command Line session). I have downsampled the data to contain 1 million paired-end reads each for two individuals (S137 and S144, from a low-altitude and high-altitude population). Let's copy these data and other reference files to your genotyping_pipeline directory in RStudio. 

```{bash}
cp ~/instructor_materials/Rena_Schweizer/2022/data_for_gatk_handson/* . 
```
## About the files

Let's see what's in the directory now. 

```{bash}
ls  
```

Output: 

<img width="996" alt="Screen Shot 2022-09-06 at 6 09 33 PM" src="https://user-images.githubusercontent.com/10552484/188761054-902bd61d-a322-4465-b8e8-0caf7ca9d362.png">

The **NW_006501067.1.fa** file is your reference genome in fasta format. The **fastq.gz** files are the raw sequence reads. One of the **vcf** files will be used as a set of "known" variants for the GATK Base Quality Score Recalibration tool, and the other **vcf** file will be used at the end of the pipeline to practice filtering variants. 

Let's finish geting organized by moving our files to our directories.

```{bash}
mv S1* raw_fastq
mv pman_n10_known* reference
mv NW_006501067.1.fa reference
mv pman_n10_congen* raw_genotypes_HC
```

## Step 1: Check sequence quality and remove adapters with trim_galore

We already ran FASTQC on two of these files during our previous session, so we know that overall the quality is high. If you'll remember, we had some small amount of adapter contamination in S144, which we can remove now. To be safe, we can also remove any adapter contamination from the other sample, S137. We're going to be defining variables in the shell then using commands with variables. 

**Remember** we can use `variable=value` to define a new shell variable, to which we can then refer by using `${variable}` in our command.  

I'll show you how you can run these commands for two samples, but think about how you might be able to put these commands into a script and iterate the analysis over thousands of samples.

**Note**: the '\' tells the command line interpreter that your command is continuing on the next line. Writing commands this way allows us to display them more clearly.

You can learn more about trim_galore here: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ 

```{bash}
LANE_NUM=6
for SAMPLE_ID in S137 S144
	do trim_galore --paired raw_fastq/${SAMPLE_ID}_L00${LANE_NUM}_R1_sub.fastq.gz \
	raw_fastq/${SAMPLE_ID}_L00${LANE_NUM}_R2_sub.fastq.gz \
	--phred33 \
	--adapter AGATCGGAAGAGC \
	--quality 20 \
	--output_dir pass_qc_fastq/
	done
```

trim_galore outputs to the terminal some useful information on the parameters and the number and sizes of sequences that are trimmed or removed. 

Task: What files does trim_galore create? 
<details>
  We can use the `ls` command to see that trim_galore makes both a file of trimmed reads for R1 and R2, as well as an output file of the summary statistics. 
</details>  

Task: What percent of reads were retained for R1 and R2 in S144? 
<details> 
  We can use the 'more' command for both report files to see that 98.9% of R1 reads were kept, and 98.1% of R2 reads were kept. 
</details> 

If we want, we can run these trimmed fastq files back through FASTQC to check that the quality is higher and the adapter sequences have been removed. 

## Step 2: Prepare reference genome

In our `reference` directory, we have downloaded the reference genome from NCBI. This is the v1 of the genome for _Peromyscus maniculatus_ and is in multiple scaffolds. We'll only be working with one scaffold (NW_006501067.1) today so that everything runs more quickly. You also have access to a set of 'known' variants for our gatk BQSR analysis later on (pman_n10_known_sites_NW_006501067.1.recode.vcf). Let's change to our reference directory and create the necessary indices that bwa, samtools, and gatk use for quicker access to the genome. This only needs to be done once per project.

```{bash}
cd reference
```

Task: How long is scaffold NW_006501067.1? How many "known" SNPs do we have in our vcf file?
<details>
  ```{bash}
	grep -v ">" NW_006501067.1.fa | wc 
	grep -v "#" pman_n10_known_sites_NW_006501067.1.recode.vcf | wc -l 
	```
</details>  
 
Answer
<details>
  Our NW_006501067.1 scaffold is 11362355 bp, and we have 57663 characterized SNPs on this scaffold. 
</details>  

```{bash}
bwa index NW_006501067.1.fa
samtools faidx NW_006501067.1.fa       
gatk CreateSequenceDictionary -R NW_006501067.1.fa
gatk IndexFeatureFile -I pman_n10_known_sites_NW_006501067.1.recode.vcf
```

You can learn more about bwa here: http://bio-bwa.sourceforge.net/bwa.shtml or in Li and Durbin 2010. 

`samtools` is another utility tool that is designed to handle SAM and BAM format files: http://www.htslib.org/doc/samtools.html

If we check the contents of our directory, we will now see multiple new index files. 

<img width="1242" alt="Screen Shot 2022-09-06 at 6 25 03 PM" src="https://user-images.githubusercontent.com/10552484/188762489-ca017d64-b91b-4175-9312-1190ec74ce4a.png">

## Step 3: Map sequence reads to reference genome

Let's return to our working directory and map our reads to the reference NW_006501067.1. 

```{bash}
cd ../
```

Today we're going to be using `bwa mem` to map, but there are lots of other options out there. Pick the software that is most appropriate for your specific project! 

We can set additional variables for the reference genome path, as well as for the number of threads to tell `bwa` to use. We will only use two threads today, but for your own work you should figure out the available resources for multi-threading. 

```{bash}
REF=`(pwd)`/reference/NW_006501067.1.fa
THREADS=2
```

Note: If you specify a COMMAND within parenthesis and the \`, the shell will interpret and run that COMMAND before it does the rest of whatever the command specifies. So, above, the "`(pwd)`" tells unix to print the working directory and thus provides the complete path for the reference genome file.


In the command below, we are telling `bwa mem` to map each of our filtered R1 and R2 to our reference contig using 2 threads, and to **redirect** the output to a SAM file. The `-M` flag is necessary for compatibility with another program we'll use. Let's run our commands now. 

```{bash}
for SAMPLE_ID in S137 S144
	do bwa mem -M  -t ${THREADS} \
	${REF} \
	pass_qc_fastq/${SAMPLE_ID}_L00${LANE_NUM}_R1_sub_val_1.fq.gz \
	pass_qc_fastq/${SAMPLE_ID}_L00${LANE_NUM}_R2_sub_val_2.fq.gz \
	> bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.sam
	done
```

This will take about 5 minutes to run. Take a quick stretch break if you need it, then, familiarize yourself with the SAM and BAM formats here: https://samtools.github.io/hts-specs/SAMv1.pdf

At the end of this step, we should have a SAM file of aligned, paired end reads for both samples. We can take a look at one of the files using `samtools view`. 

```{bash}
samtools view -h bam/S144_L006_aln-pe.sam | less
```

## Step 4: Prepare BAM files

Once we are satisfied with the mapping, we will convert the SAM files to the compressed BAM format, and then order the BAM by coordinates for downstream analysis. 

### 4a. Convert SAM to BAM

We run the command below to convert SAM to BAM
```{bash}
for SAMPLE_ID in S137 S144
	do samtools view -b -h -S bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.sam \
	> bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.bam
	done
```
**Task: While that's running, let's look at what the b, h, and S flags indicate.**
<details>
  We can run `samtools view` without any parameters to get info on the flags, or we can look at the samtools reference page: http://www.htslib.org/doc/samtools.html
 We learn that `b` specifies output in BAM, `S` specifies the input is in SAM, and `h` includes the header in the output. 
</details>

**Task: We know the compressed format saves space, but how much less space does the BAM file take up?** 
<details>
	
  We can use our handy `ls -lh` to tell us. 
	
  ```{bash}
  ls -lht bam/
  ```
  
From this we learn that the .sam files are ~543 MB, while the .bam files are only 124 MB. Imagine how much space this might save over an entire sequencing project! 
	
</details>
  
### 4b. Fix mates and fill in insert sizes, then sort BAM by coordinates. 

To remove PCR duplicates in the next step, we need to make sure the mate pair information and insert sizes are correct in our BAM using `samtools fixmate`. Software such as GATK requires the BAM file to be sorted by coordinates. Many programs that generate BAM files will automatically sort the output file, but it's good practice to double check this. We'll sort our files using `samtools sort`.

```{bash}
for SAMPLE_ID in S137 S144
	do samtools fixmate -m bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.bam bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.fix.bam
	samtools sort  -@ ${THREADS} -o bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.bam bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.fix.bam
	done
```

### 4c. Check mapping quality

As we did with FASTQC for our raw fastq reads, we can similarly check the quality of our mapped reads using Qualimap: http://qualimap.conesalab.org/. Before we run our command, we will have to change the default display of RStudio so that it outputs an .html file. This is not something you should have to do on your own computer if you run these commands in the future! 

```{bash}
unset DISPLAY
qualimap bamqc -nt ${THREADS} -bam bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.bam -outdir quality_metrics
```

In the output, we see a striking statistic that only 7% of our reads aligned to this reference! Remember that these sequence reads represent an entire exome array, but we have only mapped the reads to a small fraction of the entire genome!

### 4d. Remove PCR duplicates (Option to just mark duplicates, but why?)

Finally, we will remove PCR duplicates using `samtools markdup` and the `-r` flag. 

```{bash}
for SAMPLE_ID in S137 S144
	do samtools markdup -r bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.bam \
	bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.u.bam
	done
```

Task: As a sanity check, let's look at the sizes of the files we have made in the last few steps. Which file do we expect to be the smallest? 
<details>
  We might expect that the SAM file is the largest, and the BAM file after removing duplicates is the smallest. Of course we can check this using `ls -lh bam`
  </details>
  
## Step 5: Base Quality Score Recalibration (BQSR)

### 5a. Add read group information to BAM

BQSR relies on accurate read group information in order to discern technical issues related to quality scores in your sequence data. We can use shell variables to put together the different types of metadata, then add those to our BAM file. 

We will need to define the following: 
- RGID=String	Read Group ID 
- RGLB=String	Read Group Library
- RGPL=String	Read Group platform (e.g. illumina, solid)
- RGPU=String	Read Group platform unit (eg. run barcode of run Id + lane number) 
- RGSM=String	Read Group sample name

Now that we know the specifications that BQSR needs, we can add those to our BAM files using `picard` as implemented in gatk.

```{bash}
for SAMPLE_ID in S137 S144
	do RUN_ID=170503_100PE_HS4KA
	RGID_specify=${RUN_ID}_${LANE_NUM}
	RGLB_specify=${SAMPLE_ID}
	RGPL_specify=ILLUMINA
	RGPU_specify=${RUN_ID}_${LANE_NUM}
	RGSM_specify=${SAMPLE_ID}
	gatk AddOrReplaceReadGroups --INPUT bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.u.bam \
	--OUTPUT bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe_addRG.bam \
	--SORT_ORDER coordinate \
	--RGID ${RGID_specify} \
	--RGLB ${RGLB_specify} \
	--RGPL ${RGPL_specify} \
	--RGPU ${RGPU_specify} \
	--RGSM ${RGSM_specify} \
	--CREATE_INDEX True \
	--VALIDATION_STRINGENCY LENIENT
	done
```
We can check that our command worked by comparing two BAM files for S144. 

```{bash}
samtools view -H bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe.s.u.bam
samtools view -H bam/${SAMPLE_ID}_L00${LANE_NUM}_aln-pe_addRG.bam
```

Task: What is the difference between these two files? 

<details>
  In the header of the "addRG.bam" there is an additional line with an @RG tag:
  @RG     ID:170503_100PE_HS4KA_6 LB:S144 PL:ILLUMINA     SM:S144 PU:170503_100PE_HS4KA_6
  </details>
  
### 5b. Run GATK's BaseRecalibrator

This is the first pass of the base quality score recalibration. The `BaseRecalibrator` generates a recalibration table based on various covariates. The default covariates are read group, reported quality score, machine cycle, and nucleotide context. BQSR is most accurate at identifying technical artifacts when it is run with as many of your sequenced indivdiuals as possible. 

```{bash}
gatk BaseRecalibrator \
  -I bam/S137_L00${LANE_NUM}_aln-pe_addRG.bam \
  -I bam/S144_L00${LANE_NUM}_aln-pe_addRG.bam \
  -R reference/NW_006501067.1.fa \
  --known-sites reference/pman_n10_known_sites_NW_006501067.1.recode.vcf \
  -O bqsr/pman_recal1_data.table
```

This generates a recalibration text file with information on things such as the quality scores of each read group for mismatches, insertions and deltions, or empirical qualities for each covariate used in the dataset (e.g. context and cycle). Again, you can read more about the format of this table on the gatk website: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531 

Next, we apply this recalibration info to our original BAM file to make a new BAM with recalibrated quality scores. 

```{bash}
gatk ApplyBQSR \
  -R reference/NW_006501067.1.fa \
  -I bam/S137_L00${LANE_NUM}_aln-pe_addRG.bam \
  -I bam/S144_L00${LANE_NUM}_aln-pe_addRG.bam \
  --bqsr-recal-file bqsr/pman_recal1_data.table \
  -O bqsr/pman_recal1_output.bam
```

We'll do the above steps one more time to generate a count of covariates on our recalibrated data, then we will compare our before and after BQSR results. 

```{bash}
gatk BaseRecalibrator \
  -I bqsr/pman_recal1_output.bam \
  -R reference/NW_006501067.1.fa \
  --known-sites reference/pman_n10_known_sites_NW_006501067.1.recode.vcf \
  -O bqsr/pman_recal2_data.table
  
gatk ApplyBQSR \
  -R reference/NW_006501067.1.fa \
  -I bam/S137_L00${LANE_NUM}_aln-pe_addRG.bam \
  -I bam/S144_L00${LANE_NUM}_aln-pe_addRG.bam \
  --bqsr-recal-file bqsr/pman_recal2_data.table \
  -O bqsr/pman_recal2_output.bam  
```

We can use a gatk tool called `AnalyzeCovariates` to see how the first round of BQSR affected our data. We can look at some example plots on the gatk website to get an idea for what our before and after plots should look like: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531

```{bash}
gatk AnalyzeCovariates \
  -before bqsr/pman_recal1_data.table \
  -after bqsr/pman_recal2_data.table \
  -plots bqsr/AnalyzeCovariates_round1.pdf
```

We may choose to do additional rounds of BQSR until we are satisfied with the patterns we observe in the covariate 'before' and 'after' plots. At the end of Step 5 we should have **analysis-ready BAM files**. 

## Step 6: Run HaplotypeCaller to generate GVCF files

### 6a: Call variants per-sample using gatk's HaplotypeCaller

From the GATK website: "The HaplotypeCaller calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. HaplotypeCaller runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way."

```{bash}
for SAMPLE_ID in S137 S144
	do gatk HaplotypeCaller \
	-I  bqsr/pman_recal2_output.bam \
	--sample-name ${SAMPLE_ID} \
	-R reference/NW_006501067.1.fa \
	--emit-ref-confidence GVCF \
	-O raw_genotypes_HC/${SAMPLE_ID}.raw.snps.indels.g.vcf
	done
```
Don't be alarmed if you get some warnings output to the screen. There was a known bug with this version of GATK that has since then been fixed! This will take a few miutes to run. 

We should now have a GVCF file for each individual. We can take a look at the files, but don't be alarmed if they say there are no genotypes yet! We'll do that soon.  The GVCF mode and GVCF files can be read about in more detail here: https://software.broadinstitute.org/gatk/documentation/article.php?id=4017 

### 6b: Consolidate all GVCF files

This step consists of consolidating the contents of GVCF files across all of your samples in order to improve scalability and speed the next step, joint genotyping. We will use the `GenomicsDBImport` tool and will need to make a new temporary directory. 

```{bash}
mkdir called_genotypes_HC/temp
gatk GenomicsDBImport \
      --variant raw_genotypes_HC/S137.raw.snps.indels.g.vcf \
      --variant raw_genotypes_HC/S144.raw.snps.indels.g.vcf \
      --genomicsdb-workspace-path called_genotypes_HC/pman_database \
      -L NW_006501067.1
 ``` 

## Step 7: Genotype all GVCF files

```{bash}
gatk GenotypeGVCFs \
   -R reference/NW_006501067.1.fa \
   -V gendb://called_genotypes_HC/pman_database \
   --include-non-variant-sites \
   -O raw_genotypes_HC/pman_HC_genotypes_unfiltered.vcf.gz
```

At this point in our analysis, we should have a set of all called positions on the NW_006501067.1 scaffold of the deer mouse genome. If we look at the file more closely, though, we see some positions that were only called in one of the samples, or others that seem to have a pretty low depth of coverage. In the next step, let's do some "hard" filtering of positions. A shortcoming of working with downsampled data is that there weren't enough variable sites in our output from `GenotypeGVCFs`, so we'll be working with a different vcf file for the filtering.

## Step 8: Perform hard filtering

Filtering a set of SNPs is a delicate balance between too stringent and too lenient; you should decide what works best for your type of data and your analytical goals. For a discussion of a few aspects of filtering that you should keep in mind, see the section "Navigating the Perils of Data Filtering" in last year's ConGen review (Schweizer et al. 2021). 

For this part of the pipeline, we're going to use a different vcf file that has additional individuals and more variable site. At the beginning of this session, we moved this file and its index, to our **raw_genotypes_HC** directory. 

We will do some fairly simple filtering, but enough so that you have an idea of the tools and metrics you can use to filter your own data. 

First, let's take a look at our vcf file and notice a few things. 

```{bash}
zcat raw_genotypes_HC/pman_n10_congenSubset.vcf.gz | grep -v "#" | less
```

1. The first thing to notice is that this file has both SNPs and indels. A benefit of the GATK HaplotypeCaller algorithm is that it can determine SNPs and indels simultaneously as it does local realignment of the sequences. However, for our purposes, we do not need indels, so we'll filter them out in the next step. 
2. The HaplotypeCaller also has a number of **position-specific** and **individual-specific** default annotations that it performs, such as the strand bias (uneven amplifications of DNA strands), the allele depth (number of reads mapping to each allele), depth of coverage, genotype quality, etc. You can find out what these annotations are by looking in the header or on the GATK website: https://gatk.broadinstitute.org/hc/en-us/articles/360047216111--Tool-Documentation-Index#VariantAnnotations

There are a few options for filtering vcf files. gatk has two tools called `VariantFiltration` which filters variant calls based on INFO and/or FORMAT annotations and `SelectVariants` which selects a subset of variants from a VCF file. We can also use vcftools, which we learned about during the bioinformatics session. 

Now, we'll select only the SNPs from the file. 

```{bash}
gatk SelectVariants \
	-R reference/NW_006501067.1.fa \
	-O filtered_genotypes_HC/pman_n10_congenSubset_raw_SNPs.vcf \
	-V raw_genotypes_HC/pman_n10_congenSubset.vcf.gz \
	--select-type-to-include SNP
```

Using vcftools, we can assess aspects of our data such as the mean depth of coverage per individual or per site, or the number of missing sites per individual or number of individuals missing from each site. Once we have an idea of what our data look like, we can decide what are appropriate filters. 

```{bash}
for flag in missing-indv missing-site depth site-depth
	do vcftools --vcf filtered_genotypes_HC/pman_n10_congenSubset_raw_SNPs.vcf \
	--${flag} \
	--out quality_metrics/pman_n10_congenSubset_raw_SNPs
	done
```

The above **for loop** generates some summary statistics on missingness and depth of coverage across our samples and our SNPs. 

For example, we might choose to filter out sites below a certain genotype quality (`--minGQ`), below a minimum depth of coverage (`minDP`),  retain sites that are bi-allelic, and those sites that are called in 75% of our individuals (`max-missing`).

```{bash}
vcftools --vcf filtered_genotypes_HC/pman_n10_congenSubset_raw_SNPs.vcf \
	--minDP 10 \
	--minGQ 20 \
	--max-missing 0.75 \
	--min-alleles 2 \
	--max-alleles 2 \
	--recode \
	--out filtered_genotypes_HC/pman_n10_congenSubset_minDP10_minGQ20_biAllelic_maxMissing75
```

You might check your filters by examining genotype concordance (if you have a set of known variants, e.g. from a SNP array) with `vcftools --diff-discordance-matrix` or by measuring the transition/transversion ratio with `vcftools --TsTv`


# Other resources: 

The Broad's GATK website is an excellent resource, with very helpful tutorials, detailed manuals, and a forum for asking questions. https://gatk.broadinstitute.org/hc/en-us

# References: 
Hendricks S, Anderson EC, Antão T, Bernatchez L, Forester BR, Garner B, Hand BK, Hohenlohe PA, Kardos M, Koop B, et al. 2018. Recent advances in conservation and population genomics data analysis. Evolutionary Applications 11:1197–1211.

Li H, Durbin R. 2010. Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics 26:589.

Schweizer RM, Saarman N, Ramstad KM, Forester BR, Kelley JL, Hand BK, Malison RL, Ackiss AS, Watsa M, Nelson TC, et al. 2021. Big Data in Conservation Genomics: Boosting Skills, Hedging Bets, and Staying Current in the Field. Journal of Heredity 112:313–327.

Shafer ABA, Peart CR, Tusso S, Maayan I, Brelsford A, Wheat CW, Wolf JBW. 2016. Bioinformatic processing of RAD-seq data dramatically impacts downstream population genetic inference.Gilbert M, editor. Methods in Ecology and Evolution 8:907–917.
