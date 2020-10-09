[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)
# Aperture:  Accurate detection of structural variations and viral integrations in circulating tumor DNA using an alignment-free algorithm

Aperture is a new alignment-free SV caller designed for ctDNA sequencing with barcode adapters and multiple duplicates. Aperture applies k-mer based searching, fast intersection and breakpoint merging to detect SVs and viral integrations in high sensitivity, especially when junctions span repetitive regions, followed by a barcode based filter to ensure specificity. Aperture takes paired-end reads in fastq format as inputs and reports all SVs and viral integrations in VCF 4.2 format.  
  
If you have any trouble running Aperture, please raise an issue using the Issues tab above.  
  
[Click here to download Aperture](https://github.com/i8q8r9/Aperture/releases)

# Citation
For citing Aperture and for an overview of the Aperture algorithms, refer to our open access article:  

**Aperture: Accurate detection of structural variations and viral integrations in circulating tumor DNA using an alignment-free algorithm.**
Hongchao Liu, Guangyu Li, Xiaoyue Wang.

# Software and Hardware Requirements
## Software Requirements
To run Aperture, java 1.8 or later version must be installed in your system.
## Hardware Requirements
* **CPU** Aperture does not require or benefit from any specific modern CPU feature, but more physical cores and faster clock will significantly improve performance.
* **Memory** Typically, Aperture needs 40GB in index building and 30GB in SV calling for human genome (hg19 or hg38). The exact requirement depends on many factors including reference genome, sequencing depth, cfDNA insert size and sample quality.

# Running
Aperture takes a Aperture index and a set of ctDNA read files and outputs SV results in VCF format.  
Pre-compiled binaries are available at https://github.com/i8q8r9/Aperture/releases. 

## Building an Aperture index
Aperture needs a indexed sequence file (in FASTA and FAI format) and a corresponding common SNP database (in VCF format) to build Aperture index. If FAI file is missing, you can use `faidx` command in `samtools` to create one. Aperture outputs a set of 5 files with suffixes `.ci` `.tt` `.km` `.long.km` and `.spaced.km`. These files together constitute the index, and the original FASTA files are no longer used by Aperture once the index is built.   
  
Pre-built Aperture indexs for hg19 and hg38 are available at
### Command-line arguments
```
Usage: java -jar aperture.jar index -R <genome.fa> -V <snp.vcf> -O <out> -T <threads>
```

argument|description
---|---
-h,--help|Show help message
-O,--out <arg>|Output path
-R,--reference <arg>|Genome FASTA file with fai index
-T,--threads <arg>|Number of threads
-V,--vcf <arg>|Common SNPs database for the corresponding genome

### Example
```
java -Xmx40g -jar fusion_test/aperture12.jar index -R ~/fusion_test/ref/hg38HBV.fa -V ~/fusion_test/ref/dbsnp_common_hg38.vcf -O ~/fusion_test/ref/hg38_ap12_test -T 30
```  

## Detecting SVs and viral integrations

Aperture needs a pair of FastQ files and an Aperture index as input. The output is in compressed VCF format (.vcf.gz). Aperture supports barcode based filter to ensure specificity. So if your dataset is produced by abundant sequencing and contains barcode as unique molecular identifier, parameters including `-1BS`, `-2BS`, `-1BL`, `-2BL`, `-1S` and `-2S` should be used to specify the location of barcodes in a read.
ATCGGTTTGCCAAAGTATTGCTTATAA
                  TGCTTCTAAGTCCGTAAGTAGTGTCCG
### Command-line arguments
```
Usage: java -jar aperture.jar call  -1 <arg> -1BL <arg> -1BS <arg> -1S <arg> -2 <arg> -2BL <arg> -2BS <arg> -2S <arg> -D <arg> [-H] -I <arg> -P <arg> -T <arg>
```

Argument|Description
---|---
-1,--r1 <arg>|Path of R1.fq.gz
-1BL,--r1BarLen <arg>|Length of barcode in R1
-1BS,--r1BarStart <arg>|Barcode start index in R1 (0-based)
-1S,--r1InsStart <arg>|ctDNA fragment start index in R1 (0-based)
-2,--r2 <arg>|Path of R2.fq.gz
-2BL,--r2BarLen <arg>|Length of barcode in R2
-2BS,--r2BarStart <arg>|Barcode start index in R2 (0-based)
-2S,--r2InsStart <arg>|ctDNA fragment start index in R2 (0-based)
-D,--dir <arg>|Output path
-H,--help|Show help message
-I,--index <arg>|Path of Aperture index
-P,--project <arg>|Project name
-T,--threads <arg>|Number of threads
### Example
```
java -Xms50g -Xmx50g -jar ~/fusion_test/aperture12.jar call -1 SRR8551545_1.fastq.gz -2 SRR8551545x_2.fastq.gz -I ~/fusion_test/ref/hg38_ap12 -D SRR8551545 -P SRR8551545_ap12 -1BS 0 -2BS 0 -1BL 6 -2BL 0 -1S 0 -2S 0 -T 30
```
  
# Output interpretation
In Aperture, all SVs are described as breakends and thus all the records in Aperture VCF are identified with the tag “SYTYPE=BND” in the INFO field.  
  
Aperture VCF output follows the VCF 4.2 spec. All custom fields are described in the VCF header. 

## VCF FILTER Fields
ID | Description
------- | ---------
LOW_QUAL|Low quality call
FAKE_BP|False positive variant caused by imprecise k-mer based mapping
SMALL_EVENT|Event size is smaller than the minimum reportable size
PE_ONLY|No soft clips support this variant

## VCF INFO Fields
ID | Description
------- | ---------
SVTYPE|Type of structural variant
STRANDS|Strand orientation of the adjacency
REFQUA|K-mer mapping quality of reference junction
VARQUA|K-mer mapping quality of variant junction
REFKMER|Number of k-mers supporting reference junction in average
VARKMER|Number of k-mers supporting variant junction in average
BPSEQQUA|Quality of sequence spanning breakpoint junction
PARID|ID of partner breakend
HOMLEN|Length of base pair identical micro-homology at event breakpoints
HOMSEQ|Sequence of base pair identical micro-homology at event breakpoints	

## VCF FORMAT Fields
ID | Description
------- | ---------
GT|Genotype (Not applicable)
SR|Count of split reads supporting the breakpoint
PE|Count of paired-end reads supporting the breakpoint
REFSR|Count of split reads supporting the reference junction
VARSR|Count of split reads supporting the variant junction
BAR|Count of cfDNA molecules supporting the breakpoint
UBAR|Count of cfDNA molecules with only one read support