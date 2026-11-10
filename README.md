# CellMol_SexChrom_Acinonyx_jubatus

Name: Melon, Kris Bernadette S.

Section: B

Selected species: *Acinonyx jubatus*

Study Question: How can variations in sequencing read depth among the autosomes and sex chromosomes help determine the possible chromosomal sex of an Acinonyx jubatus (cheetah) whole-genome sequencing (WGS) sample?

**Assembly association:** GCF_027475565.1 
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_027475565.1/

**Assembly BioSample accession:** SAMN29432198 
https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN29432198

**WGS accession (Run):** SRR22273180  
https://www.ncbi.nlm.nih.gov/sra/SRX18249254

## Galaxy 

<img width="1061" height="871" alt="Screenshot 2026-11-11 000217" src="https://github.com/user-attachments/assets/0be219b0-bf63-47a9-8d21-c03a7b8edba9" />
**_Figure 1_**. Uploaded genomic FASTA file and assembly report for Acinonyx jubatus (cheetah), which provide the reference genome and supporting information for the analysis.

<img width="1061" height="871" alt="Screenshot 2026-11-11 000217 - Copy (2)" src="https://github.com/user-attachments/assets/52ee2935-65f7-4d4f-8bf6-efd9a1307791" />
**_Figure 2_**. Faster Download and Extract Reads in FASTQ format tool configured with the SRA accession number SRR22273180 to retrieve the cheetah sequencing data.

<img width="887" height="872" alt="Screenshot 2026-11-11 000325" src="https://github.com/user-attachments/assets/a2d7ecec-951f-4609-b365-78b42db4e10c" />
**_Figure 3_**. FastQC tool parameters used to assess the quality of the sequencing reads before proceeding with the downstream genomic analysis.

<img width="866" height="552" alt="Screenshot 2026-11-11 000400" src="https://github.com/user-attachments/assets/25e09bcf-8698-4edd-ad6a-0230d0c5c349" />
**_Figure 4_**. Configuration of the Unzip Collection tool, which separates the paired-end sequencing datasets into individual forward and reverse read files for subsequent processing.

<img width="868" height="848" alt="Screenshot 2026-11-11 000501" src="https://github.com/user-attachments/assets/cf3e6b28-4908-4534-aa3f-5f6fd0ed712f" />
**_Figure 5_**. BWA-MEM2 tool settings, including the selected reference genome and paired-end sequencing reads, which are used to align the sequencing data to the cheetah reference genome.

<img width="861" height="862" alt="Screenshot 2026-11-11 000625" src="https://github.com/user-attachments/assets/3fc490a2-86a3-4fde-a3a9-a3bdff045d90" />
**_Figure 6_**. Samtools Coverage tool and its parameters for calculating sequencing read coverage across the reference genome.

https://usegalaxy.org/u/krismelon/h/sex-dosage-acinonyx-jubatus

## Chromosome Groups

### Chromosome Groups and Reasons for Selection

| Chromosome Group | Accession Number | Reason for Selection |
|---|---|---|
| Autosomes (A) | NC_069380.1 | Selected to serve as a reference for comparing sequencing depth. |
| Autosomes (A) | NC_069382.1 | Chosen to help establish the average coverage of autosomes. |
| Autosomes (A) | NC_069381.1 | Used as part of the baseline for evaluating chromosome coverage. |
| Autosomes (A) | NC_069390.1 | Included to provide additional data for autosomal depth comparison. |
| Autosomes (A) | NC_069397.1 | Selected to strengthen the reference group for autosomal coverage. |
| X Chromosome | NC_069389.1 | Chosen to examine X-chromosome coverage in relation to the autosomal average. |

## R Studio Code

# Sex-chromosome dosage from whole-genome DNA read depth
# Input: Samtools coverage TSV + student chromosome labels
cov <- read.delim("chromosome_coverage.tsv", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
names(cov)[1] <- "rname"  # Accepts the '#rname' header
if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")

key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)
if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")
key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))
if (anyDuplicated(key$rname)) stop("Repeated rname in key file")
if (!all(key$group %in% c("A", "X", "Y", "Z", "W")))
  stop("Use only A, X, Y, Z or W in the group column")

x <- merge(key, cov[, c("rname", "meandepth")],
           by = "rname", all.x = TRUE, sort = FALSE)
x$meandepth <- suppressWarnings(as.numeric(x$meandepth))
if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")
if ("X" %in% x$group && "Z" %in% x$group)
  stop("Choose either an XY or a ZW species, not both")

if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
  target <- "X"; one_copy <- "XY-like"; two_copies <- "XX-like"
} else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
  target <- "Z"; one_copy <- "ZW-like"; two_copies <- "ZZ-like"
} else stop("Label exactly one X or exactly one Z chromosome")

A <- median(x$meandepth[x$group == "A"])
if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")
x$relative_depth <- x$meandepth / A
ratio <- x$relative_depth[x$group == target]
call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65) one_copy else "Ambiguous"
if (A < 1) call <- "Inconclusive: very low autosomal depth"

# Additional Y or W evidence can support, but not prove, a call
partner <- if (target == "X") "Y" else "W"
partner_ratio <- if (partner %in% x$group)
  x$relative_depth[x$group == partner][1] else NA_real_

# Export results and a simple, reproducible figure
dir.create("results", showWarnings = FALSE)
write.csv(x, "results/chromosome_relative_depth.csv", row.names = FALSE)
png("results/chromosome_depth.png", width = 1000, height = 580)
par(mar = c(8, 4, 3, 1))
barplot(x$relative_depth, names.arg = x$rname, las = 2,
        ylab = "Depth relative to median autosome",
        main = "Autosomal and sex-chromosome depth", cex.names = 0.75)
abline(h = c(0.5, 1.0), lty = 2)
dev.off()

summary <- c(paste("Median autosome depth:", round(A, 3)),
             paste(target, "/ A ratio:", round(ratio, 3)),
             paste(partner, "/ A ratio:", round(partner_ratio, 3)),
             paste("Dosage interpretation:", call),
             "Compare with BioSample sex and discuss limitations.")
writeLines(summary, "results/sex_inference.txt")
cat(paste(summary, collapse = "\n"), "\n")

<img width="1920" height="1080" alt="Screenshot (104)" src="https://github.com/user-attachments/assets/14e7191e-22c0-4aac-89cf-09b1dbb76180" />
**_Figure 7_**. Running the sex chromosome depth analysis script in RStudio and displaying the calculated sequencing depth results and dosage interpretation.

## Coverage Table

### Mean and Relative Sequencing Depth

| Accession Number | Chromosome Group | Mean Depth | Relative Depth |
|---|---|---:|---:|
| NC_069380.1 | Autosome (A) | To be calculated | To be calculated |
| NC_069382.1 | Autosome (A) | To be calculated | To be calculated |
| NC_069381.1 | Autosome (A) | To be calculated | To be calculated |
| NC_069390.1 | Autosome (A) | To be calculated | To be calculated |
| NC_069397.1 | Autosome (A) | To be calculated | To be calculated |
| NC_069389.1 | X Chromosome | To be calculated | To be calculated |

## Interpretation of the WGS Chromosome-Depth Analysis

<img width="1000" height="580" alt="chromosome_depth" src="https://github.com/user-attachments/assets/165adc76-c39b-4ac3-a520-fa37f22b8fad" />
**_Figure 8_**. Relative sequencing depth across the analyzed chromosomes of the cheetah (Acinonyx jubatus), with reference lines indicating the expected two-copy (1.0) and one-copy (0.5) coverage levels.

The bar plot shows generally consistent sequencing coverage across all analyzed chromosomes, with most values remaining close to the 1.0 baseline, which represents the expected coverage for chromosomes with two copies. Although the last bar shows a slightly lower sequencing depth compared to the others, it remains well above the 0.5 threshold. This indicates that none of the analyzed chromosomes displays a clear reduction to half the expected sequencing depth. The overall pattern suggests that the chromosomes have relatively similar coverage and that there is no obvious evidence of a chromosome with reduced copy representation. 

## Questions for the Report

**1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZWsystem?**

The species I chose is the cheetah (Acinonyx jubatus). Since cheetahs are mammals, they follow the XX/XY sex-determination system. Female cheetahs usually have two X chromosomes (XX), while males have one X and one Y chromosome (XY).

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?** 

The genome assembly accession I used is GCF_027475565.1, with the assembly name VMU_Ajub_asm_v1.0. Its BioSample accession is SAMN29432198. The WGS run accession is SRR22273180.

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

I selected five sequences as autosomes and one sequence as the X chromosome. The accession numbers are listed below:

| Chromosome Group | Accession Number |
|---|---|
| Autosome (A) | `NC_069380.1` |
| Autosome (A) | `NC_069382.1` |
| Autosome (A) | `NC_069381.1` |
| Autosome (A) | `NC_069390.1` |
| Autosome (A) | `NC_069397.1` |
| X Chromosome | `NC_069389.1` |

I identified these sequences using the NCBI genome assembly information and chromosome annotations. I also need to make sure that the accession numbers match the reference names used in my Galaxy coverage results.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

| **Measurement** | **Result** | **Interpretation** |
|---|---:|---|
| Median Autosomal Depth | 14.846 | Used as the baseline for comparing chromosome coverage. |
| X:A Ratio | 0.934 | Close to the expected value of 1.0 for an XX chromosome complement. |
| Y:A Ratio | NA | No Y-chromosome coverage result was available in the analysis. |
| Dosage Interpretation | XX-like | The result suggests a possible XX chromosome complement based on the X:A ratio. |

Caduff et al. (2024) studied how sequencing data can be used to identify possible sex-chromosome combinations. This is related to my study because I also used sequencing coverage to determine the possible sex-chromosome pattern of the cheetah. My X:A ratio was 0.934, which suggests an XX-like pattern. This result also agrees with the female sex recorded in the NCBI BioSample.However, since my analysis did not provide a Y:A ratio, I cannot completely confirm the animal's chromosome combination.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancyoruncertainty.** 

The NCBI BioSample record (SAMN29432198) identifies the cheetah as female. My sequencing coverage analysis also shows an XX-like pattern, with an X:A ratio of 0.934, which is close to the expected value of 1.0 for an XX chromosome complement. These findings are consistent with the recorded sex of the sample. However, the Y:A ratio was not available in my results, so the analysis cannot independently confirm the complete chromosome complement.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

A chromosome-level assembly alone is not enough because it only provides the reference genome and does not show the actual sequencing coverage of the sample. To determine the possible sex-chromosome complement, the WGS reads must first be aligned to the reference genome, and the sequencing depth of the sex chromosomes must be compared with the autosomes. Other factors, such as low sequencing coverage, mapping errors, repetitive sequences, and incomplete chromosome information, can also affect the results. Therefore, actual sequencing data are needed to make a more reliable interpretation.

## Limitations

This study has some limitations because it only used one cheetah WGS sample, so the results may not represent the entire species. The analysis also focused only on selected autosomes and the X chromosome, while the Y:A ratio was unavailable, making it difficult to check the Y-chromosome coverage. According to Sims et al. (2014), factors such as sequencing depth, read quality, and mapping accuracy can affect the results of genomic analysis. These factors may also affect the accuracy of the coverage values in this study. Although the X:A ratio of 0.934 suggests an XX-like chromosome pattern and agrees with the female sex recorded in the NCBI BioSample, the results cannot fully confirm the animal's chromosome combination because coverage analysis alone has its limitations.

## References

**RELATED STUDIES**

Caduff, M., Eckel, R., Leuenberger, C., & Wegmann, D. (2024). Accurate Bayesian inference of sex chromosome karyotypes and sex‐linked scaffolds from low‐depth sequencing data. Molecular Ecology Resources, 24(3), e13913. https://doi.org/10.1111/1755-0998.13913

Sims, D., Sudbery, I., Ilott, N. E., Heger, A., & Ponting, C. P. (2014). Sequencing depth and coverage: key considerations in genomic analyses. Nature Reviews Genetics, 15(2), 121–132. https://doi.org/10.1038/nrg3642

**NCBI REFERENCES**

1. National Center for Biotechnology Information. (n.d.). *Acinonyx jubatus: VMU_Ajub_asm_v1.0 (GCF_027475565.1)* [Genome assembly]. NCBI Assembly. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_027475565.1/

2. National Center for Biotechnology Information. (n.d.). *Acinonyx jubatus (SAMN29432198)* [BioSample]. NCBI BioSample. https://www.ncbi.nlm.nih.gov/biosample/SAMN29432198

3. National Center for Biotechnology Information. (n.d.). *Acinonyx jubatus (SRR22273180)* [Sequence read archive]. NCBI Sequence Read Archive. https://www.ncbi.nlm.nih.gov/sra/SRR22273180
