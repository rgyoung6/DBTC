# Dada-BLAST-Taxon Assign-Condense (DBTC)
DBTC is an R implementation of Dada2 and a BLAST approach to metabarcoding analysis.

# Description
molecular analysis of high-throughput metabarcode molecular sequence data fastq files and taxonomic assignment of unique reads (in fasta format files) using the Basic Local Alignment Search Tool (BLAST) and R to reduced obtained results (See additional documentation at [DBTC](https://github.com/rgyoung6/DBTC)). The Dada-BLAST-Taxon Assign-Condense ('DBTC') package contains the foundational ['DBTC'](https://github.com/rgyoung6/DBTC) functions run through the R command line (But also see ['DBTCShiny'](https://github.com/rgyoung6/DBTCShiny) and [DBTCShinyTutorial](https://github.com/rgyoung6/DBTCShinyTutorial)). The 'DBTC' functions have four main outcomes...

  - [Fastq](https://en.wikipedia.org/wiki/FASTQ_format) file processing using Dada in R
  - Using the Basic Local Alignment Search Tool ([BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology))), amplicon sequence variants ([ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) can be searched against local NCBI or custom libraries
  - Assign taxa to the unique reads using NCBI taxon database (obtain the database using [taxonomizr website](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html))
  - Condense the resulting ASV taxonomic assignment tables to unique taxa with the ability to combine datasets (using different sequence libraries for the same reads, or results from the same samples for different molecular regions) into a combined results table

**NOTE:** While the DBTC package has been built for the analysis of high-throughput sequencing results, the BLAST, taxonomic assignment, and taxonomic condense can be utilized with single specimen Sanger sequencing data.

Also see [DBTCShiny](https://github.com/rgyoung6/DBTCShiny) and the associated [DBTCShinyTutorial](https://github.com/rgyoung6/DBTCShinyTutorial).

# Table of Contents  
- [Installation](#installation)
- [Package Dependencies](#package-dependencies)
  * [External R Element Dependencies](#external-r-element-dependencies)
    * [NCBI BLAST+ local program to run BLAST on local databases](#ncbi-blast+-local-program-to-run-blast-on-local-databases)
    * [R package taxonomizr to establish a local NCBI taxonomy database](#r-package-taxonomizr-to-establish-a-local-ncbi-taxonomy-database)
    * [Establish a local NCBI prepared sequence database](#establish-a-local-ncbi-prepared-sequence-database)
    * [Create a custom sequence database to BLAST against](#create-a-custom-sequence-database-to-blast-against)
    * [Create a local NCBI taxonomy database to assign taxonomic identifications to BLAST results](#create-a-local-ncbi-taxonomy-database-to-assign-taxonomic-identifications-to-blast-results)
  * [R Packages Dependencies](#r-packages-dependencies)
     * [Bioconductor ShortRead and Dada2 packages](#bioconductor-shortread-and-dada2-packages)
     * [CRAN](#cran)
- [Function Descriptions](#function-descriptions)
- [Naming Convention Rules](#naming-convention-rules)
- [Package Function Details](#package-function-details)
  * [Dada Implement - dada_implement()](#dada-implement)
  * [Combine Dada Output - combine_dada_output()](#combine-dada-output)
  * [Make BLAST DB - make_BLAST_DB()](#make-blast-db)
  * [Sequence BLAST - seq_BLAST()](#sequence-blast)
  * [Taxon Assignment - taxon_assign()](#taxon-assignment)
  * [Combine Assignment Output - combine_assign_output()](#combine-assignment-output)
  * [Reduce Taxa - reduce_taxa()](#reduce-taxa)
  * [Combine Reduced Output - combine_reduced_output()](#combine-reduced-output)
- [Citation](#citation)

# Installation 

DBTC can be installed three ways.

## 1. Install from CRAN

install.packages('DBTCShiny’)

## 2. Install via GitHub
Run the following commands in your R terminal...<br/>
```
if(!require(devtools)) install.packages('devtools')
library('devtools')
devtools::install_github('rgyoung6/DBTC')
library('DBTC')
```
**Note:** the first command to install the "devtools" may not be necessary if already installed.<br/>

## 3. Install through download from GitHub
Navigate to the [DBTC](https://github.com/rgyoung6/DBTC) GitHub page. Download the files associated with this page to your local computer and place them somewhere in the main file folder named DBTC. Then run the following command pointing to that location on your local computer by replacing the HERE with the path in the below command...<br/>
```
library("DBTC", lib.loc="HERE")
```

([Back to Top](#table-of-contents))

***

# Package Dependencies
 
## External R Element Dependencies

### NCBI BLAST+ local program to run BLAST on local databases

Follow the instructions on the NCBI [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) executables page to obtain a local version of the BLAST tools. The list of the latest installation files can be found [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Note: It is best to download and install the most recent versions of the blast+ suite to your computer and place the programs in your computers [path](https://en.wikipedia.org/wiki/PATH_(variable)) so you can access the program from any folder. However, the program files for both blastn and makeblastdb have been included in the [DBTCShinyTutorial](https://github.com/rgyoung6/DBTCShinyTutorial) GitHub page for ease of use (please note that these may not be the most recent versions).

### R package taxonomizr to establish a local NCBI taxonomy database
The R package [taxonomizr website](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html) is used to establish a NCBI taxaID database (NOTE: this package is also required when using the taxon assignment elements in the DBTC pipeline).
```
install.packages('taxonomizr')
library('taxonomizr')
```

### Establish a local NCBI prepared sequence database
NCBI BLASTable databases can be established through two methods.

1. Download your desired preformatted NCBI database by using the 'update_blastdb.pl' (found in the NCBI BLAST+ local install folder). NOTE: Perl programming language needs to be installed on your local machine. Instructions can be found at [Get NCBI BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

2. You can download your desired preformatted NCBI database manually with instructions at [BLAST FTP Site](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory) and a list of the available databases at [Index of /blast/db](https://ftp.ncbi.nlm.nih.gov/blast/db/). 

### Create a custom sequence database to BLAST against
In addition to the NCBI resources, DBTC can also use custom databases. To establish these databases you will require a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the desired records with [MACER](https://github.com/rgyoung6/MACER) formatted headers. The [MACER](https://github.com/rgyoung6/MACER) R package and instructions can be found at either of the two locations:

[MACER CRAN](https://CRAN.R-project.org/package=MACER)

[MACER GitHub](https://github.com/rgyoung6/MACER) (will have the most recent version and development versions)

### Create a local NCBI taxonomy database to assign taxonomic identifications to BLAST results
In the 'Preparation' section of the [taxonomizr website](https://CRAN.R-project.org/package=taxonomizr), use the instructions and the prepareDatabase('accessionTaxa.sql', getAccessions = FALSE) taxonomizr command to establish a local taxonomy database.
```
prepareDatabase('accessionTaxa.sql', getAccessions = FALSE)
```
Note: Along with the accessionTaxa.sql two other files nodes.dmp and names.dmp files are downloaded. These two files are not necessary for downstream analyses and can be deleted.

([Back to Top](#table-of-contents))

***

## R Packages Dependencies

### Bioconductor ShortRead and Dada2 packages

The [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) package is required to run elements of the DBTC pipeline and can be obtained through Bioconductor.

The [dada2](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) package is main package to process the raw fastq files and can be obtained from Bioconductor. There is also a [dada2](https://benjjneb.github.io/dada2/) GitHub resource. 

```
if (!require('BiocManager', quietly = TRUE))
install.packages('BiocManager')
BiocManager::install('ShortRead')
BiocManager::install('dada2')
library(dada2)
```

### CRAN
Each of below CRAN packages and their dependencies are required for the DBTC package.


```
install.packages(c('ggplot2',
                   'parallel',
                   'pbapply',
                   'plyr',
                   'stats',
                   'taxonomizr',
                   'utils'))
                   
library(c('ggplot2',
          'parallel',
          'pbapply',
          'plyr',
          'stats',
          'taxonomizr',
          'utils'))
```

([Back to Top](#table-of-contents))

***

# Run DBTC

After installation of the DBTC and all of its dependencies you need to load the package using the following commands...

```
library(DBTC)
```

# Function Descriptions

## dada_implement()
The dada_implement() function takes [fastq](https://en.wikipedia.org/wiki/FASTQ_format) files as input, analyses them and produces amplicon sequence variant ([ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) files. This function requires a main directory containing folder(s) representing sequencing runs which in-turn contains [fastq](https://en.wikipedia.org/wiki/FASTQ_format) files (the location of one of the [fastq](https://en.wikipedia.org/wiki/FASTQ_format) files in one of the sequencing run folders is used as an input argument). **A run is a group of results processed at the same time on the same machine representing the same molecular methods.** All sequencing folders in the main directory need to represent data from sequencing runs that have used the same primers and protocols. Output from this function includes all processing files and final main output files in the form of [fasta](https://en.wikipedia.org/wiki/FASTA_format) files and [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables.

## combine_dada_output()
DBTC dada_implement() uses [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) output files ('YYYY_MM_DD_HH_MM_UserInputRunName_Merge' and/or 'YYYY_MM_DD_HH_MM_UserInputRunName_MergeFwdRev') and combines them into a single [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table and creates an accompanying [fasta](https://en.wikipedia.org/wiki/FASTA_format) file. This function also produces a file containing the processing information for the function. The main input argument for this function is the location of a file in a folder containing all [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables wanting to be combined. Output files are generated with the naming convention 'YYYY_MM_DD_HH_MM_combinedDada'.

## make_BLAST_DB()
This function takes a [fasta](https://en.wikipedia.org/wiki/FASTA_format) file with headers in the MACER format (Young et al. 2021) and establishes a database upon which a BLAST search can be completed. However, if a NCBI sequence database is desired, it is advisable to use, where applicable, NCBI preformatted databases and skip the make_BLAST_DB() function (https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory). The outcome of the function is a folder with a BLASTable NCBI formatted sequence database.

The MACER [fasta](https://en.wikipedia.org/wiki/FASTA_format) header format

```>GenBankAccessionOrBOLDID|GenBankAccession|Genus|species|UniqueID|Marker```

## seq_BLAST()
[Fasta](https://en.wikipedia.org/wiki/FASTA_format) file(s) are used as input along with a user selected NCBI formatted database upon which query sequences will be searched using BLAST. The outcome of the function are two files, a BLAST run file and a single file containing all of the BLAST results in tab delimited format. There are no headers in the BLAST results file but the columns (in order left to right) are: query sequence ID, search sequence ID, search taxonomic ID, query to sequence coverage, percent identity, search scientific name, search common name, query start, query end, search start, search end, e-value.

## taxon_assign()
This function takes a BLAST result file and associated [fasta](https://en.wikipedia.org/wiki/FASTA_format) files (either on their own or with accompanying [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files generated from the dada_implement() function) and collapses the multiple BLAST results into as single result for each query sequence. When an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table is present the taxonomic results will be combined with the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table.

## combine_assign_output()
The combine_assign_output() function takes a file selection and then uses all DBTC [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) taxon assign files ('_taxaAssign_YYYY_MM_DD_HHMM.tsv') in a selected directory and combines them into a single output 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' file. The files being combined should represent different samples but representing data that have all come from analysis using the same molecular methods, the same analysis arguments, and the same molecular sequence databases.

## reduce_taxa()
To reduce [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) results to unique taxa per file the reduce_taxa() function takes a file selection and then uses all '_taxaAssign_YYYY_MM_DD_HHMM.tsv' and/or 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' files in that directory. This function then reduces all [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) with the same taxonomic assignment into a single result and places these results in a '_taxaReduced_YYYY_MM_DD_HHMM.tsv' file for each of the target files in the directory.

## combine_reduced_output()
This function takes a file selection and then uses all '_taxaReduced_YYYY_MM_DD_HHMM.tsv' files in that directory and combines them into a single 'YYYY_MM_DD_HHMM_CombineTaxaReduced.txt' taxa table file with presence absence results. The files being combined should represent the same biological samples but with different molecular marker information. The output [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) can include read numbers or be reduced to presence absence results.

([Back to Top](#table-of-contents))

***

# Naming convention rules

WARNING - NO WHITESPACE!

When running DBTC functions the paths for the files selected cannot have white space! File folder locations should be as short as possible (close to the [root](https://en.wikipedia.org/wiki/Root_directory) directory) as some functions do not process long naming conventions. 

Also, special characters should be avoided (including question mark, number sign, exclamation mark). It is recommended that dashes be used for separations in naming conventions while retaining underscores for use as information delimiters (this is how DBTC functions use underscore). 

There are several key character strings used in the DBTC pipeline, the presence of these strings in file or folder names will cause errors when running DBTC functions. 
The following strings are those used in DBTC and should not be used in file or folder naming:
  - _BLAST
  - _combinedDada
  - _taxaAssign
  - _taxaAssignCombined
  - _taxaReduced
  - _CombineTaxaReduced

([Back to Top](#table-of-contents))

***

# Package Function Details

## Dada Implement
dada_implement() - Process metabarcode raw fastq files by [runs](#runs) using Dada2 (Note: molecular markers are independently analysed and combined at the end of the analysis pipeline using the [combine_reduced_output()](#combine-reduced-output) function).

### Input 
Two files are required as input for the dada_implement() function. The first are the fastq files in the appropriate folder structure (see below) and the second is a file containing the primers used for the amplification of the sequence reads (tab separated file, see [DBTCShinyTutorial](https://github.com/rgyoung6/DBTCShinyTutorial) for file examples). To select all of the fastq files simply select one file in one of the [Run](#runs) directories to point to the desired data files, as long as the file and folder structure is correct.

**Fastq File Folder Structure**

```
            Parent Directory
                  |
                  |
          -----------------
          |               |
          |               |
    Run1 Directory     Run2 Directory
    -Fastq             -Fastq
    -Fastq             -Fastq
    ...                ...
```

**Format of the primer file**

| Forward        | Reverse           | 
| :------------- |:-------------| 
| AGTGTGTAGTGATTG      | CGCATCGCTCAGACTGACTGC | 
| GAGCCCTCGATCGCT      | GGTCGATAGCTACGCGCGCATACGACT      |  
|  | GGTTCACATCGCATTCAT      |   


### Arguments
- <strong>runFolderLoc -</strong> Select a file in the one of the [run](#runs) folders with the fastq files of interest (Default NULL).
- <strong>primerFile -</strong> Select a file with the primers for this analysis (Default NULL).
- <strong>fwdIdent -</strong> Forward identifier naming string (Default '_R1_001').
- <strong>revIdent -</strong> Reverse identifier naming string (Default '_R2_001').
- <strong>bidirectional -</strong> Selection to process paired forward and reverse sequence for analysis (Default TRUE).
- <strong>unidirectional -</strong> Selection to process files independently (Default FALSE).
- <strong>printQualityPdf -</strong> Selection to process save image files showing quality metrics (Default TRUE).
- <strong>maxPrimeMis -</strong>Maximum number of mismatches allowed when pattern matching trimming the primers from the ends of the reads for the ShortRead trimLRPatterns() function (Default 2).
- <strong>fwdTrimLen -</strong> Select a forward trim length for the Dada filterAndTrim() function (Default 0).
- <strong>revTrimLen -</strong> Select a reverse trim length for the Dada filterAndTrim() function (Default 0).
- <strong>maxEEVal -</strong> Maximum number of expected errors allowed in a read for the Dada filterAndTrim() function (Default 2).
- <strong>truncQValue -</strong> Truncation value use to trim ends of reads, nucleotides with quality values less than this value will be used to trim the remainder of the reads for the Dada filterAndTrim() function (Default 2).
- <strong>truncLenValueF -</strong> Dada forward length trim value for the Dada filterAndTrim() function. This function is set to 0 when the pattern matching trim function is enabled (Default 0).
- <strong>truncLenValueR -</strong> Dada reverse length trim value for the Dada filterAndTrim() function. This function is set to 0 when the pattern matching trim function is enabled (Default 0).
- <strong>error -</strong> Percent of fastq files used to assess error rates for the Dada learnErrors() function (Default 0.1).
- <strong>nbases -</strong> The total number of bases used to assess errors for the Dada learnErrors() function (Default 1e80) NOTE: this value is set very high to get all nucleotides in the error present file subset. If the error is to be assessed using total reads and not specific fastq files then set the error to 1 and set this value to the desired number of reads.
- <strong>maxMismatchValue -</strong> Maximum number of mismatches allowed when merging two reads for the Dada mergePairs() function (Default 2).
- <strong>minOverlapValue -</strong> Minimum number of overlapping nucleotides for the forward and reverse reads for the Dada mergePairs() function (Default 12).
- <strong>trimOverhang -</strong> Trim merged reads past the start of the complimentary primer regions for the Dada mergePairs() function (Default FALSE).
- <strong>minFinalSeqLen -</strong> The minimum final desired length of the read (Default 100).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output

The output files from this function appear in four folders. See the below diagram for the saved file structure.
```            
                              Parent Directory
                                    |
                                    |
                      -----------------------------
                      |                           |
                      |                           |
                Run1 Directory                 Run2 Directory
                - Fastq                        - Fastq
                - Fastq                        - Fastq
                ...                            ...
                      |
                      |
            --------------------------------------------------------   
            |                |             |                        |          
YYYY_MM_DD_HHMM_Run1_A_Qual  |             |                        |
  -revQual.pdf               |  YYYY_MM_DD_HHMM_Run1_C_FiltQual     |
  ...                        |    -filtFwdQual.pdf                  |
                             |    -filtRevQual.pdf                  |
                             |    ...                               |
                             |                                      |
                             |                    YYYY_MM_DD_HHMM_Run1_D_Output
                             |                      -dadaSummary.txt
                             |                      -dadaSummaryTable.tsv
                             |                      -ErrorForward.pdf
                             |                      -ErrorReverse.pdf
                             |                      -MergeFwdRev.tsv
                  YYYY_MM_DD_HHMM_Run1_B_Filt       -MergeFwdRev.fas
                    -fwdFilt.fastq                  -Merge.tsv
                    -revFilt.fastq                  -Merge.fas
                    ...                             -TotalTable.tsv
                    Primer_Trim
                      -primeTrim.fastq
                      ...

```

### Interpretation
Quality pdf's in the A_Qual folder represent the quality metrics for the raw [Fastq](https://en.wikipedia.org/wiki/FASTQ_format) files (This folder may not be present if printQualityPdf is set to FALSE).
Files in the B_Filt folder represent the trimming (in the Primer_Trim folder) and trimmed and cleaned in the larger folder.
Quality pdf's in the C_FiltQual folder represent the quality metrics for the trimmed and cleaned [Fastq](https://en.wikipedia.org/wiki/FASTQ_format) files (This folder may not be present if printQualityPdf is set to FALSE).
There are numerous files in the D_Output folder. These include:

- dadaSummary.txt file which provides all of the information on the running of the dada_implement() function.
- dadaSummaryTable.tsv contains a table with summary information for the processing of the samples in the [run](#runs).
- ErrorForward.pdf and ErrorReverse.pdf provide visualizations on the assessed sequencing error for the sequencing [run](#runs).
- MergeFwdRev.tsv is the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with the data from the sequencing [run](#runs) and the MergeFwdRev.fas is a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the reads from the samples. The MergeFwdRev files include reads that were able to be merged, as well as reads that were not able to be merged.
  NOTE: The merged, forward, and reverse reads are obtained in parallel analyses and combined into a single file so MergeFwdRev files will represent triplicate molecular processing results. These files are present to see if there are reads that are not represented or poorly represented across merged and unidirectional results, perhaps indicating issues with one of the primers.
- Merge.tsv is the [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with the read data able to be merged from the sequencing [run](#runs). The companion Merge.fas is a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file with the reads from the samples in [Fasta](https://en.wikipedia.org/wiki/FASTA_format) format.
- TotalTable.tsv is an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table with all of the merged, forward, and reverse results as well as the retained results for the merged reads removed due to being suspected chimeric combinations.

### Dependencies
- dada2
- ShortRead readFastq()
- ShortRead writeFastq()
  
([Back to Top](#table-of-contents))
***

## Combine Dada Output
combine_dada_output() - Combine Dada2 [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables (YYYY_MM_DD_HHMM_FileName_MergeFwdRev.tsv OR YYYY_MM_DD_HHMM_FileName_Merge.tsv) into a single [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) output file.

### Input 
Two or more files to be combined are required as input for this function. These files need to be [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files as outputted from the [dada_inplement()](#dada-implement) and can include Merge, MergeFwdRev, or TotalTable.tsv files.

### Arguments
- <strong>fileLoc -</strong> Select a file in the file folder with [dada_inplement()](#dada-implement) results you would like to combine (YYYY_MM_DD_HHMM_FileName_MergeFwdRev OR YYYY_MM_DD_HHMM_FileName_Merge both .tsv and .fas files (Default NULL).
- <strong>minLen -</strong> The minimum final desired length of the read (Default 100).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
The output from this function includes three files.
  1. YYYY_MM_DD_HHMM_combinedDada.tsv - combined [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table
  2. YYYY_MM_DD_HHMM_combinedDada.fas - combined [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file
  3. YYYY_MM_DD_HHMM_combinedDada.txt - Summary file from the combine_dada_output()

### Interpretation
Outputted data files will come in the same [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) table format as the output [dada_inplement()](#dada-implement) [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) files. 

### Dependencies
- Base R
  
([Back to Top](#table-of-contents))
***

## Make BLAST DB
make_BLAST_DB() - Create a local database to BLAST against.

### Input 
This function takes a [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file (in [MACER](https://github.com/rgyoung6/MACER) format) and establishes a database upon which a BLAST search can be completed. The outcome of the function is a folder with an NCBI database.
- The MACER [Fasta](https://en.wikipedia.org/wiki/FASTA_format) header format - ```>UniqueID|OtherInformation|Genus|species|OtherInformation|Marker```
- An example of the header format output from the [MACER](https://github.com/rgyoung6/MACER) program is ```>GenBankAccessionOrBOLDID|GenBankAccession|Genus|species|UniqueID|Marker```
  
### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files will be used to construct a BLASTable database (Default = NULL).
- <strong>makeblastdbPath -</strong> The local path for the blast+ makeblastdbPath program (Default 'makeblastdb').
- <strong>taxaDBLoc  -</strong> The location of the NCBI [accessionTaxa.sql](#create-a-local-ncbi-taxonomy-database-to-assign-taxonomic-identifications-to-blast-results) taxonomic data base (Default NULL).
- <strong>dbName -</strong> A short 6-8 alpha character name used when building a database (Default NULL).
- <strong>minLen -</strong> The minimum sequence length used to construct the BLAST database (Default 100).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
The output from this function includes a folder with the BLAST database named according to the submitted dbName.

### Interpretation
The constructed database can then be used with the [seq_BLAST()](#sequence-blast) function.

### Dependencies
- [taxonomizr()](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html)
  
([Back to Top](#table-of-contents))

***

## Sequence BLAST
seq_BLAST() - Search [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files of unknown sequences against a BLAST formatted database.

### Input 
Provide a location for the BLAST database you would like to use by selecting a file in the target directory (This can be a built database using the [make_BLAST_DB()](#make-blast-db) function or a preformatted [NCBI BLAST database](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory)). Then provide the location of the query sequence files by indicating a file in a directory that contains the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files. Provide the path for the blast+ blastn program. Finally, provide the minimum query sequence length to BLAST (Default = 100), the maximum depth of the BLAST returned results (Default = 200), and finally the number of cores to process the function (default = 1, Windows implementation can only use a single core and will default to this value when running on Windows).

### Arguments
- <strong>databasePath -</strong> The location of a file in a directory where the desired BLAST database is located (Default NULL).
- <strong>querySeqPath -</strong> The location of a file in a directory containing all of the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files wishing to be BLASTed (Default NULL).
- <strong>blastnPath -</strong> The location of the NCBI blast+ blastn program (Default 'blastn').
- <strong>minLen -</strong> The minimum length of the sequences that will be BLASTed (Default 100).
- <strong>BLASTResults -</strong> The number of returned results, or the depth of the reported results, saved from the BLAST (Default 200).
- <strong>numCores -</strong> The number of cores used to run the function (Default 1, Windows systems can only use a single core).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
Two files are produced from this function, a BLAST run file and a BLAST results file for each of the [Fasta](https://en.wikipedia.org/wiki/FASTA_format) files in the target directory.

### Interpretation
The BLAST run file contains the command used to run the BLAST search. The BLAST results file includes all results in a tab delimited .tsv file format with the columns qseqid, sseqid, staxid, qcovs, pident, ssciname, scomname, qstart, qend, sstart, send, evalue.

### Dependencies
- Base R
  
([Back to Top](#table-of-contents))

***

## Taxon Assignment
taxon_assign() - Using BLAST results to construct a table with taxonomic assignments for each unique sequence.

### Input 
This function requires a BLAST output file and an associated [Fasta](https://en.wikipedia.org/wiki/FASTA_format) file. In addition, if present an [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) file will also be used to combine with the taxonomic results. The function will take the BLAST results and reduce the taxonomic assignment to a single result for each read. 

### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all of the paired [Fasta](https://en.wikipedia.org/wiki/FASTA_format) and BLAST (and potentially [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)) files are located (Default NULL).
- <strong>taxaDBLoc -</strong> The location of the NCBI [accessionTaxa.sql](#create-a-local-ncbi-taxonomy-database-to-assign-taxonomic-identifications-to-blast-results) taxonomic data base (Default NULL).
- <strong>numCores -</strong> The number of cores used to run the function (Default 1, Windows systems can only use a single core).
- <strong>coverage -</strong> The percent coverage used for taxonomic assignment for the above threshold results (Default 95).
- <strong>ident -</strong> The percent identity used for the taxonomic assignment for above threshold results (Default 95).
- <strong>propThres -</strong> The proportional threshold flags the final result based on the preponderance of the data. So if the threshold is set to 0.95, results will be flagged if the taxa directly below the assigned taxa has fewer than 0.95 percent of the records causing the upward taxonomic placement (Default 0.95).
- <strong>coverReportThresh -</strong> The percent coverage threshold used for reporting flags below this threshold (Default 95).
- <strong>identReportThresh -</strong> The percent identity threshold used for reporting flags below this threshold (Default 95).
- <strong>includeAllDada  -</strong> When paired Dada [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) tables are present, when set to FALSE, this will exclude records without taxonomic assignment (Default TRUE).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
A single taxonomic assignment file is created for each BLAST output file with the naming convention having the string '_taxaAssign_YYYY_MM_DD_HHMM.tsv'. In addition, there is a run file that is also generated which contains the run details with the naming convention 'YYYY_MM_DD_HHMM_taxaAssign.txt'.

### Interpretation
The number of returned BLAST results is dictated by the [seq_BLAST()](#sequence-blast) BLASTResults argument. The taxon_assign() function takes into account all returned BLAST results for each read. At each taxonomic level assignments have quality metrics in parentheses after the name. These values ("Num_Rec", "Coverage", "Identity", "Max_eVal") represent the number of records with this taxonomic placement, the minimum coverage and identity, and the maximum eValue for the reported taxa.

Column headers for the resulting taxonomic assignments include...

uniqueID, superkingdom, phylum, class, order, family, genus, species, Top_BLAST, Lowest_Single_Ran, Lowest_Single_Taxa, Lowest_Single_Rank_Above_Thres, Lowest_Single_Taxa_Above_Thres, Final_Common_Names, Final_Rank, Final_Taxa, Final_Rank_Taxa_Thres, Result_Code, Sequence, Length, Results, followed by columns of samples containing [ASV](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) values

There are three columns that deserve special explanation. 

- The Final_Rank_Taxa_Thres column contains the threshold values (Coverage, Identitiy) applied to the final rank and taxonomic values for the associated records.  
- The Results column contains the reference to the record if it was from a merged, forward or reverse analysis result.
- The Result_Code column contains flags placed on the results to better understand the quality of the resulting taxonomic assignments. Below is a list of codes:

  - SFAT(coverage, ident): Saturated filtered taxa above threshold
  - SANF(coverage, ident): Saturated non-filtered
  - BIRT(identReportThresh): Final taxa result is below the identity reporting threshold
  - BCRT(coverReportThresh): Final taxa result is below the nucleotide coverage reporting threshold
  - TBAT(propThres): Taxa Below Assigned Taxa Threshold

Note: Records with BIRT and BCRT flags should be highly scrutinized. TBAT are also concerning in that they may represent a less specific taxonomic placement due to the moving of the result to a higher taxonomic placement. The taxonomic rank directly below the final reported rank should be reviewed as there may be potential to adjust the final taxonomic assignment. SANF results should be explored and the size of the database, the trust placed in the records in the database, and the depth of the BLAST results should be considered when assessing records with this flag. Records with SFAT are among the least concerning as the BLAST results were saturated but this taxonomic assignment saturation occurred above your set quality coverage and identity threshold. Concerns with records with this result could be that the depth of the BLAST analysis was not low enough for very large databases, or that the database is not complete (taxonomic breadth) when using smaller databases.

### Dependencies
- [taxonomizr()](https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html)
- pbapply()
- parallel()
  
([Back to Top](#table-of-contents))

***

## Combine Assignment Output
combine_assign_output() - Using results from the [taxon_assign()](#taxon-assignment) function, combines all files with the string '_taxaAssign_YYYY_MM_DD_HHMM.tsv' in to a single .tsv file.

### Input 
Select a file in a folder with the taxa assigned files you would like to combine (extension '_taxaAssign_YYYY_MM_DD_HHMM.tsv'). NOTE: all '_taxaAssign_YYYY_MM_DD_HHMM.tsv' files in the folder location should originate from the same dada output file but have outputs from different BLAST sequence libraries and therefore contain the same [ASV's](https://en.wikipedia.org/wiki/Amplicon_sequence_variant).

### Arguments
- <strong>fileLoc -</strong> The location of a file in a directory where all of the '_taxaAssign_YYYY_MM_DD_HHMM.tsv' files are located.
- <strong>numCores -</strong> The number of cores used to run the function (Default 1, Windows systems can only use a single core).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
This function produces a 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' and a 'YYYY_MM_DD_HHMM_taxaAssignCombined.txt' file in the selected target directory.

### Interpretation
The interpretation of the output file for the combine_assign_output() 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' files is the same as the is the same as the [taxon_assign()](#taxon-assignment) '_taxaAssign_YYYY_MM_DD_HHMM.tsv' files.

### Dependencies
- pbapply()
- parallel()
  
([Back to Top](#table-of-contents))

***

## Reduce Taxa
reduce_taxa() - Using results from [taxon_assign()](#taxon-assignment) and/or [Combine Assignment Output](#combine-assignment-output) this function combines all reads with the same taxonomic assignment into a single result for each of the files submitted.

### Input 
This function requires a file in a directory where all '_taxaAssign_YYYY_MM_DD_HHMM.tsv' and/or 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' files in that directory will be combined. All records with the same taxonomic result will be combined. The BLAST values in parentheses ("Num_Rec", "Coverage", "Identity", "Max_eVal") are combine by the mean number of records, the mean of the minimum coverage and identity values, and the mean of the maximum eValues.

### Arguments
- <strong>fileLoc -</strong>  The location of a file in a directory where all of the '_taxaAssign_YYYY_MM_DD_HHMM.tsv' and/or 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' files are located.
- <strong>numCores -</strong> The number of cores used to run the function (Default 1, Windows systems can only use a single core).
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
This function produces a '_taxaReduced_YYYY_MM_DD_HHMM.tsv' file for every '_taxaAssign_YYYY_MM_DD_HHMM.tsv' or 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' present in the target directory.

### Interpretation
Reduced taxonomic assignment files have fewer columns in the main taxa_reduced.tsv file than the '_taxaAssign_YYYY_MM_DD_HHMM.tsv' or 'YYYY_MM_DD_HHMM_taxaAssignCombined.tsv' files as columns are collapsed. In addition, the values in the taxonomic columns in parentheses represent the average values across all of the results with the same taxonomic assignment (see [taxon_assign()](#taxon-assignment) interpretation).

The columns include, superkingdom, phylum, class, order, family, genus, species, Top_BLAST, Final_Common_Names, Final_Rank, Final_Taxa, Result_Code, RepSequence, Number_ASV, Average_ASV_Length, Number_Occurrences, Average_ASV_Per_Sample, Median_ASV_Per_Sample, Results. NOTE: The representative sequence (RepSequence) is the longest sequence for each of the collapsed taxa assigned. 

### Dependencies
- pbapply()
- parallel()

([Back to Top](#table-of-contents))

***

## Combine Reduced Output
combine_reduced_output() - This function takes '_taxaReduced_YYYY_MM_DD_HHMM.tsv' files generated from the same biological samples but representing different amplified molecular markers and collapses these data into a single file. The outcome of this process results in a presence absence matrix for all taxa and markers.

### Input 
Select a file in a folder with '_taxaReduced_YYYY_MM_DD_HHMM.tsv' files representing data for the same biological samples but representing different amplified molecular markers.

### Arguments
- <strong>fileLoc -</strong>  The location of a file in a directory where all of the '_taxaReduced_YYYY_MM_DD_HHMM.tsv' files are located.
- <strong>presenceAbsence -</strong>  A TRUE or FALSE value used to indicate if the read values should be replaced with presence/absence (1/0) data. This change is necessary when combining taxa for the same named samples across molecular markers (TRUE) but is not necessary when combining results for taxa with all unique sample names (FALSE). 
- <strong>verbose -</strong> If set to TRUE then there will be output to the R console, if FALSE then this reporting data is suppressed (Default TRUE).

### Output
Two files, a 'YYYY_MM_DD_HHMM_CombineTaxaReduced.tsv' result file and a 'YYYY_MM_DD_HHMM_CombineTaxaReduced.txt' run summary file are generated from this function. The result file contains presence/absence data in a matrix that associates the data with samples, taxa, and molecular marker. The column headers in the results file includes the following, superkingdom, phylum, class, order, family, genus, species, markers(n number of columns), samples (n number).

### Interpretation
The interpretation of the output file is the same as the [taxon_assign()](#taxon-assignment) '_taxaAssign_YYYY_MM_DD_HHMM.tsv' files.

### Dependencies
- plyr() rbind.fill

([Back to Top](#table-of-contents))

***

# Citation
Young RG, et al., Hanner RH (2024) A Scalable, Open Source, Cross Platform, MetaBarcode Analysis Method using Dada2 and BLAST. Biodiversity Data Journal (In progress)

([Back to Top](#table-of-contents))
