# Written by Rob Young at the University of Guelph in Ontario Canada, April, 2024
# ******************************************************************************
# Roxygen2 Documentation:
#' @export
#'
#' @title Make a BLAST Database
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a fasta file (in MACER format) and
#' establishes a database upon which a BLAST search can be completed.
#'
#' @details
#' The user inputs the location of a file in a directory that contains a properly formatted
#' fasta file which can be used to construct a BLASTable database. The
#' NCBI blast+ program, makeblastdb and the NCBI taxonomic database (accessionTaxa.sql) are required to
#' run this script (see readme instructions for details).
#'
#' The examples are present to display the syntax for the function.
#' These examples are not run because there are files required to run the functions,
#' in some cases multiple files are necessary and some of these are quite large. To
#' get specific examples please see https://github.com/rgyoung6/DBTCShinyTutorial/blob/main/README.md
#'
#' @examples
#' \dontrun{
#' make_BLAST_DB()
#' make_BLAST_DB(fileLoc = NULL, makeblastdbPath = "makeblastdb", taxaDBLoc = NULL,
#' inputFormat = NULL, dbName = NULL, minLen = 100)
#' }
#'
#' @param fileLoc The location of a file in a directory where all fasta files
#' will be used to construct a BLASTable database (Default NULL).
#' @param makeblastdbPath The local path for the blast+ makeblastdbPath program (Default 'makeblastdb').
#' @param taxaDBLoc The location of the NCBI taxonomic data base (Default NULL; for accessionTaxa.sql
#' see the main DBTC page for details).
#' @param dbName A short 6-8 alpha character name used when building a database (Default NULL).
#' @param minLen The minimum sequence length used to construct the BLAST database (Default 100).
#' @param verbose If set to TRUE then there will be output to the R console, if
#' FALSE then this reporting data is suppressed (Default TRUE).
#'
#' @returns
#' The output from this function includes a folder with the BLAST database named
#' according to the submitted dbName
#'
#' @references
#' <https://github.com/rgyoung6/DBTC>
#' Young, R. G., Hanner, R. H. (Submitted October 2023). Metabarcoding analysis
#' using Dada-BLAST-Taxon Assign-Condense Shiny Application (DBTCShiny). Biodiversity Data Journal.
#'
#' @note
#' WARNING - NO WHITESPACE!
#'
#' When running DBTC functions the paths for the files selected cannot have white
#' space! File folder locations should be as short as possible (close to the root
#' as some functions do not process long naming conventions.
#'
#' Also, special characters should be avoided (including question mark, number
#' sign, exclamation mark). It is recommended that dashes be used for separations
#' in naming conventions while retaining underscores for use as information
#' delimiters (this is how DBTC functions use underscore).
#'
#' There are several key character strings used in the DBTC pipeline, the presence
#' of these strings in file or folder names will cause errors when running DBTC functions.
#'
#' The following strings are those used in DBTC and should not be used in file or folder naming:
#' - _BLAST
#' - _combinedDada
#' - _taxaAssign
#' - _taxaAssignCombined
#' - _taxaReduced
#' - _CombineTaxaReduced
#'
#' @seealso
#' dada_implement()
#' combine_dada_output()
#' seq_BLAST()
#' taxon_assign()
#' combine_assign_output()
#' reduce_taxa()
#' combine_reduced_output()

##################################### make_BLAST_DB FUNCTION ##############################################################
make_BLAST_DB <- function(fileLoc = NULL, makeblastdbPath = "makeblastdb", taxaDBLoc = NULL, dbName = NULL, minLen = 100, verbose = TRUE){

  #If there are issues and I need to audit the script make this 1
  auditScript=0

  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  #Printing the start time
  if(verbose){
    print(paste0("Start time...", Sys.time()))
  }
  startTime <- paste0("Start time...", Sys.time())
  dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"))

  if(is.null(dbName)){
    if(verbose){
      print("********************************************************************************")
      print("Please rerun the function and provide a name for the outgoing database. This ")
      print("should be descriptive but short (less than 20 characters) with no special ")
      print("characters.")
      print(paste0("Current database name (dbName) is: ", dbName))
      print("********************************************************************************")
    }
  } else {

    #load in the location of the accessionTaxa.sql
    if (is.null(taxaDBLoc)){
      if(verbose){
        print(paste0("Select the NCBI taxonomic data file (accessionTaxa.sql)."))
      }
      taxaDBLoc <- file.choose()
    }

    #load in the fasta file of interest
    if (is.null(fileLoc)){
      if(verbose){
        print(paste0("Select the fasta file for use in constructing the custom database."))
      }
      fileLoc <- file.choose()
    }

    #load in the location of the makeBLASTdb file
    if (is.null(makeblastdbPath)){
      if(verbose){
        print(paste0("Select the location of the makeBLASTdb program."))
      }
      makeblastdbPath <- file.choose()
    }

    #Audit line
    if(auditScript>0){
      auditFile <- paste0(dirname(fileLoc),"/", format(Sys.time(), "%Y_%m_%d_%H%M"), "_audit.txt")
      if(verbose){
        print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"))
      }
      suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 1"), file = auditFile, append = FALSE))
    }

    if(grepl(" ",fileLoc) && grepl(" ",makeblastdbPath) && grepl(" ",taxaDBLoc)){
      if(verbose){
        print("Error: One or more of the file paths contains a space in the naming convention. Please change the naming and try again.")
      }
    } else {

      #Audit line
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 2"), file = auditFile, append = TRUE))}

      #Get the file name of the target file
      fileName <- sub("\\..*","", as.vector(basename(fileLoc)))

      #Read in the data from the target file
      seqTable <- read.delim(fileLoc, header = FALSE)

      #Audit line
      if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 3"), file = auditFile, append = TRUE))}

      #taking the read in file and changing from Fasta to tab delimited
      seqTableTemp <- data.frame(Header = (seqTable[seq(from = 1, to = nrow(seqTable), by = 2), 1]))
      seqTableTemp["Sequence"] <- seqTable[seq(from = 2, to = nrow(seqTable), by = 2), 1]
      seqTable<-seqTableTemp
      #Remove the seqTableTemp to free memory
      rm(seqTableTemp)

      if(!grepl(">", seqTable[3,1], fixed = TRUE)){
        if(verbose){
          print("********************************************************************************")
          print("The submitted fasta file is not in the single line nucleotide format which is ")
          print("needed for this script. Please correct the format and rerun this script.")
          print("********************************************************************************")
        }
      }else if (mean(sapply(seqTable[,1], function(x) sum(gregexpr("\\|", x)[[1]] >= 0))) != 5){
        if(verbose){
          print("********************************************************************************")
          print("The submitted fasta file is not in the MACER file format which is ")
          print("needed for this script. Please correct the format and rerun this script.")
          print("********************************************************************************")
        }
      }else{

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 4"), file = auditFile, append = TRUE))}

        #Create file folder to hold the database in the same location as the fasta file
        dbLoc <- paste0(dirname(fileLoc), "/", dateStamp, "_", dbName,"_Database")
        dir.create(dbLoc)

        #Set the working directory to the newly created output directory
        setwd(dbLoc)

        #Load in the MACER formatted data
        seqTable<-cbind(data.frame(do.call("rbind", strsplit(as.character(seqTable[,1]), "|", fixed = TRUE))),seqTable[,2])

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 5"), file = auditFile, append = TRUE))}

        #Get the taxonomy id's from the species names and place them into the table
        suppressWarnings(taxaIDResults<-getId(paste0(seqTable[,3], " ",seqTable[,4]),taxaDBLoc))

        #build the output format for the makedb
        seqTable<-cbind(paste(seqTable[,1],taxaIDResults, seqTable[,3],seqTable[,4],sep="|"),seqTable[,7])

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 6"), file = auditFile, append = TRUE))}

        if(max(nchar(seqTable[,1]))>50){

          #Reduce the first column if more than 50 characters
          seqTable<-cbind(substr(seqTable[,1], 1, 50),seqTable[,2])

        }

        #Remove records if the sequences are less than the supplied minimum
        seqTable<-seqTable[nchar(seqTable[,2]) > minLen,]

        #Remove duplicated records.
        seqTable<-seqTable[!duplicated(seqTable[,1]),]

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 7"), file = auditFile, append = TRUE))}

        #Write the checked fasta file to the new location to create the DB
        fileLoc <-paste0(fileName,"_DB_Create_File_",dateStamp,".fas")

        #Save the newly formatted file
        write.table(seqTable, fileLoc , row.names=FALSE, col.names=FALSE, quote = FALSE, sep="\n", append=FALSE)

        #Build the BLAST command
        BLASTMakeDBCmdString<- paste0(makeblastdbPath, " -in ", fileLoc, "  -parse_seqids -title '", dbName, "' -dbtype nucl -out ", dbName)
        if(verbose){
          print("********************************************************************************")
          print(paste0("Begin makeblastdb at time: ", Sys.time()))
          print("********************************************************************************")
        }
        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 8"), file = auditFile, append = TRUE))}

        #Build the file to run based on the operating system
        if (.Platform$OS.type == "windows"){

          #Windows command file
          blastCommandFile <- paste0(dbLoc,"/", fileName, "_MAKE_DB_BLAST_", dateStamp, ".bat")
          write(BLASTMakeDBCmdString, file = blastCommandFile, append = FALSE)

          #Run the BLAST command in a system command
          system(blastCommandFile)

        } else{

          #linux and Mac OS command file
          blastCommandFile <- paste0(dbLoc,"/", fileName, "_MAKE_DB_BLAST_", dateStamp, ".sh")

          #Initialize the file that will be run for the BLAST
          write("#!/bin/sh", file = blastCommandFile, append = FALSE)
          write("\n", file = blastCommandFile, append = TRUE)
          write(BLASTMakeDBCmdString, file = blastCommandFile, append = TRUE)

          #Run the BLAST command in a system command
          system(paste0("bash '", blastCommandFile, "'"))

        }#Closing off the checking platform if

        #Audit line
        if(auditScript>0){print(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9")); suppressWarnings(write(paste0(format(Sys.time(), "%Y_%m_%d %H:%M:%S"), " - Audit: 9"), file = auditFile, append = TRUE))}

      }#Closing the if - else if - else the format of the file is correct
    }#Closing off the check to see if the file names had spaces
  }#Closing the if checking to see if user supplied arguments were included

  if(verbose){
    print(paste0(startTime, " - End time...", Sys.time()))
  }
}#Close function

