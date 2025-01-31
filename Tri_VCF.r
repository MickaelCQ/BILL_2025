#01/29/25 Structuring and assembly of VCF files listing structural variants of interest
#01/30/25 Programming a function for read and extract columns natives of the VCFs and generate a csv
#01/30/25 Reflection and programming a function for the addition of INFO column metadata
#01/30/25 GIT update, comment translation and test set on genFile P15, P30, P50 and P65
#01/31/25 Add a cryptic column on first index with the willingness to identify for each SV the viral genFile

library(vcfR) # no std library R for import function for the VCF treatment

################################# 
# GENERIC FUNCTION
##################################

# For read and extract data in the VCF file
FunExtractVCF <- function(vcf_file) {
  vcf <- vcfR::read.vcfR(vcf_file)  # Read this file
  
  # Extract all column in VCF in dataframe (except INFO)
  nativesCol <- as.data.frame(vcf@fix)  # Extract 
  
  ################################# 
  # "INFO" FIELD METADATA SCALING
  ##################################
  
  FunExtractINFO <- function(info_col) {
    #separate the information for each line of each INFO column using the separator “;”.
    infos <- strsplit(info_col, ";")
    
    # Identifier les noms de champs uniques dans toutes les informations de la colonne INFO
    allFields <- unique(unlist
                        (lapply
                          (infos, function(x) {sapply(strsplit(x, "="), function(y) y[1])}
                          )#lapply end
                        )#lunlist end
                      )#unique end
    
    # Create a list for each element represente one line of the column INFO with the values associates with the fields
    infos_list <- lapply(infos, function(x) {
      pairs <- strsplit(x, "=")
      
      # Extract the key or the value
      values <- sapply(pairs, function(y) ifelse(
        length(y) > 1, y[2], y[1]
                                                ) #end ifelse
                      )#end sapply
      
      #Associated names with the values :
      names(values) <- sapply(pairs, function(y) y[1]) 
      
      # Fill in missing fields with NA:
      full <- setNames(rep(NA, length(allFields)), allFields)
      full[names(values)] <- values
      return(full)
    } #end function x
    
    )#end lapply
    
    # Merge all list items into a single dataframe:
    infos_df <- do.call(rbind, infos_list)
    infos_df <- as.data.frame(infos_df, stringsAsFactors = FALSE)  # Check for each column is treated as a string
    
    return(infos_df)
  }
  
  # Apply function to extract information of the column INFO:
  infoColumns <- FunExtractINFO(nativesCol$INFO)
  
  #######################################################
  # ADD A CRYCTIC COLUMN FOR AVERAGE DEPTH FOR EACH SV 
  #######################################################
  calcMoyCoverage <- function(coverage) {
    if (!is.na(coverage) && nzchar(coverage)) { # If value in cell and not empty
      val_prof <- as.numeric(unlist(strsplit(coverage, ","))) # strip ",", merge in vector, convert this items in numeric values
      if (length(val_prof) > 0) {
        return(mean(val_prof, na.rm = TRUE)) # ignore NA and calc mean
      }
    }
    return(NA) # implicit "else" to avoir empty cells
  }
  
  if ("COVERAGE" %in% names(infoColumns)) { # Conditional check existence for current SV of "COVERAGE" from the futur calc mean :
    infoColumns$MU_COVERAGE <- sapply(infoColumns$COVERAGE, calcMoyCoverage)
  } else {
    infoColumns$MU_COVERAGE <- NA
  }
  
  colonnes_reorganisees <- c(
    setdiff(names(infoColumns), c("COVERAGE", "MU_COVERAGE")),
    "COVERAGE",
    "MU_COVERAGE"
  )
  infoColumns <- infoColumns[, colonnes_reorganisees, drop = FALSE]  # Avoids columns if columns are missing
  
  ####################
  # CSV Cosmetics
  ####################
  finalData <- cbind(
    nativesCol[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")],
    infoColumns
  )
  
  #####################################################
  # A little Regex .....
  ####################################################
  # Extract relevance identification from the VCF file name (for classifying the structural variant)
  genFile <- gsub("^P([0-9]+)-.*$", "P\\1", basename(vcf_file))
  groupFile <- sub("^([^\\.]+\\.[^\\.]+).*", "\\1", basename(vcf_file))

  finalData$genFile <- genFile
  return(finalData)
}

#############################################
# successive treatments of VCFs of interest
#############################################
vcf_files <- list.files(path = "/home/mickael/Projets_GIT/BILL_2025/VCF/", pattern = "*.vcf", full.names = TRUE)
allData <- NULL  # Initialize as NULL to avoid conflicts at runtime

# Final LOOP
for (vcf_file in vcf_files) {
  vcfDATA <- FunExtractVCF(vcf_file)  # Extraire les données du fichier VCF
  
  if (is.null(allData)) {
    # If the First file processed , initialize with its data 
    allData <- vcfDATA
  } else {
    
    # Find all existing columns and accumulated data :
    allColumns <- unique(c(names(allData), names(vcfDATA)))
    
    # Add missing collumns at "allData"
    for (col in setdiff(allColumns, names(allData))) {
      allData[[col]] <- NA
    }
    
    # Add missing columns at "vcfDATA"
    for (col in setdiff(allColumns, names(vcfDATA))) {
      vcfDATA[[col]] <- NA
    }
    
    #Reorder columns and merge
    allData <- rbind(
      allData[, allColumns, drop = FALSE],
      vcfDATA[, allColumns, drop = FALSE]
    )#end rbind
  }#end else
}#end loop

# Export DATA in file csv with current date
fileName <- paste0("/home/mickael/Projets_GIT/BILL_2025/Traitement_VCF_", Sys.Date(), ".csv")
write.csv(allData, fileName, row.names = FALSE)


