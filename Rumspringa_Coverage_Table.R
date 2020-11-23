# Rumspringa_Coverage_Table.R

# Generate table where rows are samples, columns are variants, and values
#   are the allelic frequency of that variant in each sample.
# Samples should be grouped by:
#   Stat1 KO - IC
#   Stat1 KO - PO
#   WT - IC
#   WT - PO

library("data.table")
library("tidyverse")

# Read-in metadata information.  Add leading 0's to samples 1 - 9 for
#   easier merging of the sample type with the annotated VCF files.
metadata <- read.delim("../Combined_Data/Combined_Data_metadata.txt",
                       stringsAsFactors = FALSE)

sampleTypeLUT <- metadata %>%
  mutate(across(where(is.numeric), as.character),
         sampleType = paste(Genotype, Route, sep = ", "),
         Sample = ifelse(Sample %in% c("1":"9"), 
                       yes = paste0("0", Sample),
                       no = Sample)) %>% 
  select(Sample, sampleType) %>% 
  deframe()

# Read in each annotated variant file, select position, mutation type (only 
#   non-synonymous) and allelic frequency (only >= 50%)
vcfFilePath <- "../Combined_Data/variants/annotated_variants/"
vcfFileList <- list.files(path = vcfFilePath, 
                          pattern = "_variants_annotated.txt")
sampleNames <- map_chr(.x = vcfFileList,
                       .f = ~ str_extract(string = .x, pattern = "^[:digit:]+"))
names(vcfFileList) <- sampleNames

sampleVariantList <- vector(mode = "list")
emptySamples <- vector() # no variants, or no variants pass criteria

for (i in 1:length(vcfFileList)) {
  
  # samples with no variants will not have all reqCols
  reqCols <- c("Position", "Mutation Type", "Allelic Frequency (%)")
  
  currentSampleName <- names(vcfFileList)[i]
  currentSampleType <- sampleTypeLUT[currentSampleName]
  currentFile <- vcfFileList[[i]]
  
  currentIdx <- length(sampleVariantList) + 1
  
  currentVariantFile <- fread(file = paste0(vcfFilePath, currentFile),
                              header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
  
  if (all(reqCols %in% names(currentVariantFile))) {

    variants <- currentVariantFile %>% 
      select(all_of(reqCols)) %>%
      rename(position = Position, variant = `Mutation Type`, 
             allelic_frequency = `Allelic Frequency (%)`) %>% 
      filter(variant != "synonymous" & allelic_frequency >= 50) %>% 
      mutate(sample = currentSampleName,
             sampleType = currentSampleType) 
    
    # check for duplicated positions
    duplicatedPositions <- variants$position[duplicated(variants$position)]
    
    variants <- variants %>% 
      mutate(position = ifelse(position %in% duplicatedPositions,
                               yes = paste(position, variant, sep = "-"),
                               no = position)) %>% 
      select(sample, sampleType, position, allelic_frequency)
    
    if (dim(variants)[1] != 0) {
      
      sampleVariantList[[currentIdx]] <- variants
      names(sampleVariantList)[currentIdx] <- currentSampleName
      
    } else {
      emptySamples <- c(emptySamples, currentSampleName)
    }
      
  } else {
    emptySamples <- c(emptySamples, currentSampleName)
  }
}

# Bind all DFs in sampleVariantList together by position
variantsLong <- base::Reduce(f = function(df1, df2) {
  rbind(df1, df2)}, x = sampleVariantList)

variantsWide <- variantsLong %>% 
  pivot_wider(id_cols = c(sample, sampleType), names_from = position,
              values_from = allelic_frequency) %>% 
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

# arrange the numeric cols in variantsWide
numericCols <- names(variantsWide)[!(names(variantsWide)) %in% c("sample", "sampleType")]

numericColsDigit <- str_extract(numericCols, "[:digit:]+")
numericColsDigit <- as.numeric(numericColsDigit)
names(numericColsDigit) <- numericCols

numericColsArr <- numericColsDigit[order(numericColsDigit)]

columnOrder <- names(numericColsArr)

variantsWideArr <- variantsWide %>% 
  select(sample, sampleType, all_of(columnOrder))

# Add empty samples
columnNames <- c("sample", names(numericColsArr))
emptySamplesDF <- data.frame(matrix(ncol = length(columnNames),
                                    nrow = length(emptySamples)))
names(emptySamplesDF) <- columnNames

emptySamplesDF$sample <- emptySamples
emptySamplesDF <- emptySamplesDF %>% 
  mutate(across(c(names(numericColsArr)), ~ replace(., is.na(.), 0)),
         "sampleType" = sampleTypeLUT[sample]) %>% 
  select(sample, sampleType, everything())

# bind all samples, arrange by sampleType
variantsFull <- rbind(variantsWide, emptySamplesDF)

variantsFull$sampleType <- factor(variantsFull$sampleType,
                                  levels = c("Stat1-/-, IC", "Stat1-/-, PO",
                                             "WT, IC", "WT, PO"))
variantsFullArr <- variantsFull %>% 
  arrange(sampleType) %>% 
  select(sample, sampleType, all_of(columnOrder))

write.table(variantsFullArr, file = "rumspringa_variant_table.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
