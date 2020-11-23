# Rumspringa_Variant_Coverage.R

# What's the coverage for each sample at the 33 variant positions of interest?

library("data.table")
library("tidyverse")

#----- Data -----#

variants <- readLines("variants_of_interest.txt")
names(variants) <- variants

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

#----- Fxns -----#
GetCoverageAtPosition <- function(position, coverageRangeList) {
  for (pos_idx in 1:length(coverageRangeList)) {
    if (position %in% coverageRangeList[[pos_idx]]) {
      covValue <- names(coverageRangeList)[pos_idx]
      return(covValue)
      next()
    }
  }
}

#----- All BedGraph Files -----#

bgPath <- "../Combined_Data/alignment_files/"

bgFileList <- list.files(path = bgPath, 
                          pattern = "_sorted.bedGraph")

sampleNames <- map_chr(.x = bgFileList,
                       .f = ~ str_extract(string = .x, pattern = "^[:digit:]+"))
names(bgFileList) <- sampleNames

#----- Build Coverage DF/Sample -----#

variantCovDFList <- vector(mode = "list", length = length(bgFileList))

# For each file, read in bedgraph, generate coverageRangeList, then build
#   DF of sample name, sample type, position, and coverage

for (i in 1:length(bgFileList)) {
  
  currSampleName <- names(bgFileList)[i]
  
  # read in current sample's bedgraphDT
  bedgraphDT <- fread(paste0(bgPath, bgFileList[[i]]),
                      sep = "\t",
                      col.names = c("chromosome", "start", "end", "value"),
                      colClasses = c("character", "integer", "integer", "character"))
  
  # generate coverage map for current sample
  coverageRangeList <- vector(mode = "list", length = nrow(bedgraphDT))
  
  for (j in 1:nrow(bedgraphDT)) {
    coverageRangeList[[j]] <- seq(from = as.numeric(bedgraphDT[j, "start"]),
                                  to = as.numeric(bedgraphDT[j, "end"] - 1)) # -1 because end is also start of next range
    names(coverageRangeList)[j] <- bedgraphDT[j, "value"]
  }
  
  covMap <- map_chr(variants, GetCoverageAtPosition, 
                    coverageRangeList = coverageRangeList)
  
  sampleVariantDF <- data.frame("sample" = currSampleName,
                                "variant" = variants)
  sampleVariantDF <- mutate(sampleVariantDF,
                            "coverage" = covMap[variant])
  
  variantCovDFList[[i]] <- sampleVariantDF
  names(variantCovDFList)[i] <- currSampleName
  
}

variantCovMerged <- base::Reduce(f = function(df1, df2) {rbind(df1, df2)}, 
                                 x = variantCovDFList)

variantCovWide <- variantCovMerged %>% 
  mutate("sampleType" = sampleTypeLUT[sample]) %>% 
  pivot_wider(id_cols = c(sample, sampleType),
              names_from = variant,
              values_from = coverage)

# arrange rows
variantCovWide$sampleType <- factor(variantCovWide$sampleType,
                                    levels = c("Stat1-/-, IC", "Stat1-/-, PO",
                                               "WT, IC", "WT, PO"))

variantCovWide <- arrange(variantCovWide, sampleType)

write.table(variantCovWide, file = "rumspringa_variantCoverage_table.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

