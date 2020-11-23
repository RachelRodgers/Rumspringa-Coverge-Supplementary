# Rumspringa_Coverage_Plots.R

# Coverage plots as seen in VVVR app formatted for Rumspringa manuscript supplement.

library("data.table")
library("Gviz")
library("tidyverse")

options(ucscChromosomeNames=FALSE)

bedGraphFilePath <- "../Combined_Data/alignment_files/"
bedGraphFileList <- list.files(path = bedGraphFilePath, 
                               pattern = "_sorted.bedGraph")
filePathList <- paste0(bedGraphFilePath, bedGraphFileList)
sampleNames <- map_chr(.x = bedGraphFileList, 
                       .f = ~ str_extract(string = .x, pattern = "^[:digit:]+"))

names(filePathList) <- sampleNames


bedgraphDTList <- map(.x = filePathList, .f = fread,
                      col.names = c("chromosome", "start", "end", "value"))

# Generate top axis track (same for all)
gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")

# Generate the coverage track for each sample
dtrackList <- pmap(.l = list(range = bedgraphDTList),
                   .f = DataTrack,
                   genome = "ModCR6", type = "histogram", name = " ",
                   background.title = "#2C3E50", col.histogram = "gray28",
                   fontsize = 18)

for (i in 1:length(dtrackList)) {
  plotTracks(list(gtrack, dtrackList[[i]]), main = names(dtrackList)[[i]])
}


