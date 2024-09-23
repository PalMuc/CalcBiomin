library(topGO)
toplibrary(GO.db)
library(dplyr)
# set your working directory here
#setwd("<add path to repository/R_scripts/GOs")

####################################
#
# set/read files
#
###########################################################

#once the wd is defined you can point to the files using paths relative to the wd

#Molecular Function table
functionPath<-"../../inputfiles/GOterms/gene_GO_Terms/function_out_gene.tsv"

#Biological process table
processPath<-"../../inputfiles/GOterms/gene_GO_Terms/process_out_gene.tsv"

#Celullar component table
componentPath<-"../../inputfiles/GOterms/gene_GO_Terms/component_out_gene.tsv"
# List of transcript filenames
file_names <- c( "../DESeq2/Sci_body-parts_over_DEGs_osculum_vs_bodywall_p001-L2FC2.csv",
               "../wcgna/midnightblue.gene.list"
                )

#all Go terms
go_function<-readMappings(functionPath)
go_process<-readMappings(processPath)
go_component<-readMappings(componentPath)

#process each transcript list
for (file_name in file_names) {
  print (file_name)# Read in the file
  transcript_names<-read.table(file_name)

  ###########################################
  #
  # GOTerm enrichment for Molecular Function
  #
  ########################################### 
  
  module_transcripts_for_GO<-factor(as.integer(names(go_function) %in% transcript_names$V1))
  names(module_transcripts_for_GO)<-names(go_function)
  table(module_transcripts_for_GO)
  function_go_data<-new("topGOdata", ontology="MF", allGenes=module_transcripts_for_GO, annot=annFUN.gene2GO, gene2GO=go_function, nodeSize=10)
  #get the number for non trivial nodes from function_resultsFisherClassic
    function_resultsFisherClassic<-runTest(function_go_data, algorithm = "classic", statistic="fisher")
    # Extract the geneData slot
  gene_data <- function_resultsFisherClassic@geneData
  
  # Get the sigTerm (fourth element from the geneData)
  sigTerm <- gene_data[4]
 
  # Create a new filename for the output
  file_basename<-basename(file_name)
  output_filename <- paste0("function_enrichedGOs_", file_basename)
  GoTable<-GenTable(function_go_data, classicFisher=function_resultsFisherClassic,orderBy="classicFisher", ranksOf="classicFisher", topNodes=sigTerm)
  
  #replace truncated GO terms in GoTable with the complete terms from GO.db
  # Assuming df is your data frame with columns GO.ID, Term, Annotated, etc.
  go_ids <- GoTable$GO.ID

  #replace the partial Terms from TopGO with complete terms from GO.db using dplyr
  
  GoTable <- GoTable %>%
    rowwise() %>%
    mutate(Term = as.character(Term(GO.ID)))
  # Save the processed data
  write.csv(GoTable, output_filename)
  
  ###########################################
  #
  # GOTerm enrichment for Biological Process
  #
  ########################################### 
  table(module_transcripts_for_GO)
  process_go_data<-new("topGOdata", ontology="BP", allGenes=module_transcripts_for_GO, annot=annFUN.gene2GO, gene2GO=go_process, nodeSize=10)
  process_resultsFisherClassic<-runTest(process_go_data, algorithm = "classic", statistic="fisher")
  # Extract the geneData slot
  gene_data <- process_resultsFisherClassic@geneData
  # Get the sigTerm (fourth element from the geneData)
  sigTerm <- gene_data[4]
  
  GoTable<-GenTable(process_go_data, classicFisher=process_resultsFisherClassic,orderBy="classicFisher", ranksOf="classicFisher", topNodes=sigTerm)
 
  #replace truncated GO terms in GoTable with the complete terms from GO.db
  GoTable <- GoTable %>%
    rowwise() %>%
    mutate(Term = as.character(Term(GO.ID)))

   # Save the processed data
   output_filename <- paste0("process_enrichedGOs_", file_basename)
   write.csv(GoTable, output_filename)
   
   ###########################################
   #
   # GOTerm enrichment for Cellular Component
   #
   ########################################### 
   table(module_transcripts_for_GO)
   component_go_data<-new("topGOdata", ontology="CC", allGenes=module_transcripts_for_GO, annot=annFUN.gene2GO, gene2GO=go_component, nodeSize=10)
   component_resultsFisherClassic<-runTest(component_go_data, algorithm = "classic", statistic="fisher")
   # Extract the geneData slot
   gene_data <- component_resultsFisherClassic@geneData
   # Get the sigTerm (fourth element from the geneData)
   sigTerm <- gene_data[4]
   
   GoTable<-GenTable(component_go_data, classicFisher=component_resultsFisherClassic,orderBy="classicFisher", ranksOf="classicFisher", topNodes=sigTerm)
   
   #replace truncated GO terms in GoTable with the complete terms from GO.db
   GoTable <- GoTable %>%
     rowwise() %>%
     mutate(Term = as.character(Term(GO.ID)))
   
   # Save the processed data
   output_filename <- paste0("component_enrichedGOs_", file_basename)
   write.csv(GoTable, output_filename)
}



