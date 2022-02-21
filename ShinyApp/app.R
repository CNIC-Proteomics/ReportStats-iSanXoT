# import modules
library("limma")
library("logging")
library("shiny")
library("tcltk")

# set constants
basicConfig()

obj <- list(
  "jobID" = "",
  
  "filePath"="",
  "objPath"="",
  #"repData"=data.frame(),
  "colNames"=c(),
  "colTypes"=c(),
  
  "Zcols" = c(), # Index of columns with Z values
  "integrations" = c(),
  "samples" = c(),
  
  "integrationSet" = c(),
  "sampleSet" = c(),
  
  "outPath" = ""
  
  )


# define functions
readReport <- function(obj) {
  
  updatedObj = obj
  
  # read report data
  # updatedObj$repData <- read.csv(updatedObj$filePath, header=FALSE, sep="\t", skip=2)
  
  # read report metadata
  tmp <- read.csv(updatedObj$filePath, header=FALSE, sep="\t", nrows=2)
  updatedObj$colNames <- unlist(tmp[1,])
  updatedObj$colTypes <- unlist(tmp[2,])
  
  # identify column indexes containing Z values
  index <- 1
  for (i in updatedObj$colTypes) {
    
    if (i != "LEVEL" & i != "REL") {
      updatedObj$Zcols <- append(updatedObj$Zcols, index)
      updatedObj$integrations <- append(updatedObj$integrations, 
                                        updatedObj$colNames[index])
      updatedObj$samples <- append(updatedObj$samples, i)
    }
   
    index <- index + 1 
  }
  
  updatedObj$integrationSet <- unique(updatedObj$integrations)
  updatedObj$sampleSet <- unique(updatedObj$samples)
  
  # Save object for execution
  objPath <- paste(getwd(), "tmp", format(Sys.time(), "%Y%m%d%H%M%S.Rds"), sep="/")
  updatedObj["objPath"] = objPath
  saveRDS(updatedObj, file=objPath)
  
  return (updatedObj)
}


classicTTEST <- function(eset, Target, x, integration) {

  pvalues_ttest <- data.frame(row.names=1:nrow(eset))
  for (contrast in x) {
    g <- strsplit(contrast, '-')[[1]]
    g1_bool <- g[1] == Target
    g2_bool <- g[2] == Target
    pvalues <- apply(eset, 1, function(y) {
      if (sum(!is.na(y[g1_bool])) < 2 | sum(!is.na(y[g2_bool])) < 2) return (NA)
      return (t.test(x=y[g1_bool], y=y[g2_bool], alternative="two.sided", var.equal=TRUE)$p.value)
    })

    colname <- append(colnames(pvalues_ttest), paste(integration, "ttest", contrast, sep="_"))
    pvalues_ttest <- cbind(pvalues_ttest, pvalues)
    colnames(pvalues_ttest) <- colname
    
  }
  return(pvalues_ttest)
}


calculatePvalues <- function(obj) {
  
  reportData <- read.csv(obj$filePath, header=FALSE, sep="\t", skip=2)
  
  # data frame containing pvalues
  pvalues_df <- data.frame(row.names = 1:nrow(reportData))
  
  for (integration in obj$integrationSelected) {
    
    # dataframe containing working Z
    eset <- data.frame(row.names = 1:nrow(reportData))
    
    Target <- c()
    
    for (group in names(obj$sampleGroups)) {
      
      for (sample in obj$sampleGroups[[group]]) {
        
        integrationBool <- integration == obj$integrations
        sampleBool <- sample == obj$samples
        
        zColBool <- integrationBool & sampleBool
        
        zColIndex <- obj$Zcols[zColBool]
        
        eset <- cbind(eset, reportData[zColIndex])
        
        Target <- append(Target, group)
        
      }
      
    }
    
    # LIMMA
    f <- factor(Target, levels=names(obj$sampleGroups))
    design <- model.matrix(~0+f)
    colnames(design) <- names(obj$sampleGroups)
    fit <- lmFit(eset, design)
    
    x <- gsub(" vs ", "-", obj$hypTesting)
    contrast.matrix <- makeContrasts(contrasts=x, levels=names(obj$sampleGroups))
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    tmp <- fit2$p.value
    newColname <- c()
    for (i in colnames(tmp)) {
      newColname <- append(newColname, paste(integration, 'limma', i, sep="_"))
    }
    colnames(tmp) <- newColname
    pvalues_df <- cbind(pvalues_df, tmp)
    
    loginfo(paste0("LIMMA applied to ", integration), logger=obj$jobID)
    
    #TTEST
    pvalues_ttest <- classicTTEST(eset, Target, x, integration)
    pvalues_df <- cbind(pvalues_df, pvalues_ttest)
    
  }
  
  loginfo("Generating output table...", logger=obj$jobID)
  
  # Generate final table
  header <- read.csv(obj$filePath, header=FALSE, sep="\t", nrows=2)
  
  for (i in colnames(pvalues_df)) {
    header <- cbind(header, c(i, "STATS"))
  }
  
  reportData <- cbind(reportData, pvalues_df)
  #reportData <- data.frame(apply(reportData, 2, as.character))
  reportData <- data.frame(mapply('c', header, reportData))
  
  outDir <- dirname(obj$filePath)
  outFile <- paste("LIMMA", basename(obj$filePath), sep="_")
  outPath <- paste(outDir, outFile, sep="/")
  
  write.table(reportData, file = outPath, quote = F, sep = "\t", row.names = F,
              col.names = F, na="")
  
  loginfo(paste0("Output table was written: ", outPath), logger=obj$jobID)
  
  return(outPath)
}


# run server
server <- function (input, output, session) {
  
  # open report selected by the user
  observeEvent(input$fileSelect, {
    
    # set job id and create logging
    obj$jobID <- format(Sys.time(), "%Y%m%d%H%M%S")
    addHandler(writeToFile, logger=obj$jobID, file=paste(getwd(), "log", paste(obj$jobID, ".log", sep=""), sep="/"))
    loginfo("Selecting report file...", logger=obj$jobID)
    
    tmp = try(file.choose())
    
    if ('try-error' %in% class(tmp)) {
      logerror("File was not selected", logger=obj$jobID)
      return()
    }
    
    loginfo(paste0("File selected: ", tmp), logger=obj$jobID)
    obj["filePath"] = tmp
    
    session$sendCustomMessage("selectedFilePath", list("filePath" = obj$filePath, "jobID" = obj$jobID))
    
    # report processing
    loginfo("Processing report...", logger=obj$jobID)
    obj = readReport(obj)
    loginfo("Report processed", logger=obj$jobID)
    
    # send some data to front-end
    session$sendCustomMessage("integrationData", list(
      "jobID"=obj$jobID,
      "integrationSet"=obj$integrationSet,
      "sampleSet"=obj$sampleSet,
      "objPath"=obj$objPath
    ))
  })
  
  # make calculations
  observeEvent(input$execute, {
    
    obj <- readRDS(input$execute$objPath)
    obj$integrationSelected <- unlist(input$execute$integrationSelected)
    obj$sampleGroups <- input$execute$sampleGroups
    obj$hypTesting <- unlist(input$execute$hypTesting)
    
    loginfo("Executing LIMMA...", logger=obj$jobID)
    obj$outPath = calculatePvalues(obj)
    loginfo("LIMMA was executed", logger=obj$jobID)
    
    session$sendCustomMessage("execOpen", obj$outPath)
  })
  
  observeEvent(input$openRes, {
    shell.exec(dirname(input$openRes))
  })
  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui=htmlTemplate("public/index.html"), server=server)