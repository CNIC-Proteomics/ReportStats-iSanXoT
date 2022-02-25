# import modules
library("limma")
library("logging")
library("shiny")
library("tcltk")

# set constants
basicConfig()

obj <- list(
  "jobID" = "",
  
  "filePath"="", # path to report
  "objPath"="", # path to object with metadata
  
  "colNames"=c(),
  "colTypes"=c(),
  
  "integrations" = c(), # Which integrations are in the report
  "samples" = c(), # Which samples are in the report
  "Zcols" = c(), # Index of columns with Z values
  
  #"levels" = c(), # which levels are in the report
  #"Lcols" = c(), # column index of each level
  
  "integrationSet" = c(), # set of integrations
  "sampleSet" = c(), # set of samples
  
  "outPath" = "" # path to output file
  
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
    
    if (i != "LEVEL" & i != "REL" & i != "STATS" & i != "EXTRA") {
      updatedObj$Zcols <- append(updatedObj$Zcols, index)
      updatedObj$integrations <- append(updatedObj$integrations, 
                                        updatedObj$colNames[index])
      updatedObj$samples <- append(updatedObj$samples, i)
    } #else if (i == "LEVEL") {
      #updatedObj$Lcols <- append(updatedObj$Lcols, index)
      #updatedObj$levels <- append(updatedObj$levels,
      #                            updatedObj$colNames[index])
    #}
   
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

LIMMA <- function(obj, Target, x, integration, eset, type) {
  
  f <- factor(Target, levels=names(obj$sampleGroups))
  design <- model.matrix(~0+f)
  colnames(design) <- names(obj$sampleGroups)
  fit <- lmFit(eset, design)
  
  #x <- gsub(" vs ", "-", obj$hypTesting)
  contrast.matrix <- makeContrasts(contrasts=x, levels=names(obj$sampleGroups))
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  tmp <- as.data.frame(fit2$p.value)
  newColname <- c()
  for (i in colnames(tmp)) {
    newColname <- append(newColname, paste(integration, type, i, sep="_"))
  }
  colnames(tmp) <- newColname
  
  loginfo(paste0(integration, " - Prior Variance: ", fit2$s2.prior), logger=obj$jobID)
  loginfo(paste0(integration, " - Prior Degrees of Freedom: ", fit2$df.prior), logger=obj$jobID)
  
  return (tmp)
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
    
    loginfo(paste0("Performing calculations: ", integration), logger=obj$jobID)
    
    # Get vector with low level
    lowLevel <- strsplit(gsub("Z_", "", integration), "2")[[1]][1]
    lowLevelCol <- as.vector(
      unlist(reportData[obj$colTypes == "LEVEL" & obj$colNames == lowLevel][1])
    )
    pvalues_df_tmp <- data.frame(row.names = 1:nrow(reportData))
    pvalues_df_tmp$low <- lowLevelCol
    
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
    x <- gsub(" vs ", "-", obj$hypTesting)
    
    if (obj$testType['limma']) { # removing duplicates
      loginfo(paste0(integration, " - Applying LIMMA"), logger=obj$jobID)
      
      lowLevelSet_bool <- !duplicated(lowLevelCol)
      tmp <- LIMMA(obj, Target, x, integration, eset[lowLevelSet_bool,], "limma")
      
      lowLevelSet <- lowLevelCol[lowLevelSet_bool]
      tmp$low <- lowLevelSet
      
      # merge is changing the order... and we MUST preserve it
      pvalues_df_tmp$index <- 1:nrow(pvalues_df_tmp)
      pvalues_df_tmp <- merge(pvalues_df_tmp, tmp, by="low", sort=FALSE) # still changes order
      pvalues_df_tmp <- pvalues_df_tmp[order(pvalues_df_tmp$index),]
      pvalues_df_tmp <- pvalues_df_tmp[, !(names(pvalues_df_tmp) %in% c("index"))]
    }
    
    
    if (obj$testType['limmaDup']) {
      loginfo(paste0(integration, " - Applying LIMMA with duplicates"), logger=obj$jobID)
      tmp <- LIMMA(obj, Target, x, integration, eset, "limmaDup")
      pvalues_df_tmp <- cbind(pvalues_df_tmp, tmp)
    }
    
    #TTEST
    if (obj$testType['ttest']) {
      loginfo(paste0(integration, " - Calculating t-test"), logger=obj$jobID)
      pvalues_ttest <- classicTTEST(eset, Target, x, integration)
      pvalues_df_tmp <- cbind(pvalues_df_tmp, pvalues_ttest)
    }
    
    pvalues_df <- cbind(pvalues_df, pvalues_df_tmp[-1])
    
  }
  
  
  # OUTPUT TABLE
  loginfo("Generating output table...", logger=obj$jobID)

  header <- read.csv(obj$filePath, header=FALSE, sep="\t", nrows=2)
  
  for (i in colnames(pvalues_df)) {
    header <- cbind(header, c(i, "STATS"))
  }
  
  reportData <- cbind(reportData, pvalues_df)
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
    obj$testType <- unlist(input$execute$testType)
    
    loginfo("Executing LIMMA...", logger=obj$jobID)
    obj$outPath = calculatePvalues(obj)
    loginfo("LIMMA was executed", logger=obj$jobID)
    
    session$sendCustomMessage("execOpen", obj$outPath)
  })
  
  observeEvent(input$openRes, {
    shell.exec(dirname(input$openRes))
  })
  
  observeEvent(input$viewLogs, {
    shell.exec(paste(getwd(), "log", paste(input$viewLogs, ".log", sep=""), sep="/"))
  })
  
  observeEvent(input$areYouAlive, {
    cat("Bottle found\n")
  })
  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui=htmlTemplate("public/index.html"), server=server)