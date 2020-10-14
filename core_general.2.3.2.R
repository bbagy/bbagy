Go_randomForest <- function(psIN, project,data_type,cutoff, numTree,numtry=NULL,name=NULL, categorical,  numeric){
  #---- install package         ------#
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  RF_out <- file.path(sprintf("%s_%s/table/RF_out",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", RF_out)) dir.create(RF_out)
  
  print(0)
  # input
  mapping <- data.frame(sample_data(psIN))
  mapping[,categorical] <- factor(mapping[,categorical])
  mapping.na <- mapping[!is.na(mapping[,categorical]), ] 
  mapping.na[mapping.na=="#N/A"] <- "NA"
  mapping.na[,categorical] <- factor(mapping.na[,categorical])
  if (length(numeric) == 1) {
    mapping.na <- mapping[!is.na(mapping[,numeric]), ]
    sample_data(psIN) <- mapping.na
  }

  
  print(1)
  # call otutable
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN)))
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame(otu_table(psIN))
  }
  
  # get taxa 
  
  if (colnames(tax_table(psIN))[1] == "KO"){
    otu.filt[,"Path.des"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("KO.des","Path.des"),level="Path.des")
    otu.filt[,"KO.des"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("KO.des","Path.des"),level="KO.des")
    otu.filt$taxa <- paste(otu.filt[,"Path.des"], otu.filt$KO.des)
  } else {
    otu.filt[,"Phylum"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Phylum")
    otu.filt[,"Genus"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Genus")
    otu.filt[,"Species"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Species")
    otu.filt$taxa <- paste(otu.filt[,"Phylum"], otu.filt$Genus, otu.filt$Species)
  }
  

  
  
  otu.filt.sel <- otu.filt


  
  if (colnames(tax_table(psIN))[1] == "KO"){
    otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$KO.des), ]
    otu.filt.sel$Path.des  <- NULL
    otu.filt.sel$KO.des  <- NULL
  } else {
    otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$Genus), ]
    otu.filt.sel$Phylum  <- NULL
    otu.filt.sel$Genus  <- NULL
    otu.filt.sel$Species <- NULL
  }
  
  
  if (dim(otu.filt)[2] == 2){
    next
  }
  
  agg <- aggregate(as.formula(sprintf(". ~ %s" , "taxa")), otu.filt.sel, sum, na.action=na.pass)
  genera <- agg[,"taxa"]
  agg <- agg[,-1]
  
  rownames(agg) <- genera
  
  print(2)
  
  c=dim(agg)[1]
  d=dim(agg)[2]
  cat("===================\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_nonzero_counts <- apply(agg, 1, function(y) sum(length(which(y > 0))))
  #hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
  
  # cutoff by percent
  remove_rare <- function( table , cutoff_pro ) {
    row2keep <- c()
    cutoff <- ceiling( cutoff_pro * ncol(table) )  
    for ( i in 1:nrow(table) ) {
      row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
      if ( row_nonzero > cutoff ) {
        row2keep <- c( row2keep , i)
      }
    }
    return( table [ row2keep , , drop=F ])
  }
  
  otu_table_rare_removed <- remove_rare(table=agg, cutoff_pro=cutoff)
  
  print(3)
  c=dim(otu_table_rare_removed)[1]
  d=dim(otu_table_rare_removed)[2]
  cat("===================\n")
  cat("   After removed\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
  otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  
  otu_table_asinh_mean_centred <- scale(asinh(agg), center=TRUE, scale=FALSE)  
  
  set.seed(151) 
  # -----  Running model ---------#
  # table
  otu_table_scaled_Categorical <- data.frame(t(otu_table_scaled))  
  # mapping
  otu_table_scaled_Categorical[,categorical] <- mapping.na[rownames(otu_table_scaled_Categorical), categorical]  


  
  #Run RF to classify inflamed and control samples:
  #numTree = 501 # default  10001등으로 늘릴수 있다.
  
  otu_table_scaled_Categorical.na <- otu_table_scaled_Categorical[complete.cases(otu_table_scaled_Categorical), ]

  if (numtry > 0) {
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  }


  
  cat("\n")
  cat("===================\n")
  cat("    RF_classify \n")
  cat("===================\n")
  print(RF_classify) 

  #Permutation Test
  RF_classify_sig <- rf.significance( x=RF_classify, xdata=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , nperm=1000 , ntree=numTree )
  
  
  #Accuracy Estimated by Cross-validation
  fit_control <- trainControl( method = "LOOCV" ) 
  
  RF_classify_loocv <- train(otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=25 ) , trControl=fit_control )
  
  print(RF_classify_loocv$results)

  
  if (length(numeric) == 1) {
    # add continuous value
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    otu_table_scaled_Numeric <- data.frame(t(otu_table_scaled))
    otu_table_scaled_Numeric[,numeric] <- as.numeric(mapping.na[rownames(otu_table_scaled_Numeric), numeric])
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    
    
    
    if (numtry > 0) {
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
    } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
    }
    
    cat("\n")
    cat("===================\n")
    cat("    RF_regress \n")
    cat("===================\n")
    print(RF_regress)
    
    
    RF_regress_sig <- rf.significance( x=RF_regress,  xdata=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , nperm=1000 , ntree=numTree ) 
    
    RF_regress_loocv <- train( otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[, ncol(otu_table_scaled_Numeric.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=215 ) , trControl=fit_control )
    
    print(RF_regress_loocv$results)
  }
  
  
  
  
  # save data
  if (length(numeric) == 1) {
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.%s.rda", RF_out, project,categorical,numeric, numTree,name, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_agg_regress_model.%s.%s.%s.rda", RF_out, project,categorical,numeric,numTree,name,format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.rda", RF_out, project,categorical,numeric, numTree, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_agg_regress_model.%s.%s.rda", RF_out, project,categorical,numeric,numTree,format(Sys.Date(), "%y%m%d")))
      
      functionReturningTwoValues <- function() { 
        results <- list()
        results$classify <- RF_classify
        results$regress <-RF_regress
        return(results) 
      }
      cat("\n")
      print("RF$classify and RF$regress are returned.")
      
      functionReturningTwoValues()
    }
  }else{
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.rda", RF_out, project,categorical, numTree,name, format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.rda", RF_out, project,categorical, numTree, format(Sys.Date(), "%y%m%d")))
    }
    cat("\n")
    print("RF$classify are returned.")
    
    return(RF_classify)
  }

}
Go_randomForest_tab <- function(project, tab, map, cutoff, numTree,numtry=NULL,name=NULL, categorical,  numeric){
  #---- install package         ------#
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  RF_out <- file.path(sprintf("%s_%s/table/RF_out",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", RF_out)) dir.create(RF_out)

  print("input")
  # input
  mapping <- map
  
  # map + tab 정리
  sel <- intersect(rownames(mapping), colnames(tab)); head(sel)
  tab <- tab[,sel, drop=F];head(tab)
  mapping.sel <- mapping[sel,, drop=F];head(mapping)
  
  
  mapping.sel[,categorical] <- factor(mapping.sel[,categorical])
  mapping.na <- mapping.sel[!is.na(mapping.sel[,categorical]), ] 
  mapping.na[mapping.na=="#N/A"] <- "NA"
  mapping.na[,categorical] <- factor(mapping.na[,categorical])
  if (length(numeric) == 1) {
    mapping.na <- mapping[!is.na(mapping[,numeric]), ]
    sample_data(psIN) <- mapping.na
  }
  
  
  

  
  c=dim(tab)[1]
  d=dim(tab)[2]
  cat("===================\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_nonzero_counts <- apply(tab, 1, function(y) sum(length(which(y > 0))))
  #hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
  
  # cutoff by percent
  remove_rare <- function( table , cutoff_pro ) {
    row2keep <- c()
    cutoff <- ceiling( cutoff_pro * ncol(table) )  
    for ( i in 1:nrow(table) ) {
      row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
      if ( row_nonzero > cutoff ) {
        row2keep <- c( row2keep , i)
      }
    }
    return( table [ row2keep , , drop=F ])
  }
  
  otu_table_rare_removed <- remove_rare(table=tab, cutoff_pro=cutoff)
  
  print(3)
  c=dim(otu_table_rare_removed)[1]
  d=dim(otu_table_rare_removed)[2]
  cat("===================\n")
  cat("   After removed\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
  otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  
  otu_table_asinh_mean_centred <- scale(asinh(tab), center=TRUE, scale=FALSE)  
  
  set.seed(151) 
  
  print("Running model")
  # -----  Running model ---------#
  # table
  otu_table_scaled_Categorical <- data.frame(t(otu_table_scaled))  
  # mapping
  otu_table_scaled_Categorical[,categorical] <- mapping.na[rownames(otu_table_scaled_Categorical), categorical]  


  
  #Run RF to classify inflamed and control samples:
  #numTree = 501 # default  10001등으로 늘릴수 있다.
  
  otu_table_scaled_Categorical.na <- otu_table_scaled_Categorical[complete.cases(otu_table_scaled_Categorical), ]

  if (numtry > 0) {
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  }


  
  cat("\n")
  cat("===================\n")
  cat("    RF_classify \n")
  cat("===================\n")
  print(RF_classify) 

  #Permutation Test
  RF_classify_sig <- rf.significance( x=RF_classify, xdata=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , nperm=1000 , ntree=numTree )
  
  
  #Accuracy Estimated by Cross-validation
  fit_control <- trainControl( method = "LOOCV" ) 
  
  RF_classify_loocv <- train(otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=25 ) , trControl=fit_control )
  
  print(RF_classify_loocv$results)

  
  if (length(numeric) == 1) {
    # add continuous value
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    otu_table_scaled_Numeric <- data.frame(t(otu_table_scaled))
    otu_table_scaled_Numeric[,numeric] <- as.numeric(mapping.na[rownames(otu_table_scaled_Numeric), numeric])
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    
    
    
    if (numtry > 0) {
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
    } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
    }
    
    cat("\n")
    cat("===================\n")
    cat("    RF_regress \n")
    cat("===================\n")
    print(RF_regress)
    
    
    RF_regress_sig <- rf.significance( x=RF_regress,  xdata=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , nperm=1000 , ntree=numTree ) 
    
    RF_regress_loocv <- train( otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[, ncol(otu_table_scaled_Numeric.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=215 ) , trControl=fit_control )
    
    print(RF_regress_loocv$results)
  }
  
  
  
  
  # save data
  if (length(numeric) == 1) {
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_tab_classify_model.%s.%s.%s.rda", RF_out, project,categorical,numeric, numTree,name, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_tab_regress_model.%s.%s.%s.rda", RF_out, project,categorical,numeric,numTree,name,format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_tab_classify_model.%s.%s.rda", RF_out, project,categorical,numeric, numTree, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_tab_regress_model.%s.%s.rda", RF_out, project,categorical,numeric,numTree,format(Sys.Date(), "%y%m%d")))
      
      functionReturningTwoValues <- function() { 
        results <- list()
        results$classify <- RF_classify
        results$regress <-RF_regress
        return(results) 
      }
      cat("\n")
      print("RF$classify and RF$regress are returned.")
      
      functionReturningTwoValues()
    }
  }else{
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_tab_classify_model.%s.%s.rda", RF_out, project,categorical, numTree,name, format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_tab_classify_model.%s.rda", RF_out, project,categorical, numTree, format(Sys.Date(), "%y%m%d")))
    }
    cat("\n")
    print("RF$classify are returned.")
    
    return(RF_classify)
  }

}


Go_confusion <- function(model, project,legend, name, height, width){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  rf_model <- readRDS(model)
  
  # run
  conf_mat <- rf_model$confusion[, -ncol(rf_model$confusion)]
  melt_conf_mat <- reshape2::melt(conf_mat, na.rm = TRUE)
  table <- data.frame(melt_conf_mat);table
  
  plotTable <- table %>%
    mutate(goodbad = ifelse(table$Var1 == table$Var2, "good", "bad")) %>%
    group_by(Var2) %>%
    mutate(prop = value/sum(value))
  
  
  # plot
  p <- ggplot(data = plotTable, mapping = aes(x = Var2, y = Var1, fill = goodbad, alpha = prop)) +
    geom_tile() +
    geom_text(aes(label = value), vjust = .5, fontface  = "bold", alpha = 1) +
    scale_fill_manual(values = c(good = "blue", bad = "red")) +
    theme_bw() + xlab("True class") + ylab("Predicted class") + theme(legend.position=legend, plot.title = element_text(hjust = 0.5))
    #xlim(rev(levels(table$Var1))) 

  p <- p+ ggtitle(sprintf("confusionMatrix(%s) ",name)) 
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_2.confusion.%s.%s.%s.pdf",out_path, project, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_2.confusion.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}





Go_importance_plot <- function(psIN, model, project,title,aggregate, MDA, name, bySample, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  RF_model <- readRDS(model)
  tmp <- data.frame(RF_model$importance)
  
  if (aggregate == "YES" | aggregate == "Yes"| aggregate == "yes") {
    tmp$features <- rownames(tmp)
    tmp.sorted <- arrange(tmp, desc(MeanDecreaseAccuracy))
    RF.sorted.sel <- as.data.frame(subset(tmp.sorted, MeanDecreaseAccuracy > MDA))
    RF.sorted.sel$ShortName <- RF.sorted.sel$features
    
  } else if (aggregate == "NO" | aggregate == "No" | aggregate == "no") {
    RF <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF[,taxa] == "NA"
      RF[,taxa]<- as.character(RF[,taxa])
      RF[,taxa][is.na(RF[,taxa])] <- "__"
      for(i in 1:length(RF[,taxa])){
        if (RF[,taxa][i] == "s__" || RF[,taxa][i] == "g__" || RF[,taxa][i] == "f__" || RF[,taxa][i] == "o__" || RF[,taxa][i] == "c__"|| RF[,taxa][i] == "__"){
          RF[,taxa][i] <- ""
        }
      }
    }
    
    RF$ShortName <- paste(RF$Genus,"",RF$Species)
    for(taxa in c("Family", "Order", "Class","Phylum")){
      for(i in 1:length(RF[,taxa])){
        if (RF$ShortName[i] != "   "){
          next
        }      else if (RF$ShortName[i] == "   " & RF[,taxa][i] != ""){
          RF$ShortName[i] <- paste(RF[,taxa][i])
        }
      }
    }
    RF$features <- rownames(RF)
    RF.sorted <- arrange(RF, desc(MeanDecreaseAccuracy)  )
    RF.sorted.sel <- as.data.frame(subset(RF.sorted, MeanDecreaseAccuracy > MDA))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF.sorted.sel[,taxa] = NULL

    }
  }
  
  RF.sorted.sel.melt <- melt(RF.sorted.sel, id.vars=c("features","MeanDecreaseAccuracy","MeanDecreaseGini","ShortName"))
  

  if (bySample == "NO" |bySample == "No" | bySample == "no"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  } else if (bySample == "Yes" | bySample == "Yes"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value, fill=variable))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  }


  
  p<- p+ ggtitle(sprintf("RF %s (MDA>%s) ", title, MDA)) 
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_1.RF.%s.%s.(%s).%s.pdf",out_path, project, name, MDA,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_.1RF.%s.(%s).%s.pdf", out_path, project, MDA, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}






Go_importance_plot <- function(psIN, model, project,title,aggregate, MDA, name, bySample, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  RF_model <- readRDS(model)
  tmp <- data.frame(RF_model$importance)
  
  if (aggregate == "YES" | aggregate == "Yes"| aggregate == "yes") {
    tmp$features <- rownames(tmp)
    tmp.sorted <- arrange(tmp, desc(MeanDecreaseAccuracy))
    RF.sorted.sel <- as.data.frame(subset(tmp.sorted, MeanDecreaseAccuracy > MDA))
    RF.sorted.sel$ShortName <- RF.sorted.sel$features
    
  } else if (aggregate == "NO" | aggregate == "No" | aggregate == "no") {
    RF <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF[,taxa] == "NA"
      RF[,taxa]<- as.character(RF[,taxa])
      RF[,taxa][is.na(RF[,taxa])] <- "__"
      for(i in 1:length(RF[,taxa])){
        if (RF[,taxa][i] == "s__" || RF[,taxa][i] == "g__" || RF[,taxa][i] == "f__" || RF[,taxa][i] == "o__" || RF[,taxa][i] == "c__"|| RF[,taxa][i] == "__"){
          RF[,taxa][i] <- ""
        }
      }
    }
    
    RF$ShortName <- paste(RF$Genus,"",RF$Species)
    for(taxa in c("Family", "Order", "Class","Phylum")){
      for(i in 1:length(RF[,taxa])){
        if (RF$ShortName[i] != "   "){
          next
        }      else if (RF$ShortName[i] == "   " & RF[,taxa][i] != ""){
          RF$ShortName[i] <- paste(RF[,taxa][i])
        }
      }
    }
    RF$features <- rownames(RF)
    RF.sorted <- arrange(RF, desc(MeanDecreaseAccuracy)  )
    RF.sorted.sel <- as.data.frame(subset(RF.sorted, MeanDecreaseAccuracy > MDA))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF.sorted.sel[,taxa] = NULL

    }
  }
  
  RF.sorted.sel.melt <- melt(RF.sorted.sel, id.vars=c("features","MeanDecreaseAccuracy","MeanDecreaseGini","ShortName"))
  

  if (bySample == "NO" |bySample == "No" | bySample == "no"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  } else if (bySample == "Yes" | bySample == "Yes"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value, fill=variable))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  }


  
  p<- p+ ggtitle(sprintf("RF %s (MDA>%s) ", title, MDA)) 
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_1.RF.%s.%s.(%s).%s.pdf",out_path, project, name, MDA,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_.1RF.%s.(%s).%s.pdf", out_path, project, MDA, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}







Go_roc <- function(model, project,map, categorical, name, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  if (length(name) == 1) {
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.%s.pdf", out_path, name,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }else{
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.pdf", out_path,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }
  
  rf_model <- readRDS(model)
  
  
  pred <- predict(rf_model, type="prob")
  
  # map 정리 
  sel <- intersect(rownames(map), rownames(pred))
  map <- map[sel,, drop=F]
  #dim(pred)
  #dim(map.sel.sel)
  par(mfrow = c(1,length(colnames(pred))))
  for (i in 1:length(colnames(pred))){
    print(i)
    pred2 <- prediction(pred[,i], map[,categorical])
    perf <- performance(pred2, "tpr", "fpr")
    perf.auc <- performance(pred2, "auc")
    pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true = map[,categorical], stringsAsFactors=F); 
    pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
    confusion_matrix <- table(pred_df[, c("true", "predicted")])
    accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
    vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
    mccvalue <- mcc(vec.pred, vec.true)
    if (i == 1){
      color = "red"
    }else if(i == 2){
      color = "blue"
    }else if(i == 3){
      color = "green"
    }
    
    if (length(name) == 1) {
      plot(perf, main=sprintf("RF_ROC %s (%s)", colnames(pred)[i], name), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }else{
      plot(perf, main=sprintf("RF_ROC %s", colnames(pred)[i]), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }
  }
  dev.off()
}




# 200508 Ladas shoutgun 로제작
# table to cytoscape

Go_cytoscape <- function(tab, map, alpha, title){
  # package
  for(package in c("igraph", "RAM")){
    if(!package %in% installed.packages()){
      install.packages(package)
    }else{library(package, character.only = TRUE)}
  }
  if(!"RCy3" %in% installed.packages()){
    library("BiocManager")
    BiocManager::install("RCy3")
  }else{library("RCy3")}
  
  # matching map and table
  sel <- intersect(rownames(map.sel), colnames(tab)); head(sel)
  tab.sel <- tab[,sel, drop=F];dim(tab.sel)
  map.sel <- map.sel[sel,, drop=F];dim(map.sel)
  print("matching map and table")
  
  
  # make OTU table
  tab.sel$taxonomy <- rownames(tab); head(tab.sel$taxonomy )
  headers <- vector(dim(tab.sel)[1], mode="character")
  for (i in 1:dim(tab.sel)[1]) {
    headers[i] <- paste(i)
  }
  
  rownames(tab.sel) <- headers;head(rownames(tab.sel))
  # creating network files
  print("Creating network files")
  tab.sel.sel <-filter.OTU(data=list(tab.sel=tab.sel), percent=alpha)[[1]]
  tab.net<-network_data(tab.sel.sel, is.OTU=T, map.sel)
  tab.net.node<-tab.net[[1]]
  tab.net.edge<-tab.net[[2]]
  
  tab.net.edge1<-graph.data.frame(tab.net.edge,directed=F)
  ### cytoscape 로 보내기
  set.seed(1)
  cytoscapePing()
  print(sprintf("Sending network file to cytoscape, it is  named %s", title))
  return(createNetworkFromIgraph(tab.net.edge1,title))
  

}
  
# 200511 Ladas shoutgun 로제작
# table to cytoscape
# 버전이 변했지만 완전 다른 pipeline
Go_cytoscape <- function(tab, map, alpha, title){
  # package
  for(package in c("igraph", "RAM","Hmisc","Matrix")){
    if(!package %in% installed.packages()){
      install.packages(package)
    }else{library(package, character.only = TRUE)}
  }
  if(!"RCy3" %in% installed.packages()){
    library("BiocManager")
    BiocManager::install("RCy3")
  }else{library("RCy3")}
  
  # matching map and table
  sel <- intersect(rownames(map.sel), colnames(tab)); head(sel)
  tab.sel <- tab[,sel, drop=F];dim(tab.sel)
  map.sel <- map.sel[sel,, drop=F];dim(map.sel)
  print("matching map and table")
  
  
  # make OTU table
  L1 <- subset(tab.sel, grepl("d__Bacteria", rownames(tab.sel)))
  L2 <- subset(L1, grepl("p__", rownames(L1)))
  L3 <- subset(L2, grepl("c__", rownames(L2)))
  L4 <- subset(L3, grepl("o__", rownames(L3)))
  L5 <- subset(L4, grepl("f__", rownames(L4)))
  L6 <- subset(L5, grepl("g__", rownames(L5)))
  L7 <- subset(L6, grepl("s__", rownames(L6)))
  
  data <- L7[, setdiff(1:ncol(L7), grep(".Bacterial.kraken.1", colnames(L7)))]
  
  # taxa 분리
  rownames(data) <- gsub("\\|", ";", rownames(data))
  ranklist <- c("Rank1", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxlist <- lapply(rownames(data), function(x) parse_taxonomy_qiime(x))
  tax <- matrix(NA, nrow=nrow(data), ncol=length(ranklist)); colnames(tax) <- ranklist
  for (i in 1:length(taxlist)) {
    tax[i, names(taxlist[[i]])] <- taxlist[[i]]
  };head(tax)
  
  rownames(tax) <- rownames(data)
  ##  phyloseq object
  # tt <- tax_table(tax); rownames(tt) <- rownames(data)
  # ps <- phyloseq(otu_table(data, taxa_are_rows=T), tt);ps
  # sampledata <- read.csv("map/Ladas_woo_map_200504.csv",row.names=1,check.names=FALSE);head(sampledata)
  # ps <- merge_phyloseq(ps, sample_data(data.frame(sampledata)));ps
  #---   Start network  ---#  
  otu <- t(data)
  tax <- data.frame(tax)
  
  set.seed(1)
  dim(otu)
  
  otu.filter <- otu[,colSums(otu) >= alpha]
  
  print(c(ncol(otu),"versus",ncol(otu.filter)))
  otu.cor <- rcorr(as.matrix(otu.filter), type="spearman")
  otu.pval <- forceSymmetric(otu.cor$P) # Self-correlation as NA
  sel.tax <- tax[rownames(otu.pval),,drop=FALSE]
  all.equal(rownames(sel.tax), rownames(otu.pval))
  p.yes <- otu.pval < 0.05
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  
  p.yes.r <- abs(p.yes.r)>0.75 # output is logical vector
  p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
  
  adjm <- as.matrix(p.yes.rr)
  colnames(adjm) <- as.vector(sel.tax$Family)
  rownames(adjm) <- as.vector(sel.tax$Family)
  net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)
  
  edgew<-E(net.grph)$weight
  bad.vs<-V(net.grph)[degree(net.grph) == 0] 
  
  net.grph <-delete.vertices(net.grph, bad.vs)
  
  #write_graph(net.grph,"my_graph.gml",format="gml") 
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$net <- net.grph
    results$edgew <-edgew
    return(results) 
  }
  cat("\n")
  print("$net and $edgew are returned.")
  
  functionReturningTwoValues()
  
  

}
  

#' A Go_correlation
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_correlation()
#' May 22 2020



Go_correlation <- function(project, metadata, metabolicTab, abundTab, map, method, xanlgle, ncol,name, orders, height, width){
  # output
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  # map 정리 2
  sel.meta <- intersect(rownames(metabolicTab), rownames(map)); head(sel, "3")
  metabolicTab.sel <- metabolicTab[sel.meta,, drop=F];dim(metabolicTab.sel)
  
  sel.abun <- intersect(rownames(abundTab), rownames(map)); head(sel, "3")
  abundTab.sel <- abundTab[sel.abun,, drop=F];dim(abundTab.sel)
  
  
  
  x <- log((abundTab.sel+1)/(rowSums(abundTab.sel)+dim(abundTab.sel)[2]))
  x <- x[,order(colSums(x),decreasing=TRUE)]
  y <- metabolicTab.sel
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  if(length(name) == 1){
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.%s.pdf",out_path,project, method,name,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }else{
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.pdf",out_path,project, method,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }
  
  
  for (des in rownames(subset(metadata, Go_correlation=="yes"))){
    groups<-map[,des]
    #Now calculate the correlation between individual Taxa and the environmental data
    df<-NULL
    for(i in colnames(x)){
      for(j in colnames(y)){
        for(k in unique(groups)){
          a <- x[groups==k,i,drop=F]
          b <- y[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
          if(is.null(df)){
            df<-tmp  
          }
          else{
            df<-rbind(df,tmp)
          }    
        }
      }
    }
    
    df<-data.frame(row.names=NULL,df)
    colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
    df$Pvalue<-as.numeric(as.character(df$Pvalue))
    df$AdjPvalue<-rep(0,dim(df)[1])
    df$Correlation<-as.numeric(as.character(df$Correlation))
    
    adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
    adjustment<-5
    
    if(adjustment==1){
      df$AdjPvalue<-df$Pvalue
    } else if (adjustment==2){
      for(i in unique(df$Env)){
        for(j in unique(df$Type)){
          sel<-df$Env==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==3){
      for(i in unique(df$Taxa)){
        for(j in unique(df$Type)){
          sel<-df$Taxa==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==4){
      for(i in unique(df$Taxa)){
        sel<-df$Taxa==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    } else if (adjustment==5){
      for(i in unique(df$Env)){
        sel<-df$Env==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
    
    #Now we generate the labels for signifant values
    
    df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("***", "**", "*", ""))
    
    #We ignore NAs
    df<-df[complete.cases(df),]
    
    #We want to reorganize the Env data based on they appear
    #df$Env<-factor(df$Env,as.character(df$Env))
    
    #We use the function to change the labels for facet_grid in ggplot2
    Env_labeller <- function(variable,value){
      return(sel_env_label[as.character(value),"Trans"])
    }
    
    df$Env <- factor(df$Env)
    df$Type <- factor(df$Type, levels = orders)
    p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
    p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
    p<-p+theme(axis.text.x = element_text(angle = xanlgle, hjust = 1, vjust=0.5))
    p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
    p<-p+facet_wrap (~  Env, scales="free_x", ncol = ncol)
    #p<-p+facet_wrap (~  Env, ncol = 11)
    p<- p+ ggtitle(sprintf("%s %s", des, "Ladas")) 
    print(p)
  }
  dev.off()
}

Go_path <- function(project, pdf, table, path){
  # main dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  # main pdf
  if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
    out_pdf <- file.path(sprintf("%s/pdf",out)) 
    if(!file_test("-d", out_pdf)) dir.create(out_pdf)

    print("pdf is in your working dir. Use /dir$pdf/ for save.")
  }
  # main table
  if (table == "yes" | table == "Yes"|table == "YES"){
    out_tab <- file.path(sprintf("%s/table",out)) 
    if(!file_test("-d", out_tab)) dir.create(out_tab)

    print("table is in your working dir.Use /dir$tab/ for save.")
  }
  
  if(length(path) == 1){
    out_path <- file.path(sprintf("%s/%s",out,path)) 
    if(!file_test("-d", out_path)) dir.create(out_path)
    print("path is in your working dir. Use /dir$path/ for save.")
  }


  # 한개 이상 return 하기
  
  
  functionReturningTwoValues <- function() {
    dirs <- list()
    if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
      dirs$pdf <- out_pdf
    }
    if (table == "yes" | table == "Yes"|table == "YES"){
      dirs$tab <- out_tab
    }
    if(length(path) == 1){
      dirs$path <- out_path
    }

    return(dirs) 
  }

  functionReturningTwoValues ()
  
  
}
##############################################
#---- install package  by bioconductor ------#
##############################################
# 200429
# install package and reads library is combined

#source('http://bioconductor.org/biocLite.R')
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(bioconductor)
#BiocManager::install()

# version 1
# for (bioconductor in bioconductors){ BiocManager::install(sprintf("%s",bioconductor))}

# version 2 (better version)
#if(!"BiocManager" %in% installed.packages()){
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
#  library(BiocManager)
#}else{
#  library(BiocManager) 
#}

bioconductors <- c("dada2","DESeq2", "dplyr","ggpubr","ggfortify", "ggpmisc",
                   "illuminaio","msa","phyloseq","rstatix","useful")

for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor)
  }else{library(bioconductor, character.only = TRUE)}
}


#for (bioconductor in bioconductors){library(sprintf("%s",bioconductor), character.only = TRUE)}
##############################################
#----          install package         ------#
##############################################
packages <- c("ape", "car","cluster","CLME","cowplot","crayon", "caret","colorspace","e1071",
           "digest","data.table", "devtools","doParallel","ellipse", "emmeans",
           "gplots","ggplot2","grid","gridExtra","gplots","ggrepel",
           "Hmisc","huge","irlba","igraph","irr","lme4","lmerTest",
           "Matrix","magrittr","MASS","missForest","nlme","phangorn","plot3D",
           "pheatmap","pkgconfig","plyr","parallel","pscl","plotly","rfUtilities",
           "rlang","randomForest","readxl","RColorBrewer","ROCR","reshape","reshape2",
           "stringi","S4Vectors","ShortRead","tidyverse","vegan","VGAM") #"venneuler",
# version 1
#for (pack in packs){install.packages(sprintf("%s",pack))}
# version 2 (better version)
for (package in packages){
  if(!package %in% installed.packages()){
    install.packages(package)
  }else{library(package, character.only = TRUE)}
}

#for (package in packages){library(sprintf("%s",package), character.only = TRUE)}

#' A Go_filter
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords filter
#' @export
#' @examples
#' Go_filter


Go_filter <- function(psIN, project, alpha){
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  

  phylo_relabun <- transform_sample_counts(psIN, function(x) x / sum(x))
  phylo_filter <- filter_taxa(phylo_relabun, function(x) mean(x) < alpha, TRUE) #.00005
  rmtaxa <- taxa_names(phylo_filter)
  alltaxa <- taxa_names(phylo_relabun)
  myTaxa <- alltaxa[!alltaxa %in% rmtaxa]
  phylo_relabun_filtered <- prune_taxa(myTaxa,phylo_relabun)
  ps_filtered <- prune_taxa(myTaxa,psIN)

  cat("#--  Before filter  --#\n")
  print(psIN)
  cat("\n")
  cat("#--  After filter   --#\n")
  print(ps_filtered)

  prune_taxa(myTaxa,psIN)
  #saveRDS(ps_filtered, sprintf("%s/ps_filtered.%s.(%s).%s.rds", out, project, alpha,format(Sys.Date(), "%y%m%d")))
  return(ps_filtered)
  cat("\n")
  print(sprintf("ps_filtered is saved as 2_rds/ps_filtered.%s.(%s).%s.rds",  project, alpha,format(Sys.Date(), "%y%m%d")))

}

# R source file

#-- to make color table---#
# 190529
# fishtaco에서 색깔을 가져 왔다.
# Arsenic 분석을 예시로 하였고, 작동 된다.
# 200524 function으로 제작 

#-------- give color to phylum --------#

Go_color <- function(cdf, taxaName){
  cPalette <-cdf
  num_of_final_phyla = length(unique(cdf$PhylumCol)); num_of_final_phyla
  cPalette$h = 0
  cPalette$s = 0
  cPalette$v = 0
  num_Actinobacteria = length(grep("Actinobacteria",cPalette$PhylumCol));num_Actinobacteria
  num_Bacteroidetes = length(grep("Bacteroidetes", cPalette$PhylumCol));num_Bacteroidetes
  num_Firmicutes = length(grep("Firmicutes", cPalette$PhylumCol));num_Firmicutes
  num_Proteobacteria = length(grep("Proteobacteria", cPalette$PhylumCol));num_Proteobacteria
  num_Fusobacteria = length(grep("Fusobacteria", cPalette$PhylumCol));num_Fusobacteria
  num_Verrucomicrobia = length(grep("Verrucomicrobia", cPalette$PhylumCol));num_Verrucomicrobia
  num_TM7 = length(grep("TM7", cPalette$PhylumCol));num_TM7
  
  # Synergistetes
  
  number_of_other_phyla = num_of_final_phyla - ((num_Actinobacteria > 0) + (num_Bacteroidetes > 0) + (num_Firmicutes > 0) + (num_Proteobacteria > 0) + (num_Fusobacteria > 0)+ (num_Verrucomicrobia > 0)+ (num_TM7 > 0))
  
  #print(number_of_other_phyla)
  # Actinobacteria = green_pallete
  cPalette[grep("Actinobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0.4, s=seq(0.3,1,length.out=num_Actinobacteria), v=0.9)
  
  # Fusobacteria = orange_pallete
  cPalette[grep("Fusobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0.2, s=seq(0.3,1,length.out=num_Fusobacteria), v=0.9)
  
  # Bacteroidetes = purple_pallete
  cPalette[grep("Bacteroidetes", cPalette$PhylumCol), -1] = expand.grid(h=0.8, s=seq(0.3,1,length.out=num_Bacteroidetes), v=0.9)
  
  # Firmicutes = blue_pallete
  cPalette[grep("Firmicutes", cPalette$PhylumCol), -1] = expand.grid(h=0.6, s=seq(0.3,1,length.out=num_Firmicutes), v=0.9)
  
  # Proteobacteria = red_pallete
  cPalette[grep("Proteobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0, s=seq(0.3,1,length.out=num_Proteobacteria), v=0.9)
  
  # Verrucomicrobia = brown_pallete
  cPalette[grep("Verrucomicrobia", cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_Verrucomicrobia), v=1)
  
  # TM7 = yellow_pallete
  cPalette[grep("TM7", cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_TM7), v=1)
  
  #print(cPalette)
  #print(number_of_other_phyla)
  
  
  print(cPalette)
  
  
  
  # Add other and species name
  cPalette$PhylumCol <- taxaName
  other<-data.frame("[1_#Other]",0,0,0.75)
  names(other)<-c("PhylumCol", "h","s","v")
  color_table <- rbind(other, cPalette)
  
  
  ## hsv to color code ##
  ## taxa위치 변경을 용이 하게 하려면 hsv to color code를 해야 한다. 
  taxa_vs_color = cbind(color_table, apply(color_table[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)];taxa_vs_color
  colnames(taxa_vs_color) <- c("Taxa", "Color");taxa_vs_color
  class(taxa_vs_color)
  coloring <- as.character(taxa_vs_color$Color) 
  names(coloring) <- taxa_vs_color$Taxa;coloring
  
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$color_table <- color_table
    results$coloring <-coloring
    return(results) 
  }
  cat("\n")
  print("$color_table and $coloring are returned.")
  
  functionReturningTwoValues()
}
#' Make a rarefaction curve using ggplot2
#' @param physeq_object A phyloseq class object, from which abundance data are extracted
#' @param step Step Size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string. The name of the variable to map to text labels on the plot. Similar to color option but for plotting text.
#' @param color Default `NULL`. Character string. The name of the variable to map to the colors in the plot. This can be a sample variables among the set returned by sample_variables(physeq_object) or taxonomic rank, among the set returned by rank_names(physeq_object)
#' @param plot default `TRUE`. Logical. Should the graph be plotted
#' @param parallel default `FALSE`. Logical. Should rarefaction be parallelized
#' @param se default `TRUE`. Logical. Should standard errors be calculated.
#' @examples
#' good_taxon_table <- data.frame(sum.taxonomy = c("a;b;c;d;f;u", "p;q;r;s;t;u"),
#' site_1 = c(0,1), site_2 = c(10, 20))
#' good_maps <- data.frame(site = c("site_1", "site_2"),
#' season = c("wet", "dry"), host = c("oak", "sage"))
#' physeq_object <- convert_anacapa_to_phyloseq(good_taxon_table, good_maps)
#' ggrare(physeq_object, step = 20, se = TRUE)
#' @export


Go_rare <- function(physeq_object, step = 10, label = NULL, color = NULL, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  p <- p + xlim(NA, xlimit) + theme_classic()
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

calRare <- function(psdata, measures, depths, parallel=FALSE) {
  require('plyr') # ldply
  require('reshape2') # melt
  require('doParallel')
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}



Go_rare2 <- function(psIN, project, alpha_metrics, color, group, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  ps.rare <- calRare(psIN, alpha_metrics, rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
  summary(ps.rare)
  ps.summary <- ddply(ps.rare, c('Depth', 'Sample', 'Measure'), summarise, 
                      Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
  
  ps.summary$Sample <-  gsub("X", "", ps.summary$Sample); ps.summary$Sample
  ps.summ.verbose <- data.frame(merge(ps.summary, data.frame(sample_data(psIN)),by.x = 'Sample', by.y = 'row.names'))
  
  sample_variables(psIN)
  ps.summ.verbose$Measure
  
  
  # Plot
  plotlist <- list()
  for (am in alpha_metrics){
    ps.summ.verbose.sel <- subset(ps.summ.verbose, Measure == am)
    p <- ggplot(ps.summ.verbose.sel, mapping = aes_string(x = "Depth", y = "Alpha_diversity_mean", 
                                                   ymin = "Alpha_diversity_mean - Alpha_diversity_sd", 
                                                   ymax = "Alpha_diversity_mean + Alpha_diversity_sd", 
                                                   colour = color, group = group))+
      xlim(NA, xlimit) + theme_light() + geom_line() + #geom_pointrange(size = 0.1) +
      theme(legend.position="right", legend.text=element_text(size=8))+ guides(col = guide_legend(ncol = 2)) +
      ggtitle(sprintf("%s-Rarefaction curve", am )) + labs(x = "Sequence Sample Size", y = am)
    # + facet_wrap(facets = ~ StudyID, scales = 'free_y')
    #print(p)
    #plotlist[[length(plotlist)+1]] <- p 
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
  #multiplot(plotlist=plotlist, cols=cols, rows=rows)

}



#' A Go_qq
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords qq plot
#' @export
#' @examples
#' Go_qq()


Go_qq <- function(psIN, project, alpha_metrics, name, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/2_QQ.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"), project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } 
  else {
    pdf(sprintf("%s_%s/pdf/2_QQ.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"), project, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  
  # 1st adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  adiv$ShannonLn <-log(adiv$Shannon)
  # show last column name
  rev(names(adiv))[1]
  
  #----------- QQ plot and histogram -----------#
  par(mfrow = c(3,2))
  mes <- c(alpha_metrics, rev(names(adiv))[1])
  for (am in mes){
    test <- shapiro.test(adiv[,am])
    hist(adiv[,am], freq=F, xlab= am, main=sprintf("Histogram of %s (%s)", project, am ), cex.main=1) 
    lines(density(adiv[,am])) 
    rug(adiv[,am])
    # remove inf
    adiv.inf <- adiv[!is.infinite(adiv[,am]),]
    
    qqnorm(adiv.inf[,am], main=sprintf("Normal Q-Q Plot (%s p=%.2g)", "shapiro", test$p.value), cex.main=1)
    qqline(adiv.inf[,am])
    
    print(sprintf("%s %s shapiro test (p=%.2g)",am, project, test$p.value))
  }
  
  dev.off()
}




#' A Go_barchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_barchart()


Go_barchart <- function(psIN, metadata, project, taxRanks, data_type, simple = "no", x_label, facet, legend, orders, alpha, name, ncol, height, width,plotCols,  plotRows){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_taxa <- file.path(sprintf("%s_%s/table/taxa",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_taxa)) dir.create(out_taxa)

  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # logic for out file
  if (length(facet) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet,name, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  ranks <- taxRanks
  taxaname <- ranks
  # order by bdiv

  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])

    if (dim(otu.filt)[2] == 2){
      next
    }

    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    agg <- agg[,-1]


    agg <- normalizeByCols(agg)

    inds_to_grey <- which(rowMeans(agg) < alpha)
    genera[inds_to_grey] <- "[1_#Other]"
    agg[,taxaname[i]] <- genera
    
    
    #saving table
    agg_other_out <- subset(agg, agg[,taxaname[i]] != "[1_#Other]")
    write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_abundance.(%s).%s.%s.csv", out_taxa, project, alpha,taxaname[i], format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
    
    
    
    
    df <- melt(agg, variable.name="SampleID")


    # add StduyID
    df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , taxaname[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)

    #mapping.sel[df2$SampleID, "StudyID"]

    # add groups
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, mvar])
      df2[,mvar] <- mapping.sel[df2$SampleID, mvar]

      # order
      if (length(orders) >= 1) {
        df2[,mvar] <- factor(df2[,mvar], levels = orders)
      }
      else {
        df2[,mvar] <- factor(df2[,mvar])
      }
    }

    # adding facet to groups
    if (length(facet) == 1) {
      for (fa in facet){
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
      }
    }
    
    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    } 

    print(1)
    # color
    colourCount = length(unique(df2[,taxaname[i]]));colourCount
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    # pdf size height = 5, width=9

    if (legend == "bottom"){
      if (colourCount < 30) {
        col <- 4
      }
    } else if (legend == "right") {
      if (colourCount < 18) {
        col <- 1
      }
      else if (colourCount > 19 & colourCount  < 35) {
        col <- 2
      }
      else if (colourCount > 36) {
        col <- 3
      }
    }

    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    print(2)

    #facet <- "SampleType"
    #mvar <- "TreatmentGroup"
    df2[,facet] <- factor(df2[,facet], levels = orders)
    if (length(facet) == 1) {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
        }
        if (facet == mvar) {
          next
        }
        df2[,facet] <- factor(df2[,facet], levels = orders)
        
        print(3)
        
        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxaname[i], order=taxaname[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col))  + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")


        if (length(name) == 1) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
        }

        print(p)
        #plotlist[[length(plotlist)+1]] <- p
      }

    }     else if (length(facet) != "NULL"& simple == "no") {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }

        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxaname[i], order=taxaname[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    } else if (length(facet) != "NULL" & simple == "yes") {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }
        
        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxaname[i], order=taxaname[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + labs(y = "Relative abundance")
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  #multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}




#' A Go_barchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_barchart()


Go_barchart_spcies <- function(psIN, metadata, project, data_type,simple, x_label, facet, legend, orders, alpha, name, height, width,plotCols, plotRows){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_taxa <- file.path(sprintf("%s_%s/table/taxa",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_taxa)) dir.create(out_taxa)

  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # logic for out file
  if (length(facet) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_barchart_species.%s.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet,name, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_barchart_species.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_barchart_species.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_barchart_species.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  taxRanks <- c("Genus", "Species")
  ranks <- taxRanks
  taxaname <- ranks
  
  # order by bdiv
  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  
  # dada2 or nephele
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN)))
  }
  else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame(otu_table(psIN))
  }
  
  # continue
  otu.filt[,"Genus"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level="Genus")
  otu.filt[,"Species"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level="Species")
  
  otu.filt$Genus.Species <- paste(otu.filt$Genus,"",otu.filt$Species)
  otu.filt.sel <- otu.filt
  otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$Genus), ]
  otu.filt.sel$Genus  <- NULL
  otu.filt.sel$Species <- NULL
  
  if (dim(otu.filt)[2] == 2){
    next
  }
  
  agg <- aggregate(as.formula(sprintf(". ~ %s" , "Genus.Species")), otu.filt.sel, sum, na.action=na.pass)
  genera <- agg[,"Genus.Species"]
  agg <- agg[,-1]
  
  
  agg <- normalizeByCols(agg)
  
  inds_to_grey <- which(rowMeans(agg) < alpha)
  genera[inds_to_grey] <- "[1_#Other]"
  agg[,"Genus.Species"] <- genera
  
  #saving table
  agg_other_out <- subset(agg, Genus.Species != "[1_#Other]")
  write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_abundance.(%s).%s.%s.csv", out_taxa, project, alpha,"Genus.Species", format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  df <- melt(agg, variable.name="SampleID")
  
  
  # add StduyID
  df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , "Genus.Species")), df, sum)
  df2$SampleID <- as.character(df2$SampleID)
  df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
  df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)
  
  #mapping.sel[df2$SampleID, "StudyID"]
  
  # add groups
  for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
    df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, mvar])
    df2[,mvar] <- mapping.sel[df2$SampleID, mvar]
    
    # order
    if (length(orders) >= 1) {
      df2[,mvar] <- factor(df2[,mvar], levels = orders)
    }
    else {
      df2[,mvar] <- factor(df2[,mvar])
    }
  }
  
  # adding facet to groups
  if (length(facet) == 1) {
    for (fa in facet){
      df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
      df2[,fa] <- mapping.sel[df2$SampleID, fa]
    }
  }
  
  if (x_label == "SampleID"| x_label == "SampleIDfactor"){
    df2 <- df2
  } else if (length(x_label) >= 1){
    df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
  }  
  
  
  print(1)
  # color
  colourCount = length(unique(df2[,"Genus.Species"]));colourCount
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  # pdf size height = 5, width=9
  
  if (legend == "bottom"){
      col <- 4
  } else if (legend == "right") {
    if (colourCount < 18) {
      col <- 1
    }
    else if (colourCount > 19 & colourCount  < 35) {
      col <- 2
    }
    else if (colourCount > 36) {
      col <- 3
    }
  }
  
  # plot
  # df2 <- df2[order(df2$value, decreasing=T),]
  print(2)
  
  #facet <- "SampleType"
  #mvar <- "TreatmentGroup"
  df2[,facet] <- factor(df2[,facet], levels = orders)
  if (length(facet) == 1) {
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      if (length(unique(df2[,mvar])) >= 1) {
        ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
      }
      if (facet == mvar) {
        next
      }
      df2[,facet] <- factor(df2[,facet], levels = orders)
      
      print(3)
      
      p <- ggplot(df2, aes_string(x= x_label, y="value", fill="Genus.Species", order="Genus.Species")) + geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill="Species")+theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col))  + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")
      
      
      if (length(name) == 1) {
        p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
      }
      else {
        p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
      }
      
      print(p)
      #plotlist[[length(plotlist)+1]] <- p
    }
    
  }     else if (length(facet) != "NULL"& simple == "no") {
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      if (length(unique(df2[,mvar])) >= 1) {
        ncol <- length(unique(df2[,mvar]))
      }
      
      p <- ggplot(df2, aes_string(x= x_label, y="value", fill="Genus.Species", order="Genus.Species")) + geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill="Species")+theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")
      if (length(name) == 1) {
        p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
      }
      else {
        p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
      }
      #plotlist[[length(plotlist)+1]] <- p
      print(p) 
    }
  }else if (length(facet) != "NULL" & simple == "yes") {
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      if (class(ncol) == "numeric") {
        ncol <- ncol
      }else if(length(unique(df2[,mvar])) >= 1){
        ncol <- length(unique(df2[,mvar]))
      }
      
      p <- ggplot(df2, aes_string(x= x_label, y="value", fill="Genus.Species", order="Genus.Species")) + geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill="Species")+theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + labs(y = "Relative abundance")
      if (length(name) == 1) {
        p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
      }
      else {
        p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
      }
      #plotlist[[length(plotlist)+1]] <- p
      print(p)
    }
  }
  #multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}




#' A Go_overview
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_overview()


Go_overview <- function(psIN, metadata, ylabn = "", facet, Color, orders, name,xanlgle, height, width){
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  if (length(name) == 1) {
    pdf(sprintf("%s/3_taxa_overview.%s.%s.%s.%s.pdf", out_path, facet, project, name, format(Sys.Date(), "%y%m%d")), height = height, width=width)
  } else{
    pdf(sprintf("%s/3_taxa_overview.%s.%s.%s.pdf", out_path, facet, project ,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  print("Calculating relative abundance .....")
  ps3ra = transform_sample_counts(psIN, function(x){x / sum(x)})
  mphyseq <- psmelt(ps3ra)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  for (maingroup in rownames(subset(metadata, Go_overview =="yes"))) {
    mphyseq[,maingroup] <- as.character(mphyseq[,maingroup]);mphyseq[,maingroup]
    mphyseq[,maingroup][mphyseq[,maingroup]==""] <- "NA";mphyseq[,maingroup]
    mphyseq[,maingroup]<- as.factor(mphyseq[,maingroup]);mphyseq[,maingroup]
    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
    mphyseq.na <- subset(mphyseq, mphyseq[,maingroup] != "NA");mphyseq.na[,maingroup] 
    
    if (facet == "Genus") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "g__");mphyseq.na.na[,facet] 
    } else if (facet == "Family") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "f__");mphyseq.na.na[,facet] 
    } else if (facet == "Order") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "o__");mphyseq.na.na[,facet] 
    } else if (facet == "Class") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "c__");mphyseq.na.na[,facet] 
    } else if (facet == "Phylum") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "p__");mphyseq.na.na[,facet] 
    }else if (facet == "Species") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "s__");mphyseq.na.na[,facet] 
    }
    
    
    if (length(orders) >= 1) {
      mphyseq.na.na[,maingroup] <- factor(mphyseq.na.na[,maingroup], levels = orders)
    }       else {
      mphyseq.na.na[,maingroup] <- factor(mphyseq.na.na[,maingroup])
    }
    
    p<- ggplot(data = mphyseq.na.na,  mapping = aes_string(x = maingroup, y = "Abundance",color = Color, fill = Color)) +
      geom_violin(fill = NA) + theme_bw()  +#scale_colour_brewer(type="qual", palette="Set4") + #+ 
      geom_point(size = 1, alpha = 0.3, position = position_jitter(width = 0.3)) +
      theme(title=element_text(size=8), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5)) +
      facet_wrap(facets = facet) + ylab(ylabn) + scale_y_log10()
    print(p)
  }
  
  dev.off()
}

#' A Go_colbarchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_colbarchart()
#' 20200525
#' color for phylum

Go_colbarchart <- function(psIN, metadata, project, taxRanks, data_type, x_label, facet, legend, orders, alpha, name, ncol,height, width,plotCols, plotRows){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # logic for out file
  if (length(facet) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_colbarchart.%s.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet,name, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } else {
      pdf(sprintf("%s_%s/pdf/3_colbarchart.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }else {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_colbarchart.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }else {
      pdf(sprintf("%s_%s/pdf/3_colbarchart.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  ranks <- taxRanks
  taxaname <- ranks
  
  # order by bdiv
  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])
    otu.filt$PhylumCol <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks, level="Phylum")

    if (dim(otu.filt)[2] == 2){
      next
    }

    agg <- aggregate(as.formula(sprintf(". ~ %s + PhylumCol" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    PhylumCol <- agg$PhylumCol
    agg[,taxaname[i]] <- NULL
    agg$PhylumCol <- NULL

    agg <- normalizeByCols(agg)
    inds_to_grey <- which(rowMeans(agg) < sprintf("%s", alpha))
    genera[inds_to_grey] <- "[1_#Other]"
    agg[,taxaname[i]] <- genera
    agg$PhylumCol <- PhylumCol 
    
    if (taxaname[i] == "Phylum"){
      agg$Phylum <- genera
    }
    
    df <- melt(agg, variable.name="SampleID")


    # add StduyID

    df2 <- aggregate(as.formula(sprintf("value ~ %s + PhylumCol + SampleID" , taxaname[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)

    #mapping.sel[df2$SampleID, "StudyID"]

    # add groups
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, mvar])
      df2[,mvar] <- mapping.sel[df2$SampleID, mvar]

      # order
      if (length(orders) >= 1) {
        df2[,mvar] <- factor(df2[,mvar], levels = orders)
      }
      else {
        df2[,mvar] <- factor(df2[,mvar])
      }
    }
    
    # adding facet to groups
    if (length(facet) == 1) {
      for (fa in facet){
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
      }
    }
    
    


    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    }  


    print(1)
    #------------------------#
    # ---  Color table   --- #
    #------------------------#
    agg$PhylumCol <- PhylumCol 
    agg[,taxaname[i]] <- genera
    

    #-------- remove other from taxa table --------#
    TaxaTab <- agg[order(agg[,taxaname[i]] ,  decreasing = TRUE), ]
    cdf <- data.frame(subset(TaxaTab, select=c("PhylumCol", taxaname[i])))
    cdf.sel <- subset(cdf, cdf[,taxaname[i]] != "[1_#Other]");dim(cdf.sel)[1]
    
    # 몇개인지 결정후 Phylum 으로 정리
    N <- dim(cdf.sel)[1]
    cdf.sel <- cdf.sel[order(cdf.sel$PhylumCol ,  decreasing = FALSE), ]
    cdf.sel <- data.frame(as.character(cdf.sel$PhylumCol[1:N]), as.character(cdf.sel[,taxaname[i]][1:N]))
    colnames(cdf.sel) <- c("PhylumCol", taxaname[i])
    #cdf.sel[ ,c("Kingdom","Class", "Order", "Family","Genus")] <- list(NULL)
    
    cdf.sel[,taxaname[i]] <-  gsub("p__", "", gsub("c__", "", gsub("o__", "", gsub("f__", "", gsub("g__", "", gsub("s__", "", cdf.sel[,taxaname[i]]))))))
    
    
    # save species name
    taxaName <- cdf.sel[,taxaname[i]]
    cdf.sel[,taxaname[i]] <- NULL
    
    # -----  create color table   ---- #
    coltab <- Go_color(cdf=cdf.sel, taxaName=taxaName)
    
    # hsv code
    #print(coltab$color_table)
    #coltab$color_table$Phylum
    # color code
    #print(coltab$coloring)
    
    print(2)
    # pdf size height = 5, width=9
    if (legend == "bottom"){
      if (N < 30) {
        col <- 5
      }
    } else if (legend == "right") {
      if (N < 18) {
        col <- 1
      }
      else if (N > 19 & N  < 35) {
        col <- 2
      }
      else if (N > 36) {
        col <- 3
      }
    }

    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    print(3)
    level <- unique(df2[,taxaname[i]])
    #facet <- "SampleType"
    #mvar <- "TreatmentGroup"
    df2[,facet] <- factor(df2[,facet], levels = orders)
    if (length(facet) == 1) {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
        }
        if (facet == mvar) {
          next
        }
        
        df2[,facet] <- factor(df2[,facet], levels = orders)
        print(4)
        
        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxaname[i]], levels=level), order=taxaname[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col))  + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance") + labs(fill = taxaname[i])

        
        if (length(name) == 1) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
        }

        print(p)
        #plotlist[[length(plotlist)+1]] <- p
      }

    } else if (length(facet) != "NULL") {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }

        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxaname[i]], levels=level), order=taxaname[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")+ labs(fill = taxaname[i])
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, alpha))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, alpha))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  #multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}




#' A Go_adiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Alpha diversity plot
#' @export
#' @examples
#' Go_adiv


Go_adiv <- function(psIN, project, alpha_metrics){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_adiv <- file.path(sprintf("%s_%s/table/adiv",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_adiv)) dir.create(out_adiv)
  
  
  # adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  cat(sprintf("adiv table is saved in %s.\n",out_path))
  cat("                                                       \n")
  write.csv(adiv, quote = FALSE, col.names = NA, 
            file=sprintf("%s/adiv.%s.%s.csv",out_adiv,project, format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
  return(adiv)
} 
  
#' A Go_box_plot
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords basic statistic and simple plots
#' @export
#' @examples
#' Go_box_plot


Go_box_plot <- function(df, metadata, project, title, alpha_metrics,Chao1_Normality, Shannon_Normality,facet, orders, name, xanlgle=90, plot, star, height, width, plotCols, plotRows){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # out file
  if (length(facet) >= 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/4_adiv.%s.%s.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"), project, plot,facet,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s_%s/pdf/4_adiv.%s.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, plot,facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/4_adiv.%s.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,plot,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s_%s/pdf/4_adiv.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,plot,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }

  # plot
  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_box =="yes"))) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}
    
    # if(mvar == facet)
    # next
    # order
    # Na 제거
    adiv <- data.frame(df)

    adiv[,mvar] <- as.character(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]

    
    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(adiv.na)[1], dim(adiv)[1]))
    
    if (length(unique(adiv.na[,mvar])) ==1) {
      next
    }
    
    summary.adiv.na <- summary(adiv.na[,mvar])
    
    # 그룹안에 샘플이 하나면 분석이 안된다. 그것을 컨트롤 하기 위한 코드..개발중 
   # for (count.summary in 1:length(unique(adiv[,mvar]))){
    #  if (summary.adiv.na[count.summary] == 1){
     #   next
    #  }
    #}


    
    # integer control
    if (class(adiv.na[,mvar]) == "character"){
      adiv.na[,mvar] <- factor(adiv.na[,mvar])
    }
    if (class(adiv.na[,mvar]) == "integer"){
      adiv.na[,mvar] <- factor(adiv.na[,mvar])
    }
    
    if (class(adiv.na[,mvar]) == "numeric"){
      adiv.na[,mvar] <- factor(adiv.na[,mvar])
    }
    
    # make a comnination for stat
    adiv.na[,mvar] <- factor(adiv.na[,mvar])
    cbn <- combn(x = levels(adiv.na[,mvar]), m = 2)

    my_comparisons <- {}
    for(i in 1:ncol(cbn)){
      x <- cbn[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    for(am in alpha_metrics){
      adiv.na[,mvar] <- factor(adiv.na[,mvar])

      if (am == "Chao1"){
        if ( Chao1_Normality  == "no"  & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- kruskal.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "KW"
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- aov(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Anova"
          testmethod <- "t.test"
          test$p.value <- summary(test)[[1]][[1,"Pr(>F)"]]
        }

        if (Chao1_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- wilcox.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Wilcox"
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- t.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "t.test"
          testmethod <- "t.test"
        }       
      }  
      
      if (am == "Shannon") {
        if ( Shannon_Normality  == "no"  & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- kruskal.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "KW"
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- aov(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)  
          test$method <- "Anova"
          testmethod <- "t.test"
          test$p.value <- summary(test)[[1]][[1,"Pr(>F)"]]
        }
        

        if (Shannon_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- wilcox.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Wilcox"
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- t.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "t.test"
          testmethod <- "t.test"
        }
      }else{
        if ( Chao1_Normality  == "no"  & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- kruskal.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "KW"
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- aov(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Anova"
          testmethod <- "t.test"
          test$p.value <- summary(test)[[1]][[1,"Pr(>F)"]]
        }
        if (Chao1_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- wilcox.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Wilcox"
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- t.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "t.test"
          testmethod <- "t.test"
        } 
        if ( Shannon_Normality  == "no"  & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- kruskal.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "KW"
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          test <- aov(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)  
          test$method <- "Anova"
          testmethod <- "t.test"
          test$p.value <- summary(test)[[1]][[1,"Pr(>F)"]]
        }
        
        
        if (Shannon_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- wilcox.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "Wilcox"
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          test <- t.test(as.formula(sprintf("%s ~ %s", am, mvar)), adiv.na)
          test$method <- "t.test"
          testmethod <- "t.test"
        }
      }
      

      # re-order
      if (length(orders) >= 1) {
        adiv.na[,mvar] <- factor(adiv.na[,mvar], levels = orders)
      }       else {
        adiv.na[,mvar] <- factor(adiv.na[,mvar])
      }
      

      # Na 제거
      if (length(facet) >= 1) {
        for (fc in facet){
          adiv.na[adiv.na[,fc] == ""] <- "NA"
          adiv.na.sel <- adiv.na[!is.na(adiv.na[,fc]), ]
          adiv.na <- adiv.na.sel 
          # facet or not
          adiv.na[,fc] <- factor(adiv.na[,fc], levels = orders)
        }
      }

      
      
      p1 <- ggplot(adiv.na, aes_string(x=mvar, y=am, colour=mvar)) + theme_bw() + labs(y=am, x=NULL) +
         # ggtitle(sprintf("%s-%s (%s p=%.2g) ", mvar, am,test$method,test$p.value)) + 
        theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5)) +
        scale_color_brewer(palette="Dark2") + theme(legend.position="none")
        
      # significant with stars list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
      
        if (length(title) == 1) {
          p1 <- p1 + ggtitle(title)
        } else{
          p1 <- p1 + ggtitle(sprintf("%s", mvar))
        }
      
      if (star == "no") {  
      p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
      }  else if (star == "yes") {
        p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
      }
      
      # plot type
      if (plot == "box") {
        p1 = p1 + geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, alpha = 3, size = 1.5, position=position_jitter(0.2)) # alpha=0.3
      }  else if (plot == "violin") {
        p1 = p1 + geom_violin(outlier.shape = NA) + geom_jitter(shape=16, alpha = 3, size = 1.5, position=position_jitter(0.2))  #geom_boxplot(width=0.2)
      } 
      
      # facet
      if (length(facet) >= 1) {
        facetCol <- length(unique(adiv[,facet]))
        p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SampleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
        
      } else {
        p1 = p1 
      }
      plotlist[[length(plotlist)+1]] <- p1 
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
#' A Go_quantile
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Continuous value to quantile
#' @export
#' @examples
#' Go_quantile

Go_quantile <- function(df, project,metadata, alpha_metrics){
  # out dir
  #dir.create(sprintf("%s_%s", project, format(Sys.Date(), "%y%m%d")));dir.create(sprintf("%s_%s/table", project, format(Sys.Date(), "%y%m%d")))
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  # make quantile
  adiv.quant <- data_frame(df$Row.names)
  colnames(adiv.quant) <- "Row.names"
  for (quan in rownames(subset(metadata, useQuantile =="yes"))) {
    adiv <- df
    adiv[adiv==""] <- "NA"
    
    adiv.na <- adiv[!is.na(adiv[,quan]), ]; dim(adiv.na)
    
    adiv.na[,quan] <- cut(adiv.na[,quan], breaks=quantile(adiv.na[,quan]), label=c("Q1", "Q2", "Q3", "Q4"), na.rm= TRUE);adiv.na[,quan]
    
    adiv.na.sel <-  adiv.na[,quan, drop=FALSE] ; dim(adiv.na.sel)
    adiv.na.sel$Row.names <- adiv.na$Row.names
    adiv.quant <- merge(adiv.quant, adiv.na.sel, by="Row.names", all = T)
  }
  
  #adiv.quant$`rownames(adiv)` <-NULL
  if (class(alpha_metrics) == "character"){
    alpha_metrics <- c("Chao1", "Shannon")
    for(am in  alpha_metrics){
      adiv.quant[,am] <- adiv[,am] 
    }
  }
  

  
  for (mvar in rownames(subset(metadata_variables, AddQuantileTab =="yes"))) {
    adiv.quant[,mvar] <- adiv[,mvar] 
  }
  
  return(adiv.quant)
  
  #write.csv(adiv, quote = FALSE, col.names = NA, 
   #         file=sprintf("%s_%s/table/adiv.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
} 
  
#' A Go_linear
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords basic statistic and simple plots
#' @export
#' @examples
#' Go_linear


Go_linear <- function(df, metadata, project, alpha_metrics,maingroup, orders, name, height, width, plotCols, plotRows){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)

  # out file

  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/4_adiv.linear.%s.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,maingroup,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } 
  else {
    pdf(sprintf("%s_%s/pdf/4_adiv.linear.%s.%s.%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,maingroup,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  

  my.formula <- y ~ x
  my.method <- "lm"
  
  # plot
  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_linear =="yes"))) {
    # Na 제거

    adiv[,mvar] <- as.character(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv[,mvar]<- as.integer(adiv[,mvar]);adiv[,mvar]
    
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(adiv.na)[1], dim(adiv)[1]))
    if (length(unique(adiv.na[,mvar])) ==1) {
      next
    }
    summary.adiv.na <- summary(adiv.na[,mvar])
    
    # na제거 in the maingroup 
    adiv.na[,maingroup] <- as.character(adiv.na[,maingroup]);adiv.na[,maingroup]
    adiv.na[,maingroup][adiv.na[,maingroup]==""] <- "NA";adiv.na[,maingroup]
    adiv.na[,maingroup]<- as.factor(adiv.na[,maingroup]);adiv.na[,maingroup]
    adiv.na.na <- subset(adiv.na, adiv.na[,maingroup] != "NA");adiv.na.na[,maingroup]
    adiv.na.na[,maingroup] <- factor(adiv.na.na[,maingroup], levels = orders)
    
    
    
    for(i in 1:length(alpha_metrics)){
      print(alpha_metrics[i])
      p<- ggplot(adiv.na.na, aes_string(x=mvar, y=alpha_metrics[i], group= maingroup, color=maingroup, linetype = maingroup))+
        theme_classic() + geom_point(size = 0.5) + scale_colour_brewer(palette = "Set1") + 
        geom_smooth(method = my.method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5 ) + 
        ggtitle(sprintf("%s by %s", mvar, alpha_metrics[i])) + theme(title=element_text(size=10)) + labs(x = NULL)+
        theme(title=element_text(size=10),
              axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain")) +
        
        
        #stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE, size = 3) +
        stat_fit_glance(method.args = list(formula = my.formula), method = my.method, 
                        #geom = 'text', 공식이 한쪽으로 정리가 되지 않고, 라인에 수치가 붙는다.
                        aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g', 
                                            stat(r.squared), stat(p.value))),
                        parse = TRUE, size = 3)
      plotlist[[length(plotlist)+1]] <- p
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
#' A Go_barchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_clme()


Go_clme <- function(psIN, metadata, project, Descriptions, alpha_metrics, paired, node, decreasing, height,timepoint, facet,ID, orders,xangle, name, width, plotCols, plotRows){
  # Descriptions 분석 하고자 하는 variation에 subgroup
  # paired 환자나 같은 사람 ID
  # node 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # decreasing 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  if (length(facet) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_clme.%s.%s.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet,name, node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_clme.%s.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet, node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/3_clme.%s.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_clme.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  
  
  
  # adiv
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  mapping <-data.frame(sample_data(psIN))
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping)
  adiv <- merge(adiv, mapping, by="row.names"); rownames(adiv) <- adiv$SampleID

  if (length(orders) >= 1) {
    adiv[,timepoint] <- factor(adiv[,timepoint], levels = orders)
  }

  
  
  # clme
  cons <- list(order = "umbrella" , node=node, decreasing = decreasing) 
  # 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  print(cons)
  
  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_clme == "yes"))) {
    print(mvar)
    for (des in Descriptions){
      print(des)
      for (am in alpha_metrics){
        form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))
        
        clme.mod <- clme(form, data = adiv[adiv[,mvar] == des,], constraints = cons, seed = 2, nsim = 1000)
        clme.sum <- summary(clme.mod, seed=2)
        clme.globalp <- function(model) { label <- substitute(
          italic(p) == globalp,
          list(globalp <- model$p.value) )
        as.character(as.expression(format(globalp, nsmall=3))) 
        }
        
        clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))
        
        # plot
        p <- ggplot(adiv[adiv[,mvar]==des,], mapping = aes_string(x=timepoint, y=am, color=timepoint, group=paired)) + geom_line(color="grey") + geom_point(size = 1.25) + xlab(timepoint) + ylab(sprintf("%s Index\n", am)) + ggtitle(sprintf("%s-%s (%s) ", mvar, des, clme.globalp))  + scale_color_brewer(palette="Set1")+theme_bw() +theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" )
        
        if (length(ID) == 1) {
          p= p + geom_text_repel(aes_string(label = ID), size = 2)
        }
        if (length(facet) == 1) {
          p= p + facet_wrap(as.formula(sprintf("~ %s"  ,facet)))
        }
        
        plotlist[[length(plotlist)+1]] <- p
      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}


Go_regression <- function(df, metadata, project, orders, outcomes, alpha, des,name){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/regression",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # data control
  
  # fix column types
  adiv <- data.frame(df)
  
  
  for (mvar in  rownames(subset(metadata, Go_reg=="yes" | Go_regConfounder=="yes"))) {
    if (metadata[mvar, "type"] == "factor") {
      adiv[,mvar] <- factor(adiv[,mvar])
      if (!(is.na(metadata[mvar, "baseline"])) && metadata[mvar, "baseline"] != "") {
        adiv[,mvar] <- relevel(adiv[,mvar], metadata[mvar, "baseline"])
      }
    } else if (metadata[mvar, "type"] == "numeric") {
      adiv[,mvar] <- as.numeric(as.character(adiv[,mvar]))
    } else if (metadata[mvar, "type"] == "date") {
      adiv[,mvar] <- as.Date(sprintf("%06d", adiv[,mvar]), format="%m%d%y")
      adiv[,mvar] <- factor(as.character(adiv[,mvar]), levels=as.character(unique(sort(adiv[,mvar]))))
    }
  }
  
  
  #-------------------------------------------------------------------#
  #--------------    regression model     -------------#
  #-------------------------------------------------------------------#
  
  for (outcome in outcomes){
    res <- {}
    for (mvar in rownames(subset(metadata, Go_reg =="yes"))) {
      
      if (outcome == mvar | outcome == "Chao1" & mvar == "Shannon" | outcome == "Shannon" & mvar == "Chao1") {
        print(sprintf("Stop function bacause out was %s and mvar was %s", outcome, mvar))
        next
      }
      # NA 제거
      adiv[,mvar] <- as.character(adiv[,mvar]);adiv[,mvar]
      adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
      #adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]
      # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
      adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
      
      print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                    mvar, dim(adiv.na)[1], dim(adiv)[1]))
      
      if (length(unique(adiv.na[,mvar])) ==1) {
        next
      }
      
      # column 정리
      if (metadata[mvar, "type"] == "factor") {
        adiv[,mvar] <- factor(adiv[,mvar])
        if (metadata[mvar, "baseline"] != "") {
          adiv[,mvar] <- relevel(adiv[,mvar], metadata[mvar, "baseline"])
        }
      } else if (metadata[mvar, "type"] == "numeric") {
        adiv[,mvar] <- as.numeric(as.character(adiv[,mvar]))
      } else if (metadata[mvar, "type"] == "date") {
        adiv[,mvar] <- as.Date(sprintf("%06d", adiv[,mvar]), format="%m%d%y")
        adiv[,mvar] <- factor(as.character(adiv[,mvar]), levels=as.character(unique(sort(adiv[,mvar]))))
      }
      
      if (length(rownames(subset(metadata, Go_regConfounder =="yes"))) >= 1){
        
        regConfounder <- rownames(subset(metadata, Go_regConfounder =="yes"))[mvar != rownames(subset(metadata, Go_regConfounder =="yes"))]
        form <- as.formula(sprintf("%s ~ %s + %s", outcome, mvar, paste(setdiff(regConfounder, "SampleType"), collapse="+")))
        print(form)
        print(1)
      } else{
        form <- as.formula(sprintf("%s ~ %s", outcome, mvar))
        print(form)
        print(3)
      }
      
      if (class(adiv[,outcome]) == "numeric"){
        mod <- lm(form, adiv)  # lm or glm or lmer
        m <- "lm"
      } else if (class(adiv[,outcome]) == "factor"){
        mod <- glm(form, adiv[adiv[,outcome] %in% levels(adiv[,outcome]),],  family = binomial(link='logit'))
        m <- "glm"
      }
      
      summary(mod)
      coef <- as.data.frame(summary(mod)$coefficients)
      coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
      colnames(coef) <- c("Estimate", "SE", "t", "pval")
      coef$outcome <- outcome
      coef$mvar <- mvar
      coef$model <- m
      
      if (length(rownames(subset(metadata, Go_regConfounder =="yes"))) >= 1){
        coef$multi <- paste(setdiff(rownames(subset(metadata, Go_regConfounder=="yes")), "SampleType"), collapse="+")
        type <-"multi"
      }else{
        type <-"uni"
      }
      
      
      
      res <- rbind(res, coef)
      
    }
    
    if (length(des) == 1) {
      res$des <- des
    }
    res$padj <- p.adjust(res$pval, method="fdr")
    #res <- res[order(res$time_point),]
    res$comp <- factor(rownames(res), levels=rownames(res))
    res$dir <- ifelse(res$pval < alpha, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
    
    if (length(des) == 1) {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.%s.csv",out_table, project,outcome, des,name, type, format(Sys.Date(), "%y%m%d"),sep="/"))
      }else{
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.csv",out_table, project,outcome, des,type, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }else{  
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.csv",out_table, project,outcome, name, type, format(Sys.Date(), "%y%m%d"),sep="/"))
      } else{ 
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.csv",out_table, project,outcome,type,format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  }
}






#' A Go_bdiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity ordination plot
#' @export
#' @examples
#' Go_bdiv()

Go_bdiv <- function(psIN, metadata, project, orders, distance_metrics, plot,shapes,ellipse, ID, facet, des, name, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)

  # out file
  if (length(des) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1) {
      pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1 & length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }


  # out file
    plotlist <- list()
    mapping.sel <- data.frame(sample_data(psIN))

    for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
      if (length(facet) >= 1){
        if (facet == mvar){
          next
        }
      } else {}
      
      if (length(shapes) >= 1){
        if (shapes == mvar){
          next
        }
      } else {}
      
        for(distance_metric in distance_metrics){
          #na remove
          mapping.sel <- data.frame(sample_data(psIN))
          mapping.sel[mapping.sel==""] <- "NA"
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          na.count <- length(mapping.sel.na)
          psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))

          if (length(des) == 1) {
            print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                          des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          } else{
            print(sprintf("##-- %s (total without NA: %s/%s) --##",
                          mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          }
          
          if (class(mapping.sel.na.rem[,mvar]) == "integer"){
            mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
            sample_data(psIN.na) <- mapping.sel.na.rem
          }
          #print(1)
          ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
          plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
            ordi = ordinate(psIN.na, method=i, distance=distance_metric)
            plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
          }, psIN.na, distance_metric)
          #print(2)
          names(plist) <- ord_meths

          pdataframe = ldply(plist, function(x){
            df = x$data[, 1:2]
            colnames(df) = c("Axis_1", "Axis_2")
            return(cbind(df, x$data))
          })
          names(pdataframe)[1] = "method"
          
          # re-order
          #if (length(orders) >= 1 & length(shapes) >= 1) {
          #  pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)
          #  pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
          #}       else {
          #  pdataframe[,mvar] <- factor(pdataframe[,mvar])
          #}

          #print(3)

          if (length(shapes) == 1) {
            pdataframe[,shapes] <- factor(pdataframe[,shapes])
            p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar, shape=shapes)) #, fill=PerioStatus_3))
            #print(5)
            p = p + geom_point(size=1.5, alpha = 3)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) #+ geom_polygon() alpha = 0.5
            p = p + scale_shape_manual(values=c(16,1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) # + scale_shape_manual(values=c(16,1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) # open(1), cross(10), closed(2)
            
           # print(4)
          } else{
            p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
            #print(5)
            p = p + geom_point(size=1.5, alpha = 3)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
          }

          p = p + facet_wrap(~ method, scales="free") + theme_bw() 
          p = p + scale_fill_brewer(type="qual", palette="Dark2")
          
          if (length(ID) == 1) {
            p = p + geom_text_repel(aes_string(label = ID), size = 2)
          } else {
            p = p 
          }
          
          if (ellipse == "yes" | ellipse == "Yes" ) {
            p = p + stat_ellipse(type = "norm", linetype = 2) 
          } else if (ellipse == "no" | ellipse == "No" ){
            p = p 
          }

          if (length(orders) >= 1) {
            p = p + scale_colour_brewer(type="qual", palette="Dark2", breaks=orders)
          }
          else {
            p = p + scale_colour_brewer(type="qual", palette="Dark2")
          }
          
          if (length(facet) == 1) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }
          else {
            p = p
          }
          
          #plotlist[[length(plotlist)+1]] <- p
          print(p)
        }
      }
    dev.off()
  }

    else {
      if (length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1) {
        pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1 & length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else {
        pdf(sprintf("%s_%s/pdf/5_ordi.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }

      plotlist <- list()
      for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
        if (length(facet) >= 1){
          if (facet == mvar){
            next
          }
        } else {}
        
        if (length(shapes) >= 1){
          if (shapes == mvar){
            next
          }
        } else {}
        
      
          for(distance_metric in distance_metrics){
            #na remove
            mapping.sel <- data.frame(sample_data(psIN))
            mapping.sel[mapping.sel==""] <- "NA"
            mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
            na.count <- length(mapping.sel.na)
            psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
            mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))


            if (length(des) == 1) {
              print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                            des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            } else{
              print(sprintf("##-- %s (total without NA: %s/%s) --##",
                            mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            }


            if (class(mapping.sel.na.rem[,mvar]) == "integer"){
              mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
              sample_data(psIN.na) <- mapping.sel.na.rem
            }

            ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
            plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
              ordi = ordinate(psIN.na, method=i, distance=distance_metric)
              plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
            }, psIN.na, distance_metric)

            names(plist) <- ord_meths

            pdataframe = ldply(plist, function(x){
              df = x$data[, 1:2]
              colnames(df) = c("Axis_1", "Axis_2")
              return(cbind(df, x$data))
            })
            names(pdataframe)[1] = "method"
            
            pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)
            if (length(shapes) == 1) {
              pdataframe[,shapes] <- factor(pdataframe[,shapes])
              p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar, shape=shapes)) #, fill=PerioStatus_3))
            }

            else{
              p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
            }

            p = p + geom_point(size=1.5, alpha = 3)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) #+ geom_polygon()
            p = p + facet_wrap(~ method, scales="free")+theme_bw() + scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) # open(1), cross(10), closed(2)
            p = p + scale_fill_brewer(type="qual", palette="Dark2")

            if (length(ID) == 1) {
              p = p + geom_text_repel(aes_string(label = ID), size = 2)
            } else {
              p = p 
            }
            
            if (ellipse == "yes" | ellipse == "Yes" ) {
              p = p + stat_ellipse(type = "norm", linetype = 2) 
            } else if (ellipse == "no" | ellipse == "No" ){
              p = p 
            }

            if (length(orders) >= 1) {
              p = p + scale_colour_brewer(type="qual", palette="Dark2", breaks=orders)
            }
            else {
              p = p + scale_colour_brewer(type="qual", palette="Dark2")
            }
            if (length(facet) == 1) {
              ncol <- length(unique(mapping.sel.na.rem[,facet]))
              p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
            }
            else {
              p = p
            }

            #plotlist[[length(plotlist)+1]] <- p
            print(p)
          }
      }
      dev.off()
    }
}
#' A Go_bdiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Go_bdiv()

Go_dist <- function(psIN, project, distance_metrics){
  # out dir
  #out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out)) dir.create(out)
  #out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_path)) dir.create(out_path)
  #out_perm <- file.path(sprintf("%s_%s/table/perm ",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  # run distance
  dm <- list()
  for (distance_metric in distance_metrics) {
    dm[[length(dm)+1]] <- phyloseq::distance(psIN, method=distance_metric)
  }
  
  names(dm) <- distance_metrics
  class(dm)
  
  return(dm)
}
  #' A Go_perm
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Mar 07 2020
#' adjsted 기능을 추가 하였다.
#' 분석 할때 마다 수치가 조금 변하는 것을 수정 하였다.set.seed(1)
#' dm를 따로 분리 하여 시간을 단축 하였고, dm를 다른 방법으로 분석 할수 있게 되었다.
#' Go_perm()


Go_perm <- function(psIN, metadata, project, distance,distance_metrics, adjust,des, name){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s/table",out)) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_perm <- file.path(sprintf("%s/perm",out_path)) 
  if(!file_test("-d", out_perm)) dir.create(out_perm)
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  # Run
  if (length(des) == 1) {
    # Uni
    print(sprintf("#--- Running Paired-PERMANOVA (%s) ---#", des))
  }
  else {
    print("#--- Running Paired-PERMANOVA  ---#")
  }
  
  set.seed(1)
  mapping.sel <-data.frame(sample_data(psIN))
  res.pair <-{}
  for (mvar in rownames(subset(metadata, Go_perm =="yes"))) {
    if (length(unique(as.character(mapping.sel[,mvar]))) == 1) {
      next
    }
    for (distance_metric in distance_metrics) {
      mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
      if (length(unique(mapping.sel.na[,mvar])) == 1)
        next
      
      psIN.sel <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)

      if (length(adjust) >= 1) {
        if(mvar == adjust )
          next
        
      pair.ado <- pairwise.adonis.ad(x=as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar], adjust=adjust, map=mapping.sel.na, mvar=mvar)
      
      
      } else{
        pair.ado <- pairwise.adonis(as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar],map=mapping.sel.na, mvar=mvar)
      }
      
      tmp <- as.data.frame(pair.ado)
      tmp$distance_metric <- distance_metric
      tmp$mvar <- mvar
      tmp$adjusted <- paste(setdiff(adjust, "SampleType"), collapse="+")
      res.pair <- rbind(res.pair, tmp)
    }
  }
  
  
  if (length(adjust) >= 1) {
    if (length(des) == 1) {
      if (length(name) == 1) {
        print(1)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.%s.csv",out_perm, project, des, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.csv",out_perm, project, des, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
    else{
      if (length(name) == 1) {
        print(2)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.csv",out_perm, project,name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.csv",out_perm, project, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  } else{
    if (length(des) == 1) {
      if (length(name) == 1) {
        print(3)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.%s.csv",out_perm, project, des, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.csv",out_perm, project, des, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
    
    else{
      if (length(name) == 1) {
        print(4)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.csv",out_perm, project, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.csv",out_perm, project, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  }

  return(res.pair)
}




Go_dist_plot <- function(psIN, project, distance_metrics, distance, group,orders, name,height,width,plot = TRUE) {

  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  if (plot == TRUE) {
    if (length(name) == 1) {
      pdf(sprintf("%s/5_distplot.%s.%s.%s.pdf", out_path, project, name,format(Sys.Date(), "%y%m%d")),height = height, width=width)
    } else{
      pdf(sprintf("%s/5_distplot.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")),height = height, width=width)
    }
  
  # calc distances
  for (dist in distance_metrics){
    wu <- as.dist(distance[[dist]])
    wu.m <- melt(as.matrix(wu))
    
    # remove self-comparisons
    wu.m <- wu.m %>%
      filter(as.character(Var1) != as.character(Var2)) %>%
      mutate_if(is.factor, as.character)
    
    
    wu.m  = wu.m  %>%
      rowwise() %>%      # for each row
      mutate(Samples = paste(sort(c(Var1, Var2)), collapse = "-")) %>%  # sort the teams alphabetically and then combine them separating with -
      ungroup()
    
    wu.m.sel  = distinct(wu.m, Samples, .keep_all=T)
    
    
    # get sample data (S4 error OK and expected)
    mapping <- data.frame(sample_data(psIN))
    mapping$ID <- as.character(rownames(mapping))
    mapping[,group] <- as.character(mapping[,group])

    
    sd <- mapping %>%
      select("ID", group) %>%
      mutate_if(is.factor,as.character)
    
    # combined distances with sample data
    # sample1
    colnames(sd) <- c("Var1", "Type1")
    wu.m.sel$Var1 <- factor(wu.m.sel$Var1)
    sd$Var1 <- factor(sd$Var1)
    wu.sd <- left_join(wu.m.sel, sd, by = "Var1")
    # sample2
    wu.sd$Var2 <- as.factor(wu.sd$Var2)
    colnames(sd) <- c("Var2", "Type2")
    sd$Var2 <- factor(sd$Var2)
    wu.sd <- left_join(wu.sd, sd, by = "Var2")
    
    wu.sd$Type3 <- ifelse(wu.sd$Type1 == wu.sd$Type2, wu.sd$Type1,"across_group")
    wu.sd$Type3 <- factor(wu.sd$Type3)
    # make a combination for stat
    wu.sd.sel <- subset(wu.sd, Type3 != "across_group")
    wu.sd.sel$Type3 <- factor(wu.sd.sel$Type3)

    # cbn <- combn(x = levels(wu.sd$Type3), m = 2)
    baseline <- "across_group"
    cbn <-{}
    for (x in levels(wu.sd.sel$Type3)){
      cbn <- cbind(cbn, c(baseline, x))
    }
    
    cbn.sel<-cbn[, !duplicated(t(cbn))]
    
    my_comparisons <- {}
    for(i in 1:ncol(cbn.sel)){
      x <- cbn.sel[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    
    # plot
    if (length(orders) >= 1) {
      wu.sd$Type2 <- factor(wu.sd$Type2, levels = orders)
    }       else {
      wu.sd$Type2 <- factor(wu.sd$Type2)
    }
    
    p <- ggplot(wu.sd, aes(x = Type3, y = value, colour=Type3)) + theme_classic()+
      geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, alpha = 0.3,position=position_jitter(0.2)) +
      scale_color_brewer(palette="Dark2") +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      #facet_wrap(~ Type2, scales = "free_x") +
      ylab(dist) + xlab(NULL) +  theme(legend.position="none")+
      stat_compare_means(method= "wilcox.test", label = "p.format", comparisons = my_comparisons, size = 2.5)
    
    if (length(name) == 1) {
      p<- p+ ggtitle(sprintf("%s_%s",dist,name))
    } else{
      p <- p+ ggtitle(sprintf("%s",dist))
    }
    print(p)
  }
    dev.off()
    
  } else {
    return(wu.sd)
  }
}
#' A Go_lmem
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run LMEM
#' @export
#' @examples
#' Go_lmem()

## thresholds for association/permutation tests
# nsamps_threshold <- 0.01 fraction of relabund to call a sample positive
#filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nperm <- 100000


Go_lmem <- function(psIN, metadata, StudyID, project, nsamps_threshold, filt_threshold, taxRanks, data_type, des, name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/lmem",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  ranks <- taxRanks
  taxaname <- ranks
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)

  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    else if (data_type == "other" | data_type == "Other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])
    print(1)

    if (dim(otu.filt)[2] == 2){
      next
    }
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    rownames(agg) <- genera
    dim(agg)
    ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    agg <- agg[intersect(ftk,ftk),]
    # control data set after filter
    if (dim(agg)[1] == 0)
      next

    agg[,taxaname[i]] <- rownames(agg)

    # metatable에서 useForNB는 여러개의 yes가 가능 하지만, useAsConfounder 는 그렇지 않다.
    ## baseline등을 관리 하려면 다음이 필요하다.
    mapping <- data.frame(sample_data(psIN))
    sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
    metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
    mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)

    dim(mapping.sel)


    print(2)
    #--------------    lmer    -------------#
    res <- {}
    for (f in agg[,taxaname[i]]) {
      # clean bacteria name
      if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
        next
      }

      df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)

      df$StudyID <- mapping.sel[df$SampleID, StudyID]
      
      for (cvar in rownames(subset(metadata, Go_lmemConfounder =="yes"))) {
        df[, cvar] <- mapping.sel[df$SampleID, cvar]
      }
      for (mvar in rownames(subset(metadata, Go_lmem =="yes"))) {

        # na remove
        mapping <- data.frame(sample_data(psIN))
        mapping[mapping==""] <- "NA"
        mapping.na <- mapping[!is.na(mapping[,mvar]), ]
        na.count <- length(mapping.na)
        if (length(unique(mapping.na[,mvar])) == 1)
          next

        
        #------------ fix column types------------#
        if (metadata[mvar, "type"] == "factor") {
          mapping.na[,mvar] <- factor(mapping.na[,mvar])
          if (length(unique(mapping.na[,mvar])) ==1 ){
            next
          }
          if (metadata[mvar, "baseline"] != "") {
            mapping.na[,mvar] <- relevel(mapping.na[,mvar], metadata[mvar, "baseline"])
          }
        } else if (metadata[mvar, "type"] == "numeric") {
          mapping.na[,mvar] <- factor(mapping.na[,mvar])
        }
          
        print(3)
        
        
        # na count
        if (length(des) == 1) {
          print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                        des,mvar, dim(mapping.na)[1], dim(mapping)[1]))
        } else{
          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.na)[1], dim(mapping)[1]))
        }
        print(4)

        df[,mvar] <- mapping.na[df$SampleID, mvar]
        form <- as.formula(sprintf("value ~ %s + %s + (1 | StudyID)", mvar, paste(rownames(subset(metadata, Go_lmemConfounder=="yes")), collapse="-")))
        # form <- as.formula(sprintf("value ~ %s +  (1 | StudyID)", mvar))

        mod <- lmer(form, data=df,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
        ## lmer에서 control=은 "number of levels of each grouping~~" 오류가 있을때만 사용한다.
        ##
        # mod2 <- lmer(form, data=df)

        coef <- summary(mod)$coefficients
        coef <- coef[grep(mvar, rownames(coef)),,drop=F]

        res <- rbind(res, cbind(f, mvar, rownames(coef), coef))
        dim(res)
      }
    }


    #-- create table --#
    res <- as.data.frame(res)
    colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue")
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$Estimate <- as.numeric(as.character(res$Estimate))
    res$SE <- as.numeric(as.character(res$SE))
    res$padj <- p.adjust(res$pvalue, method="fdr")
    res <- res[order(res$pvalue),]
    if (length(des) == 1) {
      res$des <- des
    }

    if (length(des) == 1) {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE, col.names = NA, file=sprintf("%s_%s/table/lmem/%s.%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], name, des, project, format(Sys.Date(), "%y%m%d"), sep="/"))
      }
      else {
        write.csv(res, quote = FALSE,col.names = NA, sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i],des, project, format(Sys.Date(), "%y%m%d"), sep="/"))
      }
    }
    else {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], name, project,format(Sys.Date(), "%y%m%d"), sep="/"))
      }
      else{
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], project,format(Sys.Date(), "%y%m%d"), sep="/"))
      }
    }
  }
}
#' A Go_deseq2
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run Deseq2
#' @export
#' @examples
#' Go_deseq2()


#des = des
#alpha = 0.05
Go_deseq2_fs <- function(psIN, metadata, project, adjust, des, name, alpha){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  # map 정리
  mapping <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)


  dim(mapping.sel)
  # start
  res <- {}
  for (mvar in rownames(subset(metadata.sel, useForDESeq2 =="yes"))) {
    if (length(unique(mapping.sel[, mvar])) == 1) {
      next
    }

    #na remove
    mapping.sel <- data.frame(sample_data(psIN))
    mapping.sel[mapping.sel==""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    na.count <- length(mapping.sel.na)
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
    if (length(unique(mapping.sel.na.rem[,mvar])) == 1 )
      next

    if (length(des) == 1) {
      print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                    des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    } else{
      print(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
    }

    if (length(mapping.sel.na.rem[,mvar]) < 4){
      next
      print(sprintf("%s is removed because length(%s) less than 4", mvar, length(mapping.sel.na.rem[,mvar])))
    }



    # integer control
    if (class(mapping.sel.na.rem[,mvar]) == "character"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }

    #-- DESeq2 for phyloseq --#
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    if (length(adjust) == 1) {
      dds = phyloseq_to_deseq2(psIN.na, as.formula(sprintf("~ %s + %s", mvar, adjust)))
      
      print(sprintf("~ %s + %s", mvar, adjust))
    }
    else {
      dds = phyloseq_to_deseq2(psIN.na, as.formula(sprintf("~ %s", mvar)))
      print(sprintf("~ %s", mvar))
    }

    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
    dds = estimateDispersions(dds)
    vst = getVarianceStabilizedData(dds)

    dds = DESeq(dds, fitType="local")
    resultsNames(dds)

    # go back to the default order
    mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    
    for (smvar in levels(mapping.sel[,mvar])) {
      if(smvar == metadata.sel[mvar, "baseline"] | smvar == "" )
        next
      print("pass1")
      basline <- metadata.sel[mvar, "baseline"]
      
      # basline이  존재 하는지 확인 # 200202 ladas 분석중 다음부분이 오류가 발생 하여 일단 사용중지 하였다
      #if (basline != levels(mapping.sel[,mvar]))
       #next
      
        
      print("pass2")
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))

      tmp$dir <- ifelse(tmp$padj < alpha, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      tmp$mvar <- mvar
      tmp$basline<-basline
      tmp$smvar <- smvar
      if (length(des) == 1) {
        tmp$des <- des
      }
      
      
      #-- give taxa name --#
      res = cbind(as(tmp, "data.frame"), as(tax_table(ps.sel)[rownames(tmp), ], "matrix"))
      print("pass3")


      
      #res$ShortName <- paste(res$Phylum,res$Family," ",res$Genus," ",res$Species)
      res$ShortName <- paste(res$Genus," ",res$Species)
      print("pass4")

      #--- give simple name to res---#
      headers <- vector(dim(res)[2], mode="character")
      for (i in 1:dim(res)[1]) {
        headers[i] <- paste("ASV", i, sep="_")
      }
      res$taxa <- headers
      print("pass5")
      
      #-- create table --#
      res <- as.data.frame(res)
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")

      if (length(des) == 1) {
        if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar, name, des, project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        } else {
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.csv",project,  format(Sys.Date(), "%y%m%d"),smvar,mvar,des,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      } else {
        if (length(name) == 1) {
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar,name,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        } else{
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      }
    }
  }
}

#' A Go_deseq2
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run Deseq2
#' @export
#' @examples
#' Go_deseq2()


#des = des
#alpha = 0.05
Go_deseq2 <- function(psIN, metadata, project, adjust, des, name, alpha){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  # map 정리
  mapping <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)


  dim(mapping.sel)
  # start
  res <- {}
  for (mvar in rownames(subset(metadata.sel, Go_deseq2 =="yes"))) {
    if (length(unique(mapping.sel[, mvar])) == 1) {
      next
    }

    #na remove
    mapping.sel <- data.frame(sample_data(psIN))
    mapping.sel[mapping.sel==""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    na.count <- length(mapping.sel.na)
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
    if (length(unique(mapping.sel.na.rem[,mvar])) == 1 )
      next

    if (length(des) == 1) {
      print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                    des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    } else{
      print(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
    }

    if (length(mapping.sel.na.rem[,mvar]) < 4){
      next
      print(sprintf("%s is removed because length(%s) less than 4", mvar, length(mapping.sel.na.rem[,mvar])))
    }



    # integer control
    if (class(mapping.sel.na.rem[,mvar]) == "character"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }

    #-- DESeq2 for phyloseq --#
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    if (length(adjust) >= 1) {
      form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
      print(form)
      dds = phyloseq_to_deseq2(psIN.na, form)
      
    }
    else {
      dds = phyloseq_to_deseq2(psIN.na, as.formula(sprintf("~ %s", mvar)))
      print(sprintf("~ %s", mvar))
    }

    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
    dds = estimateDispersions(dds)
    vst = getVarianceStabilizedData(dds)

    dds = DESeq(dds, fitType="local")
    resultsNames(dds)

    # go back to the default order
    mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    
    for (smvar in levels(mapping.sel[,mvar])) {
      if(smvar == metadata.sel[mvar, "baseline"] | smvar == "" )
        next
      print("pass1")
      basline <- metadata.sel[mvar, "baseline"]
      
      # basline이  존재 하는지 확인 # 200202 ladas 분석중 다음부분이 오류가 발생 하여 일단 사용중지 하였다
      #if (basline != levels(mapping.sel[,mvar]))
       #next
      
        
      print("pass2")
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))

      tmp$dir <- ifelse(tmp$padj < alpha, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      tmp$mvar <- mvar
      tmp$basline<-basline
      tmp$smvar <- smvar
      if (length(des) == 1) {
        tmp$des <- des
      }

      #-- give taxa name --#
      res <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
      print("pass3")
      #res$TaxaName <- paste("p__",res$Phylum,";c__",res$Class,";o__",res$Order,";f__",res$Family,";g__",res$Genus,";s__",res$Species)
      # removing no name and NA
      #for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      for(taxa in c("Phylum","Class","Order","Family","Genus","Species")){
        res[,taxa] == "NA"
        res[,taxa]<- as.character(res[,taxa])
        res[,taxa][is.na(res[,taxa])] <- "__"
        for(i in 1:length(res[,taxa])){
          if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
            res[,taxa][i] <- ""
          }
        }
      }
      print("pass4")
      res$TaxaName <- paste(res$Phylum,"",res$Class,"",res$Order,"",res$Family,"",res$Genus,"",res$Species)
      #res$ShortName <- paste(res$Phylum,res$Family," ",res$Genus," ",res$Species)
      res$ShortName <- paste(res$Genus,"",res$Species)

      # use last taxa name
      for(taxa in c("Family", "Order", "Class","Phylum")){
        for(i in 1:length(res[,taxa])){
          if (res$ShortName[i] != "  "){
            next
          }      else if (res$ShortName[i] == "  " & res[,taxa][i] != ""){
            res$ShortName[i] <- paste(res[,taxa][i])
          }
        }
      }


      #--- give simple name to res---#
      headers <- vector(dim(res)[2], mode="character")
      for (i in 1:dim(res)[1]) {
        headers[i] <- paste("ASV", i, sep="_")
      }
      res$taxa <- headers
      print("pass5")
      #-- create table --#
      res <- as.data.frame(res)
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")

      if (length(des) == 1) {
        if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar,  des, name,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        } else {
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.csv",project,  format(Sys.Date(), "%y%m%d"),smvar,mvar,des,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      } else {
        if (length(name) == 1) {
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar,name,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        } else{
          write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/deseq2/%s.%s.%s.%s.csv",project, format(Sys.Date(), "%y%m%d"),smvar,mvar,project, "Deseq2",format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      }
    }
  }
}

#' A Go_deseq2_volc
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Volcano plot
#' @export
#' @examples
#' Go_deseq2_volc()

dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_deseq2_volc <- function(project, file_path, alpha,beta, name, height, width){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=sprintf(".%s.Deseq2.csv", project));filenames
  sample.names <- sapply(strsplit(filenames, sprintf(".%s.Deseq2.csv", project)), `[`, 1);sample.names
  
  
  print(path)
  print(filenames)
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/6_deseq2.volcano.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } 
  else {
    pdf(sprintf("%s_%s/pdf/6_deseq2.volcano.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, beta, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }

  for (sn in 1:length(sample.names)) {
    file <- list.files(path, pattern = sprintf("^%s%s", sample.names[sn], sprintf(".%s.Deseq2.csv", project)), recursive = T, full.names =T);file
    df <- read.csv(file, row.names=NULL ,check.names=FALSE)
    # remove NA
    df[df==""] <- "NA"
    df$dir <- ifelse(df$padj < alpha & abs(df$log2FoldChange) > beta, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$dir), ]

    basline <- unique(df$basline)
    smvar <- unique(df$smvar)

    dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
    p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue), colour=dir)) +theme_bw() +
      geom_point(alpha=1, size=2) + scale_color_manual(values=dircolors) +
      geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < alpha & abs(df.na$log2FoldChange) > beta, as.character(ShortName),'')), size=2.5) +xlab("log2 fold change") + ylab("-log10 p-value")+
      geom_vline(xintercept = -beta,col = "blue", linetype = "dotted", size = 1) +
      geom_vline(xintercept = beta,col = "red", linetype = "dotted", size = 1) + theme(legend.position = "right", plot.title = element_text(size=10)) 
    
      p1 <- p1 + ggtitle(sprintf("%s %s vs %s (p < %s,cutoff=%s) ", sample.names[sn], basline, smvar,  alpha,beta)) 

    print(p1)
  } 
  dev.off()
}

#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
#dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
dircolors <- c("#4f86f7", "#e10000", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_deseq2_fore <- function(project,file_path, alpha, beta,font, name, height, width){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=sprintf(".%s.Deseq2.csv", project));filenames
  sample.names <- sapply(strsplit(filenames, sprintf(".%s.Deseq2.csv", project)), `[`, 1);sample.names
  
  
  print(path)
  print(filenames)
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/7_deseq2.forest.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } 
  else {
    pdf(sprintf("%s_%s/pdf/7_deseq2.forest.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha,beta, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }

  for (sn in sample.names) {
    file <- list.files(path, pattern = sprintf("^%s%s", sn, sprintf(".%s.Deseq2.csv", project)), full.names =T)
    df <- read.csv(file,row.names=NULL ,check.names=FALSE)

    df.sel <- df
    resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > beta))
   if (dim(resSig)[1] == 0){
     next
   }
    
    resSig$smvar <- factor(resSig$smvar)
    lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
    
    p1 <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=log2FoldChange, color=dir)) + 
      geom_point() + geom_hline(yintercept=0) + coord_flip() + theme_bw() + #theme_classic() +
      geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2)  +  scale_color_manual(values=dircolors) + ylim(c(-lims, lims))+scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s__%s__%s",as.character(resSig$Phylum),as.character(resSig$Family), as.character(resSig$ShortName))) +  xlab("Taxa") + ylab("log2FoldChange")+theme(text = element_text(size=font), plot.title = element_text(hjust=1))
    
      p1<- p1+ ggtitle(sprintf("%s (p < %s, cutoff=%s) ", sn, alpha,beta)) 

    print(p1)
  } 
  dev.off()
}

#' A Go_deseq2_heat
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Heatmap
#' @export
#' @examples
#' Go_deseq2_heat()

Go_deseq2_heat <- function(df, project, data_type, facet,groupby,font, alpha,beta, orders, name, height, width){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  

  
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/8_deseq2.heatmap.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1) {
    pdf(sprintf("%s_%s/pdf/8_deseq2.heatmap.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1 & length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/8_deseq2.heatmap.%s.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/8_deseq2.heatmap.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha,beta, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  resSig <- as.data.frame(subset(df, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
  resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > beta))
  #print("c")
  #if (length(unique(resSig$smvar)) >=2 ){
  if (dim(resSig)[1] >= 1) {
    # re-order
    if (length(orders) >= 1) {
      if (groupby == "smvar"){
        resSig.top$smvar <- factor(resSig.top$smvar)
        if (length(unique(resSig.top$smvar)) <= 1) 
          next
        resSig.top$smvar <- factor(resSig.top$smvar, levels = orders)
      } else{
        resSig.top$des <- factor(resSig.top$des, levels = orders)
      }
    }       else {
      if (groupby == "smvar"){
        if (length(unique(resSig.top$smvar)) <= 1) 
          next
        resSig.top$smvar <- factor(resSig.top$smvar)
      } else{
        if (length(unique(resSig.top$des)) <= 1) 
          next
        resSig.top$des <- factor(resSig.top$des)
      }
    }
    
    print(1)
    if (groupby == "smvar"){
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=smvar, color=smvar)) + theme_classic()+ coord_flip() + #x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
        geom_tile(aes(fill = log2FoldChange), colour = "white") + 
        labs(y = "Comparison Group") +labs(x = NULL) +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "red")+
        ggtitle(sprintf("%s baseline %s vs %s (p < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig$basline), "All groups",  alpha,beta))  + 
        theme(plot.title = element_text(hjust = 0.5))+ #0.5
        theme(legend.position= "right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) 
    }  else {
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=des, color=des)) + theme_classic()+ coord_flip() + #x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
        geom_tile(aes(fill = log2FoldChange), colour = "white") + 
        labs(y = "Comparison Group") +labs(x = NULL) +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "red")+
        ggtitle(sprintf("%s baseline %s vs %s (p < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig$basline), "All groups",  alpha,beta))  + 
        theme(plot.title = element_text(hjust = 0.5))+ #0.5
        theme(legend.position= "right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) 
    }
    
    
    print(2)
    if (data_type == "dada2" | data_type == "DADA2") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$Phylum, resSig$ShortName)))
    } else if (data_type == "Other" | data_type == "other") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$Phylum, resSig$Species)))
    }
    
    
    
    
    print(3)
    if (groupby == "smvar"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"smvar"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "smvar", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  smvar, scales="free_x", ncol = 10)
      }
    }else if (groupby == "des"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"des"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "des", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  des, scales="free_x", ncol = 10)
      }
    }
    #print(4)
    #plotlist[[length(plotlist)+1]] <- p
    p3 = p2 + theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + theme(text = element_text(size=font), plot.title = element_text(hjust=1))
   # print(p3)
  }else{
    next
  }

  
  p4 <- ggplotGrob(p3)
  id <- which(p4$layout$name == "title")
  p4$layout[id, c("l","r")] <- c(1, ncol(p4))
  #grid.newpage()
  grid.draw(p4)
  dev.off()
}
# 200505 ladas genefamily 로제작
# table to deseq2

Go_deseq2_tab <- function(project, tab, map, metadata, name,alpha, height, width){
  
  # package
  if(!'colorspace' %in% installed.packages()){
    install.packages('colorspace')
  }else{
    library('colorspace')
  }
  
  meta <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # map 정리 2
  sel <- intersect(rownames(meta), colnames(map)); head(sel, "3")
  meta.sel <- meta[sel,, drop=F];head(meta.sel)
  map.sel <- map[rownames(map), sel, drop=F];head(map.sel)
  
  # 중요 match map and gene table
  # remove NA row and 반올림 하여 정수 만들기 (integer)
  tab <- tab*10000
  tab.ceiling <- ceiling(tab[-c(99),])
  gene3 <- tab.ceiling+1
  
  
  
  sel <- intersect(rownames(map.sel), colnames(gene3)); head(sel)
  gene4 <- gene3[,sel, drop=F];head(gene4)
  map.sel.sel <- map.sel[sel,, drop=F];head(map.sel.sel)
  
  print(sprintf("%s %s","table", dim(gene4)))
  print(sprintf("%s %s","map", dim(map.sel.sel)))
  
  res <- {}
 
  
  dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
  for (mvar in rownames(subset(meta, Go_deseq2=="yes"))) {

    print(sprintf("Analyzong for %s",mvar))
    
    # map 정리 2
    if (meta[mvar, "type"] == "factor") {
      map[,mvar] <- factor(map[,mvar])
      if (!(is.na(meta[mvar, "baseline"])) && meta[mvar, "baseline"] != "") {
        map[,mvar] <- relevel(map[,mvar], meta[mvar, "baseline"])
      }
    } else if (meta[mvar, "type"] == "numeric") {
      map[,mvar] <- as.numeric(as.character(map.sel[,mvar]))
    } else if (meta[mvar, "type"] == "date") {
      map[,mvar] <- as.Date(sprintf("%06d", map.sel[,mvar]), format="%m%d%y")
      map[,mvar] <- factor(as.character(map[,mvar]), levels=as.character(unique(sort(map.sel[,mvar]))))
    }
    
    # run deseq2
    dds <- DESeqDataSetFromMatrix(countData = gene4,
                                  colData = map.sel.sel,
                                  design= as.formula(sprintf("~ %s", mvar)))
    dds <- DESeq(dds)
    tmp <- results(dds)
    
    print(1)
    for (smvar in levels(map[,mvar])) {
      if(smvar == meta.sel[mvar, "baseline"] | smvar == "" )
        next
      basline <- meta.sel[mvar, "baseline"]
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))
      tmp$dir <- ifelse(tmp$padj < 0.05, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      tmp$mvar <- mvar
      res <- rbind(res, tmp)
      resSig <- as.data.frame(subset(tmp, padj<0.05)); resSig <- resSig[order(resSig$log2FoldChange),]
      # 다음이 왜 있는지 모르겠네. 멈주니까 잠시 지우자 
      # resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa)
      
      # p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + theme(axis.text.y=element_blank())
      #print(p)
      # forest plot of significant results (FDR adjusted for this variable only)
      
      print(2)
      
      lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
      p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange, color=dir)) + geom_point() + geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
      
      
      if(length(name) == 1){
        pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.%s.pdf",out_path, project, mvar,name, format(Sys.Date(), "%y%m%d")), height = height, width=width)
        print(p)
      } else{
        pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.pdf",out_path, project, mvar,format(Sys.Date(), "%y%m%d")), height = height, width=width)
        print(p)
      }
    }
    dev.off()
    
    # heatmap
    print(3)
    ntd <- normTransform(dds)
    assay.ntd <- assay(ntd)
    
    sub.res <- subset(res, padj < alpha)
    sub.res.uniq <- sub.res[!duplicated(sub.res$taxa),]
    # 이름 정리
    sub.res.uniq$taxa <- gsub(".*:_", "", sub.res.uniq$taxa);head(sub.res.uniq$taxa)
    
    assay.ntd.sel <- assay.ntd[rownames(sub.res.uniq),]
    rownames(assay.ntd.sel) <- sub.res.uniq$taxa
    
    dim(assay.ntd.sel)
    
    
    
    ####  my_colours 아 진짜 내가 이걸 해내는 구만.. mvar별로 다른 색 입히기.
    #display.brewer.pal(6, "Set2")
    #display.brewer.pal(8, "Set3")
    mvar <- rownames(subset(meta, Go_deseq2=="yes"))
    my_colours <- list()
    for(i in 1:length(mvar)){
      if (i==1){
        cols <- colorRampPalette(brewer.pal(2, "Set1"))
      } else if (i==2){
        cols <- colorRampPalette(brewer.pal(5, "Paired"))
      }else if (i==3){
        cols <- colorRampPalette(brewer.pal(8, "Set3"))
      }
      
      colours <- cols(length(unique(map[,mvar[i]])))
      names(colours) <- unique(map[,mvar[i]])
      my_colour <- list(x = colours)
      names(my_colour) <- mvar[i]
      my_colours <- append(my_colours,my_colour)
    }
    
    
    df <- as.data.frame(colData(dds)[,rownames(subset(meta, Go_deseq2=="yes"))])
    
    if (length(rownames(subset(meta, Go_deseq2=="yes"))) ==1){
      colnames(df) <- mvar
    }
    
    rownames(df) <- colnames(assay.ntd.sel)
    
    # 한개 이상 return 하기
  
  # multiple list 만들기

  }
  
  functionReturningTwoValues <- function() {
    results <- list()
    results$ntd[[mvar]] <- assay.ntd.sel
    results$cols[[mvar]] <- my_colours
    results$df[[mvar]] <- df
    return(results) 
  }
  cat("\n")
  print("$ntd, $df and $cols are returned.")
  functionReturningTwoValues()
}
#' A Go_deseq2_heat
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Heatmap
#' @export
#' @examples
#' Go_deseq2_heat()

Go_mergeTab <- function(pattern, file_path){
  
   # add input files
  path <- file_path
 
  
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  filenames <- list.files(path, pattern=pattern);filenames
  
  
  cat(sprintf("Files location: %s\n",path))
  cat("=======================================================================\n")
  cat("Merged files:\n")
  cat(sprintf("%s\n",filenames))

  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  return(df)
}
# 200505 ladas genefamily 로제작
# table to deseq2

Go_deseq2_ntd <- function(psIN, project,metadata, name, nsamps_threshold,filt_threshold, taxRanks,data_type,alpha, height, width){
  
  # package
  if(!'colorspace' %in% installed.packages()){
    install.packages('colorspace')
  }else{
    library('colorspace')
  }
  
  ranks <- taxRanks
  taxaname <- ranks
  meta <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # map 정리 2
  map <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(meta), colnames(map)); head(sel, "3")
  meta.sel <- meta[sel,, drop=F];head(meta.sel)
  map.sel <- map[rownames(map), sel, drop=F];head(map.sel)
  

  
  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    else if (data_type == "other" | data_type == "Other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    
    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])
    print(1)
    
    if (dim(otu.filt)[2] == 2){
      next
    }
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    rownames(agg) <- genera
    dim(agg)
    ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    agg <- agg[intersect(ftk,ftk),]
    agg <- agg*10000
    # control data set after filter
    if (dim(agg)[1] == 0)
      next
    
    #agg[,taxaname[i]] <- rownames(agg)
    
    print("complete agg")
    # 중요 match map and gene table
    # remove NA row and 반올림 하여 정수 만들기 (integer)
    tab <- agg
    tab.ceiling <- ceiling(tab[-c(99),])
    gene3 <- tab.ceiling+1
    #gene3 <- tab
    sel <- intersect(rownames(map.sel), colnames(gene3)); head(sel)
    gene4 <- gene3[,sel, drop=F];head(gene4)
    map.sel.sel <- map.sel[sel,, drop=F];head(map.sel.sel)
    
    print(sprintf("%s %s","table", dim(gene4)))
    print(sprintf("%s %s","map", dim(map.sel.sel)))
    
    res <- {}
    
    
    dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
    for (mvar in rownames(subset(meta, Go_deseq2=="yes"))) {
      
      print(sprintf("Analyzing for %s",mvar))
      
      # map 정리 2
      if (meta[mvar, "type"] == "factor") {
        map[,mvar] <- factor(map[,mvar])
        if (!(is.na(meta[mvar, "baseline"])) && meta[mvar, "baseline"] != "") {
          map[,mvar] <- relevel(map[,mvar], meta[mvar, "baseline"])
        }
      } else if (meta[mvar, "type"] == "numeric") {
        map[,mvar] <- as.numeric(as.character(map.sel[,mvar]))
      } else if (meta[mvar, "type"] == "date") {
        map[,mvar] <- as.Date(sprintf("%06d", map.sel[,mvar]), format="%m%d%y")
        map[,mvar] <- factor(as.character(map[,mvar]), levels=as.character(unique(sort(map.sel[,mvar]))))
      }
      
      # run deseq2
      dds <- DESeqDataSetFromMatrix(countData = gene4,
                                    colData = map.sel.sel,
                                    design= as.formula(sprintf("~ %s", mvar)))
      dds <- DESeq(dds)
      tmp <- results(dds)
      
      print("run deseq2")
      for (smvar in levels(map[,mvar])) {
        if(smvar == meta.sel[mvar, "baseline"] | smvar == "" )
          next
        basline <- meta.sel[mvar, "baseline"]
        tmp <- results(dds, contrast = c(mvar, smvar, basline))
        
        tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
          tmp <- unlist(strsplit(x, ";"))
          tmp[length(tmp)]
        }))
        tmp$dir <- ifelse(tmp$padj < 0.05, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
        tmp$mvar <- mvar
        res <- rbind(res, tmp)
        resSig <- as.data.frame(subset(tmp, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
        # 다음이 왜 있는지 모르겠네. 멈주니까 잠시 지우자 
        # resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa)
        
        # p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + theme(axis.text.y=element_blank())
        #print(p)
        # forest plot of significant results (FDR adjusted for this variable only)
        
        print(2)
        
        lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
        p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange, color=dir)) + geom_point() + geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
        
        
        if(length(name) == 1){
          pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.%s.pdf",out_path, project, mvar,name, format(Sys.Date(), "%y%m%d")), height = height, width=width)
          print(p)
        } else{
          pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.pdf",out_path, project, mvar,format(Sys.Date(), "%y%m%d")), height = height, width=width)
          print(p)
        }
      }
      dev.off()
      
      # heatmap
      print(3)
      ntd <- normTransform(dds)
      assay.ntd <- assay(ntd)
      
      sub.res <- subset(res, padj < alpha)
      sub.res.uniq <- sub.res[!duplicated(sub.res$taxa),]
      # 이름 정리
      sub.res.uniq$taxa <- gsub(".*:_", "", sub.res.uniq$taxa);head(sub.res.uniq$taxa)
      
      assay.ntd.sel <- assay.ntd[rownames(sub.res.uniq),]
      rownames(assay.ntd.sel) <- sub.res.uniq$taxa
      
      dim(assay.ntd.sel)
      
      
      
      ####  my_colours 아 진짜 내가 이걸 해내는 구만.. mvar별로 다른 색 입히기.
      #display.brewer.pal(6, "Set2")
      #display.brewer.pal(8, "Set3")
      mvar <- rownames(subset(meta, Go_deseq2=="yes"))
      my_colours <- list()
      for(i in 1:length(mvar)){
        if (i==1){
          cols <- colorRampPalette(brewer.pal(2, "Set1"))
        } else if (i==2){
          cols <- colorRampPalette(brewer.pal(5, "Paired"))
        }else if (i==3){
          cols <- colorRampPalette(brewer.pal(8, "Set3"))
        }
        
        colours <- cols(length(unique(map[,mvar[i]])))
        names(colours) <- unique(map[,mvar[i]])
        my_colour <- list(x = colours)
        names(my_colour) <- mvar[i]
        my_colours <- append(my_colours,my_colour)
      }
      
      
      df <- as.data.frame(colData(dds)[,rownames(subset(meta, Go_deseq2=="yes"))])

      if (length(rownames(subset(meta, Go_deseq2=="yes"))) ==1){
        colnames(df) <- mvar
      }
      rownames(df) <- colnames(assay.ntd.sel)
      # 한개 이상 return 하기
      # multiple list 만들기
    }
  }
  functionReturningTwoValues <- function() {
    results <- list()
    results$ntd <- assay.ntd.sel
    results$cols <- my_colours
    results$df <- df
    return(results) 
  }
  cat("\n")
  print("$ntd, $df and $cols are returned.")
  functionReturningTwoValues()
}

# alternative
#functionReturningTwoValues <- function() {
#  results <- list()
#  results$ntd[[mvar]] <- assay.ntd.sel
#  results$cols[[mvar]] <- my_colours
#  results$df[[mvar]] <- df
#  return(results) 
#}
#' A Go_zinb
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run LMEM
#' @export
#' @examples
#' Go_lmem()

## thresholds for association/permutation tests
# nsamps_threshold <- 0.01 fraction of relabund to call a sample positive
#filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing



nperm <- 100000
filt_threshold= 0.1 

Go_zinb <- function(psIN, metadata, StudyID, project, ranks, alpha, nsamps_threshold, taxRanks, data_type,name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/zinb",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  ranks <- taxRanks
  taxaname <- ranks
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)

  # compute GMPR normalization
  sizeFactors <- GMPR(as.matrix(as.data.frame(otu_table(psIN))), intersect.no=9)
  # fix column types
  mapping <- data.frame(sample_data(psIN))
  for (mvar in  rownames(subset(metadata, Go_zinb=="yes" | Go_zinbConfounder=="yes"))) {
    if (metadata[mvar, "type"] == "factor") {
      mapping[,mvar] <- factor(mapping[,mvar])
      if (!(is.na(metadata[mvar, "baseline"])) && metadata[mvar, "baseline"] != "") {
        mapping[,mvar] <- relevel(mapping[,mvar], metadata[mvar, "baseline"])
      }
    } else if (metadata[mvar, "type"] == "numeric") {
      mapping[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
    } else if (metadata[mvar, "type"] == "date") {
      mapping[,mvar] <- as.Date(sprintf("%06d", mapping.sel[,mvar]), format="%m%d%y")
      mapping[,mvar] <- factor(as.character(mapping[,mvar]), levels=as.character(unique(sort(mapping.sel[,mvar]))))
    }
  }
  

  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    else if (data_type == "other" | data_type == "Other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])

    if (dim(otu.filt)[2] == 2){
      next
    }
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    rownames(agg) <- genera
    dim(agg)
    ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    agg <- agg[intersect(ftk,ftk),]
    # control data set after filter
    if (dim(agg)[1] == 0)
      next

    agg[,taxaname[i]] <- rownames(agg)

    # metatable에서 useForNB는 여러개의 yes가 가능 하지만, useAsConfounder 는 그렇지 않다.
    ## baseline등을 관리 하려면 다음이 필요하다.
    mapping <- data.frame(sample_data(psIN))
    sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
    metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
    mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)

    dim(mapping.sel)

    #-------------------------------------------------------------------#
    #--------------    ZINB+GMPR, emmeans for estimates    -------------#
    #-------------------------------------------------------------------#
    res <- {}; res.interaction <- {}
    for (f in agg[,taxaname[i]]) {
      # clean bacteria name
      if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
        next
      }
      
      df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)

      nsamps_detected <- length(which(df$value>=nsamps_threshold))
      
      for (con in rownames(subset(metadata, Go_zinb=="yes"))) {
        df[, con] <- mapping.sel[df$SampleID, con]
      }
      for (var in rownames(subset(metadata, Go_zinbConfounder=="yes"))) {
        df[, var] <- mapping.sel[df$SampleID, var]
      }
      
      df$size_factor <- log(sizeFactors[df$SampleID])
      for (mvar in rownames(subset(metadata, Go_zinb=="yes"))) {
        tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s * %s + offset(size_factor) | 1",mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*"))), data = df, dist = "negbin", EM = F, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
        
        if (class(tt) == "zeroinfl") {
          coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
          coef <- coef[grep(".", rownames(coef)),,drop=F] #:
          colnames(coef) <- c("Estimate", "SE", "t", "pval")
          res.interaction <- rbind(res.interaction, cbind(f, nsamps_detected, mvar, "ZINB", rownames(coef), coef))
          
          if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
            form <- as.formula(sprintf("pairwise ~ %s | %s", mvar, mvar))
            emm <- emmeans(m, form, adjust="none")
            coef <- as.data.frame(emm$contrasts)
            coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
            res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
          }else{
            form <- as.formula(sprintf("pairwise ~ %s | %s", mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")))
            emm <- emmeans(m, form, adjust="none")
            coef <- as.data.frame(emm$contrasts)
            coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
            #coef$confounder <- paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")
            res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
          }
          
        } else if (class(tt) == "try-error") {
          tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s * %s + offset(size_factor)",mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*"))), data = df), silent=T)
          
          if (class(tt)[1] == "negbin") {
            coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
            coef <- coef[grep(".", rownames(coef)),,drop=F]#:
            colnames(coef) <- c("Estimate", "SE", "t", "pval")
            res.interaction <- rbind(res.interaction, cbind(f, nsamps_detected, mvar, "NB", rownames(coef), coef))
            
            print(1)
            
            if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
              form <- as.formula(sprintf("pairwise ~ %s | %s", mvar, mvar))
              emm <- emmeans(m, form, adjust="none")
              coef <- as.data.frame(emm$contrasts)
              coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
              res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
            }else{
              form <- as.formula(sprintf("pairwise ~ %s | %s", mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")))
              emm <- emmeans(m, form, adjust="none")
              coef <- as.data.frame(emm$contrasts)
              coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
              #coef$confounder <- paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")
              res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
              }
          }
        }
      }
    }
    print(1)
    
    #form <- as.formula(sprintf("pairwise ~ %s", mvar))
    #emm <- emmeans(m, form, adjust="none")
    
    #-- tidy up results and create table --#
    if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "Estimate", "SE", "df", "t", "pval")
    } else if (mvar == paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")){
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "Estimate", "SE", "df", "t", "pval")
    }    else{
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "condition", "Estimate", "SE", "df", "t", "pval")
    }


    res$padj <- p.adjust(res$pval, method="fdr")
    res <- res[order(res$pval, decreasing=F),]
    print(2)
    
    res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
    res.interaction <- as.data.frame(res.interaction)
    colnames(res.interaction) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "coefficient", "Estimate", "SE", "t", "pval"); rownames(res.interaction) <- {}
    res.interaction$Estimate <- as.numeric(as.character(res.interaction$Estimate))
    res.interaction$SE <- as.numeric(as.character(res.interaction$SE))
    res.interaction$pval <- as.numeric(as.character(res.interaction$pval))
    res.interaction$padj <- p.adjust(res.interaction$pval, method="fdr")
    res.interaction$sigstr <- ifelse(res.interaction$padj < alpha,"*", "")
  }
  
  
  
  
  if (length(name) == 1) {
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.%s.zinb.res.%s.csv",out_table, name, project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
    write.csv(res.interaction, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.%s.zinb.res.interaction.%s.csv",out_table, name, project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
  }else{
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.zinb.res.%s.csv",out_table,  project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
    write.csv(res.interaction, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.zinb.res.interaction.%s.csv", out_table,  project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
  }
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$res <- res
    results$interaction <-res.interaction
    return(results) 
  }
  cat("\n")
  print("$res and $interaction are returned.")
  functionReturningTwoValues()
}

#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
#dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
dircolors <- c("#4f86f7", "#e10000", "grey"); names(dircolors) <- c("down", "up", "NS")
Go_lmem_fore <- function(project,file_path, alpha, pattern, name, order, height, width){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  

  print(path)
  print(filenames)

  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/9_lmem.forest.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/9_lmem.forest.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }

  
  
  for (sn in sample.names) {
    file <- list.files(path, pattern = sprintf("%s%s", sn, sprintf("%s", pattern)), full.names =T)
    df <- read.csv(file,row.names=NULL ,check.names=FALSE)
    df.sel <- df
    resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$Estimate),]
    # resSig$smvar <- factor(resSig$smvar)
    print(1)
    if (dim(resSig)[1] == 0)
      next

    print(2)
    resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)==1, "up", "down"), "NS")
    print(3)
    
    # 중복 이름 처리 하기
    headers <- vector(dim(resSig)[1], mode="character")
    
    for (i in 1:dim(resSig)[1]) {
      headers[i] <- paste("ASV", i, sep="_")
    }
    resSig$ASV <- headers
    
    
    for (plot in unique(resSig$metadata)){
      resSig.sel <- subset(resSig, metadata == plot)
      if (length(unique(resSig.sel$coefficient)) ==1 ){
        
        lims <- max(abs(resSig.sel$Estimate) + abs(resSig.sel$SE))*1.0
        p1 <- ggplot(resSig.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
          geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
          geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
          scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel$ASV), labels = resSig.sel$taxa)
        
        p1 <- p1+ ggtitle(sprintf("LMEM-%s (%s p < %s) ",plot, sn, alpha)) +
          labs(y = "Estimate") +labs(x = NULL)
        
        print(p1)
      } else if (length(unique(resSig.sel$coefficient)) >=1){
        for (cof in unique(resSig.sel$coefficient)){
          resSig.sel.sel <- subset(resSig.sel, coefficient == cof)
          lims <- max(abs(resSig.sel.sel$Estimate) + abs(resSig.sel.sel$SE))*1.0
          p1 <- ggplot(resSig.sel.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
            geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
            geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
            scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel.sel$ASV), labels = resSig.sel.sel$taxa)
          
          p1 <- p1+ ggtitle(sprintf("LMEM-%s (%s p < %s) ",cof, sn, alpha)) +
            labs(y = "Estimate") +labs(x = NULL)
          print(p1)
        }
      }
    }
  }
  dev.off()
}





#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_lmem_heat <- function(project,file_path, alpha, pattern, facet, name, orders, height, width){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  
  
  print(path)
  print(sample.names)
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1 & length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  
  df.sel <- df
  resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$Estimate),]
  # resSig$smvar <- factor(resSig$smvar)
  
  
  resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)==1, "up", "down"), "NS")
  print(1)
  for (plot in unique(resSig$metadata)){
    resSig.sel <- subset(resSig, metadata == plot)
    print(2)
    if (length(unique(resSig.sel$coefficient)) >=1 ){
      resSig.sel$comparison <-  gsub(sprintf("%s", unique(resSig.sel$metadata)) ,"" , resSig.sel$coefficient)
      resSig.sel$comparison <- factor(resSig.sel$comparison , levels = orders)
      resSig.sel$stars <- cut(resSig.sel$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
      print(3)
      if (length(facet) == 1) {
        resSig.sel[,facet] <- factor(resSig.sel[,facet] , levels = orders)
      }
      print(4)
      p <- ggplot(resSig.sel, aes(x=reorder(taxa,Estimate), y=comparison, color=comparison)) + #, alpha=padj)) +
        theme_classic()+ coord_flip() + geom_tile(aes(fill = Estimate), colour = "white") + 
        geom_text(aes(label=stars), color = "black") + scale_fill_gradient2(low = "darkblue", mid = "white", high = "red") + 
        ggtitle(sprintf("LMEM All comparison group (%s p < %s) ", plot, alpha)) +
        theme(plot.title = element_text(hjust = 0.5))+ #0.5
        theme(legend.position= "right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + labs(y = "comparison group") +labs(x = NULL)
      p = p + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
      
      
      if (length(facet) == 1) {
        print(5)
        ncol <- length(unique(resSig.sel[,facet]))*length(unique(resSig.sel[,"comparison"]))
        p = p + facet_wrap(as.formula(sprintf("~ %s+%s", facet, "comparison")), scales="free_x", ncol = ncol)
      }
      else {
        p = p + facet_wrap(~  comparison, scales="free_x", ncol = 10) 
      }
      
      #plotlist[[length(plotlist)+1]] <- p
      print(p)
    }
  }  
  dev.off()
}
#' A Go_lmem
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run LMEM
#' @export
#' @examples
#' Go_lmem()

## thresholds for association/permutation tests
# nsamps_threshold <- 0.01 fraction of relabund to call a sample positive
#filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nperm <- 100000


Go_lmem_species <- function(psIN, metadata, StudyID, project, adjust, nsamps_threshold, filt_threshold,  data_type, des, name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/lmem",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  taxRanks <- c("Genus", "Species")
  ranks <- taxRanks
  taxaname <- ranks
  
  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
    # dada2 or nephele
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN)))
  }
  else if (data_type == "Nephele" | data_type == "nephele") {
    otu.filt <- as.data.frame(otu_table(psIN))
  }
  else if (data_type == "other" | data_type == "Other") {
    otu.filt <- as.data.frame(otu_table(psIN))
  }
  
  # continue
  otu.filt[,"Genus"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level="Genus")
  otu.filt[,"Species"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level="Species")
  otu.filt$Genus.Species <- paste(otu.filt$Genus," ",otu.filt$Species)
  otu.filt.sel <- otu.filt
  otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$Genus), ]
  otu.filt.sel$Genus  <- NULL
  otu.filt.sel$Species <- NULL
  
  
  if (dim(otu.filt)[2] == 2){
    next
  }
  agg <- aggregate(as.formula(sprintf(". ~ %s" , "Genus.Species")), otu.filt.sel, sum, na.action=na.pass)
  genera <- agg[,"Genus.Species"]
  agg <- agg[,-1]
  
  agg <- normalizeByCols(agg)
  
  rownames(agg) <- genera
  dim(agg)
  ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
  
  agg <- agg[intersect(ftk,ftk),]
  # control data set after filter
  if (dim(agg)[1] == 0)
    next
  
  agg[,"Genus.Species"] <- rownames(agg)
  
  # metatable에서 useForNB는 여러개의 yes가 가능 하지만, useAsConfounder 는 그렇지 않다.
  ## baseline등을 관리 하려면 다음이 필요하다.
  mapping <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)
  
  dim(mapping.sel)
  
  
  
  #--------------    lmer    -------------#
  res <- {}
  for (f in agg[,"Genus.Species"]) {
    # clean bacteria name
    if (f == "g__   s__" || f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
      next
    }
    
    
    df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
    
    df$StudyID <- mapping.sel[df$SampleID, StudyID]
    
    for (cvar in rownames(subset(metadata, useAsConfounder=="yes"))) {
      df[, cvar] <- mapping.sel[df$SampleID, cvar]
    }
    for (mvar in rownames(subset(metadata, useForLmem=="yes"))) {
      
      # na remove
      mapping <- data.frame(sample_data(psIN))
      mapping[mapping==""] <- "NA"
      mapping.na <- mapping[!is.na(mapping[,mvar]), ]
      na.count <- length(mapping.na)
      if (length(unique(mapping.na[,mvar])) == 1)
        next
      
      
      #------------ fix column types------------#
      if (metadata[mvar, "type"] == "factor") {
        mapping.na[,mvar] <- factor(mapping.na[,mvar])
        if (length(unique(mapping.na[,mvar])) ==1 ){
          next
        }
        if (metadata[mvar, "baseline"] != "") {
          mapping.na[,mvar] <- relevel(mapping.na[,mvar], metadata[mvar, "baseline"])
        }
      } else if (metadata[mvar, "type"] == "numeric") {
        mapping.na[,mvar] <- factor(mapping.na[,mvar])
      }
      
      # na count
      if (length(des) == 1) {
        print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                      des,mvar, dim(mapping.na)[1], dim(mapping)[1]))
      } else{
        print(sprintf("##-- %s (total without NA: %s/%s) --##",
                      mvar, dim(mapping.na)[1], dim(mapping)[1]))
      }
      
      
      df[,mvar] <- mapping.na[df$SampleID, mvar]
      # 공식에 -1을 추가 하여 intercept 를 삭제 하였다. 20200406
      
      if (length(adjust) >= 1) {
        form <-as.formula(sprintf("value ~ %s + %s + (1 | StudyID)", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
        print(form)
        
      }
      else {
        form <- as.formula(sprintf("value ~ %s + (1 | StudyID)", mvar))
        print(form)
      }
      
      # form <- as.formula(sprintf("value ~ %s +  (1 | StudyID)", mvar))
      
      mod <- lmer(form, data=df,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      
      ## lmer에서 control=은 "number of levels of each grouping~~" 오류가 있을때만 사용한다.
      ##
      # mod2 <- lmer(form, data=df)
      
      coef <- summary(mod)$coefficients
      # remove intercept근데 작동을 안할때가 있어서 일단 불활성 하였다.
      coef <- coef[grep(mvar, rownames(coef)),,drop=F]
      
      res <- rbind(res, cbind(f, mvar, rownames(coef), coef))
      dim(res)
    }
  }

  
  #-- create table --#
  res <- as.data.frame(res)
  colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue")
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$Estimate <- as.numeric(as.character(res$Estimate))
  res$SE <- as.numeric(as.character(res$SE))
  res$padj <- p.adjust(res$pvalue, method="fdr")
  res <- res[order(res$pvalue),]
  if (length(des) == 1) {
    res$des <- des
  }
  
  if (length(des) == 1) {
    if (length(name) == 1) {
      write.csv(res, quote = FALSE, col.names = NA, file=sprintf("%s_%s/table/lmem/%s.%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"), des, project, "Species", name, format(Sys.Date(), "%y%m%d"), sep="/"))
    }
    else {
      write.csv(res, quote = FALSE,col.names = NA, sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"), des, project, "Species", format(Sys.Date(), "%y%m%d"), sep="/"))
    }
  }
  else {
    if (length(name) == 1) {
      write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"), project,"Species", name, format(Sys.Date(), "%y%m%d"), sep="/"))
    }
    else{
      write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  project,"Species", format(Sys.Date(), "%y%m%d"), sep="/"))
    }
  }
    
}
#' A Go_rpkm
#'
#' getting RPKM table
#' @param love getting RPKM table
#' @keywords rpkm
#' @export
#' @examples
#' Go_rpkm


Go_rpkmToPs <- function(project, speciesTab, genoemSizeTab, taxaTab) {
  # out path
  tem_path <- file.path(sprintf("RPKM_tem_%s", format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", tem_path)) dir.create(tem_path)
  rds_path <- file.path(sprintf("RPKM_rds_%s", format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", rds_path)) dir.create(rds_path)
  

  bacTab <- speciesTab
  
  #-- get simple species and genus names --#  species 이름이 길다.
  bacTab$Species <- rownames(bacTab)
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(bacTab$Species),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab$Genus <- cleanSpecies$X1
  bacTab$Species <- cleanSpecies$Species
  bacTab$SubSpecies <- rownames(bacTab)
  
  # read genome size table
  gSizeTab <- genoemSizeTab
  colnames(gSizeTab) <- c("AccessionVersion", "TaxId", "Completeness", "GenomeLen", "Title")
  cleanSpecies_gSizeTab <- data.frame(do.call('rbind', strsplit(as.character(gSizeTab$Title),' ',fixed=TRUE)))
  cleanSpecies_gSizeTab$Species <- paste(cleanSpecies_gSizeTab$X1, cleanSpecies_gSizeTab$X2, sep="_")
  cleanSpecies_gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Species <- cleanSpecies_gSizeTab$Species
  
  
  # create new data frame and exrtact genome size for calculation
  genomeFinal <- data.frame(cbind(bacTab$Species, bacTab$Genus)) # 두번째 col은 이름이 안따라온다. 다음 명령어를 사용해야 한다.
  colnames(genomeFinal) <- c("Species","Genus")
  genomeFinal$Genus <- bacTab$Genus
  genomeFinal$GenomeLen <- gSizeTab$GenomeLen[match(bacTab$Genus, gSizeTab$Genus)] # index match
  
  
  write.csv(genomeFinal, quote = FALSE, #col.names = NA, #row.names = FALSE,
            file=sprintf("%s/%s.genomeFinal2.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  bacTab$Genus <- NULL
  bacTab$Species <- NULL
  bacTab$SubSpecies <- NULL
  ## RPKM
  geneLength <- 1
  bacNor <- data.frame(sapply(bacTab , function(column) 10^9 * column / geneLength / sum(column)))
  rownames(bacNor) <- rownames(bacTab); head(rownames(bacNor))
  
  #brac.rpkmS <- bracNor*100/genomeFinal$SlenS
  #brac.rpkmS <- ceiling(brac.rpkmS[-c(99),]) # 반올림
  #brac.rpkmS[is.na(brac.rpkmS)] <- 0 # remove NA
  
  bacTab.rpkmG <- bacNor*100/genomeFinal$GenomeLen
  bacTab.rpkmG <- ceiling(bacTab.rpkmG[-c(99),]) # 반올림
  bacTab.rpkmG[is.na(bacTab.rpkmG)] <- 0 # remove NA
  
  ## add LCA
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(rownames(bacTab.rpkmG)),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab.rpkmG$Species <- cleanSpecies$Species
  bacTab.rpkmG$Genus <- cleanSpecies$X1
  
  Taxa <- taxaTab
  
  ranks <- c("Rank1", "Phylum", "Class", "Order", "Family")
  for(taxa in ranks){
    bacTab.rpkmG[,taxa]  <- Taxa[,taxa][match(bacTab.rpkmG$Genus, Taxa$Genus)]
  }
  
  bacTab.rpkmG.sel <- subset(bacTab.rpkmG, Rank1 == "Bacteria")
  colnames(bacTab.rpkmG.sel) <-  gsub("X", "", colnames(bacTab.rpkmG.sel));head(colnames(bacTab.rpkmG.sel))
  
  ##########################
  #---      Make ps     ---#
  ##########################
  #------------- tax table -------------#
  headers <- rownames(bacTab.rpkmG.sel);head(headers)
  
  tax <- data.frame(bacTab.rpkmG.sel$Rank1, bacTab.rpkmG.sel$Phylum, bacTab.rpkmG.sel$Class, bacTab.rpkmG.sel$Order, bacTab.rpkmG.sel$Family, bacTab.rpkmG.sel$Genus, bacTab.rpkmG.sel$Species)
  rownames(tax) <- headers
  colnames(tax) <- c("Kingdom","Phylum","Class", "Order", "Family", "Genus", "Species") ;tax
  
  #------------- otu table -------------#
  for(i in colnames(tax)){
    bacTab.rpkmG.sel$Rank1 <- NULL
    bacTab.rpkmG.sel[,i] <- NULL
  }
  
  
  otu <- bacTab.rpkmG.sel
  rownames(otu) <- headers
  otu <- otu[rowSums(otu[, -1])>0, ]
  
  write.csv(otu, quote = FALSE, col.names = NA, #row.names = FALSE,
            file=sprintf("%s/otu_%s.genome_size.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  #--- create phyloseq file ---#
  tax<-as.matrix(tax)
  otu<-as.matrix(otu)
  
  OTU <- otu_table(otu, taxa_are_rows = TRUE);head(OTU)
  TAX <- tax_table(tax);head(TAX)
  ps1 <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), tax_table(TAX))
  saveRDS(ps1, sprintf("%s/ps_rpkmG.%s.%s.rds",rds_path, project, format(Sys.Date(), "%y%m%d")))
  
  return(ps1)
}




#' A Go_pickTaxa
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity ordination plot
#' @export
#' @examples
#' Go_pickTaxa()

Go_pickTaxa <- function(psIN, project, TaxaLevel, data_type, bestPick){
  # ------------- Aggregate
  cat("#--  Getting Taxa  --#\n")
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN))) # for dada2 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame((otu_table(psIN))) # not for dada2 
  }
  
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN),taxRanks=ranks, level= TaxaLevel)
  agg <- aggregate(. ~ Genus, otu.filt, sum, na.action=na.pass) #add na.action=na.pass if have error "no rows to aggregate"
  genera <- agg$Genus
  agg <- agg[,-1]
  rownames(agg) <- genera
  agg_t <-t(agg)
  dim(agg)
  
  
  # ------------- biplot Genus 선택 
  

  #-----------------Parameters----------------#
  cmethod<-"spearman" 
  #Correlation method to use: pearson, spearman, kendall
  fmethod<-"bray" 
  #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  vmethod<-"bray" 
  #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  nmethod<-"bray" 
  #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao jaccard
  
  best5 <- bestPick #100
  
  if (data_type == "dada2" | data_type == "DADA2") {
    res.bv.step.biobio <- bv.step(wisconsin(agg_t), wisconsin(agg_t), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    res.bv.step.biobio <- bv.step(wisconsin(agg), wisconsin(agg), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  }
  
  taxaNames<-colnames(agg_t);taxaNames
  bestTaxaFit <- ""
  for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl))){
    bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
  }
  
  bestTaxaFit <- data.frame(bestTaxaFit);bestTaxaFit 
  colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"
  
  bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")));bio.keep
  
  taxaNames <- colnames(agg_t);taxaNames
  
  cat("\n Use this numbers for biplot;  \n")

  
  print(bio.keep)
  for (taxa in bio.keep){
    cat(sprintf("%s : %s \n", taxa, taxaNames[taxa]))
  }
  return(taxaNames)
}
#' A Go_biplot
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity ordination plot
#' @export
#' @examples
#' Go_biplot()

Go_biplot <- function(psIN, metadata, project, orders, distance_metrics, data_type, biplot, shapes, TaxaLevel, ID, facet, des, name, height, width){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)
  
  cat("#--  Getting Taxa  --#")
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN))) # for dada2 
  }
  else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame((otu_table(psIN))) # not for dada2 
  }
  

  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN),taxRanks=ranks, level= TaxaLevel)
  agg <- aggregate(. ~ Genus, otu.filt, sum, na.action=na.pass) #add na.action=na.pass if have error "no rows to aggregate"
  genera <- agg$Genus
  agg <- agg[,-1]
  rownames(agg) <- genera
  
  agg_t <-t(agg)
  
  dim(agg)
  
  
  
  # out file
  if (length(des) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1 & length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }


  # out file
    mapping.sel <- data.frame(sample_data(psIN))

    for (mvar in rownames(subset(metadata, useForBeta =="yes"))) {
        for(distance_metric in distance_metrics){
          #na remove
          mapping.sel <- data.frame(sample_data(psIN))
          mapping.sel[mapping.sel==""] <- "NA"
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          na.count <- length(mapping.sel.na)
          psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))

          if (length(des) == 1) {
            print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                          des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          } else{
            print(sprintf("##-- %s (total without NA: %s/%s) --##",
                          mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          }


          if (class(mapping.sel.na.rem[,mvar]) == "integer"){
            mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
            sample_data(psIN.na) <- mapping.sel.na.rem
          }

          
          
          ordi <- ordinate(psIN , method = "NMDS", distance = distance_metric)
          #Get site information
          df <- scores(ordi,display=c("sites"));head(df)
          
          # Add grouping information
          df <- merge(df, mapping.sel.na.rem, by="row.names")
          
          #Get the vectors for bioenv.fit
          bio.fit <- envfit(ordi, agg_t[,biplot,drop=F],  perm = 999); head(bio.fit)
          df_biofit <- scores(bio.fit, display=c("vectors"))
          df_biofit <- df_biofit*vegan:::ordiArrowMul(df_biofit)
          df_biofit <- as.data.frame(df_biofit);df_biofit
          

          if (length(shapes) == 1) {
            p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar, shape=shapes), size=2) + 
              stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                           inherit.aes = TRUE) 
          }
          else{
            p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar), size=2) +
              stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                           inherit.aes = TRUE) 
          }
          
          p = p + theme_bw() +  ggtitle(sprintf("%s - %s - %s - %s", project, mvar, "NMDS", distance_metric)) #+ geom_polygon()
          # open(1), cross(10), closed(2)
          p = p + scale_fill_brewer(type="qual", palette="Set1")#scale_colour_manual(values=colours) 
          p = p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                         arrow = arrow(length = unit(0.2, "cm")), color="#808080",alpha=0.7)+
            geom_text(data=as.data.frame(df_biofit*1.1), aes(NMDS1, NMDS2, label = rownames(df_biofit)), 
                      color="#808080",alpha=0.7, size = 3)


          if (length(ID) == 1) {
            p = p +  geom_text_repel(aes_string(label = ID), size = 2)
          }

          else {
            p = p 
          }

          if (length(orders) >= 1) {
            p = p + scale_colour_brewer(type="qual", palette="Set1", breaks=orders) + 
              scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) 
          }
          else {
            p = p + scale_colour_brewer(type="qual", palette="Set1") + 
              scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16)) 
          }
          
          if (length(facet) == 1) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }
          else {
            p = p
          }
          print(p)
        }
      }
    dev.off()
  }

    else {
      if (length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1 & length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      
      for (mvar in rownames(subset(metadata, useForBeta =="yes"))) {
          for(distance_metric in distance_metrics){
            #na remove
            mapping.sel <- data.frame(sample_data(psIN))
            mapping.sel[mapping.sel==""] <- "NA"
            mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
            na.count <- length(mapping.sel.na)
            psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
            mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))


            if (length(des) == 1) {
              print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                            des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            } else{
              print(sprintf("##-- %s (total without NA: %s/%s) --##",
                            mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            }


            if (class(mapping.sel.na.rem[,mvar]) == "integer"){
              mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
              sample_data(psIN.na) <- mapping.sel.na.rem
            }

            ordi <- ordinate(psIN , method = "NMDS", distance = distance_metric)
            #Get site information
            df <- scores(ordi,display=c("sites"));head(df)
            
            # Add grouping information
            df <- merge(df, mapping.sel.na.rem, by="row.names")
            #Get the vectors for bioenv.fit
            bio.fit <- envfit(ordi, agg_t[,biplot,drop=F],  perm = 999); head(bio.fit)
            df_biofit <- scores(bio.fit, display=c("vectors"))
            df_biofit <- df_biofit*vegan:::ordiArrowMul(df_biofit)
            df_biofit <- as.data.frame(df_biofit);df_biofit
            
            
            if (length(shapes) == 1) {
              p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar, shape=shapes), size=2) + 
                stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                             inherit.aes = TRUE) 
            }
            else{
              p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar), size=2) +
                stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                             inherit.aes = TRUE) 
            }
            
            p = p + theme_bw() +  ggtitle(sprintf("%s - %s - %s - %s", project, mvar, "NMDS", distance_metric)) #+ geom_polygon()
            # open(1), cross(10), closed(2)
            p = p + scale_fill_brewer(type="qual", palette="Set1")#scale_colour_manual(values=colours) 
            p = p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                               arrow = arrow(length = unit(0.2, "cm")), color="#808080",alpha=0.7)+
              geom_text(data=as.data.frame(df_biofit*1.1), aes(NMDS1, NMDS2, label = rownames(df_biofit)), 
                        color="#808080",alpha=0.7, size = 3)
            
            
            if (length(ID) == 1) {
              p = p +  geom_text_repel(aes_string(label = ID), size = 2)
            }
            
            else {
              p = p 
            }
            
            if (length(orders) >= 1) {
              p = p + scale_colour_brewer(type="qual", palette="Set1", breaks=orders) + 
                scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) 
            }
            else {
              p = p + scale_colour_brewer(type="qual", palette="Set1") + 
                scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16)) 
            }
            
            if (length(facet) == 1) {
              ncol <- length(unique(mapping.sel.na.rem[,facet]))
              p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
            }
            else {
              p = p
            }
            print(p)
          }
      }
      dev.off()
  }
}
removeRows <- function(rowNum, data) {
    newData <- data[-rowNum, , drop = FALSE]
    rownames(newData) <- NULL
    newData
}




getTaxonomy <- function(otus, tax_tab, level, taxRanks,na_str = c( "NA")) {
    ranks <- taxRanks
    sel <- ranks[1:match(level, ranks)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
    retval[inds!=match(level, ranks)] <- paste(na_str[1], retval[inds!=match(level, ranks)], sep="_")
    return(retval)
}


getTaxonomyAll <- function(otus, tax_tab, level, na_str = c( "NA")) {
    rank <- rank
    sel <- rank[1:match(level, rank)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, rank[inds])]
    retval[inds!=match(level, rank)] <- paste(na_str[1], retval[inds!=match(level, rank)], sep="_")
    return(retval)
}



# X is indicator matrix of predictions, Y is indicator matrix of truth
# columns are classes, rows are samples
mcc <- function(preds=NULL, actuals=NULL, x=NULL, y=NULL) {
    # if preds and actuals are provided, x and y will be ignored
    if (!is.null(preds)) {
        nclasses <- length(union(preds, actuals))
        x <- matrix(0, nrow=length(preds), ncol=nclasses)
        y <- matrix(0, nrow=length(actuals), ncol=nclasses)
        x[cbind(1:nrow(x), preds+1)] <- 1
        y[cbind(1:nrow(y), actuals+1)] <- 1
    }
    if (!all(dim(x) == dim(y))) {
        stop("X and Y must have the same dimensions")
    }
    
    cov_biased <- function(x, y) {
        sum(sapply(1:ncol(x), function(k) {
            cov(x[,k], y[,k]) # unbiased estimate with (n-1) denominator as opposed to (n), but cancels out anyways so identical result
        }))
    }
    numerator <- cov_biased(x,y)
    denominator <- sqrt(cov_biased(x,x) * cov_biased(y,y))
    numerator / denominator
}


multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    i = 1
    while (i < numPlots) {
        numToPlot <- min(numPlots-i+1, cols*rows)
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
        if (numToPlot==1) {
            print(plots[[i]])
        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            # Make each plot, in the correct location
            for (j in i:(i+numToPlot-1)) {
                # Get the i,j matrix positions of the regions that contain this subplot
                matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
                print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col))
            }
        }
        i <- i+numToPlot
    }
}

normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
    tab <- table(fvec)
    retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
    #    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
    retval <- factor(retval, levels=newlevels)
    if (originalLevelsAsNames) {
        names(retval) <- fvec
    }
    return(retval)
}






# ============================================================
# Tutorial on plotting significant taxa and environmental variables on an NMDS plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(vegan)
library(ggplot2)
library(grid)


bv.step <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
max.rho=0.95,
min.delta.rho=0.001,
random.selection=TRUE,
prop.selected.var=0.2,
num.restarts=10,
var.always.include=NULL,
var.exclude=NULL,
output.best=10
){
    
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
    require(vegan)
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    
    fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
    
    #an initial removal phase
    var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
    full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
    var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
    RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
    for(i in 1:dim(var.comb)[2]){
        var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
        RES$rho[i] <- temp$estimate
    }
    delta.rho <- RES$rho - full.cor
    exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
    
    if(random.selection){
        num.restarts=num.restarts
        prop.selected.var=prop.selected.var
        prob<-rep(1,ncol(var.mat))
        if(prop.selected.var< 1){
            prob[exclude]<-0
        }
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    } else {
        num.restarts=1
        prop.selected.var=1
        prob<-rep(1,ncol(var.mat))
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    }
    
    RES_TOT <- c()
    for(i in 1:num.restarts){
        step=1
        RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
        attr(RES$step.dir, "levels") <- c("F","B")
        best.comb <- which.max(RES$rho)
        best.rho <- RES$rho[best.comb]
        delta.rho <- Inf
        selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
        while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
            #forward step
            step.dir="F"
            step=step+1
            var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
            if(RES$n.var[best.comb] == 0){
                var.comb.incl<-1:length(var.comb)
            } else {
                var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                temp <- NA*1:length(var.comb)
                for(j in 1:length(temp)){
                    temp[j] <- all(var.keep %in% var.comb[[j]])
                }
                var.comb.incl <- which(temp==1)
            }
            
            RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
            for(f in 1:length(var.comb.incl)){
                var.incl <- var.comb[[var.comb.incl[f]]]
                var.incl <- var.incl[order(var.incl)]
                var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                RES.f$var.incl[f] <- paste(var.incl, collapse=",")
                RES.f$rho[f] <- temp$estimate
            }
            
            last.F <- max(which(RES$step.dir=="F"))
            RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
            best.comb <- which.max(RES$rho)
            delta.rho <- RES$rho[best.comb] - best.rho
            best.rho <- RES$rho[best.comb]
            
            if(best.comb == step){
                while(best.comb == step & RES$n.var[best.comb] > 1){
                    #backward step
                    step.dir="B"
                    step <- step+1
                    var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                    var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
                    RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
                    for(b in 1:length(var.comb)){
                        var.incl <- var.comb[[b]]
                        var.incl <- var.incl[order(var.incl)]
                        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                        RES.b$var.incl[b] <- paste(var.incl, collapse=",")
                        RES.b$rho[b] <- temp$estimate
                    }
                    RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
                    best.comb <- which.max(RES$rho)
                    best.rho<- RES$rho[best.comb]
                }
            } else {
                break()
            }
            
        }
        
        RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
        print(paste(round((i/num.restarts)*100,3), "% finished"))
    }
    
    RES_TOT <- unique(RES_TOT[,3:5])
    
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
    } else {
        order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
    }
    rownames(order.by.best)<-NULL
    
    order.by.i.comb <- c()
    for(i in 1:length(selected.var)){
        f1 <- which(RES_TOT$n.var==i)
        f2 <- which.max(RES_TOT$rho[f1])
        order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
    }
    rownames(order.by.i.comb)<-NULL
    
    if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
    out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
    )
    out
    
}

bio.env <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
output.best=10,
var.max=ncol(var.mat)
){
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(var.max > dim(var.mat)[2]){stop("var.max cannot be larger than the number of variables (columns) in var.mat")}
    
    require(vegan)
    
    combn.sum <- sum(factorial(ncol(var.mat))/(factorial(1:var.max)*factorial(ncol(var.mat)-1:var.max)))
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    fix.dist <- vegdist(fix.mat, method=fix.dist.method)
    RES_TOT <- c()
    best.i.comb <- c()
    iter <- 0
    for(i in 1:var.max){
        var.comb <- combn(1:ncol(var.mat), i, simplify=FALSE)
        RES <- data.frame(var.incl=rep(NA, length(var.comb)), n.var=i, rho=0)
        for(f in 1:length(var.comb)){
            iter <- iter+1
            var.dist <- vegdist(as.matrix(var.mat[,var.comb[[f]]]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES$var.incl[f] <- paste(var.comb[[f]], collapse=",")
            RES$rho[f] <- temp$estimate
            if(iter %% 100 == 0){print(paste(round(iter/combn.sum*100, 3), "% finished"))}
        }
        
        order.rho <- order(RES$rho, decreasing=TRUE)
        best.i.comb <- c(best.i.comb, RES$var.incl[order.rho[1]])
        if(length(order.rho) > output.best){
            RES_TOT <- rbind(RES_TOT, RES[order.rho[1:output.best],])
        } else {
            RES_TOT <- rbind(RES_TOT, RES)
        }
    }
    rownames(RES_TOT)<-NULL
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)[1:output.best]
    } else {
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)
    }
    OBB <- RES_TOT[order.by.best,]
    rownames(OBB) <- NULL
    
    order.by.i.comb <- match(best.i.comb, RES_TOT$var.incl)
    OBC <- RES_TOT[order.by.i.comb,]
    rownames(OBC) <- NULL
    
    out <- list(
    order.by.best=OBB,
    order.by.i.comb=OBC,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(OBB$var.incl[1], ",")))], collapse=",") ,
    best.model.rho=OBB$rho[1]
    )
    out
}



heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
distfun = dist, hclustfun = hclust, dendrogram = c("both",
"row", "column", "none"), reorderfun = function(d, w) reorder(d,
w), symm = FALSE, scale = c("none", "row", "column"),
na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
col = "heat.colors", colsep, rowsep, sepcolor = "white",
sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
na.color = par("bg"), trace = c("column", "row", "both",
"none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks),
linecol = tracecol, margins = c(5,5,5,5), ColSideColors, RowSideColors, side.height.fraction=0.3,
cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,
NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
key = TRUE, keysize = 1.5, density.info = c("histogram",
"density", "none"), denscol = tracecol, symkey = any(x <
0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL,
key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL,
key.par = list(), main = NULL, xlab = NULL, ylab = NULL,
lmat = NULL, lhei = NULL, lwid = NULL, ColSideColorsSize = 1, RowSideColorsSize = 1, extrafun = NULL, ...)
{
    library(gtools)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
    "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
    "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
    else if (all(Colv == "Rowv"))
    Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 4)
    stop("`margins' must be a numeric vector of length 4")
    if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
            dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
            dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
        nr))
        stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        browser()
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd >
        nc))
        stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
    (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
    (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
    1) {
        if (missing(col) || is.function(col))
        breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
    col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
            stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
            stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[4]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            #            par(mar = c(0.5, 0, 0, margins[4]))
            par(mar = c(0.5, margins[2], 0, margins[4]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
    par(mar = margins)
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
    retval$rowDendrogram <- ddr
    if (exists("ddc"))
    retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
    }
    if (is.null(srtCol))
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
    offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
    padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol))
            adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
            tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
            strheight("M"), labels = labCol, adj = adjCol,
            cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        par(mar = c(margins[1L], 0, 0, margins[4L]))
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
        tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2,
            line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
            y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
            srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[4] - 1.25)
    if (!missing(add.expr))
    eval(substitute(add.expr))
    if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
    xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
    1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
    1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
    1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
    col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol,
                lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
    col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE,
        yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[4]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i",
        leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab))
        mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab))
        mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title))
        mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0)
        do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row")
            key.xlab <- "Row Z-Score"
            else if (scale == "column")
            key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) *
                0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Density"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Count"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (is.null(key.title))
        title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun))
    extrafun()
    invisible(retval)
}


normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}





data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}



pairwise.adonis <- function(x,factors, map,mvar, sim.function = 'vegdist', sim.method = 'bray', adjust,p.adjust.m ='bonferroni',reduce=NULL)
{
  set.seed(1)
  co <- combn(unique(as.character(map[,mvar])),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    
    #run
    
    map.pair <- subset(map, map[,mvar] %in% c(co[1,elem],co[2,elem]))
    
    if (count(map.pair[,mvar])[1,2] <=2 | count(map.pair[,mvar])[2,2] <=2){
      next
    }
    form <- as.formula(sprintf("x1 ~ %s", mvar))
    print(form)
    ad <- adonis(form, data= map.pair,permutations=999, by="margin")
    
    
    #원본
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$aov.tab[1,1])
    SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
    F.Model <- c(F.Model,ad$aov.tab[1,4]);
    R2 <- c(R2,ad$aov.tab[1,5]);
    p.value <- c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  sig[p.adjusted <= 0.0001] <-'****'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}


### Method summary
summary.pwadonis = function(object, ...) {
    cat("Result of pairwise.adonis:\n")
    cat("\n")
    print(object, ...)
    cat("\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

pairwise.adonis.ad <- function(x,factors, map,mvar, sim.function = 'vegdist', sim.method = 'bray', adjust,p.adjust.m ='bonferroni',reduce=NULL)
{
  set.seed(1)
  co <- combn(unique(as.character(map[,mvar])),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    #test

    map.pair <- subset(map, map[,mvar] %in% c(co[1,elem],co[2,elem]))
  
    if (count(map.pair[,mvar])[1,2] <=2 | count(map.pair[,mvar])[2,2] <=2){
      next
    }
    form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
    print(form)
    ad <- adonis(form, data= map.pair,permutations=999, by="margin")
    
    
    #원본
    #ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],permutations = 999)
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$aov.tab[1,1])
    SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
    F.Model <- c(F.Model,ad$aov.tab[1,4]);
    R2 <- c(R2,ad$aov.tab[1,5]);
    p.value <- c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  sig[p.adjusted <= 0.0001] <-'****'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

# Title: Geometric Mean of Pairwise Ratios (GMPR) for Microbiome Sequencing data normalization
# Version: 0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Date: 2017/02/07
# Description: The function calculates the normalizing factors for microbiome sequencing data or, more generally, zeroinflated sequencing data.
# The size factors can be used as offsets in count-based regression models or as divisors to produce normalized data


require(matrixStats)

GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n',
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}
#------------------------------------------------------#
#---------      [Core general] analysis       ---------#
#------------------------------------------------------#
library(crayon)
#functions <- c("Go_filter()","Go_qq()","Go_adiv()", "Go_quantile()","Go_box_plot()","Go_barchart()", "Go_barchartAll()","Go_bdiv()","Go_perm()", "Go_deseq2()", "Go_deseq2_heat()","Go_deseq2_fore()","Go_deseq2_volc()")


cat(blue("#--------------------------------------------------------------# \n"))
cat(blue("#------    Quick statistics and visualization tools      ------# \n"))
cat(blue("#--------------------------------------------------------------# \n"))
cat(red("                                  Version: core_general.2.3.2 \n"))
cat("                                              Write by Heekuk \n")
cat(yellow("Automated install for the required package and read libraries.\n"))
cat(blue("#--------------------------------------------------------------# \n"))


#for (func in functions){
#  cat(sprintf("    Function is activated for %s\n", func))
#}

#cat("                                                       \n")
#cat(yellow("What's new: \n"))
#cat(yellow("core_general.2.0.1 (Jan 03,2020): Go_quantile() is added  and bug fixed.\n"))
#cat(yellow("core_general.2.0.0 (Jan 02,2020): Able to use numeric levels as factors \n"))
#cat(yellow("core_general.1.9.2 (Dec 30,2019): NA is ignored all the analysis \n"))
#cat(yellow("core_general.1.9.1 (Dec 27,2019): Go_deseq2_volc() is added for volcano plot.\n"))
#cat(yellow("core_general.1.9.0 (Dec 13,2019): add box or violin in plot of Go_adiv().\n"))
#cat(yellow("core_general.1.8.1 (Dec 05,2019): Number of Multi plot added in Go_adiv().\n"))
#cat(yellow("core_general.1.8.0 (Dec 02,2019): Go_filter() is added and bug fixed.\n"))
#cat(yellow("core_general.1.7.2 (Nov 30,2019): bug fixed.\n"))
#cat(yellow("core_general.1.7.1 (Nov 29,2019): Go_barchart() is added variation taxa ranks (taxaRanks =).\n"))
#cat(yellow("core_general.1.7.0 (Nov 26,2019): Go_barchart_dada() is merged to Go_barchart() (data_type =).\n"))
#cat(yellow("core_general.1.6.0 (Nov 25,2019): Go_barchart_dada() is added for x label with certain variation (x_label =).\n"))
#cat(yellow("core_general.1.5.0 (Nov 22,2019): Go_deseq2() is added for analysis with certain variation (addvari = ).\n"))
#cat(yellow("core_general.1.4.1 (Nov 19,2019): bug fixed.\n"))
#cat(yellow("core_general.1.4.0 (Nov 14,2019): out multi in Go_perm().\n"))
#cat(yellow("core_general.1.3.0 (Nov 11,2019): shapes and ID variations (for labelling each dot) are added in Go_bdiv().\n"))
