#!/usr/bin/Rscript 

library(dplyr)
library(tibble)
library(data.table)
library(caret)
library(pROC)

library(survcomp)
library(gtools)
library(survival)
library(survminer)
library(gtools)
library(caretEnsemble)
library(cowplot)

# rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

grp <- c("S1", "S2")
grp.col <- c("#6C326C", "#77A2D1")

library(optparse)
# get options
option_list = list(
  make_option(c("-p", "--phen"), type="character", default=".",
        help="phenotype information [default= %default]", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=".",
        help="profile file [default= %default]", metavar="character"),
  make_option(c("-a", "--algorithm"), type="character", default=".",
        help="Machine learning algorithm [default= %default]", metavar="character"),
  make_option(c("-k", "--kind"), type="character", default=".",
        help="the kind of run model [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out",
        help="output directory name [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="out",
        help="output name [default= %default]", metavar="character")
);
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

# input parameters
phen <- fread(opt$phen)
profile <- fread(opt$file)
algorithm <- opt$algorithm
kind <- opt$kind
out <- opt$out 
name <- opt$name


###############################################################
###      Funtions of building and assessing model       #######
###############################################################
### ROC/AUC
ROC_fun <- function(model_fit=model){
  
  # model_fit=model
  
  fit <- model_fit$model
  tData <- model_fit$tData
  
  rocobj <- roc(tData$Cluster, predict(fit, newdata = tData, type = "prob")[, grp[1]])
  auc <- round(auc(tData$Cluster, predict(fit, newdata = tData, type = "prob")[, grp[1]]),4)
  plotroc <- tibble(tpr = rocobj$sensitivities,
                     fpr = 1 - rocobj$specificities)
  plroc <- ggplot(data = plotroc, aes(x=fpr, y=tpr))+
              geom_path(color="black", size = 1)+
              geom_abline(intercept = 0, slope = 1, color="grey", size = 1, linetype=2)+
              labs(x = "False Positive Rate (1 - Specificity)",
                   y = "True Positive Rate (Sensivity or Recall)")+
              annotate("text",x = .75, y = .25,label=paste("AUC =", auc),
                       size = 5, family="serif")+
              coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
              theme_bw()+
              theme(panel.background = element_rect(fill = 'transparent'),
                    axis.ticks.length = unit(0.4, "lines"),
                    axis.ticks = element_line(color='black'),
                    axis.line = element_line(size=.5, colour = "black"),
                    axis.title = element_text(colour='black', size=12,face = "bold"),
                    axis.text = element_text(colour='black',size=10,face = "bold"),
                    text = element_text(size=8, color="black", family="serif"))
  res <- list(rocobj=rocobj, auc=auc, pl=plroc)
  return(res)  
}

### Evaluation
Evaluation_fun <- function(model_fit=all_model_svm){
  
  # model_fit=all_model_svm
  
  fit <- model_fit$model
  tData <- model_fit$tData
  tData_all <- model_fit$tData_all
  
  pred <- predict(fit, newdata=tData)
  
  # test Data
  testData_new <- data.frame(SampleID=rownames(tData),
                             Predict_Cluster=pred) %>%
    inner_join(tData_all %>% rownames_to_column("SampleID"),
               by = "SampleID") %>%
    dplyr::select(SampleID, Cluster, Predict_Cluster, OS, OS.Time) %>%
    na.omit()

  # Concordance index
  Concordance_index_fun <- function(dataset=testData_new, 
                                    group_name="Cluster"){
    # dataset=testData_new
    # group_name="Cluster"
    
    require(survcomp)
    colnames(dataset)[which(colnames(dataset) == group_name)] <- "Group"
    
    ci_res <- with(dataset, concordance.index(x=Group, surv.time=OS.Time, surv.event=OS, 
                      method="noether"))  
    
    res <- data.frame(Type=group_name,
                       C_index=signif(ci_res$c.index, 3),
                       pvalue=ci_res$p.value,
                       se=signif(ci_res$se, 3),
                       odd=with(ci_res, paste0("[", signif(lower, 3), ", ", signif(upper, 3), "]")))
    
    return(res)
  }
  Cluster_ci <- Concordance_index_fun(dataset=testData_new, group_name="Cluster")
  Predict_Cluster_ci <- Concordance_index_fun(dataset=testData_new, group_name="Predict_Cluster")
  ci_res <- rbind(Cluster_ci, Predict_Cluster_ci)
  #DT::datatable(ci_res)

  # Log-rank p-value of Cox-PH regression
  Log_rank_fun <- function(dataset=testData_new, 
                           group_name="Cluster"){
    # dataset=testData_new
    # group_name="Cluster"
    
    require(survival)
    require(survminer)
    colnames(dataset)[which(colnames(dataset) == group_name)] <- "Group"
    cox <- coxph(Surv(OS.Time, OS) ~ Group, data = dataset) 
    tmp <- summary(cox)
    tmp.wald <- data.frame(t(tmp$waldtest)) %>%
          setNames(c("Wald_test", "Wald_df", "Wald_pvlaue"))
    tmp.lg <- data.frame(t(tmp$logtest)) %>%
          setNames(c("lg_rank", "lg_rank_df", "lg_rank_pvlaue"))
    tmp.total <- cbind(tmp.wald, tmp.lg) 
      
    pvalue <- paste(paste0("Log-Rank P=", signif(tmp.lg$lg_rank_pvlaue, 3)), 
                    paste0("Cox P=", signif(tmp.wald$Wald_pvlaue, 3)), sep = "\n")
    factors <- levels(dataset$Group)
    surv_fit <- survfit(Surv(OS.Time, OS) ~ Group, data = dataset)
    
    sur_pl <- ggsurvplot(surv_fit,
               data = dataset,  
               surv.median.line = "hv",
               add.all = TRUE,
               palette = "aaas",
               risk.table = TRUE, 
               xlab = "Follow up time(Years)", 
               legend = c(0.8, 0.2), 
               legend.title = "",
               legend.labs = c("all", factors), 
               break.x.by = 2,
               font.legend = c(10, "italic"),
               ggtheme = theme_bw())
    sur_pl$plot <- sur_pl$plot + 
               annotate("text", x=3, y=0.2, label=pvalue)
      
    pl <- cowplot::plot_grid(sur_pl$plot, sur_pl$table, 
                             ncol = 1, align = "hv", 
                             rel_heights = c(2, 1))
      
    pval <- data.frame(Type=group_name,
                             tmp.total)
    res <- list(plot=pl, pvalue=pval)
    return(res)
  }
  
  Cluster_surv <- Log_rank_fun(dataset=testData_new, group_name="Cluster")
  if(length(unique(testData_new$Predict_Cluster)) == 1){
    Predict_Cluster_surv <- list()
    Predict_Cluster_surv$pvalue <- data.frame(Type="Predict_Cluster",
                                              Wald_test=NA,
                                              Wald_df=NA,
                                              Wald_pvlaue=NA,
                                              lg_rank=NA,
                                              lg_rank_df=NA,
                                              lg_rank_pvlaue=NA)
    Predict_Cluster_surv$plot <- NA
  }else{
    Predict_Cluster_surv <- Log_rank_fun(dataset=testData_new, group_name="Predict_Cluster")    
  }
  surv_res <- rbind(Cluster_surv$pvalue, Predict_Cluster_surv$pvalue)
  #DT::datatable(surv_res)

  surv_pl <- cowplot::plot_grid(Cluster_surv$plot, Predict_Cluster_surv$plot,
                     ncol = 2, align = "hv", labels = c("Actual", "Predict"))
    
  # Brier score
  Brier_score_fun <- function(dataset=testData_new, 
                              group_name="Cluster"){
    # dataset=testData_new
    # group_name="Cluster"
    
    require(survcomp)
    colnames(dataset)[which(colnames(dataset) == group_name)] <- "Group"
    
    dat <- dataset %>% dplyr::select(c("Group", "OS", "OS.Time")) %>%
      setNames(c("score", "event", "time"))
    bsc <- sbrier.score2proba(data.tr=dat, data.ts=dat, method="cox")  
    
    res <- list(res_dat=dat, res_bsc=bsc)
    return(res)
  }
  
  Cluster_bsc <- Brier_score_fun(dataset=testData_new, group_name="Cluster")
  Predict_Cluster_bsc <- Brier_score_fun(dataset=testData_new, group_name="Predict_Cluster")
  
  # if(!all(Cluster_bsc$res_dat$time == Predict_Cluster_bsc$res_dat$time)) {
  #   stop("the two vector of BSCs must be computed for the same points in time!") 
  # }
  # ibsc.comp(bsc1=Cluster_bsc$res_bsc$bsc, 
  #           bsc2=Predict_Cluster_bsc$res_bsc$bsc, time=Cluster_bsc$res_dat$time)  
  
  bsc_res <- data.frame(Type=c("Actual", "Predict"),
                        Brier_score=c(Cluster_bsc$res_bsc$bsc.integrated,
                                      Predict_Cluster_bsc$res_bsc$bsc.integrated))
  
  res <- list(ci_res=ci_res, surv_res=surv_res, bsc_res=bsc_res)
  return(res)
}


# building model
build_model_fun <- function(dataset=profile,
                            algorithm="svmRadial"){
  
  # dataset=profile
  # algorithm="svmRadial"
  
  
  set.seed(123)
  mdat <- dataset %>% mutate(Cluster=factor(Cluster, levels = grp)) %>% 
    column_to_rownames("V1")
  
  # split data into 5 folds based on the cluster
  require(sampling)
  split_folds <- list()
  sample_size <- rev(as.numeric(table(mdat$Cluster))/5)
  for(i in 1:5){
    if(i == 1){
      mdat_test <- mdat
    }else{
      mdat_test <- mdat_remain
    }
    picked <- sampling::strata(mdat_test, 
                               stratanames = ("Cluster"),
                               size = sample_size,
                               method = "srswor")
    
    mdat_remain <- mdat_test[-picked$ID_unit, ]
    # print(table(mdat_test[picked$ID_unit, ]$Cluster))
    split_folds[[i]] <- mdat_test[picked$ID_unit, ]
  }
  
  # randomly choose 60% data as train test and 40% as test data
  require(gtools)
  split_folds_combination <- combinations(length(split_folds), 2,
                                          c(1:length(split_folds)))
  ROC_list <- list()
  Evaluation_list <- list()
  rocobj_list <- list()
  auc_list <- list()  
  for(i in 1:nrow(split_folds_combination)){
    print(i)
    testData_all <- do.call(rbind, split_folds[split_folds_combination[i, ]])
    trainData_all <- do.call(rbind, split_folds[-split_folds_combination[i, ]])
    trainData <- trainData_all %>% dplyr::select(colnames(phen)[2], colnames(mdat)[-c(1:11)]) 
    testData <- testData_all %>% dplyr::select(colnames(phen)[2], colnames(mdat)[-c(1:11)])  

    # set parameters for model
    myControl <- trainControl(method = "repeatedcv", 
                             number = 10,
                             repeats = 3,
                             search = "grid",
                             classProbs = TRUE,
                             savePredictions = TRUE,
                             verboseIter = FALSE)  
    
    metrics <- "Accuracy"
    
    if(algorithm == "all"){
      algorithmList  <- c("svmRadial", "rf", "gbm", "LogitBoost", "nb")
      set.seed(123)
      require(caretEnsemble)
      fit <-  caretList(Cluster ~ ., data = trainData, 
                        methodList = algorithmList , 
                        trControl = myControl, 
                        metric = metrics,
                        verbose = FALSE) 
    }else{
      if(algorithm == "svmRadial"){        # Support Vector Machines with Radial Basis Function Kernel
        myGrid <- expand.grid(sigma = seq(0.01, 0.1, 0.02),
                              C = seq(0, 2, length = 10))
        # new_algorithm <- getModelInfo()$lssvmLinear
        # new_algorithm$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
        #   kernlab::lssvm(x = as.matrix(x), y = y,
        #                  tau = param$tau)}
        new_algorithm <- "svmRadial"
      }else if(algorithm == "rf"){         # random forest 
        myGrid <- expand.grid(.mtry=c(sqrt(ncol(trainData)-1)))
        new_algorithm <- "rf"   
      }else if(algorithm == "LogitBoost"){ # Boosted Logistic Regression
        myGrid <- expand.grid(nIter  = 500)
        new_algorithm <- "LogitBoost"     
      }else if(algorithm == "gbm"){        # Stochastic Gradient Boosting
        myGrid <- expand.grid(interaction.depth = c(1, 5, 9),  
                              n.trees = seq(1, 30, 4)*50,
                              shrinkage = 0.1,
                              n.minobsinnode = 20) 
        new_algorithm <- "gbm" 
      }else if(algorithm == "nb"){      # naive_bayes
        myGrid <- expand.grid(fL=c(0,0.5,1.0), 
                              usekernel = TRUE, 
                              adjust=c(0,0.5,1.0)) 
        new_algorithm <- "nb"     
      }
      set.seed(123)
      fit <- train(Cluster ~ ., 
                   data = trainData, 
                   method = new_algorithm, 
                   trControl = myControl, 
                   tuneGrid = myGrid,
                   metric = metrics,
                   verbose = FALSE) 
    pred <- predict(fit, newdata=testData)
    print(confusionMatrix(pred, testData$Cluster))
    }
    model_res <- list(tData_all=testData_all, tData=testData, model=fit)
    
    # ROC 
    ROC_list[[i]] <- ROC_fun(model_fit = model_res)
    rocobj_list[[i]] <- ROC_list[[i]]$rocobj
    auc_list[[i]] <- ROC_list[[i]]$auc
    # Evalution index
    Evaluation_list[[i]] <- Evaluation_fun(model_fit = model_res)
    
  }
  
  res <- list(ROC=ROC_list, Eva=Evaluation_list,
              rocobj=rocobj_list, auc=auc_list)
  return(res)
}

# trainning model using 4-omics layer data, but test model by using 4-omics layer data and each single omic layer of data respectively.
build_model_fun_v2 <- function(dataset=all_norm,
                               algorithm="svmRadial"){
  
  # dataset=all_norm
  # algorithm="svmRadial"
  
  set.seed(123)
  mdat <- dataset %>% mutate(Cluster=factor(Cluster, levels = grp)) %>% 
    column_to_rownames("V1")
  
  # split data into 5 folds based on the cluster
  require(sampling)
  split_folds <- list()
  sample_size <- rev(as.numeric(table(mdat$Cluster))/5)
  for(i in 1:5){
    if(i == 1){
      mdat_test <- mdat
    }else{
      mdat_test <- mdat_remain
    }
    picked <- sampling::strata(mdat_test, 
                               stratanames = ("Cluster"),
                               size = sample_size,
                               method = "srswor")
    
    mdat_remain <- mdat_test[-picked$ID_unit, ]
    # print(table(mdat_test[picked$ID_unit, ]$Cluster))
    split_folds[[i]] <- mdat_test[picked$ID_unit, ]
  }
  
  # randomly choose 60% data as train test and 40% as test data
  require(gtools)
  split_folds_combination <- combinations(length(split_folds), 2,
                                          c(1:length(split_folds)))
  ROC_list_all <- list()
  ROC_list_mRNA <- list()
  ROC_list_CNV <- list()
  ROC_list_Methylation <- list()
  ROC_list_RPPA <- list()
  
  Evaluation_list_all <- list()
  Evaluation_list_mRNA <- list()
  Evaluation_list_CNV <- list()
  Evaluation_list_Methylation <- list()
  Evaluation_list_RPPA <- list()
  
  rocobj_list_all <- list()
  rocobj_list_mRNA <- list()
  rocobj_list_CNV <- list()
  rocobj_list_Methylation <- list()
  rocobj_list_RPPA <- list()

  auc_list_all <- list()
  auc_list_mRNA <- list()
  auc_list_CNV <- list()
  auc_list_Methylation <- list()
  auc_list_RPPA <- list()    
  
  for(i in 1:nrow(split_folds_combination)){
    print(i)
    testData_all <- do.call(rbind, split_folds[split_folds_combination[i, ]])
    trainData_all <- do.call(rbind, split_folds[-split_folds_combination[i, ]])
    trainData <- trainData_all %>% dplyr::select(colnames(phen)[2], colnames(mdat)[-c(1:11)]) 
    
    testData <- testData_all %>% dplyr::select(colnames(phen)[2], colnames(mdat)[-c(1:11)]) 
    
    testData_convert <- function(tag="geneExp"){
      # tag="geneExp"
      dat_tag <- testData %>% dplyr::select(c("Cluster", ends_with(tag)))
      dat_no_tag <- testData %>% dplyr::select(-colnames(dat_tag))
      
      dat_no_tag_0 <- data.frame(matrix(0, nrow = nrow(dat_no_tag), ncol = ncol(dat_no_tag)))
      rownames(dat_no_tag_0) <- rownames(dat_no_tag)
      colnames(dat_no_tag_0) <- colnames(dat_no_tag)
      
      res <- inner_join(dat_tag %>% rownames_to_column("tmp"),
                        dat_no_tag_0 %>% rownames_to_column("tmp"),
                        by = "tmp") %>%
        column_to_rownames("tmp") %>%
        dplyr::select(colnames(testData))
      return(res)
    }
    
    testData_mRNA <- testData_convert(tag="geneExp")
    testData_CNV <- testData_convert(tag="copyNumber")
    testData_Methylation <- testData_convert(tag="methylation")
    testData_RPPA <- testData_convert(tag="protein_RPPA")

    # set parameters for model
    myControl <- trainControl(method = "repeatedcv", 
                             number = 5,
                             repeats = 3,
                             search = "grid",
                             classProbs = TRUE,
                             savePredictions = TRUE,
                             verboseIter = FALSE)  
    
    metrics <- "Accuracy"
    
    if(algorithm == "all"){
      algorithmList  <- c("svmRadial", "rf", "gbm", "LogitBoost", "nb")
      set.seed(123)
      require(caretEnsemble)
      fit <-  caretList(Cluster ~ ., data = trainData, 
                        methodList = algorithmList , 
                        trControl = myControl, 
                        metric = metrics,
                        verbose = FALSE) 
    }else{
      if(algorithm == "svmRadial"){        # Support Vector Machines with Radial Basis Function Kernel
        myGrid <- expand.grid(sigma = seq(0.01, 0.1, 0.02),
                              C = seq(0, 2, length = 10))
        # new_algorithm <- getModelInfo()$lssvmLinear
        # new_algorithm$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
        #   kernlab::lssvm(x = as.matrix(x), y = y,
        #                  tau = param$tau)}
        new_algorithm <- "svmRadial"
      }else if(algorithm == "rf"){         # random forest 
        myGrid <- expand.grid(.mtry=c(sqrt(ncol(trainData)-1)))
        new_algorithm <- "rf"   
      }else if(algorithm == "LogitBoost"){ # Boosted Logistic Regression
        myGrid <- expand.grid(nIter  = 500)
        new_algorithm <- "LogitBoost"     
      }else if(algorithm == "gbm"){        # Stochastic Gradient Boosting
        myGrid <- expand.grid(interaction.depth = c(1, 5, 9),  
                              n.trees = seq(1, 30, 4)*50,
                              shrinkage = 0.1,
                              n.minobsinnode = 20) 
        new_algorithm <- "gbm" 
      }else if(algorithm == "nb"){      # naive_bayes
        myGrid <- expand.grid(fL=c(0,0.5,1.0), 
                              usekernel = TRUE, 
                              adjust=c(0,0.5,1.0)) 
        new_algorithm <- "nb"     
      }
      set.seed(123)
      fit <- train(Cluster ~ ., 
                   data = trainData, 
                   method = new_algorithm, 
                   trControl = myControl, 
                   tuneGrid = myGrid,
                   metric = metrics,
                   verbose = FALSE) 
    pred <- predict(fit, newdata=testData)
    print(confusionMatrix(pred, testData$Cluster))
    }
    model_res <- list(tData_all=testData_all, tData=testData, model=fit)
    model_res_mRNA <- list(tData_all=testData_all, tData=testData_mRNA, model=fit)
    model_res_CNV <- list(tData_all=testData_all, tData=testData_CNV, model=fit)
    model_res_Methylation <- list(tData_all=testData_all, tData=testData_Methylation, model=fit)
    model_res_RPPA <- list(tData_all=testData_all, tData=testData_RPPA, model=fit)
    
    # ROC 
    ROC_list_all[[i]] <- ROC_fun(model_fit = model_res)
    rocobj_list_all[[i]] <- ROC_list_all[[i]]$rocobj
    auc_list_all[[i]] <- ROC_list_all[[i]]$auc
    
    ROC_list_mRNA[[i]] <- ROC_fun(model_fit = model_res_mRNA)
    rocobj_list_mRNA[[i]] <- ROC_list_mRNA[[i]]$rocobj
    auc_list_mRNA[[i]] <- ROC_list_mRNA[[i]]$auc

    ROC_list_CNV[[i]] <- ROC_fun(model_fit = model_res_CNV)
    rocobj_list_CNV[[i]] <- ROC_list_CNV[[i]]$rocobj
    auc_list_CNV[[i]] <- ROC_list_CNV[[i]]$auc
    
    ROC_list_Methylation[[i]] <- ROC_fun(model_fit = model_res_Methylation)
    rocobj_list_Methylation[[i]] <- ROC_list_Methylation[[i]]$rocobj
    auc_list_Methylation[[i]] <- ROC_list_Methylation[[i]]$auc
    
    ROC_list_RPPA[[i]] <- ROC_fun(model_fit = model_res_RPPA)
    rocobj_list_RPPA[[i]] <- ROC_list_RPPA[[i]]$rocobj
    auc_list_RPPA[[i]] <- ROC_list_RPPA[[i]]$auc    
    
    # Evalution index
    Evaluation_list_all[[i]] <- Evaluation_fun(model_fit = model_res)
    Evaluation_list_mRNA[[i]] <- Evaluation_fun(model_fit = model_res_mRNA)
    Evaluation_list_CNV[[i]] <- Evaluation_fun(model_fit = model_res_CNV)
    Evaluation_list_Methylation[[i]] <- Evaluation_fun(model_fit = model_res_Methylation)
    Evaluation_list_RPPA[[i]] <- Evaluation_fun(model_fit = model_res_RPPA)
  }
  
  res <- list(ROC_all=ROC_list_all,
              ROC_mRNA=ROC_list_mRNA,
              ROC_CNV=ROC_list_CNV,
              ROC_Methylation=ROC_list_Methylation,
              ROC_RPPA=ROC_list_RPPA,
              
              Eva_all=Evaluation_list_all,
              Eva_mRNA=Evaluation_list_mRNA,
              Eva_CNV=Evaluation_list_CNV,
              Eva_Methylation=Evaluation_list_Methylation,
              Eva_RPPA=Evaluation_list_RPPA, 
              
              rocobj_all=rocobj_list_all,
              rocobj_mRNA=rocobj_list_mRNA,
              rocobj_CNV=rocobj_list_CNV,
              rocobj_Methylation=rocobj_list_Methylation,
              rocobj_RPPA=rocobj_list_RPPA,
              
              auc_all=auc_list_all,
              auc_mRNA=auc_list_mRNA,
              auc_CNV=auc_list_CNV,
              auc_Methylation=auc_list_Methylation,
              auc_RPPA=auc_list_RPPA)
  return(res)
}

### Summary ROC
summary_ROC <- function(dat_list=model){
  
  # dat_list=model
  
  rocboj_list <- dat_list$rocobj
  auc_list <- dat_list$auc
  
  colors_value <- c("#A6CEE3", "#1F78B4", "#08306B", "#B2DF8A", "#006D2C", 
              "#8E0152", "#DE77AE", "#CAB2D6", "#6A3D9A", "#FB9A99")
  
  pl <- ggroc(rocboj_list, linetype = 1, size = 1, alpha = 1, legacy.axes = T)+
          geom_abline(intercept = 0, slope = 1, color="grey", size = 1, linetype=1)+
          labs(x = "False Positive Rate (1 - Specificity)",
               y = "True Positive Rate (Sensivity or Recall)")+
          geom_label(aes(x = .75, y = .50, label=paste0("AUC = ", auc_list[[1]], "(i=1)")),
                   size = 5, family="serif", color=colors_value[1])+
          geom_label(aes(x = .75, y = .45, label=paste0("AUC = ", auc_list[[2]], "(i=2)")),
                   size = 5, family="serif", color=colors_value[2])+
          geom_label(aes(x = .75, y = .40, label=paste0("AUC = ", auc_list[[3]], "(i=3)")),
                   size = 5, family="serif", color=colors_value[3])+
          geom_label(aes(x = .75, y = .35, label=paste0("AUC = ", auc_list[[4]], "(i=4)")),
                   size = 5, family="serif", color=colors_value[4])+
          geom_label(aes(x = .75, y = .30, label=paste0("AUC = ", auc_list[[5]], "(i=5)")),
                   size = 5, family="serif", color=colors_value[5])+
          geom_label(aes(x = .75, y = .25, label=paste0("AUC = ", auc_list[[6]], "(i=6)")),
                   size = 5, family="serif", color=colors_value[6])+
          geom_label(aes(x = .75, y = .20, label=paste0("AUC = ", auc_list[[7]], "(i=7)")),
                   size = 5, family="serif", color=colors_value[7])+
          geom_label(aes(x = .75, y = .15, label=paste0("AUC = ", auc_list[[8]], "(i=8)")),
                   size = 5, family="serif", color=colors_value[8])+    
          geom_label(aes(x = .75, y = .10, label=paste0("AUC = ", auc_list[[9]], "(i=9)")),
                   size = 5, family="serif", color=colors_value[9])+
          geom_label(aes(x = .75, y = .05, label=paste0("AUC = ", auc_list[[10]], "(i=10)")),
                   size = 5, family="serif", color=colors_value[10])+            
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
          scale_colour_manual(values = colors_value)+
          guides(color=F)+
          theme_bw()+
          theme(panel.background = element_rect(fill = 'transparent'),
                axis.ticks.length = unit(0.4, "lines"), 
                axis.ticks = element_line(color='black'),
                axis.line = element_line(size=.5, colour = "black"),
                axis.title = element_text(colour='black', size=12,face = "bold"),
                axis.text = element_text(colour='black',size=10),
                text = element_text(size=8, color="black", family="serif"))
  return(pl)
}

### Summary Evalution
summary_Evalution <- function(dat_list=model){
  
  #dat_list=model
  
  Evalution_list <- dat_list$Eva
  ci_res <- data.frame()
  surv_res <- data.frame()
  bsc_res <- data.frame()
  for(i in 1:length(Evalution_list)){
    ci_res <- rbind(ci_res, 
                    Evalution_list[[i]]$ci_res %>% mutate(Times=i))
    surv_res <- rbind(surv_res, 
                    Evalution_list[[i]]$surv_res %>% mutate(Times=i))
    bsc_res <- rbind(bsc_res, 
                    Evalution_list[[i]]$bsc_res %>% mutate(Times=i))    
  }
  
  # mean 
  ci_final <- ci_res %>% group_by(Type) %>%
    summarise(C_index_mean=mean(C_index),
              C_index_sd=sd(C_index),
              pvalue_mean=mean(pvalue))
  surv_final <- surv_res %>% group_by(Type) %>%
    summarise(Wald_pvlaue_mean=mean(Wald_pvlaue),
              lg_rank_pvlaue_mean=sd(lg_rank_pvlaue)) 
  bsc_final <- bsc_res %>% group_by(Type) %>%
    summarise(Brier_score_mean=mean(Brier_score),
              Brier_score_sd=sd(Brier_score))  
  
  res <- list(ci=ci_final, surv=surv_final, bsc=bsc_final)
  
  return(res)
}

###############################################################
#########      Run functions and save results        ##########
###############################################################
if(kind == "single"){
  result <- build_model_fun(dataset=profile, algorithm = algorithm)
}else if(kind == "all"){
  result <- build_model_fun_v2(dataset=profile, algorithm = algorithm)
}

if(!dir.exists(out)){
  dir.create(out)
}

filename <- paste0(out, "/", name, "_", kind,".RDS")
saveRDS(result, file=filename)
