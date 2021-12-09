# 4-omics data 
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a svmRadial -k single -o ./Evaluation -n all_model_svm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a rf -k single -o ./Evaluation -n all_model_rf
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a LogitBoost -k single -o ./Evaluation -n all_model_lgb
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a gbm -k single -o ./Evaluation -n all_model_gbm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a nb -k single -o ./Evaluation -n all_model_nb

# copyNumber 
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/copyNumber_top_feature_norm.csv -a svmRadial -k single -o ./Evaluation -n copyNumber_model_svm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/copyNumber_top_feature_norm.csv -a rf -k single -o ./Evaluation -n copyNumber_model_rf
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/copyNumber_top_feature_norm.csv -a LogitBoost -k single -o ./Evaluation -n copyNumber_model_lgb
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/copyNumber_top_feature_norm.csv -a gbm -k single -o ./Evaluation -n copyNumber_model_gbm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/copyNumber_top_feature_norm.csv -a nb -k single -o ./Evaluation -n copyNumber_model_nb

# geneExp 
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/geneExp_top_feature_norm.csv -a svmRadial -k single -o ./Evaluation -n geneExp_model_svm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/geneExp_top_feature_norm.csv-a rf -k single -o ./Evaluation -n geneExp_model_rf
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/geneExp_top_feature_norm.csv -a LogitBoost -k single -o ./Evaluation -n geneExp_model_lgb
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/geneExp_top_feature_norm.csv -a gbm -k single -o ./Evaluation -n geneExp_model_gbm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/geneExp_top_feature_norm.csv -a nb -k single -o ./Evaluation -n geneExp_model_nb

# DNA methylation 
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/methylation_top_feature_norm.csv -a svmRadial -k single -o ./Evaluation -n methylation_model_svm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/methylation_top_feature_norm.csv -a rf -k single -o ./Evaluation -n methylation_model_rf
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/methylation_top_feature_norm.csv -a LogitBoost -k single -o ./Evaluation -n methylation_model_lgb
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/methylation_top_feature_norm.csv -a gbm -k single -o ./Evaluation -n methylation_model_gbm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/methylation_top_feature_norm.csv -a nb -k single -o ./Evaluation -n methylation_model_nb

# protein_RPPA
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/protein_RPPA_top_feature_norm.csv -a svmRadial -k single -o ./Evaluation -n protein_RPPA_model_svm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/protein_RPPA_top_feature_norm.csv -a rf -k single -o ./Evaluation -n protein_RPPA_model_rf
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/protein_RPPA_top_feature_norm.csv -a LogitBoost -k single -o ./Evaluation -n protein_RPPA_model_lgb
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/protein_RPPA_top_feature_norm.csv -a gbm -k single -o ./Evaluation -n protein_RPPA_model_gbm
#Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/protein_RPPA_top_feature_norm.csv -a nb -k single -o ./Evaluation -n protein_RPPA_model_nb


# 4-omics data (test on single omic layer data)
Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a svmRadial -k all -o ./Evaluation -n all_model_svm
Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a rf -k all -o ./Evaluation -n all_model_rf
Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a LogitBoost -k all -o ./Evaluation -n all_model_lgb
Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a gbm -k all -o ./Evaluation -n all_model_gbm
Rscript 02Classification_building_v2.R -p phenotype_cluster.csv -f /disk/user/zouhua/project/KIRC/Assay/06.Classification/classification/all_top_feature_norm.csv -a nb -k all -o ./Evaluation -n all_model_nb
