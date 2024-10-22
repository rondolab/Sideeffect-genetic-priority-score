library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(caret)
library(parallel)
library(logistf)
							      

                    
#genetic features used which comprise the SE-GPS
# geneticpredictors=c('clinicalvariant', 'hgmd', 'omim','geneburden_ot', 'geneburden_ravar','singlevar_gb', 'singleva_ravar','eqtl_phenotype','locus2gene')
geneticpredictors=c('clinicalvariant', 'geneburden', 'singlevar','eqtl_phenotype','locus2gene')

#genetic features used which comprise the GPS-direction of effect
geneticpredictors_doe=c('clinicalvariant_doe','geneburden_doe','singlevar_doe','eqtl_phenotype_doe','locus2gene_doe')

#phecode categories used as covariates in regression models
covariates=c('categorycongenital_anomalies','categorydermatologic','categorydigestive','categoryendocrine_metabolic','categorygenitourinary','categoryhematopoietic','categorymental_disorders','categorymusculoskeletal','categoryneoplasm', 'categoryneurological','categoryrespiratory','categorysense_organs','categorysymptoms')
#Make clinicalvariant predictor
clincicalpred=c('clinvar', 'hgmd', 'omim')
dataset$clinicalvariant=0
dataset$clinicalvariant=0
dataset$clinicalvariant[rowSums(dataset[, clincicalpred]) == 1]<-1
dataset$clinicalvariant[rowSums(dataset[, clincicalpred]) == 2]<-2
dataset$clinicalvariant[rowSums(dataset[, clincicalpred]) == 3]<-3
#Make geneburden predictor
geneburden=c('geneburden_ot', 'geneburden_ravar')
dataset$geneburden=0
dataset$geneburden[rowSums(dataset[, geneburden]) != 0]<-1
#Make singlevariant predictor
singlevar=c('singlevar_gb', 'singlevar_ravar')
dataset$singlevar=0
dataset$singlevar[rowSums(dataset[, singlevar]) != 0]<-1
#SE no main indication outcome
dataset$senomi=ifelse(dataset$se==1 & dataset$mi.phase4!=1,1,0)


## Make Clinical Variant pre
 dataset$clinicalvariant_DOE=0
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == 1]<-1
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == 2]<-2
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == 3]<-3
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == 0]<-0

 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == -1]<--1
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == -2]<--2
 dataset$clinicalvariant_DOE[rowSums(dataset[, clincicalpred_DOE]) == -3]<--3

geneburden_doe=c('OTgeneburden_doe','ravar_genburden')
dataset$geneburden_doe=0
dataset$geneburden_doe[rowSums(dataset[, geneburden_doe]) != 0]<-1

singlevar_doe=c('ravar_snp_predicted','genebass_predicted')
dataset$singlevar_doe=0
dataset$singlevar_doe[rowSums(dataset[, singlevar_doe]) > 0]<-1
dataset$singlevar_doe[rowSums(dataset[, singlevar_doe]) < 0]<--1


# Analysis 1a-1e: Code required to run five-fold cross-validation (CV) model to construct the side effect genetic priority score (SE-GPS) using the Open Targets dataset and then applied to OnSIDES and the all genes dataset (19,422 protein-coding genes and 470 drug side side effects).
# Analysis 1a: Get weights for each of the 5-CV Open target datasets 
## load specific libraries for regression model
library(Matrix)
library(Rcpp)
library(StanHeaders)
library(rstan)
library(lme4)
library(brms)

Mixedmodelreg_weights<-mclapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  OT_dataset=fread(paste0('OT_drugdataset_80_CV', CVsample, '.txt'),data.table=F) #80% training Open target dataset 
  #genetic features + phecode category covariates
  Predictors=paste(paste(geneticpredictors,collapse='+'),'+', paste(covariates,collapse='+'),collapse='+')
   ##Run mixed regression model across each CV training dataset to get weights, using drug name as a random effect, and a severity score to weight the side effect outcome
  mod_se_mixed <- glmer(as.formula(paste0(outcome, " ~ ", paste(Predictors), " +(1 | drugname)")), family = "binomial", data = OT_dataset,control = glmerControl(optimizer = "bobyqa"), weights=severity_score)
  results <- cbind.data.frame(beta = summary(mod_se_mixed)$coefficient[, c(1)], lowerCI = summary(mod_se_mixed)$coefficient[, c(1)] - 1.96 * summary(mod_se_mixed)$coefficient[, c(2)], upperCI = summary(mod_se_mixed)$coefficient[, c(1)] + 1.96 * summary(mod_se_mixed)$coefficient[, c(2)], P.val = summary(mod_se_mixed)$coefficient[, c(4)])
  write.table(results, paste0('Mixedmodel_weights_Opentargets_',CVsample,'.txt'), sep='\t',quote=F)
}, mc.cores=10)

#Analysis 1b: Use weights from Mixedmodel and create score using remaining 20% of data for each CV

Genescore_sum<-lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){

  Mixedmodel<-fread(paste0('Mixedmodel_weights_Opentargets_',CVsample,'.txt'),data.table=F) #weights from CV sample
  Mixedmodel=Mixedmodel %>% select(V1, beta) %>% rename(predictor=V1)
  Mixedmodel_weights=as.data.frame(t(Mixedmodel$beta))
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F) # OT remaining 20% test set 
  #extract gene, parentterm and genetic predictor columns
  OT_dataset_20_gene_phenotype=OT_dataset_20[, grepl(paste0('\\bgene\\b|parentterm|',  paste(geneticpredictors, collapse='|')), colnames(OT_dataset_20))] %>% distinct() 
  ## multiply betas from firth with values from 20% dataset
  Genescores_beta_genetic_values=mapply("*", OT_dataset_20_gene_phenotype[intersect(names(OT_dataset_20_gene_phenotype), names(Mixedmodel_weights))],
    Mixedmodel_weights[intersect(names(OT_dataset_20_gene_phenotype), names(Mixedmodel_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'), sep='\t', row.names=F, quote=F )

})

#Analysis 1c: Choose the CV sample with the max OR for validation set
max_filetype<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(samplecv){   
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  genescorefile=fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'), data.table=F)
  # combine genescore sum file with drug and mi data for each test set 
  Dataset_genescores=merge(OT_dataset_20[c('drugname','gene','parentterm','category','se','mi','senomi')] ,genescorefile, by=c('gene', 'parentterm') )
  #run logistic model 
  model1 =glm(as.formula(paste0('senomi ~ genescoresum + category ')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=samplecv,OR=exp(summary(model2)$coefficient[2,1]),lowerCI=exp(summary(model2)$coefficient[2,1]-(1.96* summary(model2)$coefficient[2,2])),upperCI=exp(summary(model2)$coefficient[2,1]+(1.96* summary(model2)$coefficient[2,2])),P.val=summary(model2)$coefficient[2,4]))
  return(mod_output)
  }))

## take CV with the max OR
max_filetype<-max_filetype[order(max_filetype$OR, decreasing=T),]
OT_CV=as.character(max_filetype$CV[1])

#Analysis 1d: Combine 5-CV 20% datasets to create one dataset

#scores with predictors 
Combine_genescores<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'),data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','parentterm','category','se','mi','senomi')], Genescore_sumfile_20 ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores1<-do.call(rbind.fill,Combine_genescores)
write.table(Combine_genescores1, gzfile(paste0('All_genescoresum_across_all_predictors_opentargets.txt.gz')), sep='\t',quote=F,row.names=F)

#scores without predictors 
Combine_genescores_nopredictors<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'.txt'),data.table=F)
  Genescore_sumfile_20_nopred<-Genescore_sumfile_20 %>% distinct(gene,parentterm,genescoresum )
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','parentterm','category','se','mi','senomi')], Genescore_sumfile_20_nopred ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_nopredictors1<-do.call(rbind.fill,Combine_genescores_nopredictors)
write.table(Combine_genescores_nopredictors1, gzfile(paste0('All_genescoresum_opentargets.txt.gz')), sep='\t',quote=F,row.names=F)

#Analysis 1e: Calculate GPS in OnSIDES dataset and all genes dataset using the OT mixedmodel weights from Analysis 1c which gave the max OR 

samplecv=OT_CV
Mixedmodel<-fread(paste0('Mixedmodel_weights_Opentargets_',samplecv,'.txt'),data.table=F)
Mixedmodel=Mixedmodel %>% select(V1, beta) %>% rename(predictor=V1)

Genescore_sum<-lapply(c('OnSIDES','Allgenes'), function(valdataset){
  Validation_dataset=fread(paste0(valdataset,'_drugdataset.txt'),data.table=F) #OnSIDES/allgenes dataset 
  Mixedmodel_weights=as.data.frame(t(Mixedmodel$beta))
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|',  paste(geneticpredictors, collapse='|')), colnames(Validation_dataset))] %>% distinct()
  ## multiply betas from firth with values from OnSIDES/Allgenes dataset
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Mixedmodel_weights))],
    Mixedmodel_weights[intersect(names(Validation_dataset_gene_phenotype), names(Mixedmodel_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('All_genescoresum_across_all_predictors_',valdataset,'.txt'), sep='\t', row.names=F, quote=F )
    #for OnSIDES add GPS sum to the OnSIDES-drug dataset with mi
  if(valdataset=='OnSIDES'){
  Genescores_beta_genetic_genept=Validation_dataset %>% distinct(gene, parentterm, genescoresum) 
  genescorefile_drugs=inner_join(Validation_dataset[c('drugname','gene','parentterm','category','se','mi','senomi')],Genescores_beta_genetic_genept, by=c('gene', 'parentterm'))
  write.table(genescorefile_drugs, paste0('All_genescoresum_across_drugs_',valdataset,'.txt'), sep='\t', row.names=F, quote=F )
  }
})

#Analysis 2a-2e: Create SE-GPS-DOE
###For score GOF annotated as -1, LOF annotated as +1 and missing/no evidence/neutral annotated as 0
### when calculating weights GOF/LOF recorded as 1

#Analysis 2a: Get weights for each of the 5-CV Open target datasets 
Mixedmodelreg_weights_doe<-mclapply(c(paste0('CVsample',rep(1:5))), function(CVsample){

  OT_dataset_doe=fread(paste0('OT_drugdataset_80_CV', CVsample, '_doe.txt'),data.table=F) #80% training Open target dataset. Restricted to drugs with activator/inhibitor mechanism 
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='GOF'] <- 1
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='LOF'] <- 1
  OT_dataset_doe[8:15][OT_dataset_doe[8:15]=='Neutral'] <-0
 
  #genetic features + phecode category covariates
  Predictors_doe=paste(paste(geneticpredictors_doe,collapse='+'), '+', paste(covariates,collapse='+'),collapse='+' )
   ##Run mixedmodel regression across each CV training dataset to get weights - don't take DOE into account here
  mod_se_mixed <- glmer(as.formula(paste0(outcome, " ~ ", paste(Predictors), " +(1 | drugname)")), family = "binomial", data = OT_dataset_doe,control = glmerControl(optimizer = "bobyqa"), weights=severity_score)
  results <- cbind.data.frame(beta = summary(mod_se_mixed)$coefficient[, c(1)], lowerCI = summary(mod_se_mixed)$coefficient[, c(1)] - 1.96 * summary(mod_se_mixed)$coefficient[, c(2)], upperCI = summary(mod_se_mixed)$coefficient[, c(1)] + 1.96 * summary(mod_se_mixed)$coefficient[, c(2)], P.val = summary(mod_se_mixed)$coefficient[, c(4)])
   write.table(results, paste0('Mixedmodel_weights_Opentargets_',CVsample,'_doe.txt'), sep='\t',quote=F)
}, mc.cores=10)

#Analysis 2b: Use weights from Mixedmodel and create score using remaining 20% of data for each CV
Genescore_sum<-lapply(c(paste0('CVsample',rep(1:5))), function(CVsample){
  Mixedmodel<-fread(paste0('Mixedmodel_weights_Opentargets_',CVsample,'_doe.txt'),data.table=F) #weights from CV sample
  Mixedmodel=Mixedmodel %>% select(V1, beta) %>% rename(predictor=V1)
  Mixedmodel_weights=as.data.frame(t(Mixedmodel$beta))
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F) # OT remaining 20% test set 
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='GOF'] <- -1 
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='LOF'] <-1
  OT_dataset_20[8:15][OT_dataset_20[8:15]=='Neutral'] <-0
  #extract gene, parentterm and genetic predictor columns
  OT_dataset_20_gene_phenotype=OT_dataset_20[, grepl(paste0('\\bgene\\b|parentterm|moa|',  paste(geneticpredictors_doe, collapse='|')), colnames(OT_dataset_20))] %>% distinct() 
  ## multiply betas from firth with values from 20% dataset
  Genescores_beta_genetic_values=mapply("*", OT_dataset_20_gene_phenotype[intersect(names(OT_dataset_20_gene_phenotype), names(Mixedmodel_weights))],
    Mixedmodel_weights[intersect(names(OT_dataset_20_gene_phenotype), names(Mixedmodel_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm","moa")])
  write.table(Genescores_beta_genetic, paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'), sep='\t', row.names=F, quote=F )

})

#Analysis 2c: Choose the CV sample with the max OR for validation set

max_filetype_doe<-do.call(rbind,lapply(c(paste0('CVsample',seq(1:5))), function(samplecv){   
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  genescorefile=fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'), data.table=F)
  # combine genescore sum file with drug and mi data for each test set 
  Dataset_genescores=merge(OT_dataset_20[c('drugname','gene','parentterm','category','moa','se','mi','senomi')] ,genescorefile, by=c('gene', 'parentterm') )
  #run logistic model 
  model1 =glm(as.formula(paste0('senomi ~ abs(genescoresum) + category ')), data=Dataset_genescores,family = 'binomial')
  mod_output<-rbind(cbind.data.frame(CV=samplecv,OR=exp(summary(model2)$coefficient[2,1]),lowerCI=exp(summary(model2)$coefficient[2,1]-(1.96* summary(model2)$coefficient[2,2])),upperCI=exp(summary(model2)$coefficient[2,1]+(1.96* summary(model2)$coefficient[2,2])),P.val=summary(model2)$coefficient[2,4]))
  return(mod_output)
  }))

## take CV with the max OR
max_filetype_doe<-max_filetype_doe[order(max_filetype_doe$OR, decreasing=T),]
OT_CV_doe=as.character(max_filetype_doe$CV[1])

#Analysis 2d: Combine 5-CV 20% datasets to create one dataset

#scores with predictors 
Combine_genescores_doe<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'),data.table=F)
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','parentterm','category','moa','se','mi','senomi')], Genescore_sumfile_20 ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_doe1<-do.call(rbind.fill,Combine_genescores_doe)
write.table(Combine_genescores_doe1, gzfile(paste0('All_genescoresum_across_all_predictors_opentargets_doe.txt.gz')), sep='\t',quote=F,row.names=F)

#scores without predictors 
Combine_genescores_nopredictors_doe<-lapply(paste0('CVsample', seq(1:5)), function(CVsample){
  OT_dataset_20=fread(paste0('OT_drugdataset_20_CV', CVsample, '_doe.txt'),data.table=F)
  Genescore_sumfile_20<-fread(paste0('Genescore_sum_across_predictor_Opentargets_',CVsample,'_doe.txt'),data.table=F)
  Genescore_sumfile_20_nopred<-Genescore_sumfile_20 %>% distinct(gene,parentterm,genescoresum )
  Genescore_sumfile_20_drugs<-inner_join(OT_dataset_20[c('drugname','gene','parentterm','category','moa','se','mi','senomi')], Genescore_sumfile_20_nopred ) %>% distinct()
  Genescore_sumfile_20_drugs$CV=CV
  return(Genescore_sumfile_20_drugs)
})
Combine_genescores_nopredictors_doe1<-do.call(rbind.fill,Combine_genescores_nopredictors_doe)
write.table(Combine_genescores_nopredictors_doe1, gzfile(paste0('All_genescoresum_opentargets_doe.txt.gz')), sep='\t',quote=F,row.names=F)

#Analysis 2e: Calculate GPS-D in OnSIDES dataset and all genes dataset using the OT firth weights from Analysis 2c which gave the max OR 

samplecv=OT_CV_doe
Mixedmodel_doe<-fread(paste0('Mixedmodel_weights_Opentargets_',samplecv,'_doe.txt'),data.table=F)
Mixedmodel_doe=Mixedmodel_doe %>% select(V1, beta) %>% rename(predictor=V1)

Genescore_sum<-lapply(c('OnSIDES','Allgenes'), function(valdataset){
  Validation_dataset=fread(paste0(valdataset,'_drugdataset_doe.txt'),data.table=F) #OnSIDES/allgenes dataset 
  Validation_dataset[8:15][Validation_dataset[8:15]=='GOF'] <- -1 
  Validation_dataset[8:15][Validation_dataset[8:15]=='LOF'] <-1
  Validation_dataset[8:15][Validation_dataset[8:15]=='Neutral'] <-0
  Mixedmodel_weights=as.data.frame(t(Mixedmodel_doe$beta))
  Validation_dataset_gene_phenotype=Validation_dataset[, grepl(paste0('\\bgene\\b|parentterm|moa|',  paste(geneticpredictors_doe, collapse='|')), colnames(Validation_dataset))] %>% distinct()
  ## multiply betas from firth with values from OnSIDES/Allgenes dataset
  Genescores_beta_genetic_values=mapply("*", Validation_dataset_gene_phenotype[intersect(names(Validation_dataset_gene_phenotype), names(Mixedmodel_weights))],
    Mixedmodel_weights[intersect(names(Validation_dataset_gene_phenotype), names(Mixedmodel_weights))])
  Genescores_beta_genetic=as.data.frame(cbind(gene_parentterm,Genescores_beta_genetic_values))
  #sum all predictor values for each gene - phenotype row
  Genescores_beta_genetic$genescoresum=rowSums(Genescores_beta_genetic[, !names(Genescores_beta_genetic) %in% c("gene", "parentterm")])
  write.table(Genescores_beta_genetic, paste0('All_genescoresum_across_all_predictors_',valdataset,'_doe.txt'), sep='\t', row.names=F, quote=F )
    #for OnSIDES add GPS sum to the OnSIDES-drug dataset with mi
  if(valdataset=='OnSIDES'){
  Genescores_beta_genetic_genept=Validation_dataset %>% distinct(gene, parentterm, genescoresum) 
  genescorefile_drugs=inner_join(Validation_dataset[c('drugname','gene', 'parentterm','category','moa','se','mi','senomi')],Genescores_beta_genetic_genept, by=c('gene', 'parentterm'))
  write.table(genescorefile_drugs, paste0('All_genescoresum_across_drugs_',valdataset,'_doe.txt'), sep='\t', row.names=F, quote=F )
  }
})
