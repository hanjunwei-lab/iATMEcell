
##' @title Constructing the cox risk regression model with cell's marker genes
##' @description This function is used to perform regression analysis and build risk regression models.
##' @param cellname A cell whose marker genes will be used to perform regression analysis. The format of the entered cell name should refer to the cell information we provide.
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param clinical A dataframe with three columns which are "sample" (sample id),"status" (survival status of samples, "0" represents live and "1" represents dead) and "time" (survival time of samples).
##' @param marker A character vector composed with marker genes. If you does not want to use the marker genes provided by us, you can specify the marker genes you need with this parameter.
##' @param method This parameter specifies the method of regression analysis. "method=cox": Only univariate regression analysis was performed, and the model was constructed with coefficients. "method=lasso"(default): The significant variables in univariate analysis were used for LASSO regression analysis, and the coefficients of LASSO analysis were used to construct the model.
##' @param p.cutoff Statistical significance threshold for regression analysis (default: 0.05).
##' @return A list with two dataframes which are riskscores of samples and result of cox regression analysis respectively.
##' @importFrom survival coxph
##' @importFrom survival Surv
##' @importFrom stats as.formula
##' @importFrom glmnet cv.glmnet
##' @importFrom stats coef
##' @usage RiskRegressModel(cellname,ExpData,clinical,marker=NULL,method = "lasso",p.cutoff=0.05)
##' @details In the default method, users can specify a cell, and then the function will perform cox regression analysis on expression value of marker genes of the cell and survival data. Statistical significant genes will be selected for further lasso regression analysis. Finally, the lasso regression coefficients were used to weight the gene expression values to calculate the risk score for samples.
##' @export
##' @examples
##' library(survival)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,
##'   clinical=clinicaldata,method = "lasso",p.cutoff=0.05)


RiskRegressModel<-function(cellname,ExpData,clinical,marker=NULL,method = "lasso",p.cutoff=0.05){
  if(is.null(marker)){
    TMEcellinfo<-GetExampleSet('TMEcellinfo')
    markergenes<-as.data.frame(TMEcellinfo)[which(TMEcellinfo[,1]== cellname),4]
    markergenes<-as.character(markergenes)
    markergenes<-unlist(strsplit(markergenes,','))
  }else{
    markergenes<-marker
  }

  congenes<-intersect(markergenes,rownames(ExpData))
  markergep<-ExpData[congenes,]
  markergep<-t(markergep)

  samplename<-intersect(clinical[,"sample"],rownames(markergep))
  markergep<-markergep[samplename,]
  clinicaldata<-clinical[samplename,]

  markergep<-markergep[match(clinicaldata[,"sample"],rownames(markergep)),]
  survframe<-cbind(clinicaldata,markergep)

  covariates <- colnames(survframe[,4:length(survframe[1,])])

  univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(time, status) ~',"`",x,"`",sep = "")))

  univ_models <- lapply( univ_formulas, function(x){coxph(x, survframe)})

  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           #wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2)
                           HR <-signif(x$coef[2], digits=2)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           #HR <- paste0(HR, " (",
                           #            HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR,HR.confint.lower,HR.confint.upper, p.value)
                           names(res)<-c("beta", "HR", "HR.95L", "HR.95H","p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })

  result_cox <- data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  result_cox[,'beta']<-as.numeric(as.character(result_cox[,'beta']))
  result_cox[,'p.value']<-as.numeric(as.character(result_cox[,'p.value']))
  result_cox<-result_cox[which(result_cox[,'p.value']<p.cutoff),]
  sig.genes<-rownames(result_cox)
  sig.survframe<-cbind(survframe[,1:3],survframe[,which(colnames(survframe)%in%sig.genes)])

  if(method == "lasso"){
    x=as.matrix(sig.survframe[,-c(1,2,3)])

    y <- data.matrix(Surv(sig.survframe$time,sig.survframe$status))
    set.seed(4)
    lasso_basal <-cv.glmnet(x,y,family = "cox",alpha = 1,type.measure = 'deviance')

    coeff<-coef(lasso_basal,s=lasso_basal$lambda.min)
    coeff1<-as.vector(coeff)
    names(coeff1)<-rownames(coeff)
    coeff1<-coeff1[which(coeff1!=0)]
    final_character1<-rownames(coeff)[which(as.numeric(coeff)!=0)]
    sig.survframe=cbind(sig.survframe[,c(1,2,3)],
                        sig.survframe[,colnames(sig.survframe)%in%final_character1])

    riskscore<-c()
    for (i in 1:dim(sig.survframe)[1]) {
      rs<-sum(as.numeric(sig.survframe[i,4:dim(sig.survframe)[2]])*coeff1)
      riskscore<-c(riskscore,rs)
    }

    risk.result<-sig.survframe[,1:3]
    risk.result$riskscore<-riskscore

    result.list<-list()
    result.list[[1]]<-risk.result
    result_cox<-result_cox[match(final_character1,rownames(result_cox)),]
    result_cox$lasso.coef<-coeff1
    result_cox<-result_cox[order(result_cox$HR),]
    result.list[[2]]<-result_cox
  }

  if(method == "cox"){
    riskscore<-c()
    for (i in 1:dim(sig.survframe)[1]) {
      rs<-sum(as.numeric(sig.survframe[i,4:dim(sig.survframe)[2]])*result_cox[match(colnames(sig.survframe)[4:dim(sig.survframe)[2]],rownames(result_cox)),'beta'])
      riskscore<-c(riskscore,rs)
    }

    risk.result<-sig.survframe[,1:3]
    risk.result$riskscore<-riskscore

    result.list<-list()
    result.list[[1]]<-risk.result
    result_cox<-result_cox[order(result_cox$HR),]
    result.list[[2]]<-result_cox
  }


  names(result.list)<-c('riskscore.frame','cox.result')
  return(result.list)
}
