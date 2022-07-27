##' @title Draw a Kaplan-Meier curve.
##' @description The function "plotKMcurve" is used to draw the Kaplan-Meier curve according to the riskscore of samples from function "RiskRegressModel".
##' @param Regress.list The result of function "RiskRegressModel".
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param risk.table TRUE or FALSE specifying whether to show or not the risk table. Default is TRUE.
##' @param labs A character vector for specifying legend labels.
##' @param title legend title
##' @param line.col Vector of colors for specifying the color of curve.
##' @return Kaplan-Meier curve
##' @importFrom survival survfit
##' @importFrom survminer ggsurvplot
##' @usage plotKMcurve(Regress.list,ExpData,risk.table=TRUE,labs=c("High risk","Low risk"),
##'                    title="Group",line.col=c("#FFAA2C", "#2CBADA"))
##' @export
##' @examples
##' library(survival)
##' library(survminer)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,clinical=clinicaldata,
##'        p.cutoff=0.05)
##' plotKMcurve(Regress.list=R.result,ExpData=GEP)


plotKMcurve<-function(Regress.list,ExpData,risk.table=TRUE,labs=c("High risk","Low risk"),title="Group",
                      line.col=c("#FFAA2C", "#2CBADA")){
  riskscore.frame<-Regress.list$riskscore.frame
  riskscore.frame[which(riskscore.frame[,'riskscore']>=median(riskscore.frame[,'riskscore'])),'riskscore']='high'
  riskscore.frame[which(riskscore.frame[,'riskscore']!='high'),'riskscore']='low'

  ggsurvplot(survfit(Surv(time, status) ~ riskscore, data = riskscore.frame),
             data = riskscore.frame,pval = T,risk.table = risk.table,
             legend.labs=labs,
             legend.title=title,
             palette=line.col,
             censor.size=3)
}

