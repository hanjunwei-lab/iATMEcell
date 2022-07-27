
##' @title Draw a heat map.
##' @description The function "plotHeatmap" is used to draw a heat map of marker genes.
##' @param Regress.list The result of function "RiskRegressModel".
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param cut.off Samples will be grouped according to this threshold. If not specified, the sample will be grouped according to the median risk score.
##' @param p.cutoff Statistical significance threshold of cox regression analysis, based on which to determine the genes used to draw the graph.
##' @param cluster.rows Boolean values determining if rows should be clustered or hclust object.
##' @param cluster.cols Boolean values determining if columns should be clustered or hclust object.
##' @param bk A numeric vector that covers the range of values. Users could adjust color depth through this parameter.
##' @param show.rownames Boolean specifying if row names are be shown.
##' @param show.colnames Boolean specifying if column names are be shown.
##' @param ann_colors Vector of colors for specifying the color of column annotation.
##' @param col Vector of colors used in heat map.
##' @return A heat map
##' @importFrom pheatmap pheatmap
##' @importFrom grDevices colorRampPalette
##' @usage plotHeatmap(Regress.list,ExpData,cut.off=NULL,p.cutoff=0.05,cluster.rows=F,
##' cluster.cols=F,bk=c(-2.4,2.3),show.rownames=T,show.colnames=F,
##' ann_colors=c("#FFAA2C","#2CBADA"),col=c("#2A95FF","#FF1C1C"))
##' @export
##' @examples
##' library(pheatmap)
##' library(survival)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,clinical=clinicaldata,
##'       p.cutoff=0.05)
##' plotHeatmap(Regress.list=R.result,ExpData=GEP,p.cutoff=0.05)



plotHeatmap<-function(Regress.list,ExpData,cut.off=NULL,p.cutoff=0.05,cluster.rows=F,
         cluster.cols=F,bk=c(-2.4,2.3),show.rownames=T,show.colnames=F,
         ann_colors=c("#FFAA2C","#2CBADA"),col=c("#2A95FF","#FF1C1C")){
  GEP<-ExpData
  cox.result<-Regress.list$cox.result
  riskscore.frame<-Regress.list$riskscore.frame
  sig.genes<-rownames(cox.result[which(cox.result[,"p.value"]<p.cutoff),])
  GEP <- GEP[match(sig.genes,rownames(GEP)), ]

  if(is.null(cut.off)==T){
    riskscore.frame[which(riskscore.frame[,'riskscore']>=median(riskscore.frame[,'riskscore'])),'riskscore']='High.risk'
    riskscore.frame[which(riskscore.frame[,'riskscore']!='High.risk'),'riskscore']='Low.risk'
  }else{
    riskscore.frame[which(riskscore.frame[,'riskscore']>= cut.off),'riskscore']='High.risk'
    riskscore.frame[which(riskscore.frame[,'riskscore']!='High.risk'),'riskscore']='Low.risk'
  }

  colann<-riskscore.frame[,c(1,4)]
  colann.high<-colann[which(colann[,2]=='High.risk'),]
  colann.low<-colann[which(colann[,2]=='Low.risk'),]
  colann<-rbind(colann.high,colann.low)
  colann.1<-as.data.frame(colann[,-1])
  rownames(colann.1)<-rownames(colann)
  colnames(colann.1)<-'Group'
  GEP<-GEP[,match(rownames(colann.1),colnames(GEP))]

  breaks<-c(seq(bk[1],-0.1,by=0.1),seq(0,bk[2],by=0.1))
  ann_colors=list(Group=c(High.risk=ann_colors[1],Low.risk=ann_colors[2]))
  pheatmap(GEP,
           scale = 'row',
           color = colorRampPalette(c(col[1], "#FFFFFF", col[2]))(50),
           breaks = breaks,
           cluster_rows=cluster.rows,
           cluster_cols=cluster.cols,
           annotation_col =colann.1,
           annotation_colors = ann_colors,
           show_rownames=show.rownames,
           show_colnames=show.colnames
           #main=""
  )

}
