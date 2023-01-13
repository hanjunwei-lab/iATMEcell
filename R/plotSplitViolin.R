


geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL,
                               trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,GeomSplitViolin=NULL) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position,
        show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale,
                                                                            draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  #datac <- rename(datac, c("mean" = measurevar))
  colnames(datac)[4]<-measurevar
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

##' @title Draw a split violin plot.
##' @description The function "plotSplitViolin" is used to draw a split violin plot of gene expression.
##' @param Regress.list The result of function "RiskRegressModel".
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param gene.name A gene symbol in inputted gene expression profile.
##' @param method A character string indicating which method to be used for comparing means. The default method is "t.test". Other three methods are "wilcox.test", "anova" and "kruskal.test".
##' @param compare.label A character string specifying label type. Allowed values include "p.signif" (shows the significance levels), "p.format" (shows the formatted p value).
##' @param col Vector of colors used to specify the color of different groups.
##' @param x.ceiling The maximum value of the y axis.
##' @param y.lab Setting the title of the y-axis.
##' @param x.lab Setting the title of the x-axis.
##' @param title Setting the title
##' @return A split violin plot
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 layer
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 geom_errorbar
##' @importFrom ggplot2 scale_fill_manual
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 position_dodge
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 Geom
##' @importFrom ggplot2 Stat
##' @importFrom ggplot2 Position
##' @importFrom ggplot2 Scale
##' @importFrom stats sd
##' @importFrom stats qt
##' @importFrom grid grobTree
##' @importFrom reshape2 melt
##' @importFrom plyr ddply
##' @importFrom plyr arrange
##' @importFrom ggpubr stat_compare_means
##' @usage plotSplitViolin(Regress.list,ExpData,gene.name,method="t.test",
##' compare.label="p.signif",col=c("#E69F00", "#56B4E9"),x.ceiling=15,
##' y.lab="Gene Expression",x.lab=NULL,title=NULL)
##' @export
##' @examples
##' library(ggplot2)
##' library(reshape2)
##' library(plyr)
##' library(ggpubr)
##' library(survival)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,clinical=clinicaldata,
##'         p.cutoff=0.05)
##' plotSplitViolin(Regress.list=R.result,ExpData=GEP,gene.name="PDCD1")



plotSplitViolin<-function(Regress.list,ExpData,gene.name,method="t.test",compare.label="p.signif",
                          col=c("#E69F00", "#56B4E9"),x.ceiling=15,
                          y.lab="Gene Expression",x.lab=NULL,title=NULL){

  GeomSplit<-GetExampleSet("GeomSplitViolin")

  riskscore.frame<-Regress.list$riskscore.frame
  riskscore.frame[which(riskscore.frame[,'riskscore']>=median(riskscore.frame[,'riskscore'])),'riskscore']='High risk'
  riskscore.frame[which(riskscore.frame[,'riskscore']!='High risk'),'riskscore']='Low risk'
  riskscore.frame<-riskscore.frame[,c(1,4)]
  GeneEXP<-as.data.frame(ExpData)[,match(rownames(riskscore.frame),colnames(ExpData))]
  GeneEXP<-GeneEXP[which(rownames(GeneEXP)==gene.name),]
  GeneEXP.new<-data.frame(t(GeneEXP))
  GeneEXP.new$sample = rownames(GeneEXP.new)
  GeneEXP.new <- merge(GeneEXP.new, riskscore.frame,by.x = "sample",by.y = 0)
  GeneEXP.new = melt(GeneEXP.new)
  colnames(GeneEXP.new) = c("sample","subject","group", "gene","expression")
  GeneEXP.new<-GeneEXP.new[which(GeneEXP.new[,5]<=x.ceiling),]

  Data_summary <- summarySE(GeneEXP.new, measurevar="expression", groupvars=c("group","gene"))

  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"),
                   axis.text = element_text(size= 12,color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(hjust = 1 ),
                   #panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  )

  ggplot(GeneEXP.new,aes(x= gene ,y= expression,fill= group))+
    geom_split_violin(trim= F,color="white",scale = "area",GeomSplitViolin=GeomSplit) +
    geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
               position=position_dodge(0.5),size= 2)+
    geom_errorbar(data = Data_summary,aes(ymin = expression-ci,
                                          ymax= expression+ci),
                  width= 0.05,
                  position= position_dodge(0.5),
                  color="black",
                  alpha = 0.8,
                  size= 0.5) +
    scale_fill_manual(values = c(col[1], col[2]))+
    labs(y=y.lab,x=x.lab,title = title) +
    theme_bw()+ mytheme +
    stat_compare_means(aes(group = group),
                       label = compare.label,
                       method = method,
                       label.y = max(Data_summary$expression),
                       hide.ns = T)
}




