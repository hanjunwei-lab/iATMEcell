##' @title Draw a forest plot.
##' @description The function "plotforest" is used to draw a forest plot according to the result of cox analysis from function "RiskRegressModel".
##' @param Regress.list The result of function "RiskRegressModel".
##' @param p.cutoff Statistical significance threshold of cox regression analysis, based on which to determine the genes used to draw the graph.
##' @param g.pos A number to control the position of the graph element in forestplot.
##' @param b.size A number to control the box size.
##' @param col Vector of colors including three color code which are corresponding to box, box line and reference line.
##' @param lwd.zero A number to control the thickness of the reference line.
##' @param lwd.ci A number to control the thickness of the box line.
##' @param x.lab Setting the title.
##' @return A forest plot
##' @importFrom forestplot forestplot
##' @importFrom forestplot fpColors
##' @importFrom forestplot fpTxtGp
##' @importFrom grid gpar
##' @usage plotforest(Regress.list,p.cutoff=0.05,g.pos=2,b.size=3,
##'         col=c("#FE0101","#1C61B6","#A4A4A4"),
##'         lwd.zero=2,lwd.ci=3,x.lab="Hazard Ratio Plot")
##' @export
##' @examples
##' library(forestplot)
##' library(survival)
##' #Obtain input data
##' GEP<-GetExampleSet('GEP')
##' clinicaldata<-GetExampleSet('clinicaldata')
##' #Run the function
##' R.result<-RiskRegressModel(cellname='NK cells',ExpData=GEP,clinical=clinicaldata,
##'      p.cutoff=0.05)
##' plotforest(Regress.list=R.result,p.cutoff=0.05)



plotforest<-function(Regress.list,p.cutoff=0.05,g.pos=2,b.size=3,col=c("#FE0101","#1C61B6","#A4A4A4"),
                     lwd.zero=2,lwd.ci=3,x.lab="Hazard Ratio Plot"){
  cox.result=Regress.list$cox.result
  cox.result=cox.result[which(cox.result[,'p.value']<p.cutoff),]
  HR=cox.result[,2:4]
  hr=HR[,"HR"]
  hrLow=HR[,"HR.95L"]
  hrHigh=HR[,"HR.95H"]
  beta<-cox.result[,"beta"]
  pVal=cox.result[,"p.value"]

  tabletext <- list(c(NA, rownames(HR)),
                    append("HR (95% CI for HR)",paste(hr,paste("(",hrLow,"-",hrHigh,")",sep = ''),sep = ' ')),
                    append("Cox.Beta", beta),
                    append("P.value", pVal))   #定义图片文字


  boxsize<-c(NA,(cox.result$HR)/b.size)


  forestplot(tabletext,
             rbind(rep(NA, 3), HR),
             graph.pos=g.pos,#设置箱线图的位置
             #调节箱、线和参考线的颜色
             col=fpColors(box = col[1],lines = col[2],zero = col[3]),
             xlog=F,
             zero=1,#修改了T
             lwd.zero=lwd.zero,#修改参考线的粗细
             lwd.ci=lwd.ci,#修改箱线的粗细
             boxsize=boxsize,#修改箱子的大小
             xlab=x.lab,#设置标题
             ci.vertices=T,#是否添加箱线两边的竖线
             ci.vertices.height=0.2,#箱线两边竖线的长短
             hrzl_lines=list("2"=gpar(lwd=2,col="gray")),
             txt_gp=fpTxtGp(label = gpar(cex=1),#设置字号大小
                            ticks=gpar(cex=0.8),#设置刻度字号大小
                            xlab=gpar(cex = 1.5),)#设置x轴小标题字号大小
  )
}



