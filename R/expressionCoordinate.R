.plotExpressionCoordinates2Db <- function(expr, xCoor, yCoor, textVar, textCol="black",pch=21,
  bgPar="lightgrey",sizeRange=c(0.4,10),textRange=c(1,3),xlab=NULL,ylab=" ",minIs0=TRUE,...){
 q = signif(quantile(expr),2)
 q = paste(paste(names(q),q,sep=": "),collapse="  |  ")
 if(is.null(xlab)) xlab = q
 sizes = expr;   if(!minIs0) sizes = sizes-min(sizes);  sizeText = sizes
 sizes = sizes*(sizeRange[2]-sizeRange[1])/max(sizes);  sizes = sizes+sizeRange[1]
 sizeText = sizeText*(textRange[2]-textRange[1])/max(sizeText);  sizeText = sizeText+textRange[1]
 exprColors = .numbers_to_colors(expr, c("white", "red"))
 plot(xCoor,yCoor,col=bgPar,xlab=xlab,ylab=ylab,...)
 for (i in 1:length(expr))
  points(xCoor[i],yCoor[i],col=exprColors[i],bg=exprColors[i],pch=pch,cex=sizes[i])
 text(xCoor,yCoor,textVar,cex=sizeText,col=textCol)
}