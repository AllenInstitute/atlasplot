#    Copyright (C) 2017 Allen Institute
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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