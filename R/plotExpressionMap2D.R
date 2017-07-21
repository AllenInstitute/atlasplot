.plotExpressionMap2D <- function(expr, xVar, yVar, transformation=NULL, bgPar="lightblue",
sizeRange=c(0.4,6),sizeText=2,xlab=NULL,ylab=" ",pch=21,sizeLabel=1,
minIs0=FALSE,sampleLabel=names(expr),textShift=0.3,isXSubVar=TRUE,
ylines=TRUE,xlines=TRUE,colLine="grey",equalScale=TRUE,returnXY=FALSE,collapse=NULL,...){
# Function to get the XY coordinates, shapes, & sizes from yVar and xVar, then make a 2D plot. The 
#  function assumes that samples are ordered within xVar/yVar in some biologically reasonable way.
# "..." is for additional plotting parameters.  Use cex=___ for point size.
# Set minIs0=TRUE if you want to scale the minimum dot size to 0 rather than to minumum(expr)
# Typically, xVar = lobe and yVar = layer for these data, and sampleLabel = lobe subregions

    
    
    # sampleLabel same as ontology except when multiple subregions with same name
    # apply a transformation if supplied
    if (!is.null(transformation)) { 
        expr <- match.fun(transformation)(expr) 
    }

    # create copy of expr
    expr2 <- expr

    # create title for xlabel if none provided
    q = signif(quantile(expr),2)
    q = paste(paste(names(q),q,sep=": "),collapse="  |  ")
    if(is.null(xlab)) { 
        xlab = q 
    }

    # select the unique factor levels of xVar... 
    if(class(xVar)=="factor") {
        xVars = levels(xVar)
    } else {
        xVars = unique(xVar)
    }
    
    # ...and yVar
    if(class(yVar)=="factor") {
        yVars = levels(yVar)
    } else {
        yVars = unique(yVar)
    }

    # reverse the order of yVars    
    yVars = yVars[length(yVars):1]
    
    # create  copy of factor label before coercing to character
    sampleLabelF = sampleLabel  # names of expr (column names)
    sampleLabel = as.character(sampleLabel)

    # set x and y pos to NULL and instantiate a counter
    xPos  <- yPos <- NULL
    count=0
    
    # load the WGCNA library
    # library(WGCNA)
    # detach("package:WGCNA")
    
    # places a position based on which of the catagories of yVars it fall in; 
    # Ex: If yVar[i] is the ninth factor in yVars it's `9` - yVars order set lines 32-44
    # yPos has length of yVar
    for (i in 1:length(yVar)) {
        yPos = c(yPos, which(is.element(yVars,yVar[i])))
    }
    
    # see yPos above
    for (i in 1:length(xVar)) {
        xPos = c(xPos, which(is.element(xVars,xVar[i])))
    }
    
    # equal scale means all of the subplots created have the same scale
    # create a zero vector; has to do if there's only one catagory
    if(equalScale) xPos = xPos-xPos
    
    # crazy loop
    for (i in 1:length(xVars)){ 
        if(isXSubVar){
            # if X is the subvariable(?) then for each one grab the ones associated
            # with it; for instance grab all of the frontal positions
            # unique values of things that belong to that xVars (subcats of f)
            subs   = unique(sampleLabel[is.element(xVar,xVars[i])])
            if(class(sampleLabelF)=="factor") {
                subs = intersect(levels(sampleLabelF),subs)
            }

            # values are physical positions are the plotting
            values = ((1:length(subs))/(length(subs)+1)-0.5)

            if(equalScale) {
                values = count + 1:length(subs) # if the scales are same positions same
            }
            names(values) = subs  # assigns the unique subs to the nambes?
        }

        # for each of the xVars were gunna do a thing!!!
        # create the xPos from values. If there are subcatagoires we'll just use the 
        # values we computed before and add them to xPos; else we'll create them based
        # on the values we keep for i and j
        for (j in 1:length(yVars)){
            # keep ones that are both the right x cat and y cat (i and j)
            keep = (is.element(xVar,xVars[i]))&(is.element(yVar,yVars[j]))

            if(!isXSubVar) {
                # if xVars are not subvariables then we need to decide thir positions
                xPos[keep] = xPos[keep]+((1:sum(keep))/(sum(keep)+1)-0.5)
            } else {
                xPos[keep] = xPos[keep]+values[sampleLabel[keep]]
            }
            count = max(xPos)+1
        }
    } 
#------------------------------NO COMMENTS BELOW----------------------------------------#
    if (!is.null(collapse)) {
        combine = paste(xVar,yVar,sampleLabel)
        for (comb in unique(combine)) {
            expr[combine==comb] = match.fun(collapse)(expr2[combine==comb])
        }
    }

    if(equalScale) {
        xPos = 0.5+length(xVars)*xPos/max(xPos)
    }

    if(is.null(collapse)) {
        expr = expr2
    }

    sizes = expr
    if(!minIs0) {
        sizes = sizes-min(sizes)
    }

    if(!isXSubVar) { # This will not work unless a subregion variable is entered s
        equalScale = FALSE
    } 
    sizes = sizes*(sizeRange[2]-sizeRange[1])/max(sizes)
    sizes = sizes+sizeRange[1]
    
    # set exprcolors
    exprColors = WGCNA::numbers2colors(expr,signed=FALSE)
    
    par(bg=bgPar)
    
    plot(0,0,col=bgPar,xlim=c(0,max(xPos)),ylim=c(0,max(yPos)+1),axes=FALSE,xlab=xlab,ylab=ylab,...)
    
    for (i in 1:length(expr)) {
      points(xPos[i],yPos[i],col=exprColors[i],bg=exprColors[i],pch=pch,cex=sizes[i])
    }
    
    for (i in 1:length(yVars)) {
        text(0,i,yVars[i],cex=sizeText,col="black")
    }
    
    if(!is.null(sampleLabel)) {
        for (i in 1:length(xPos)) {
            text(xPos[i], ifelse(isXSubVar,0,yPos[i]-textShift),sampleLabel[i],
                 cex=sizeLabel,srt=90,col="black")
        }
    }

    if(ylines) {
        for (i in 0:length(yVars)) abline(h=i+0.5,col=colLine)
    }
    
    if(xlines) {
        abline(v=min(xPos)-0.25,col=colLine); abline(v=max(xPos)+0.25,col=colLine)
        for (i in 1:(length(xVars)-1)) {
            abline(v=(0.5*(max(xPos[xVar==xVars[i]])+min(xPos[xVar==xVars[i+1]]))),
                   col=colLine)
        }
     }

    headerPos = rep(0,length(xVars))
    for (i in 1:length(xVars)) {
       headerPos[i] = 0.5*(max(xPos[xVar==xVars[i]])+min(xPos[xVar==xVars[i]]))
    }
    
    for (i in 1:length(xVars)) {
        text(headerPos[i],max(yPos)+1,xVars[i],cex=sizeText,col="black")
    }
    
    if (returnXY) {
        return(list(x=xPos,y=yPos,expr=expr))
    }
}