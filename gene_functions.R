# functions to make visualizations of unsupervised clustering methods
# Jessica Minnier
# Updated: October 19, 2016
# 
# Adapted from DESeq2 vignette - MI Love, S Anders, W Huber
# 
require(ggplot2)
require(RColorBrewer)
require(ggrepel)
require(ggthemes)

# annotation_row is a data.frame

#' gene_pheatmap
#'
#' @param exprdat 
#' @param sampleid 
#' @param annotation_row 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' emat = matrix(rnorm(100),25,4)
#' samp = c("H1","H2","C1","C2")
#' ann = data.frame("fac1"=rep(c("M","N"),2),"fac2"=c("A","B","C","C"))
#' gene_pheatmap(emat,samp,ann)
#' gene_pheatmap(emat,samp,ann, dist="correlation")
#' gene_pheatmap(emat,samp,ann, dist="correlation",display_numbers=TRUE)
gene_pheatmap <- function(exprdat,sampleid,annotation_row=NULL, dist="euclidean",
                          display_numbers=FALSE,number_color="grey30",
                          cor_method = "pearson") {
  if(dist=="euclidean") {
    sampleDists <- dist(t(exprdat))
    sampleDistMatrix <- as.matrix(sampleDists)
  }else if(dist=="correlation") {
    sampleDists <- cor(exprdat,use = "pairwise.complete.obs",method=cor_method)
    sampleDistMatrix <- as.matrix(sampleDists)
    sampleDists = as.dist(1-sampleDists) # for distance
  }else{ stop("'dist' must be 'euclidean' or 'correlation'") }
  
  rownames(sampleDistMatrix) <- sampleid
  if(!is.null(annotation_row)) rownames(annotation_row) <- sampleid
  colnames(sampleDistMatrix) <- NULL
  colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  if(dist=="correlation") colors <- rev(colors)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     annotation_row=annotation_row,
                     col=colors,
                     display_numbers=display_numbers,
                     number_color=number_color)
}

#' gene_pheatmap
#'
#' @param exprdat 
#' @param sampleid 
#' @param annotation_row 
#'
#' @return
#' @export
#'
gene_heatmap <- function(exprdat,sampleid,annotation_row=NULL,
                         display_numbers=FALSE,number_color="grey30",mypallette = "RdBu",
                         ...) {
  exprdat_t = t(exprdat)
  if(!is.null(annotation_row)) rownames(annotation_row) <- sampleid
  colnames(exprdat) <- NULL
  #colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "PuOr")) )(255)
  
  paletteLength <- 100
  myColor  <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, mypallette)))(paletteLength)
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(exprdat), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(exprdat)/paletteLength, max(exprdat), length.out=floor(paletteLength/2)))
  
  pheatmap::pheatmap(exprdat_t,
                     annotation_row=annotation_row,
                     color=myColor, breaks=myBreaks,
                     display_numbers=display_numbers,
                     number_color=number_color,...)
}

# add filter on columns (samples)
# gene x sample matrix like in array data
# 
#' gene_pcaplot
#'
#' @param exprdat 
#' @param sampleid 
#' @param groupdat 
#' @param colorfactor 
#' @param shapefactor 
#' @param plot_sampleids 
#' @param pcnum 
#' @param plottitle 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' emat = cbind(matrix(rnorm(400,5),ncol=2),matrix(rnorm(400,-3),ncol=2))
#' samp = c("H1","H2","C1","C2")
#' grdat = data.frame("fac1"=rep(c("M","N"),2),"fac2"=c("A","B","C","C"))
#' gene_pcaplot(exprdat=emat,sampleid=samp,groupdat=grdat,colorfactor="fac1")
#' gene_pcaplot(exprdat=emat,sampleid=samp,
#'    groupdat=grdat,colorfactor="fac1",
#'    shapefactr="fac2")
gene_pcaplot <- function(exprdat,sampleid,groupdat=NULL,colorfactor=NULL,shapefactor=NULL,
                         plot_sampleids=TRUE, pcnum=1:2, plottitle = "PCA Plot") {
  #adapted from DESeq2:::plotPCA.DESeqTransform
  if((!is.null(groupdat))&(!is.null(colorfactor))) {
    if(!colorfactor%in%colnames(groupdat)) stop("colorfactor must be a column name of groupdat")
  }
  if((!is.null(groupdat))&(!is.null(shapefactor))) {
    if(!shapefactor%in%colnames(groupdat)) stop("shapefactor must be a column name of groupdat")
    if(plot_sampleids) {
      warning("plot_sampleids=TRUE and shapefactor is specified, changing plot_sampleids=FALSE to use shapefactor")
      plot_sampleids=FALSE 
    }
    
  }
  pca <- prcomp(t(exprdat))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if(is.null(groupdat)) groupdat = data.frame("group"=rep(1,ncol(exprdat)))
  intgroup = colnames(groupdat)
  allgroup <- if (length(intgroup) > 1) {
    factor(apply(groupdat, 1, paste, collapse = ":"))
  }else{allgroup <- intgroup}
  d <- data.frame(PC1 = pca$x[, pcnum[1]], PC2 = pca$x[, pcnum[2]], uniquegroup = allgroup, 
                  groupdat, name = sampleid)
  percentVar <- round(100 * percentVar)
  if(is.null(colorfactor)) {d$color=as.factor(1)}else{
    colnames(d)[colnames(d)==colorfactor] <- "color"}
  if(is.null(shapefactor)) {d$shape=as.factor(1)}else{
    colnames(d)[colnames(d)==shapefactor] <- "shape"
  }
  if(identical(shapefactor,colorfactor)) {d$shape = d$color}
  p <- ggplot(d, aes(PC1, PC2, color=color, shape=shape, size=3)) 
  if(plot_sampleids) {
    p <- p + geom_text(aes(label=name,size=10))
  }else{
    p <- p + geom_point()
  }
  
  if(!is.null(colorfactor)) {
    p <- p + guides(color=guide_legend(title=colorfactor))
  }else {
    p <- p + guides(color = "none")
  }
  if(!is.null(shapefactor)) {
    p <- p + guides(shape=guide_legend(title=shapefactor))
  }else{
    p <- p + guides(shape = "none")
  }
  p <- p + guides(size= "none") + theme_bw() + 
    xlab(paste0("PC",pcnum[1],": ",percentVar[pcnum[1]],"% variance")) +
    ylab(paste0("PC",pcnum[2],": ",percentVar[pcnum[2]],"% variance")) + ggtitle(plottitle)
  
  return(p)
}


#gene_pcaplot(intensdat_combat,sampleid,sampdat[,c("group","batch")],colorfactor="group",shapefactor="batch")
#

gene_volcanoplot_old  <- function(resdat,log2fc.cut=2,padj.cut = 0.05,pval.cut = 0.05,
                                  pval.name="pval",padj.name="padj",log2fc.name="log2fc",
                                  main="Volcano Plot",...) {
  
  check_colnames = !(c(pval.name, padj.name, log2fc.name)%in%colnames(resdat))
  if(sum(check_colnames)>0) stop(
    paste("The following argument values are not column names of resdat:\n",
          paste0(c("pval.name", "padj.name", "log2fc.name")[check_colnames],collapse=" ") )
  )
  
  resdat = resdat[,match(c(pval.name,padj.name,log2fc.name),colnames(resdat))]
  colnames(resdat) = c("pval","padj","log2fc")
  if(is.null(xlim)) {xlim = range(resdat$log2fc)}
  if(is.null(ylim)) {ylim = range(-log10(resdat$pval))}
  with(resdat, plot(log2fc, -log10(pval), pch=20, main=main,...))
  legend("topright",pch=20,legend = c(paste0("p<",pval.cut),
                                      paste0("padj<",padj.cut),
                                      paste0("|log2fc|>",log2fc.cut),
                                      paste0("p<",pval.cut," & |log2fc|>",log2fc.cut),
                                      paste0("padj<",padj.cut," & |log2fc|>",log2fc.cut)),
         col=c("grey","red","orange","purple","green"))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  # grey if p < .05 but not adj significant, purple if log2fc >1 and p<.05
  with(subset(resdat, pval<pval.cut ), points(log2fc, -log10(pval), pch=20, col="grey"))
  with(subset(resdat, padj<padj.cut ), points(log2fc, -log10(pval), pch=20, col="red"))
  with(subset(resdat, abs(log2fc)>log2fc.cut), points(log2fc, -log10(pval), pch=20, col="orange"))
  with(subset(resdat, pval<pval.cut & abs(log2fc)>log2fc.cut), points(log2fc, -log10(pval), pch=20, col="purple"))
  with(subset(resdat, padj<padj.cut & abs(log2fc)>log2fc.cut), points(log2fc, -log10(pval), pch=20, col="green"))
}




gene_volcanoplot <- function(resdat,
                             geneids = NULL,
                             log2fc.cut = 2,
                             padj.cut = 0.05,
                             pval.cut = 0.05,
                             pval.name = "pval",
                             padj.name = "padj",
                             log2fc.name = "log2fc",
                             id.name = "geneid",
                             main = "Volcano Plot",
                             sel_genes = NULL,
                             interactive = FALSE,
                             ...) {
  
  
  check_colnames = !(c(pval.name, padj.name, log2fc.name,id.name)%in%colnames(resdat))
  if(sum(check_colnames)>0) stop(
    paste("The following argument values are not column names of resdat:\n",
          paste0(c("pval.name", "padj.name", "log2fc.name","id.name")[check_colnames],collapse=" ") )
  )
  
  res = resdat[,match(c(pval.name,padj.name,log2fc.name,id.name),colnames(resdat))]
  colnames(res) = c("pval","padj","log2fc","unique_id")
  if(is.null(xlim)) {xlim = range(res$log2fc)}
  if(is.null(ylim)) {ylim = range(-log10(res$pval))}
  
  res$color="None"
  res$color[which((abs(res$log2fc)>log2fc.cut)*(res$pval<pval.cut)==1)] = 
    paste0("pval","<",pval.cut," & abs(logfc)>",log2fc.cut)
  res$color[which((abs(res$log2fc)<log2fc.cut)*(res$pval<pval.cut)==1)] =  
    paste0("pval","<",pval.cut, " & abs(logfc)<",log2fc.cut)
  res$color[which((abs(res$log2fc)>log2fc.cut)*(res$padj<padj.cut)==1)] = 
    paste0("padj","<",padj.cut," & abs(logfc)>",log2fc.cut)
  res$color[which((abs(res$log2fc)<log2fc.cut)*(res$padj<padj.cut)==1)] = 
    paste0("padj","<",padj.cut, " & abs(logfc)<",log2fc.cut)
  # if pval.cut is high and only have genes < padj.cut, padj.cut dominates, is this the best way, or should it
  # be intersection?
  
  # levels of color will be a subset of all_levels, but we want the color to match the all_levels
  tmplevels = levels(as.factor(res$color))
  all_levels = c("None",
                 paste0("pval","<",pval.cut, " & abs(logfc)<",log2fc.cut), # only exists if pval != pval.name 
                 paste0("padj","<",padj.cut," & abs(logfc)<",log2fc.cut), 
                 paste0("pval","<",pval.cut, " & abs(logfc)>",log2fc.cut), # only exists if pval != pval.name 
                 paste0("padj","<",padj.cut, " & abs(logfc)>",log2fc.cut))
  res$color = factor(
    res$color,
    levels = intersect(all_levels,
                       tmplevels
    ))
  tmplevels = levels(res$color)
  
  # add selected genes
  shapedata = data.frame()
  if(!is.null(sel_genes)) {
    tmpind = sapply(unlist(sel_genes),function(k) match(k,res$unique_id))
    tmpind = unique(unlist(tmpind))
    shapedata <- res[tmpind,]
  }
  
  
  p <- ggplot(res,aes(x=log2fc,y=-log10(pval),color=color,text=unique_id))+
    geom_point(shape=19,fill="black")
  p <- p + scale_color_manual(values=
                                c("grey40","grey60","green3","grey70","red2")[match(tmplevels,all_levels)],
                              limits=levels(res$color),
                              name="Significance")
  
  if(nrow(shapedata)>0) {
    p <- p + geom_point(data=shapedata,fill="orange",shape=23,size=3,color="grey40") +
      geom_text_repel(data=shapedata,aes(label=unique_id),color="black",
                      size = 5,
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines"))+
      #scale_size_manual(values = c(1,3))+
      guides(size=FALSE,shape=FALSE,fill=FALSE)
  }
  
  p <- p + theme_base() + ggtitle(main)
  
  if(!interactive){p}else{
    
    p <- p + theme(plot.margin = unit(c(2,2,2,2), "cm"))
    
    
    g <- plotly_build(p)
    
    
    #Match order of text to proper gene order
    newtext =  paste("Gene ID:",res$unique_id,"<br>",
                     "Comparison",res$test,"<br>",
                     "log2fc",signif(res$log2fc,3),"<br>",
                     "pval",signif(res$pval,3),"<br>",
                     "padj",signif(res$padj,3))
    
    for(ii in 1:length(g$x$data)) {
      tmpid = do.call(rbind,strsplit(g$x$data[[ii]]$text,"<br>"))[,4]
      g$x$data[[ii]]$text <- newtext[match(tmpid,res$unique_id)]
    }
    
    g
  }
  
}


gene_dge_table <- function(resdat,
                           log2fc.cut=2,
                           padj.cut = 0.05,
                           pval.cut = 0.05,
                           pval.name="pval",
                           padj.name="padj",
                           log2fc.name="log2fc",
                           gene.name="geneid") {
  
  check_colnames = !(c(pval.name, padj.name, log2fc.name, gene.name)%in%colnames(resdat))
  if(sum(check_colnames)>0) stop(
    paste("The following argument values are not column names of resdat:\n",
          paste0(c("pval.name", "padj.name", "log2fc.name", "gene.name")[check_colnames],collapse=" ") )
  )
  
  outres = resdat[,match(c(gene.name,pval.name,padj.name,log2fc.name),colnames(resdat))]
  colnames(outres) = c(gene.name,"pval","padj","log2fc")
  
  
  outres = outres%>%mutate(direction=ifelse(log2fc<0,"down","up"),
                           significance_padj = ifelse(padj<padj.cut,paste0("padj<",padj.cut),paste0("padj>",padj.cut)),
                           significance_p = ifelse(pval<pval.cut,paste0("pval<",pval.cut),paste0("pval>",pval.cut)),
                           log2fccut = ifelse(abs(log2fc)>log2fc.cut,paste0("abs(log2fc)>",log2fc.cut),paste0("abs(log2fc)<",log2fc.cut)))
  
  tmp1 = outres%>%group_by(significance_padj,significance_p,log2fccut,direction)%>%summarize("number"=n())
  
  tmp2 = with(outres,table(log2fccut,direction,significance_p,significance_padj))
  list(tmp1,tmp2,outres)
  
}


# from julja's normalization code

density_plots <- function(mat,col.labels=NULL,colvec=NULL,ltyvec=NULL,
                          xlab="",ylab="",title="") {
  if(is.null(col.labels)) col.labels = colnames(mat)
  if(is.null(col.labels)) col.labels = 1:ncol(mat)
  if(is.null(colvec)) colvec = rep(1,ncol(mat))
  if(is.null(ltyvec)) ltyvec = rep(1,ncol(mat))
  dl = NULL; yl = c(0,0); xl = range( mat )
  for( i in 1:ncol( mat ) ){
    dl = c(dl,list(density(mat[,i],from=xl[1],to=xl[2])))
    yl = range(yl,dl[[i]]$y)
  }
  plot(x=0,y=0,type="n",xlim=xl,ylim=yl,xlab=xlab,ylab=ylab,main=title)
  for(i in 1:ncol( mat ) ){
    lines(dl[[i]]$x,dl[[i]]$y,
          col=colvec[i],
          lty=ltyvec[i])
  }
  legend(x="topright", legend=col.labels, 
         col=colvec, lty=ltyvec, lwd=2, cex=.8)
}

scatter_compare <- function(xmat,ymat,colvec,col.labels,
                            xlab="",ylab="",main="",...) {
  xmat = as.matrix(xmat)
  ymat = as.matrix(ymat)
  plot(xmat,ymat,
       pch='.',
       xlab=xlab,ylab=ylab,main=main,
       col=matrix(colvec, ncol=length(colvec), nrow=nrow(xmat), byrow=T),...)
  legend(x="bottomright", legend=col.labels, fill = colvec,cex=.8)
}

# work in progress
gene_dotplots <- function(dat,
                          selectgenes,
                          title="") {
  #dat in in long format with columns genename, y, SampleID, group
  tmpdat = dat%>%filter(genename%in%selectgenes)
  ggplot(tmpdat,aes(SampleID,y,color=group))+
    geom_boxplot(position="dodge")+
    facet_wrap(~genename)+
    xlab("")+theme_classic()+ggtitle(title)
}
