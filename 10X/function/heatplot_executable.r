# https://github.com/YuLab-SMU/enrichplot/blob/master/R/heatplot.R

##' @rdname heatplot
##' @exportMethod heatplot
setMethod("heatplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, foldChange = NULL) {
            heatplot.enrichResult(x, showCategory, foldChange)
          })

##' @rdname heatplot
##' @exportMethod heatplot
setMethod("heatplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, foldChange = NULL) {
            heatplot.enrichResult(x, showCategory, foldChange)
          })



##' @rdname heatplot
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @author Guangchuang Yu

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)
  
  if(x@readable) {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}

list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}


update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}



extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- DOSE::geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

#axis.text.x = element_text(size = 15 ... : x axis label size
#axis.text.y = element_text(size = 15 ... : y axis label size

heatplot.enrichResult <- function(x, showCategory=30, foldChange=NULL) {
  n <- update_n(x, showCategory)
  geneSets <- extract_geneSets(x, n)
  
  foldChange <- fc_readable(x, foldChange)
  d <- list2df(geneSets)
  
  if (!is.null(foldChange)) {
    d$foldChange <- foldChange[as.character(d[,2])]
    ## palette <- fc_palette(d$foldChange)
    p <- ggplot(d, aes_(~Gene, ~categoryID)) +
      geom_tile(aes_(fill = ~foldChange), color = "white") +
      scale_fill_continuous(low="blue", high="red", name = "fold change")
    ## scale_fill_gradientn(name = "fold change", colors = palette)
    
  } else {
    p <- ggplot(d, aes_(~Gene, ~categoryID)) + geom_tile(color = 'white')
  }
  p + xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1), 
          #axis.text.x = element_text(size = 15, angle = 60, hjust = 1), 
          axis.text.y = element_text(size = 15, angle = 0))
}

### excute 
#edox1 <- setReadable(ekk, 'org.Hs.eg.db', 'ENTREZID') 
heatplot.enrichResult(edox1, showCategory=30, foldChange = geneList)
