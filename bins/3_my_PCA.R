plotPCA.sny <- function (object, intgroup = "condition", ntop = 500, returnData=TRUE){
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d12 <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  d13 <- data.frame(PC1=pca$x[,1], PC3=pca$x[,3], group=group, intgroup.df, name=colnames(object))
  d23 <- data.frame(PC2=pca$x[,2], PC3=pca$x[,3], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d12, "percentVar") <- percentVar[1:2]
    attr(d13, "percentVar") <- percentVar[c(1,3)]
    attr(d23, "percentVar") <- percentVar[2:3]
    d <- list(d12=d12, d13=d13, d23=d23)
    return(d)
  }
}
