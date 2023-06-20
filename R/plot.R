maxColors = 5

.repClusterLinetype<-function(clusters,rep=F,jump=1){
  n = length(unlist(clusters))
  resVec = integer(n)
  if(rep){
    for (i in 1:length(clusters)) {
      resVec[clusters[[i]]]<-rep((i+jump),length(clusters[[i]]))
    }
  }else{
    for (i in 1:length(clusters)) {
      resVec[clusters[[i]]]<-c((1+jump):(jump+length(clusters[[i]])))
    }
  }
  return(resVec)
}

.plotPmfs <- function(data, clusterResult, centroid, ylim, title, lineSize, int_hue_pal,colors){
  # create the dataframe for plotting.
  epmfs_df = NULL
  pmf_levels = NULL
  for (i in 1:length(data)) {
    epmf <- epmf(data[[i]])
    epmf$pmfs <- rep(paste0("pmf-",toString(i)), nrow(epmf))
    pmf_levels = c(pmf_levels, paste0("pmf-",toString(i)))
    epmf$groups <- rep(i, nrow(epmf))
    epmf$clusters <- rep(-1, nrow(epmf))
    epmfs_df = rbind(epmfs_df, epmf)
  }
  sort(unique(epmfs_df$category))->cateNames
  if(is.null(ylim)){
    ylim = c(0, max(epmfs_df$epmf))
  }

  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- ggplot2::position_dodge(.7) # move them .05 to the left and right
  epmfs_df$pmfs = factor(epmfs_df$pmfs,levels = pmf_levels)
  epmfs_df$category = factor(epmfs_df$category)

  if(is.null(clusterResult)){
    ggplot2::ggplot(epmfs_df, ggplot2::aes(color=pmfs, y=epmf, x=category, ymin = 0,
                                           ymax=epmf, shape=pmfs,
                         linetype=pmfs)) +
      ggplot2::scale_shape_manual(values=1:nlevels(epmfs_df$pmfs))+
      ggplot2::geom_linerange(position = pd, size=lineSize, na.rm=TRUE) +
      ggplot2::geom_point(position=pd, size=1.6, stroke = 1.0, na.rm=TRUE)+
      ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
      ggplot2::ylim(ylim) +
      ggplot2::labs(title=title, y = "", x="")
  }else{
    centroids_df=NULL
    legend_pmfs_ids = rep(0,length(data))
    for (j in 1:length(clusterResult$clusters)) {
      ct <- data.frame(category = cateNames[(clusterResult$centroids[[j]][[1]]+1)],
                       epmf=clusterResult$centroids[[j]][[2]])
      ct$clusters <- rep(j,nrow(ct))
      epmfs_df$clusters[epmfs_df$groups %in% clusterResult$clusters[[j]]]<-j
      centroids_df = rbind(centroids_df,ct)
      legend_pmfs_ids[clusterResult$clusters[[j]]]<-j
    }
    epmfs_df$clusters = factor(epmfs_df$clusters)
    #centroids_df$category=factor(centroids_df$category)
    centroids_df$clusters=factor(centroids_df$clusters)

    legend_pmfs_colors=NULL
    if(is.null(int_hue_pal)){
      if(is.null(colors)){
        if(length(clusterResult$clusters)>maxColors){
          color_values = scales::hue_pal()(length(clusterResult$clusters))
        }else{
          color_values = scales::hue_pal()(maxColors)[c(1:length(clusterResult$clusters))]
        }
      }else{
        color_values = colors
      }
    }else{
      color_values = scales::hue_pal()(int_hue_pal)
    }

    for (e in legend_pmfs_ids) {
      legend_pmfs_colors = c(legend_pmfs_colors, color_values[e])
    }
    if(centroid){
      plt = ggplot2::ggplot(data=epmfs_df) +
        ggplot2::scale_color_manual(breaks = c(1:(length(clusterResult$clusters))), values=color_values)+
        ggplot2::scale_shape_manual(values=1:nlevels(epmfs_df$pmfs))+
        ggplot2::geom_linerange(data=epmfs_df, mapping=ggplot2::aes(color=clusters, y=epmf, x=category, ymin = 0,
                                                                    ymax=epmf, linetype=pmfs),position = pd, size=lineSize, na.rm=TRUE) +
        ggplot2::geom_point(data=epmfs_df, mapping=ggplot2::aes(color=clusters, y=epmf, x=category, shape=pmfs),
                            position=pd, size=1.6, stroke = 1.0, na.rm=TRUE)+
        ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(color = legend_pmfs_colors)))+
        ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
        ggplot2::geom_point(data=centroids_df,mapping=ggplot2::aes(color=clusters, y=epmf, x=category),
                            size=2.5, na.rm=TRUE)+
        ggplot2::geom_line(data=centroids_df,mapping=ggplot2::aes(color=clusters, y=epmf, x=category,
                                                                  group = clusters), size=0.3, na.rm=TRUE)+
        ggplot2::guides(color=ggplot2::guide_legend(title="cluster mean: g"))+
        #ggplot2::guides(color = "none")+
        ggplot2::ylim(ylim)+
        ggplot2::labs(title=title, y = "", x="")
    }else{
      plt = ggplot2::ggplot(data=epmfs_df) +
        ggplot2::scale_color_manual(breaks = c(1:(length(clusterResult$clusters))), values=color_values)+
        ggplot2::scale_shape_manual(values=1:nlevels(epmfs_df$pmfs))+
        ggplot2::geom_linerange(data=epmfs_df, mapping=ggplot2::aes(color=clusters, y=epmf, x=category, ymin = 0,
                                                                    ymax=epmf, linetype=pmfs),position = pd, size=lineSize, na.rm=TRUE) +
        ggplot2::geom_point(data=epmfs_df, mapping=ggplot2::aes(color=clusters, y=epmf, x=category, shape=pmfs),
                            position=pd, size=1.6, stroke = 1.0, na.rm=TRUE)+
        ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(color = legend_pmfs_colors)))+
        ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
        ggplot2::ylim(ylim)+
        ggplot2::labs(title=title, y = "", x="")
    }
    return(plt)
  }
}

.plotPdfs <- function(data, clusterResult, centroid, centroid.colors, rep, jump, xlim, ylim, title,
                      lineSize, lineTypes,centroidLineSize, int_hue_pal,colors){
  num_pdfs = length(data)
  min_data=Inf
  max_data=-Inf
  for (e in data) {
    min(e)->l
    max(e)->r
    if(min_data>l){
      min_data=l
    }
    if(max_data<r){
      max_data=r
    }
  }
  epdfs_df = NULL
  pdf_levels = NULL
  for (i in 1:length(data)) {
    epdf = as.data.frame(epdf(c(data[[i]],min_data,max_data),scale = 0))
    epdf$pdfs = rep(paste0("pdf-",toString(i)), nrow(epdf))
    pdf_levels = c(pdf_levels, paste0("pdf-",toString(i)))
    epdf$groups <- rep(i, nrow(epdf))
    epdf$clusters <- rep(-1, nrow(epdf))
    epdfs_df = rbind(epdfs_df, epdf)
  }
  if(is.null(xlim)){
    xlim = c(min(epdfs_df$positions), max(epdfs_df$positions))
  }
  if(is.null(ylim)){
    ylim = c(0, max(epdfs_df$probabilities))
  }

  epdfs_df$pdfs = factor(epdfs_df$pdfs, levels = pdf_levels)

  if(is.null(clusterResult)){
    ggplot2::ggplot(data=epdfs_df, ggplot2::aes(y=probabilities, x=positions,
                                                color=pdfs, linetype=pdfs)) +
      ggplot2::geom_line(size=lineSize, na.rm=TRUE)+
      ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
      ggplot2::xlim(xlim)+
      ggplot2::ylim(ylim)+
      ggplot2::labs(title=title, y = "", x="")

  }else{
    centroids_df=NULL
    legend_pdfs_ids = rep(0,length(data))
    for (j in 1:length(clusterResult$clusters)) {
      cluster_data_list=list()
      for (e in clusterResult$clusters[[j]]) {
        cluster_data_list = append(cluster_data_list,list(c(data[[e]],min_data,max_data)))
      }
      ct <- as.data.frame(epdf.cluster(cluster_data_list,scale = 0))
      ct$clusters <- rep(j,nrow(ct))
      epdfs_df$clusters[epdfs_df$groups %in% clusterResult$clusters[[j]]]<-j
      centroids_df = rbind(centroids_df,ct)
      legend_pdfs_ids[clusterResult$clusters[[j]]]<-j
    }
    epdfs_df$clusters = factor(epdfs_df$clusters)
    centroids_df$clusters=factor(centroids_df$clusters)

    if(is.null(int_hue_pal)){ # if we don't use int_hue_pal
      if(is.null(colors)){
        if(length(clusterResult$clusters)>maxColors){
          color_values = scales::hue_pal()(length(clusterResult$clusters))
        }else{
          color_values = scales::hue_pal()(maxColors)[c(1:length(clusterResult$clusters))]
        }
      }else{
        color_values = colors
      }
    }else{ # if we use int_hue_pal
      color_values = scales::hue_pal()(int_hue_pal)
    }

    if(is.null(centroid.colors)){
      clusterColors = color_values
      centroids_df$col = centroids_df$clusters
    }else{
      clusterColors = centroid.colors
      centroids_df$col=0
      for (i in 1:length(clusterResult$clusters)) {
        centroids_df$col[centroids_df$clusters==i]<-(i+length(clusterResult$clusters))
      }
      centroids_df$col<-as.factor(centroids_df$col)
    }

    legend_pdfs_colors=NULL
    for (e in legend_pdfs_ids) {
      legend_pdfs_colors = c(legend_pdfs_colors, color_values[e])
    }
    if(is.null(lineTypes)){
      lineTypes = .repClusterLinetype(clusterResult$clusters,rep,jump)
    }
    if(centroid){
      plt = ggplot2::ggplot(data=epdfs_df) +
        ggplot2::scale_color_manual(breaks = c(1:(2*length(clusterResult$clusters))), values=c(color_values,centroid.colors))+
        ggplot2::geom_line(data=epdfs_df, ggplot2::aes(y=probabilities, x=positions,
                                                       color=clusters, linetype=pdfs), size=lineSize, na.rm=TRUE)+
        ggplot2::scale_linetype_manual(values=lineTypes)+
        ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(color = legend_pdfs_colors,size=0.8)))+
        ggplot2::geom_line(data=centroids_df, mapping=ggplot2::aes(y=probabilities, x=positions,
                                                                   color=col), size=centroidLineSize, na.rm=TRUE) +
        ggplot2::guides(color=ggplot2::guide_legend(override.aes = list(title="cluster mean: g",
                                                                        color = c(clusterColors,rep(0,length(clusterColors))))))+
        # ggplot2::guides(color = "none")+
        ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
        ggplot2::xlim(xlim)+
        ggplot2::ylim(ylim)+
        ggplot2::labs(title=title, y = "", x="")
    }else{
      plt = ggplot2::ggplot(data=epdfs_df) +
        ggplot2::scale_color_manual(breaks = c(1:(length(clusterResult$clusters))), values=color_values)+
        ggplot2::geom_line(data=epdfs_df, ggplot2::aes(y=probabilities, x=positions,
                                                       color=clusters, linetype=pdfs), size=lineSize, na.rm=TRUE)+
        ggplot2::scale_linetype_manual(values=lineTypes)+
        ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(color = legend_pdfs_colors,size=0.8)))+
        ggplot2::theme(legend.key = ggplot2::element_rect(fill = "transparent"))+
        ggplot2::xlim(xlim)+
        ggplot2::ylim(ylim)+
        ggplot2::labs(title=title, y = "", x="")
    }
    return(plt)
  }
}

#' Plot a set of distributions
#'
#' Plot distributions that are estimated from the datasets stored in \code{data}.
#'
#' @param data A list of vectors. Each vector stores a dataset that represents a distribution.
#' @param distribution_type A char variable to indicate the types of distributions, "pdf": continuous or "pmf": discrete.
#' @param clusterResult A list stores the best clustering result which includes the clusters, ams and centroids (check ?KMeansCluster).
#' @param centroid A bool variable controlling plot the centroid distribution or not.
#' @param xlim A vector of two values represent the limit of horizontal axis.
#' @param ylim A vector of two values represent the limit of vertical axis.
#' @param lineSize A real value to set the size of line.
#'
#' @details Plot distributions that are estimated from the datasets stored in \code{data}.
#'
#' @examples
#'
#'# Discrete distributions
#' set.seed(123)
#' data1 = rbinom(n=50,size=20,prob=0.3)
#' data2 = rbinom(n=30,size=10,prob=0.5)
#' data3 = rpois(35,4)
#' data4 = rpois(45,9)
#' data5 = rpois(15,8)
#' ds = list(data1,data2,data3,data4,data5)
#' plot.pmdfs(data=ds,distribution_type="pmf")
#'
#' @export
#'
plot.pmdfs <- function(data, distribution_type, clusterResult=NULL,centroid=TRUE,
                       centroid.colors=NULL, rep=T, jump=1, xlim=NULL,ylim=NULL, title="",
                       lineSize=0.8, lineTypes=NULL, centroidLineSize=1.5, int_hue_pal=NULL, colors=NULL,...){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if(!(distribution_type=="pdf" | distribution_type=="pmf")){
    return("Please specify type to pmf or pdf.")
  }
  if(is.data.frame(data) | is.matrix(data) | is.array(data)){
    pmdfs <- .toList(data)
  }else if(!is.list(data)){
    return("Acceptable data types include vector, list, data.frame, matrix or array.")
  }else{
    pmdfs <- data
  }
  if(distribution_type=="pdf" & (!is.numeric(data[[1]]))){
    return("Data type should be numeric.")
  }
  if(distribution_type=="pmf"){
    .plotPmfs(data=pmdfs, clusterResult=clusterResult, centroid, ylim, title, lineSize, int_hue_pal,colors)
  }else{
    .plotPdfs(data=pmdfs, clusterResult=clusterResult, centroid,centroid.colors, rep, jump,
              xlim, ylim, title, lineSize, lineTypes, centroidLineSize, int_hue_pal,colors)
  }
}

