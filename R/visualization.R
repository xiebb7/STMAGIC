#----cornplot basic----
#' Basic Cornplot
#'
#' Plot information on Cornplot, five stage is support now (E6.5, E6.75, E6.75_rep, E7.0, E7.25, E7.5, E7.5_suo).
#'
#' @param coord Coordinate object from spatial object.
#' @param data A vector of plot data on each spot.
#' @param title The title plot on the head of cornplot.
#' @param stage The support stage of cornplot.
#' @param cols The color of cornplot.
#' @param shapes The shape of spot (1: "circle" or 2: "gaussian").
#' @param limits NULL to use the default scale range, or a numeric vector of length two providing limits of the scale.
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @import ggplot2
#'
#'
#'
corn_plot = function(coord,
                     data,
                     title = NULL,
                     stage,
                     cols = c("azure2", "green",'yellow1','yellow','#FF3333',"#FF0000"),
                     shapes = 1,
                     limits = NULL){

  if(!shapes %in% 1:2){
    stop('Please input the right shape, 1 or 2!')
  }

  if(is.null(title)){
    title = stage
  }

  #---plot data---
  data = data.frame(exp = data, location = names(data))
  plot_data = merge(coord@coordinate, data, by = "location",all = T)
  plot_data$shape = 'circle'
  plot_data$shape[grep('M', plot_data$location)] = 'diamond'
  #---plot data---


  #---axis label---
  x_axis_label = data.frame(text_x = coord@coordinate[which(coord@coordinate$V2 == coord@text_height - 1), 1],
                            text_y = rep(coord@text_height, length(coord@text)),
                            text = coord@text)
  y_axis_label = data.frame(text_x = rep(coord@tick_start_of_layers, coord@text_height - 1),
                            text_y = 1:(coord@text_height -1),
                            text = 1:(coord@text_height -1))
  #---axis label---

  #---tick---
  tick_x = c(rep(coord@coordinate[which(coord@coordinate$V2 == coord@text_height - 1), 1], each = 2),
             rep(c(coord@tick_start_of_layers + coord@shift, coord@tick_start_of_layers + coord@shift + 0.15), coord@text_height-1))
  tick_y = c(rep(c(coord@text_height - 0.4,coord@text_height - 0.3), length(coord@text)),
             rep(1:(coord@text_height -1), each = 2))
  group = rep(rep(1:(length(tick_x)/2), each = 2))
  tick_df = data.frame(tick_x, tick_y, group = group)
  #---tick---


  #---border---
  broder_x = c(min(coord@coordinate[, 1]) - 0.8, min(coord@coordinate[, 1]) - 0.8, max(coord@coordinate[, 1]) + 0.8, max(coord@coordinate[, 1]) + 0.8, min(coord@coordinate[, 1]) - 0.8)
  broder_y = c(min(coord@coordinate[, 2]) - 0.6, max(coord@coordinate[, 2]) + 0.6, max(coord@coordinate[, 2]) + 0.6, min(coord@coordinate[, 2]) - 0.6, min(coord@coordinate[, 2]) - 0.6)
  broder = data.frame(broder_x, broder_y)
  #---border---

  #---plot---
  p = ggplot(plot_data, aes(x1, V2))

  if(shapes == 1){
    p = p + geom_point(color = "transparent", aes(fill = exp, shape = as.factor(shape), size = as.factor(shape)))
  }
  if(shapes == 2){
    df_polygon = data.frame(polygon_x = c(f_x(plot_data$x1)),
                            polygon_y = c(f_y(plot_data$x1, plot_data$V2, plot_data$exp)),
                            type = rep(plot_data$location, each = 100),
                            exp = rep(plot_data$exp, each = 100))
    p = p + geom_polygon(data = df_polygon, aes(x = polygon_x, y = polygon_y, group = type, fill = exp))
  }

  p = p +
    scale_shape_manual(values=c(21, 23)) +
    scale_size_manual(values=c(coord@size, coord@size - 2)) +
    guides(shape='none', size='none') +
    theme_void() +
    xlim(coord@x) +
    ylim(coord@y) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = coord@title_height),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.position = coord@legend_position) +
    labs(fill = NULL) +
    geom_text(data = x_axis_label, aes(text_x, text_y, label = text), size = coord@text_size + 1, color = 'black') +
    geom_text(data = y_axis_label, aes(text_x, text_y, label = text), size = coord@text_size + 2, color = 'black') +
    geom_path(data = broder, aes(broder_x, broder_y)) +
    geom_line(data = tick_df, aes(tick_x, tick_y, group = group)) +
    scale_fill_gradientn(colors = colorRampPalette(cols)(100))
  #---plot---


  if(all(plot_data$exp == 0, na.rm = T)){
    p = p + scale_fill_gradientn(colors = 'azure2')
  }

  if(!is.null(limits) & length(limits) == 2){
    p = p + scale_fill_gradientn(colors = colorRampPalette(cols)(100), limits = limits)
  }

  return(p)

}

#----cornplot basic----

#----CornPlot----
#' CornPlot
#'
#' CornPlot of multi-omics on cornplot.
#'
#' @param object A spatial object.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a cornplot.
#'
#' @export
#'
#' @rdname CornPlot
#' @export CornPlot
#'
CornPlot = function(object, ...){
  UseMethod(generic = 'CornPlot', object = object)
}

#' @param key A gene or peak in our dataset.
#' @param omics The omics to plot: 'snRNASpatial', 'snATACSpatial', 'Celltype_spatial', 'GEOseq'.
#' @param slot_name The slot name of data.
#' @param shapes The shape of spot (1: "circle" or 2: "gaussian").
#' @param return_data Return cornplot data, default False.
#' @param smooth Smooth cornplot data, default False.
#' @param sigmas parameter in gaussian kernel, default:0.5.
#' @param title The title plot on the head of cornplot.
#' @param grid_smooth Set smooth method: 'all', 'col', 'row', default 'col'.
#'
#' @rdname CornPlot
#' @export
#' @method CornPlot spatial
#'
#' @importFrom methods slot
#'
CornPlot.spatial = function(object, key, omics = 'snRNASpatial', slot_name = 'exp_mat', shapes = 'circle', return_data = F, smooth = F, sigmas = 0.5, title = NULL, grid_smooth = 'col'){

  if(is.null(title)){
    title = key
  }

  data = slot(slot(object, name = omics), name = slot_name)

  if(!key %in% rownames(data)){
    stop('Please input the gene or peak in our dataset!')
  }

  key_data = unlist(data[key, ])
  names(key_data) = sapply(colnames(data), function(x) unlist(strsplit(x, '_'))[2])

  if(smooth){
    key_data = SmoothData(object, key,
                          omics = omics,
                          slot_name = slot_name,
                          sigmas = sigmas,
                          grid_smooth = grid_smooth)
  }

  p = corn_plot(object@Coordinate, key_data, title = title, object@stage, shapes = shapes)

  if(return_data){
    result = list()
    result$plot = p
    result$data = key_data
    return(result)
  }else{
    return(p)
  }

}
#----CornPlot----

#----CornPlot keysets----
#' CornPlotKeysets
#'
#' Keysets of CornPlot of multi-omics on cornplot.
#'
#' @param object A spatial object.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a cornplot.
#'
#' @export
#'
#' @rdname CornPlotKeysets
#' @export CornPlotKeysets
#'
#'
#'
CornPlotKeysets = function(object, ...){
  UseMethod(generic = 'CornPlotKeysets', object = object)
}

#' @param key A gene or peak in our dataset.
#' @param omics The omics to plot: 'snRNASpatial', 'snATACSpatial', 'Celltype_spatial', 'GEOseq'.
#' @param slot_name The slot name of data.
#' @param shapes The shape of spot (1: "circle" or 2: "gaussian").
#' @param return_data Return cornplot data, default False.
#' @param smooth Smooth cornplot data, default False.
#' @param sigmas parameter in gaussian kernel, default:0.5.
#' @param title The title plot on the head of cornplot.
#' @param grid_smooth Set smooth method: 'all', 'col', 'row', default 'col'.
#'
#' @rdname CornPlotKeysets
#' @export
#' @method CornPlotKeysets spatial
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom methods slot
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#'
#'
CornPlotKeysets.spatial = function(object, key, omics = 'snRNASpatial', slot_name = 'exp_mat', shapes = 'circle', return_data = F, smooth = F, sigmas = 0.5, title = NULL, grid_smooth = 'col'){

  if(omics == 'snATACSpatial'){

    if(class(key) == 'data.frame'){
      colnames(key)[1:3] = c('chr','start','end')
      key_grange = makeGRangesFromDataFrame(key, keep.extra.columns = TRUE)
      peak_name = rownames(slot(slot(object, name = omics), name = slot_name))
      peak_df = t(sapply(peak_name, function(x) unlist(strsplit(x, '-'))))
      colnames(peak_df)[1:3] = c('chr','start','end')
      peak_grange = makeGRangesFromDataFrame(as.data.frame(peak_df), keep.extra.columns = TRUE)
      use_key = peak_name[unique(findOverlaps(peak_grange, key_grange)@from)]
    }

    if(class(key) == 'character'){
      keys = rownames(slot(slot(object, name = omics), name = slot_name))
      use_key = intersect(key, keys)
    }

  } else if(omics == 'snRNASpatial' | omics == 'GEOseq'){
    keys = rownames(slot(slot(object, name = omics), name = slot_name))
    use_key = intersect(key, keys)
  }

  if(length(use_key) == 0){
    stop('Please input the geneset or peakset in our dataset!')
  }

  if(class(key) == 'data.frame'){
    print(paste0('The length of input keysets is ', nrow(key), ', ', length(use_key), ' are in our dataset.'))
  }

  if(class(key) == 'character'){
    print(paste0('The length of input keysets is ', length(key), ', ', length(use_key), ' are in our dataset.'))
  }

  cells_rankings = AUCell_buildRankings(as.matrix(slot(slot(object, name = omics), name = slot_name)), plotStats=FALSE, splitByBlocks=TRUE)
  cells_AUC = AUCell_calcAUC(list(signal = use_key), cells_rankings, aucMaxRank = floor(nrow(cells_rankings)), verbose = F)
  key_data = c(t(getAUC(cells_AUC)))
  names(key_data) = sapply(colnames(slot(slot(object, name = omics), name = slot_name)), function(x){unlist(strsplit(x,'_'))[2]}, USE.NAMES = F)

  if(smooth){
    key_data = SmoothData(object, use_key,
                          omics = omics,
                          slot_name = slot_name,
                          sigmas = sigmas,
                          grid_smooth = grid_smooth)
  }

  p = corn_plot(object@Coordinate, key_data, title = title, object@stage, shapes = shapes)

  if(return_data){
    result = list()
    result$plot = p
    result$data = key_data
    return(result)
  }else{
    return(p)
  }

}
#----CornPlot keysets--------

#----CornPlot MultiStage----
#' Plot Multi stage CornPlot
#'
#' Plot Multi stage CornPlot.
#'
#' @param object A list of spatial object.
#' @param key A gene or peak in our dataset.
#' @param omics The omics to plot: 'snRNASpatial', 'snATACSpatial', 'Celltype_spatial', 'GEOseq'.
#' @param slot_name The slot name of data.
#' @param title The title of cornplot, the default is gene sysmbol.
#' @param shapes The shape of spot (1: "circle" or 2: "gaussian").
#' @param scale Unify expression into the same range, default False
#' @param return_data Return cornplot data, default False.
#' @param smooth Smooth cornplot data, default False.
#' @param sigmas parameter in gaussian kernel, default:0.5.
#' @param grid_smooth Set smooth method: 'all', 'col', 'row', default 'col'.
#'
#' @rdname CornPlotMultiStage
#' @export
#'
#' @importFrom methods slot
#'
#'
CornPlotMultiStage = function(object, key, omics = 'snRNASpatial', title = NULL, shapes = 1, scale = F, slot_name = 'exp_mat', return_data = F, smooth = F, sigmas = 0.5, grid_smooth = 'col'){

  if(length(object) == 1 & class(object) == 'spatial'){
    p = CornPlot(object, key, omics = omics, shapes = shapes, smooth = smooth, sigmas = sigmas, grid_smooth = grid_smooth, return_data = return_data)
  }

  if(length(object) > 1 & all(sapply(object, class) == 'spatial')){

    keys = unique(unlist(lapply(object, function(x) rownames(slot(slot(x, name = omics), name = slot_name)))))

    if(!key %in% keys){
      stop('Please input the gene or peak in our dataset!')
    }

    key_data = lapply(object, function(x){
      if(key %in% rownames(slot(slot(x, name = omics), name = slot_name))){
        exp = unlist(slot(slot(x, name = omics), name = slot_name)[key, ])
      }else(
        exp = rep(0, ncol(slot(slot(x, name = omics), name = slot_name)))
      )
      names(exp) = sapply(colnames(slot(slot(x, name = omics), name = slot_name)), function(x) unlist(strsplit(x, '_'))[2])
      return(exp)
    })

    if(smooth){

      key_data = lapply(object, function(x){
        SmoothData(x, key,
                   omics = omics,
                   slot_name = slot_name,
                   sigmas = sigmas,
                   grid_smooth = grid_smooth)
      })

    }

    plot_list = list()

    if(scale == T){

      limits = c(min(unlist(key_data), na.rm = T), max(unlist(key_data), na.rm = T))

      for(i in 1:length(object)){

        if(!is.null(title) & length(title) == length(object)){
          plot_list[[i]] = corn_plot(object[[i]]@Coordinate, key_data[[i]], title = title[i], object[[i]]@stage, shapes = shapes, limits = limits)
        }else{
          plot_list[[i]] = corn_plot(object[[i]]@Coordinate, key_data[[i]], object[[i]]@stage, shapes = shapes, limits = limits)
        }

      }

    }else{

      for(i in 1:length(object)){

        if(!is.null(title) & length(title) == length(object)){
          plot_list[[i]] = corn_plot(object[[i]]@Coordinate, key_data[[i]], title = title[i], object[[i]]@stage, shapes = shapes)
        }else{
          plot_list[[i]] = corn_plot(object[[i]]@Coordinate, key_data[[i]], object[[i]]@stage, shapes = shapes)
        }

      }

    }

    if(return_data){
      result = list()
      result$plot = plot_list
      result$data = key_data
      return(result)
    }else{
      return(plot_list)
    }
  }

}
#----CornPlot MultiStage--------

#----CornPlot MultiStage keysets----
#' Plot Multi stage CornPlot
#'
#' Plot Multi stage CornPlot.
#'
#' @param object A list of spatial object.
#' @param key A vector of geneset or peakset.
#' @param omics The omics to plot: snRNASpatial, snATACSpatial, Celltype_spatial, GEOseq.
#' @param slot_name The slot name of data.
#' @param title The title of cornplot, the default is gene sysmbol.
#' @param smooth Smooth cornplot data, default False.
#' @param sigmas parameter in gaussian kernel, default:0.5.
#' @param grid_smooth Set smooth method: 'all', 'col', 'row', default 'col'.
#' @param shapes The shape of spot (1: "circle" or 2: "gaussian").
#' @param return_data Return cornplot data, default False.
#'
#' @rdname CornPlotMultiStageKeySets
#' @export
#'
CornPlotMultiStageKeySets = function(object, key, omics = 'snRNASpatial', title = NULL, shapes = 1, slot_name = 'exp_mat', return_data = F, smooth = F, sigmas = 0.5, grid_smooth = 'col'){

  if(length(object) == 1 & class(object) == 'spatial'){
    p = CornPlotKeysets(object, key, omics = omics, shapes = shapes, smooth = smooth, sigmas = sigmas, grid_smooth = grid_smooth, return_data = return_data, title = title)
  }

  if(length(object) > 1 & all(sapply(object, class) == 'spatial')){
    p = list()
    for(i in 1:length(object)){
      p[[i]] = CornPlotKeysets(object[[i]], key, omics = omics, shapes = shapes, smooth = smooth, sigmas = sigmas, grid_smooth = grid_smooth, return_data = return_data, title = title[i])
    }
  }

  return(p)

}
#----CornPlot MultiStage keysets--------

