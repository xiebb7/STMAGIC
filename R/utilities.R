f_x = function(x){
  sapply(x, function(xx){
    points = seq(xx - 0.5, xx + 0.5, length.out = 100)
  })
}

f_y = function(x, y, z){
  z_scale = scales::rescale(z, to = c(0, 0.8))
  result = c()
  for(i in 1:length(x)){
    points = seq(x[i] - 0.5, x[i] + 0.5, length.out = 100)
    result = c(result, z_scale[i] - z_scale[i]*scales::rescale(points, to = c(-1, 1))^2 + y[i] - 0.3)
  }
  return(result)
}

gaussian_kernel = function(distance, sigmas = 0.51){
  return(exp(-(distance^2) / (2 * sigmas^2)))
}

#----Spatial diff feature----
#' SpatialFeature
#'
#' Spatial Feature of multi-omics on cornplot.
#'
#' @param object A spatial object.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a feature list.
#'
#' @export
#'
#' @rdname SpatialFeature
#' @export SpatialFeature
#'
#'
#'
SpatialFeature = function(object, ...){
  UseMethod(generic = 'SpatialFeature', object = object)
}

#' @param key A vector of spot name.
#' @param omics The omics to plot: snRNA_spatial, snATAC_spatial, Celltype_spatial, GEOseq.
#' @param slot_name The slot name of data.
#' @param top_feature Number of top ranked features returned in the result.
#'
#' @rdname SpatialFeature
#' @export
#' @method SpatialFeature spatial
#'
#' @importFrom methods slot
#' @importFrom proxyC simil
#' @importFrom Matrix rowSums
#'
#'
SpatialFeature.spatial = function(object, key, omics = 'snRNA_spatial', slot_name = 'exp_mat', top_feature = 100){

  omics_data = slot(slot(object, name = omics), name = slot_name)
  key_data = matrix(0, nrow = 1, ncol = ncol(omics_data))
  colnames(key_data) = colnames(omics_data)
  key = paste(unlist(strsplit(object@stage, '_'))[1], key, sep = '_')
  key_data[1, match(key, colnames(key_data))] = 1

  #----code from COSG https://github.com/genecell/COSGR/blob/main/R/cosg.R----

  cosine_sim=simil(omics_data, key_data, method = "cosine",drop0=TRUE)
  pos_nonzero = cosine_sim != 0
  pos_nonzero = which(as.matrix(pos_nonzero),arr.ind = TRUE)
  genexlambda = cosine_sim*cosine_sim
  e_power2_sum = rowSums(genexlambda)
  genexlambda[pos_nonzero] = genexlambda[pos_nonzero]/(replicate(ncol(genexlambda),e_power2_sum)[as.matrix(pos_nonzero)])
  genexlambda = genexlambda*cosine_sim
  result = genexlambda[order(genexlambda[,1], decreasing = T),][1:top_feature]
  return(result)

  #----code from COSG https://github.com/genecell/COSGR/blob/main/R/cosg.R----

}

#----Spatial diff feature----


#----smooth cornplot----
#' SmoothData
#'
#' Smooth cornplot data.
#'
#' @param object A spatial object.
#' @param ... Arguments passed to other methods.
#'
#' @return Returns a cornplot.
#'
#' @export
#'
#' @rdname SmoothData
#' @export SmoothData
#'
SmoothData = function(object, ...){
  UseMethod(generic = 'SmoothData', object = object)
}

#' @param key A vector of gene or peak in our dataset.
#' @param omics The omics to plot: snRNA_spatial, snATAC_spatial, Celltype_spatial, GEOseq.
#' @param slot_name The slot name of data.
#' @param sigmas parameter in gaussian kernel, default:0.5.
#' @param grid_smooth Set smooth method: all, column, row, default all.
#'
#' @rdname CornPlot
#' @export
#' @method CornPlot spatial
#'
#' @importFrom methods slot
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @importFrom stats dist
#'
#'
SmoothData.spatial = function(object, key, omics = 'snRNA_spatial', slot_name = 'exp_mat', sigmas = 0.5, grid_smooth = 'col'){

  data = slot(slot(object, name = omics), name = slot_name)

  if(length(key) == 1){

    if(!key %in% rownames(data)){
      stop('Please input the gene or peak in our dataset!')
    }

    key_data = unlist(data[key, ])
    names(key_data) = sapply(colnames(data), function(x) unlist(strsplit(x, '_'))[2])

  }else{

    if(length(intersect(key, rownames(data))) == 0){
      stop('Please input the gene or peak in our dataset!')
    }

    cells_rankings = AUCell_buildRankings(as.matrix(slot(slot(object, name = omics), name = slot_name)), plotStats=FALSE, splitByBlocks=TRUE)
    cells_AUC = AUCell_calcAUC(list(signal = key), cells_rankings, aucMaxRank = floor(nrow(cells_rankings)), verbose = F)
    key_data = c(t(getAUC(cells_AUC)))
    names(key_data) = sapply(colnames(slot(slot(object, name = omics), name = slot_name)), function(x){unlist(strsplit(x,'_'))[2]}, USE.NAMES = F)

  }


  if(grid_smooth == 'all'){
    distance = as.matrix(dist(object@Coordinate@coordinate[,1:2]))
  }

  if(grid_smooth == 'col'){
    grid = matrix(10000, nrow = nrow(object@Coordinate@coordinate),
                  ncol = nrow(object@Coordinate@coordinate),
                  dimnames = list(object@Coordinate@coordinate[,3], object@Coordinate@coordinate[,3]))
    layers = sapply(colnames(grid), function(x) gsub('\\d', '', x))
    for(i in 1:nrow(grid)){
      layer = gsub('\\d','',rownames(grid)[i])
      ind = which(layers %in% layer)
      if(layer == 'EAP'){
        ind = c(which(layers %in% layer), which(layers %in% 'EA'), which(layers %in% 'EP'))
      }
      if(layer == 'AP'){
        ind = c(which(layers %in% layer), which(layers %in% 'A'), which(layers %in% 'P'))
      }
      grid[i, ind] = 1
    }
    grid[lower.tri(grid)] = t(grid)[lower.tri(t(grid))]
    distance = as.matrix(dist(object@Coordinate@coordinate[,1:2])) * grid
  }

  if(grid_smooth == 'row'){
    grid = matrix(10000, nrow = nrow(object@Coordinate@coordinate),
                  ncol = nrow(object@Coordinate@coordinate),
                  dimnames = list(object@Coordinate@coordinate[,3], object@Coordinate@coordinate[,3]))
    heights = sapply(colnames(grid), function(x) gsub('\\D', '', x))
    for(i in 1:nrow(grid)){
      height = gsub('\\D','',rownames(grid)[i])
      ind = which(heights %in% height)
      grid[i, ind] = 1
    }
    distance = as.matrix(dist(object@Coordinate@coordinate[,1:2])) * grid
  }

  rownames(distance) = colnames(distance) = object@Coordinate@coordinate[,3]
  gaussian_kernel_weight = matrix(sapply(unlist(distance), gaussian_kernel, sigmas), nrow = nrow(distance), ncol = ncol(distance), dimnames = list(rownames(distance), colnames(distance)))

  key_data_smooth = c()
  for(i in 1:length(key_data)){
    ind = match(names(key_data)[i], rownames(gaussian_kernel_weight))
    inds = match(names(key_data), colnames(gaussian_kernel_weight))
    key_data_smooth[i] = key_data %*% gaussian_kernel_weight[ind, inds]
  }
  names(key_data_smooth) = names(key_data)
  key_data = key_data_smooth

  return(key_data)



}


#----smooth cornplot----
