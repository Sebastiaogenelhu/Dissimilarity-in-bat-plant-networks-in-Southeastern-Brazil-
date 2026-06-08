ggord.default <- function(obs, vecs, axes = c('1', '2'), grp_in = NULL, cols = NULL, facet = FALSE, nfac = NULL,
                          addpts = NULL, obslab = FALSE, ptslab = FALSE, ellipse = TRUE, ellipse_pro = 0.95, poly = TRUE,
                          polylntyp = 'solid', hull = FALSE, arrow = 0.4, labcol = 'black', veccol = 'black', vectyp = 'solid',
                          veclsz = 0.5, ext = 1.2, repel = FALSE, vec_ext = 1, vec_lab = NULL, size = 4, sizelab = NULL,
                          addsize = size/2, addcol = 'blue', addpch = 19, txt = 4, alpha = 1, alpha_el = 0.4, xlims = NULL,
                          ylims = NULL, var_sub = NULL, coord_fix = TRUE, parse = TRUE, grp_title = 'Groups', force = 1,
                          max.overlaps = 10, exp = c(0, 0), ...){
  
  # extend vectors by scale
  vecs <- vecs * vec_ext
  
  # tweaks to vecs for plotting
  # create vecs label  from vecs for labels
  names(vecs) <- c('one', 'two')
  vecs <- vecs[, na.omit(names(vecs))]
  
  vecs_lab <- ext * vecs
  if(is.null(vec_lab)) vecs_lab$labs <- as.character(row.names(vecs_lab))
  else{
    vecs_lab$labs <- vec_lab[row.names(vecs_lab)]
  }
  vecs$lab <- row.names(vecs)
  
  # add groups if grp_in is not null
  if(!is.null(grp_in)) obs$Groups <- factor(grp_in)
  
  # remove vectors for easier viz
  if(!is.null(var_sub)){
    
    var_sub <- paste(var_sub, collapse = '|')
    vecs <- vecs[grepl(var_sub, vecs$lab), ]
    vecs_lab <- vecs_lab[grepl(var_sub, vecs_lab$lab), ]
    
  }
  
  # add size to obs
  if(length(size) > 1){
    
    if(length(size) != nrow(obs))
      stop('size must have length equal to 1 or ', nrow(obs))
    
    obs$size <- size
    
    if(is.null(sizelab))
      sizelab <- 'Size'
    
  }
  
  ## plots
  
  # individual points
  nms <- names(obs)[1:2]
  names(obs)[1:2] <- c('one', 'two')
  obs$lab <- row.names(obs)
  p <- ggplot(obs, aes_string(x = 'one', y = 'two')) +
    scale_x_continuous(name = nms[1], limits = xlims, expand = c(exp, exp)) +
    scale_y_continuous(name = nms[2], limits = ylims, expand = c(exp, exp)) +
    theme_bw()
  
  # map size as aesthetic if provided
  if(length(size) > 1){
    
    # observations as points or text, colour if groups provided
    if(obslab){
      if(!is.null(obs$Groups))
        p <- p + geom_text(aes_string(group = 'Groups', fill = 'Groups', colour = 'Groups', label = 'lab', size = 'size'), alpha = alpha, parse = parse)
      else
        p <- p + geom_text(aes_string(size = 'size'), label = row.names(obs), alpha = alpha, parse = parse)
    } else {
      if(!is.null(obs$Groups))
        p <- p + geom_point(aes_string(group = 'Groups', fill = 'Groups', colour = 'Groups', shape = 'Groups', size = 'size'), alpha = alpha) +
          scale_shape_manual('Groups', values = rep(16, length = length(obs$Groups)))
      else
        p <- p + geom_point(aes_string(size = 'size'), alpha = alpha)
    }
    
    # change size legend title if provided
    p <- p + guides(size=guide_legend(title = sizelab))
    
  } else {
    
    # observations as points or text, colour if groups provided
    if(obslab){
      if(!is.null(obs$Groups))
        p <- p + geom_text(aes_string(group = 'Groups', fill = 'Groups', colour = 'Groups', label = 'lab'), size = size, alpha = alpha, parse = parse)
      else
        p <- p + geom_text(label = row.names(obs), size = size, alpha = alpha, parse = parse)
    } else {
      if(!is.null(obs$Groups))
        p <- p + geom_point(aes_string(group = 'Groups', fill = 'Groups', colour = 'Groups', shape = 'Groups'), size = size, alpha = alpha) +
          scale_shape_manual('Groups', values = rep(16, length = length(obs$Groups)))
      else
        p <- p + geom_point(size = size, alpha = alpha)
    }
    
  }
  
  # add species scores if addpts not null, for triplot
  if(!is.null(addpts)){
    
    addpts$lab <- paste0('italic(', row.names(addpts), ')')
    nms <- names(addpts)[1:2]
    names(addpts)[1:2] <- c('one', 'two')
    
    # pts as text labels if TRUE
    if(ptslab){
      p <- p +
        geom_text(data = addpts, aes_string(x = 'one', y = 'two', label = 'lab'),
                  size = addsize, col = addcol, alpha = alpha, parse = parse)
    } else {
      p <- p +
        geom_point(data = addpts, aes_string(x = 'one', y = 'two'),
                   size = addsize, col = addcol, alpha = alpha, shape = addpch)
    }
  }
  
  # fixed coordinates if TRUE
  if(coord_fix)
    p <- p + coord_fixed()
  
  # concentration ellipse if there are groups, from ggbiplot
  if(!is.null(obs$Groups) & ellipse) {
    
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- plyr::ddply(obs, 'Groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$one, x$two))
      mu <- c(mean(x$one), mean(x$two))
      ed <- sqrt(qchisq(ellipse_pro, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    })
    
    if(nrow(ell) == 0 & ellipse)
      stop('Insufficient observations for confidence ellipses, set ellipse = F to plot')
    
    names(ell)[2:3] <- c('one', 'two')
    
    # get convex hull for ell object, this is a hack to make it work with geom_polygon
    ell <- plyr::ddply(ell, 'Groups', function(x) x[chull(x$one, x$two), ])
    
    if(poly){
      
      p <- p + geom_polygon(data = ell, aes(group = Groups, fill = Groups), alpha = alpha_el)
      
    } else {
      
      # setup line type as grp_in or standard
      if(identical(polylntyp, obs$Groups))
        p <- p + geom_polygon(data = ell, aes_string(color = 'Groups', group = 'Groups', linetype = 'Groups'), fill = NA, alpha = alpha)
      else
        p <- p + geom_polygon(data = ell, aes_string(color = 'Groups', group = 'Groups'), linetype = polylntyp, fill = NA, alpha = alpha)
      
    }
    
  }
  
  # add convex hull if true
  if(hull){
    
    if(!is.null(obs$Groups)){
      
      # get convex hull
      chulls <- plyr::ddply(obs, 'Groups', function(x) x[chull(x$one, x$two), ])
      
      if(poly){
        
        p <- p + geom_polygon(data = chulls, aes(group = Groups, fill = Groups), alpha = alpha_el)
        
      } else {
        
        # setup line type as grp_in or standard
        if(identical(polylntyp, obs$Groups))
          p <- p + geom_polygon(data = chulls, aes(group = Groups, colour = Groups, linetype = Groups), fill = NA, alpha = alpha)
        else
          p <- p + geom_polygon(data = chulls, aes(group = Groups, colour = Groups), linetype = polylntyp, fill = NA, alpha = alpha)
        
      }
      
    } else {
      
      chulls <- obs[chull(obs$one, obs$two), ]
      
      if(poly){
        
        p <- p + geom_polygon(data = chulls, alpha = alpha_el)
        
      } else {
        
        p <- p + geom_polygon(data = chulls, linetype = polylntyp, alpha = alpha, fill = NA)
        
      }
      
    }
    
  }
  
  # set colors if provided
  if(!is.null(cols) & !is.null(obs$Groups)){
    if(length(cols) != length(unique(obs$Groups)))
      stop('col vector must have length equal to unique values in grp_in')
    p <- p +
      scale_colour_manual('Groups', values = cols) +
      scale_fill_manual('Groups', values = cols)
  }
  
  # add vectors
  if(!is.null(arrow))
    p <- p + geom_segment(
      data = vecs,
      aes_string(x = 0, y = 0, xend = 'one', yend = 'two'),
      arrow = grid::arrow(length = grid::unit(arrow, "cm")),
      colour= veccol, linetype = vectyp, size = veclsz
    )
  
  # facet if true and groups
  nlabs <- 1
  if(facet & !is.null(obs$Groups)){
    
    p <- p +
      facet_wrap(~Groups, ncol = nfac)
    
    # for labels if facetting
    nlabs <- length(unique(obs$Groups))
    
  }
  
  # add labels
  if(!is.null(txt)){
    
    # repel overlapping labels
    if(repel){
      
      p <- p + ggrepel::geom_text_repel(data = vecs_lab, aes_string(x = 'one', y = 'two'),
                                        label = rep(unlist(lapply(vecs_lab$labs, function(x) as.character(as.expression(x)))), nlabs),
                                        size = txt,
                                        parse = parse,
                                        point.padding = NA,
                                        color = labcol,
                                        force = force,
                                        max.overlaps = max.overlaps
                                        
      )
      
    } else {
      
      p <- p + geom_text(data = vecs_lab, aes_string(x = 'one', y = 'two'),
                         label = rep(unlist(lapply(vecs_lab$labs, function(x) as.character(as.expression(x)))), nlabs),
                         size = txt,
                         parse = parse,
                         color = labcol
      )
      
    }
    
  }
  
  # set legend titles to all scales
  p <- p +
    guides(
      fill = guide_legend(title = grp_title),
      colour = guide_legend(title = grp_title),
      shape = guide_legend(title = grp_title)
    )
  
  return(p)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord PCA
ggord.PCA <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- ord_in$eig[as.numeric(axes), 'percentage of variance']
  axes <- paste0('Dim.', axes)
  obs <- data.frame(ord_in$ind$coord)
  obs <- obs[, names(obs) %in% axes]
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$var$coord)
  vecs <- vecs[, names(vecs) %in% axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord MCA
ggord.MCA <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- ord_in$eig[as.numeric(axes), 'percentage of variance']
  axes <- paste0('Dim.', axes)
  obs <- data.frame(ord_in$ind$coord)
  obs <- obs[, names(obs) %in% axes]
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$var$coord)
  vecs <- vecs[, names(vecs) %in% axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord mca
ggord.mca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  obs <- data.frame(ord_in$rs[, as.numeric(axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$cs[, as.numeric(axes)])
  exp_var <-  100 * ord_in$d/(ord_in$p - 1)
  exp_var <- exp_var[as.numeric(axes)]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord acm
ggord.acm <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$co[, paste0('Comp', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord prcomp
ggord.prcomp <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('PC', axes)
  obs <- data.frame(ord_in$x[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$rotation[, axes])
  exp_var <- 100 * summary(ord_in)$importance[2, ][axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord princomp
ggord.princomp <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('Comp.', axes)
  obs <- data.frame(ord_in$scores[, axes])
  obs$Groups <- grp_in
  vecs <- loadings(ord_in)
  dims <- dim(vecs)
  rownms <- row.names(loadings(ord_in))
  colnms <- colnames(loadings(ord_in))
  vecs <- matrix(vecs, nrow = dims[1], ncol = dims[2])
  vecs <- data.frame(vecs, row.names = rownms)
  names(vecs) <- colnms
  vecs <- vecs[, axes]
  exp_var <- 100 * (ord_in$sdev^2)/sum(ord_in$sdev^2)
  exp_var <- exp_var[axes]
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord metaMDS
ggord.metaMDS <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('MDS', axes)
  obs <- data.frame(ord_in$points[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$species[, axes])
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord lda
ggord.lda <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- 100 * ord_in$svd^2 / sum(ord_in$svd^2)
  exp_var <- exp_var[as.numeric(axes)]
  axes <- paste0('LD', axes)
  obs <- data.frame(predict(ord_in)$x[, axes])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$scaling[, axes])
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord pca
ggord.pca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
  exp_var <- exp_var[as.numeric(axes)]
  ax_nms <- paste0('Axis', axes)
  obs <- data.frame(ord_in$li[, ax_nms])
  obs$Groups <- grp_in
  
  # vector locations, uses scaling from factoextra
  vc_nms <- paste0('Comp', axes)
  vecs <- data.frame(ord_in$co[, vc_nms])
  r <- min((max(obs[, ax_nms[1]]) - min(obs[, ax_nms[1]])/(max(vecs[, vc_nms[1]]) -
                                                             min(vecs[, vc_nms[1]]))), (max(obs[, ax_nms[2]]) - min(obs[, ax_nms[2]])/(max(vecs[,
                                                                                                                                                vc_nms[2]]) - min(vecs[, vc_nms[2]]))))
  vecs[, vc_nms] <- vecs[, vc_nms] * r * 0.7
  
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord coa
ggord.coa <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$co[, paste0('Comp', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord ca
ggord.ca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # axis contribution
  exp_var <- with(ord_in, 100 * sv^2/sum(sv^2))
  exp_var <- exp_var[as.numeric(axes)]
  
  # data to plot
  sv <- ord_in$sv[as.numeric(axes)]
  colmass <- ord_in$colmass[as.numeric(axes)]
  axes <- paste0('Dim', axes)
  obs <- ord_in$rowcoord[, axes]
  
  # scale the observations and vectors correctly
  obs <- as.data.frame(t(apply(obs, 1, "*", sv)))
  obs$Groups <- grp_in
  vecs <- ord_in$colcoord[, axes]
  vecs <- as.data.frame(t(apply(vecs, 1, "*", sv)))
  
  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord ppca
ggord.ppca <- function(ord_in, grp_in = NULL, axes = NULL, ...){
  
  # data to plot
  obs <- as.data.frame(ord_in$li)
  vecs <- as.data.frame(ord_in$c1)
  
  if(is.null(axes)){
    axes <- gsub('^PC', '',  names(obs))
  } else {
    chk <- any(!names(obs) %in% paste0('PC', axes))
    if(chk) stop('axes not found in input')
  }
  
  # axis contribution
  exp_var <- summary(ord_in, printres = FALSE)$ppca$var[as.numeric(axes)]
  exp_var <- 100 * exp_var
  
  # scale the observations and vectors correctly
  obs$Groups <- grp_in
  
  # make a nice label for the axes
  axes <- paste0('PC', axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord rda
ggord.rda <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('RDA', axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  addpts <- data.frame(ord_in$CCA$v[, axes])
  
  # vectors for constraining matrix
  constr <- data.frame(ord_in$CCA$biplot[, axes])
  
  # explained variance of constrained axes
  exp_var <- summary(ord_in)$concont$importance[2, axes]
  
  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(100 * exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs = constr, axes, addpts = addpts, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord capscale
ggord.capscale <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('CAP', axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  addpts <- data.frame(ord_in$CCA$v[, axes])
  
  # vectors for constraining matrix
  constr <- data.frame(ord_in$CCA$biplot[, axes])
  
  # explained variance of constrained axes
  exp_var <- summary(ord_in)$concont$importance[2, axes]
  
  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(100 * exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs = constr, axes, addpts = addpts, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord dbrda
ggord.dbrda <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('dbRDA', axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  
  # vectors for constraining matrix
  constr <- data.frame(ord_in$CCA$biplot[, axes])
  
  # explained variance of constrained axes
  exp_var <- summary(ord_in)$concont$importance[2, axes]
  
  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(100 * exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs = constr, axes, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord cca
ggord.cca <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # data to plot
  axes <- paste0('CCA', axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  addpts <- data.frame(ord_in$CCA$v[, axes])
  
  # vectors for constraining matrix
  constr <- data.frame(ord_in$CCA$biplot[, axes])
  
  # explained variance of constrained axes
  exp_var <- summary(ord_in)$concont$importance[2, axes]
  
  # make a nice label for the axes
  axes <- paste0(axes, ' (', round(100 * exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs = constr, axes, addpts = addpts, ...)
  
}

#' @rdname ggord
#'
#' @export
#'
#' @method ggord dpcoa
ggord.dpcoa <- function(ord_in, grp_in = NULL, axes = c('1', '2'), ...){
  
  # ord_in$dls = coordinates of the species
  # ord_in$li  = coordinates of the communities
  # ord_in$c1  = scores of principal axes of the species
  
  # data to plot
  exp_var <- 100 * ord_in$eig / sum(ord_in$eig)
  exp_var <- exp_var[as.numeric(axes)]
  obs <- data.frame(ord_in$li[, paste0('Axis', axes)])
  obs$Groups <- grp_in
  vecs <- data.frame(ord_in$dls[, paste0('CS', axes)])
  axes <- paste0('Axis', axes)
  axes <- paste0(axes, ' (', round(exp_var, 2), '%)')
  names(obs)[1:2] <- axes
  
  ggord.default(obs, vecs, axes, ...)
  
}

