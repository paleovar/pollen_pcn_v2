# larger and/or reused plotting functions for networks on maps


#' plot a palaeo climate network on a world map
#' 
#' @rdname acer.networks.plot
#'
#' @param graph \code{\link[tidygraph:tbl_graph]{tbl_graph}}. graph object created by function \code{\link{create_network_from_corr_list}}
#' @param projection character. map projection, one of c('fixed', 'robinson')
#' @param bathy logical. should bathymetry data be included?
#' @param bg raster. additional background layer. bathy value is ignored if a background layer is passed
#' @param split_pos_neg logical. should positive and negative correlation link be split in the plot?
#' @param split_sgn logical. should positive and negative signum links be split in the plot? Only used in case a column diff_sgn is given in graph-edges list (see \code{\link{create_network_from_corr_list}}).
#' @param straight logical. should links be plotted straight, thus without curvature?
#' @param forced_bundle logical. should links be forced to bundles? - Just try out which options looks better. Has effect only if graph has hierarchial structure (see \code{\link{create_network_from_corr_list}}).
#' @param as_list logical. should facetting be done internally or should a distinct plot be returned for each facet? Helpful for additional plot combination via e.g. \code{\link[cowplot:cowplot-package]{cowplot}}.
#' @param zoom numeric. vector of length four giving zoom window coordinates in the following way c(lowerleftx, lowerlefty, upperrightx, upperrighty)
#' @param node_degree logical. should the node degree be included in the plot?
#' @param color_strength logical. map colour strength to link strength? If FALSE would just map colour to the sign of the link
#' @param split_regions_by_sign logical. should different regions be assigned for positive (incl zero) and negative (excl zero) links? Might enhance readability when plotted. Use same value as in create_network_from_corr_list()
#' @param plot_title logical. add a title to the plot showing the number of nodes that are used.
#' @param save_plot list. of shape list(activate = logical, dir = character)
#' @param ... passed to underlying shapefile_map_plot_
#'
#' @return ggraph object (the plot)
#' @export
#'
#' @examples
plot_network_spatial <- function(graph, projection, bathy = FALSE, bg = NULL, split_pos_neg = FALSE, split_sgn = FALSE, ndg_lims = NULL, leg_opt = 'single',
                                 straight = FALSE, forced_bundle = TRUE, as_list = FALSE, split_regions_by_sign = F,
                                 zoom = NULL, node_degree = FALSE, link_width = 0.4, color_strength = TRUE, plot_title = FALSE, save_plot, ...) {
  if (!(leg_opt %in% c('single','cmb'))) stop('unknown leg_opt given to plot_network_spatial')
  if ('corr_bundle' %in% names(graph)) {corr_bundle <- graph$corr_bundle; bundled <- TRUE}else {bundled <- FALSE}
  if ('occurs' %in% names(graph)) {node_occurences <- graph$occurs; occurs <- TRUE}else {occurs <- FALSE}
  isorigcb <- FALSE
  if ('orig_corr_bundle' %in% names(graph) & node_degree) {orig_corr_bdl <- graph$orig_corr_bundle; isorigcb <- TRUE} else if (node_degree) {isorigcb <- FALSE; warning('no original corr bundle provided to plot_network_spatial for node degree plotting')}
  graph <- graph$graph %>% activate(nodes) %>% arrange(node_id)
  if (!is.null(bg)) shp_map <- shapefile_map_plot_bg_ else shp_map <- shapefile_map_plot_
  zoom_orig <- zoom
  if (!(projection %in% names(GLOBAL_CRS))) {
    stop('unknown projection given to plot_network_spatial')
  }
  else if (!(projection == 'fixed')) {
    sites <- as.data.frame(graph %>% activate(nodes))
    sites_transformed <- as_tibble(project(cbind(sites$long, sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(long_new = V1, lat_new = V2) %>% 
      bind_cols(select(sites, node_id))
    
    graph <- graph %>% 
      activate(nodes) %>% 
      full_join(sites_transformed,
                by = 'node_id') %>% 
      select(-lat, -long) %>% 
      rename(long = long_new, lat = lat_new)
    if (!is.null(zoom)) {
      yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
      zoom <- c(zoom[1], 0, 0, zoom[2], zoom[3], 0, 0, zoom[4])
      zoom <- project(matrix(c(unlist(zoom)), nrow = 4, byrow = T), proj = GLOBAL_CRS[[projection]]) #;zoom <- project(matrix(c(unlist(zoom)), nrow = 2, byrow = T), proj = GLOBAL_CRS[[projection]])
      zoom <- c(zoom[1,1], zoom[2,2], zoom[3,1], zoom[4,2])
      zoom <- matrix(zoom, nrow = 2, byrow = T)
    } else {
      yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
    }
    
    cat(paste('transformed network to projection:', projection, '\n'))
  }
  else {cat('using default fixed projection\n')}
  
  if (length(unique(activate(graph, edges) %>% as_tibble %>% .$frac_abs_strongest)) == 1) {graph <- activate(graph, edges) %>% select(-frac_abs_strongest)}
  if (split_pos_neg) {graph <- activate(graph, edges) %>% mutate(corr_sgn = if_else(corr >= 0, 1, -1)) %>% mutate(corr_sgn = factor(as.character(corr_sgn), levels = c('1', '-1')))}
  if (split_sgn) {graph <- activate(graph, edges) %>% mutate(diff_sgn = factor(as.character(diff_sgn), levels = c('1', '-1')))}
  
  if (!is.null(ndg_lims)) {
    ndg_breaks <- seq(from = ndg_lims[1], to = round(ndg_lims[2],-1)+11, by = 10)
    ndg_breaks <- c(ndg_breaks[1], round(ndg_breaks[2:length(ndg_breaks)], digits = -1))
  }
  
  graph <- activate(graph, nodes) %>% 
    rename(x  = long, y = lat)
  
  if (occurs) {node_occurences <- inner_join(node_occurences, activate(graph, nodes) %>% as_tibble(), by = 'node_id')}
  
  if (bundled) {
    corr_bundle <- filter(corr_bundle, from != to)
    plot <- lapply(unique(arrange(corr_bundle, window)$window), function (w, cb, g, no, z) {
      no <- filter(no, window == w) %>% {if (node_degree & isorigcb) {left_join(., compute_nw_node_degree(orig_corr_bdl, id_type = 'node_id'), by = c('node_id', 'window')) %>%
          mutate(node_degree = if_else(is.na(node_degree) | node_degree <= 0, 0, node_degree))} else {.}} %>% select(-window)
      cb <- filter(cb, window == w) %>% select(-window)# %>% mutate(corr = if_else(corr >= 0, 'pos', 'neg'))
      if (!color_strength) {
        cb <- cb %>% 
          mutate(corr = if_else(corr >= 0, 'pos', 'neg')); cb$corr <- factor(cb$corr, levels = c('pos', 'neg'))
          cbp <- filter(cb, corr == 'pos'); cbn <- filter(cb, corr == 'neg')
          pathsp <- as.matrix(select(cbp, from, to)) %>% split(., row(select(cbp, from, to))) %>% unname() %>% lapply(., as.integer)
          pathsn <- as.matrix(select(cbn, from, to)) %>% split(., row(select(cbn, from, to))) %>% unname() %>% lapply(., as.integer)
      } else {
        paths <- as.matrix(select(cb, from, to)) %>% split(., row(select(cb, from, to))) %>% unname() %>% lapply(., as.integer)
        bincorr <- function(c) {
          bins <- base::cut(c, breaks = c(-Inf, -0.4, 0, 0.4, Inf), 
                            labels = 1:4) %>% #labels = c('c < -0.4', '-0.4 < c < 0', '0 <= c < 0.4', 'c > 0.4'))
            as.character()
          #levels(bins) <- c('c < -0.4', '-0.4 < c < 0', '0 <= c < 0.4', 'c > 0.4')
          bins
        }
      }
      #pal <- colorRampPalette(c('blue3', 'white', 'red2'))
      pal <- tibble(corr = as.character(1:4),#c('c < -0.4', '-0.4 < c < 0', '0 <= c < 0.4', 'c > 0.4'), 
                    col = c('#0f26a6', '#13f0f0', '#f0e913', '#cc0000')) # c('blue1', 'red1', 'blue3', 'red3'))
      if (color_strength) {
        # based on REGION_LOCS_pn
        cbn <- cb %>% filter(corr < 0)
        cbp <- cb %>% filter(corr >= 0)
        if (split_regions_by_sign) {
          no <- no %>% filter(node_id <= length(no$node_id)/2)# %>% mutate(node_degree = if_else(node_degree > 40, 40, node_degree)) 
        }
        plot <- ggraph(g, layout = 'manual', x = activate(g, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$x,
                       y = activate(g, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$y) %>% 
          shp_map(ggobj = ., projection = projection, bg = bg, bathy = bathy, lgm_ice = F, no_bg_clr = c('grey55', 'grey98'), bg_clr = rev(c('#96c1e3', '#e8e8d1', '#fcfcfc')),
                  grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)), zoom = zoom_orig,
                  latlabs = 'outside')
        if (nrow(cbn) > 0) {
          plot <- plot + 
            geom_conn_bundle(#color = 'blue3',
              mapping = aes(colour = bincorr(corr)),#'pos'),
              data = {if (forced_bundle) {get_con(from = as.integer(cbn$from), to = as.integer(cbn$to), frac_abs_strongest = cbn$frac_abs_strongest, corr = cbn$corr)}
                else {get_con(from = as.integer(cbn$from), to = as.integer(cbn$to), paths = paths, frac_abs_strongest = cbn$frac_abs_strongest, corr = cbn$corr)}},
              width = I(link_width), # 0.15
              alpha = I(0.5), tension = 0.85)
        }
        if (nrow(cbp) > 0) {
          plot <- plot + 
            geom_conn_bundle(#color = 'blue3',
              mapping = aes(colour = bincorr(corr)),#'pos'),
              data = {if (forced_bundle) {get_con(from = as.integer(cbp$from), to = as.integer(cbp$to), frac_abs_strongest = cbp$frac_abs_strongest, corr = cbp$corr)}
                else {get_con(from = as.integer(cbp$from), to = as.integer(cbp$to), paths = paths, frac_abs_strongest = cbp$frac_abs_strongest, corr = cbp$corr)}},
              width = I(link_width), # 0.15
              alpha = I(0.5), tension = 0.85)
        }
        plot <- plot +
          scale_edge_color_manual(values = pal %>% filter(corr %in% unique(bincorr(cb$corr))) %>% .$col, labels = c('c < -0.4', expression('-0.4'<='c < 0'), expression('0'<='c < 0.4'), 'c > 0.4'),
                                  guide = global_legend(title = 'correlation\nstrength', override.aes = list(edge_width = 7, alpha = I(1)), ncol = 2, byrow = F))# +
          #scale_edge_color_gradient2(low = 'blue3', high = 'red2', mid = 'grey90', midpoint = 0, guide = global_legend(title = 'correlation', override.aes = list(edge_width = 5))) +
          #new_scale('color')
      } else {
        if (split_regions_by_sign) {
          no <- no %>% filter(node_id <= length(no$node_id)/2)# %>% mutate(node_degree = if_else(node_degree > 40, 40, node_degree)) 
        }
        plot <- ggraph(g, layout = 'manual', x = activate(g, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$x,
                       y = activate(g, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$y) %>% 
          shp_map(ggobj = ., projection = projection, bg = bg, bathy = bathy, yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip, lgm_ice = FALSE,
                  no_bg_clr = c('grey55', 'grey98'), bg_clr = rev(c('#96c1e3', '#e8e8d1', '#fcfcfc')),
                  grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)), zoom = zoom_orig,
                  latlabs = 'outside') +
          geom_conn_bundle(#color = 'blue3',
            mapping = aes(colour = corr),#'pos'),
            data = {if (forced_bundle) {get_con(from = as.integer(cbp$from), to = as.integer(cbp$to), frac_abs_strongest = cbp$frac_abs_strongest, corr = cbp$corr)}
              else {get_con(from = as.integer(cbp$from), to = as.integer(cbp$to), paths = pathsp, frac_abs_strongest = cbp$frac_abs_strongest, corr = cbp$corr)}},
            width = I(link_width), # 0.15
            alpha = I(0.4), tension = 0.85)
        if (nrow(cbn) != 0) {
          plot <- plot + 
            geom_conn_bundle(#color = 'red2',
              mapping = aes(colour = corr),#'neg'),
              data = {if (forced_bundle) {get_con(from = as.integer(cbn$from), to = as.integer(cbn$to), frac_abs_strongest = cbn$frac_abs_strongest, corr = cbn$corr)}
                else {get_con(from = as.integer(cbn$from), to = as.integer(cbn$to), paths = pathsn, frac_abs_strongest = cbn$frac_abs_strongest, corr = cbn$corr)}},
              width = I(link_width), # 0.15
              alpha = I(0.4), tension = 0.845) + 
            scale_edge_color_manual(values = c('#0f26a6', '#cc0000'), #c('blue3', 'red2'), 
                                    guide = global_legend(title = 'correlation\nstrength', override.aes = list(edge_width = 5))) #+ 
            #new_scale('color')
        }
        else {
          plot <- plot + 
            scale_edge_color_manual(values = c('red2'), guide = FALSE) #+ new_scale('color')
        }
      }
      
      if (!is.null(no)) {
        if (node_degree) {
          plot <- plot + 
            geom_point(data = no, mapping = aes(x = x, y = y, shape = occurs, colour = occurs, size = node_degree), fill = GLOBAL_GREY_LIGHT_ALPHA_LOW) + 
            scale_size(breaks = {if (!is.null(ndg_lims)) {ndg_breaks} else {NULL}}, #c(1, 2.5, 5, 7.5), }
                       #labels = c(1, 2.5, 5, 7.5), 
                       limits = {if (!is.null(ndg_lims)) {ndg_lims} else {NULL}},
                       range = NDG_RANGE,#0.75 * c(2, max(no$node_degree)), 
                       guide = {if(leg_opt == 'single') {global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 3)} else {FALSE}})
            #scale_size(breaks = c(10,20,30,40), labels = c('10','20','30',expression("">=40)), limits = c(0,60), range = NDG_RANGE,
            #           guide = global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 3))
        } else {
          plot <- plot + 
            geom_point(data = no, mapping = aes(x = x, y = y, shape = occurs, colour = occurs), fill = GLOBAL_GREY_LIGHT_ALPHA_LOW, size = I(3))
        }
        plot <- plot + 
          scale_shape_manual(values = c(23, 21), guide = global_legend('samples in\ntime window', override.aes = list(size = 6))) + 
          scale_color_manual(values = c(GLOBAL_BLUE_DARK, 'black'), guide = global_legend('samples in\ntime window', override.aes = list(size = 6)))
      }
      else {
        plot <- plot + 
          geom_node_point(data = activate(g, nodes) %>% filter(!hierarchial), mapping = aes(x = x, y = y))#!!!
      }
      plot <- plot + #{if (!is.null(z)) {coord_equal(xlim = z[,1], ylim = z[,2])}} + 
        theme(panel.border = element_blank(), strip.background.y = element_blank(), strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, vjust = 1, angle = 180),
              legend.direction = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}}, legend.box = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}})
      return(plot)
    }, cb = corr_bundle, g = graph, no = {if(occurs) {node_occurences}else {NULL}}, z = zoom)
  }
  else {
    graph <- activate(graph, edges) %>% filter(from != to)
    
    plot <- lapply(unique(activate(graph, edges) %>% as_tibble() %>% arrange(window) %>% .$window), function(w, g, no, z) {
      no <- filter(no, window == w) %>% {if (node_degree) {left_join(., compute_nw_node_degree(corr_bundle, id_type = 'node_id'), by = c('node_id', 'window')) %>% mutate(node_degree = if_else(is.na(node_degree), 0, node_degree))} else {.}} %>% select(-window)
      g <- activate(g, edges) %>% filter(window == w)
      
      plot <- shp_map(ggobj = ggraph(g, layout = 'nicely'), projection = projection, bg = bg, bathy = bathy,  zoom = zoom_orig,
                      yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip, grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)))
      if (!is.null(no)) {
        if (node_degree) {
          plot <- plot + 
            geom_point(data = no, mapping = aes(x = x, y = y, shape = occurs, color = occurs, size = node_degree), fill = GLOBAL_GREY_LIGHT_ALPHA_LOW) + 
            scale_size(breaks = ndg_breaks, #c(1, 2.5, 5, 7.5), 
                       #labels = c(1, 2.5, 5, 7.5), 
                       limits = ndg_lims,
                       range = NDG_RANGE,#0.75 * c(2, max(no$node_degree)), 
                       guide = {if(leg_opt == 'single') {global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 3)} else {FALSE}})
          #scale_size(breaks = c(10,20,30,40), labels = c('10','20','30',expression("">=40)), limits = c(0,60), range = NDG_RANGE,
          #           guide = global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 3))
        } else {
          plot <- plot + 
            geom_point(data = no, mapping = aes(x = x, y = y, shape = occurs, color = occurs), fill = GLOBAL_GREY_LIGHT_ALPHA_LOW, size = I(2))
        }
        plot <- plot + 
          scale_shape_manual(values = c(23, 21), guide = global_legend('samples in\n time window', order = 2, override.aes = list(size = 6))) + 
          scale_color_manual(values = c(GLOBAL_BLUE_DARK, 'black'), guide = global_legend('samples in\n time window', order = 2, override.aes = list(size = 6)))
      }
      else {plot <- plot + geom_node_point(data = activate(g, nodes), mapping = aes(x = x, y = y))}
      
      plot <- plot + {if (straight) {geom_edge_fan2(mapping = aes(x = x, y = y, color = corr), alpha = I(0.4), width = I(0.7))}else {geom_edge_hive2(mapping = aes(x = x, y = y, color = corr), alpha = I(0.4), width = I(0.7))}} + 
        scale_edge_color_gradient2(low = 'blue3', mid = GLOBAL_GREY_MEDIUM, high = 'red2', limits = c(-1, 1), guide = {if(!is.null(zoom)) {global_edge_colorbar_horiz('correlation\nstrength', order = 1)}else {global_edge_colorbar('correlation\nstrength', order = 1)}}) + 
        #{if (!is.null(z)) {coord_equal(xlim = z[,1], ylim = z[,2])}} + 
        theme(panel.border = element_blank(), strip.background.y = element_blank(), strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, vjust = 1, angle = 180),
              legend.direction = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}}, legend.box = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}})
    }, g = graph, no = {if(occurs) {node_occurences}else {NULL}}, z = zoom)
  }
  
  if (bundled) {
    if (occurs) {nms <- filter(node_occurences, occurs == TRUE) %>% group_by(window) %>% summarise(noccurs = n()) %>% rowwise() %>% mutate(nms = paste0(window, ' - # ', noccurs)) %>% arrange(window) %>% .$nms} else {nms <- unique(arrange(corr_bundle, window)$window)}
    plot <- setNames(plot, nms)
  } else {
    if (occurs) {nms <- filter(node_occurences, occurs) %>% group_by(window) %>% summarise(noccurs = n()) %>% rowwise() %>% mutate(nms = paste0(window, ' - # ', noccurs)) %>% arrange(window) %>% .$nms} else {nms <- unique(activate(graph, edges) %>% as_tibble() %>% arrange(window) %>% .$window)}
    plot <- setNames(plot, nms)
  }
  
  if (!as_list & is.null(zoom) & plot_title) {
    plot <- ggdraw() + 
      draw_plot(plot_grid(plotlist = lapply(plot, function(x) return(x + theme(legend.position = 'none', axis.title = element_blank()))), ncol = 1, labels = names(plot),
                          label_x = 0.322, label_y = 0.98, hjust = 0), 0, 0, 1, 1) + 
      draw_plot(get_legend(plot[[1]] + theme(legend.position = c(0.5, 0.1), legend.justification = 'bottom')), 0.67, 0, 0.1, 0.1)
  }
  else if (!as_list & !is.null(zoom) & plot_title) {
    ly <- abs(zoom[1,]) + abs(zoom[2,]); ly <- (ly[2] / (3 *ly[1]))
    plot <- ggdraw() + 
      draw_plot(plot_grid(plotlist = lapply(plot, function(x) return(x + theme(legend.position = 'none', axis.title = element_blank()))), nrow = 1, labels = names(plot),
                          label_x = 0.085, label_y = 1, hjust = 0), 0, 0.15, 1, 0.85) + 
      draw_plot(get_legend(plot[[1]] + theme(legend.position = c(0, 0.1), legend.justification = 'right')), 0.98, 0.1, 0.1, 0.1)
  } else if (length(plot) == 1) {
    plot <- plot[[1]] + theme(axis.title = element_blank())
  }
  return(plot)
}



plot_bubble_network_spatial <- function(graph, width_lims, bg = NULL, leg_opt = 'single',
                                        zoom = c(-180, -65, 180, 65), site_data = NULL, ndg_lims = NULL, ...) {
  projection <- 'robinson'
  if (!is.null(bg)) shp_map <- shapefile_map_plot_bg_ else shp_map <- shapefile_map_plot_
  zoom_orig <- zoom
  #if (!(projection %in% names(GLOBAL_CRS))) stop('unknown projection given to plot_bubble_network_spatial')
  if (!(leg_opt %in% c('single','cmb'))) stop('unknown leg_opt given to plot_bubble_network_spatial')
  
  reg_sites <- as.data.frame(graph %>% activate(nodes))
  reg_sites_transformed <- as_tibble(project(cbind(reg_sites$long, reg_sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
    rename(x = V1, y = V2) %>% 
    bind_cols(select(reg_sites, node_id))
  
  if (!is.null(site_data)) {
    sites_transformed <- as_tibble(project(cbind(site_data$long, site_data$lat), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(x = V1, y = V2) %>% 
      bind_cols(site_data %>% select(-lat,-long))
    sites_transformed$occurs <- factor(sites_transformed$occurs, levels = c(TRUE,FALSE))
  }
  
  graph <- graph %>% 
    activate(nodes) %>% 
    full_join(reg_sites_transformed,
              by = 'node_id') %>% 
    select(-lat, -long) %E>% 
    mutate(adj = case_when(adj == 1 ~ '+',
                           adj == -1 ~ '-')) %>%
    activate(edges) %>% 
    mutate(adj = factor(adj, levels = c('+','-')))
  
  #graph <- graph %>% activate(edges) %>%  #%>% View()
  
  if (!is.null(zoom)) {
    yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
    zoom <- c(zoom[1], 0, 0, zoom[2], zoom[3], 0, 0, zoom[4])
    zoom <- project(matrix(c(unlist(zoom)), nrow = 4, byrow = T), proj = GLOBAL_CRS[[projection]])
    zoom <- c(zoom[1,1], zoom[2,2], zoom[3,1], zoom[4,2])
    zoom <- matrix(zoom, nrow = 2, byrow = T)
  } else {
    yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
  }
  cat(paste('transformed network to projection:', projection, '\n'))
  
  width_lims <- c(width_lims[1],round(width_lims[2],-1))
  nlink_breaks <- seq(from = width_lims[1], to = width_lims[2]+50, by = 50)
  nlink_breaks <- c(nlink_breaks[1], round(nlink_breaks[2:length(nlink_breaks)], digits = -1))
  width_lims <- c(width_lims[1],nlink_breaks[length(nlink_breaks)])
  
  ndg_breaks <- seq(from = ndg_lims[1], to = round(ndg_lims[2],-1)+10, by = 10)
  ndg_breaks <- c(ndg_breaks[1], round(ndg_breaks[2:length(ndg_breaks)], digits = -1))
  
  
  
  plot <- ggraph(graph, layout = 'manual', x = activate(graph, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$x,
                 y = activate(graph, nodes) %>% as_tibble() %>% arrange(node_id) %>% .$y) %>% 
    shp_map(ggobj = ., projection = projection, bg = bg, bathy = FALSE, lgm_ice = F, no_bg_clr = c('grey55', 'grey98'), bg_clr = rev(c('#96c1e3', '#e8e8d1', '#fcfcfc')),
            grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)), zoom = zoom_orig,
            latlabs = 'outside') + 
    geom_point(data = sites_transformed, 
               mapping = aes(x = x, y = y, shape = occurs, color = occurs, size = ndg), 
               fill = GLOBAL_GREY_LIGHT_ALPHA_LOW) + 
    scale_size(breaks = ndg_breaks, #c(1, 2.5, 5, 7.5), 
               #labels = c(1, 2.5, 5, 7.5), 
               limits = ndg_lims,
               range = NDG_RANGE,#0.75 * c(2, max(no$node_degree)), 
               guide = global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 4)) + #{if(leg_opt == 'single') {global_legend_horiz('node degree', override.aes = list(shape = 23, color = GLOBAL_BLUE_DARK), order = 4)} else {FALSE}}) + 
    scale_shape_manual(values = c(23, 21), 
                       guide = {if(leg_opt == 'single') {guide_legend('samples in\ntime window', title.position = 'top', override.aes = list(size = 6), order = 3)} else {FALSE}}) + 
    scale_color_manual(values = c(GLOBAL_BLUE_DARK, 'black'), 
                       guide = {if(leg_opt == 'single') {guide_legend('samples in\ntime window', title.position = 'top', order = 3)} else {FALSE}}) + 
    geom_edge_parallel0(mapping = aes(x = x, y = y, from = from, to = to, edge_color = adj, edge_width = nlinks),
                        alpha = 0.75, sep = unit(3,'mm')) + 
    scale_edge_color_manual(values = rev(GLOBAL_PN_COLORS),
                            guide = guide_legend(title = 'correlation\nsign', 
                                                 override.aes = list(edge_width = 7, alpha = 1),
                                                 title.position = 'top',
                                                 order = 1)) +
    scale_edge_width_continuous(limits = width_lims, breaks = nlink_breaks,
                                guide = guide_legend(title = '# links', 
                                                     title.position = 'top',
                                                     order = 2, nrow = 3)) + 
    theme(panel.border = element_blank(), strip.background.y = element_blank(), strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, vjust = 1, angle = 180),
          legend.direction = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}}, legend.box = {if(!is.null(zoom)) {'horizontal'}else {'vertical'}},
          legend.background = element_blank())
  
  return(plot)
}



