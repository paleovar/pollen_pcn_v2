# larger and/or reused plotting functions for maps

# pre-load shapefiles if available
if (file.exists(file.path(DIR_CACHE,'pre_transformed_shapefiles.RData'))) load(file.path(DIR_CACHE,'pre_transformed_shapefiles.RData'))

load_natural_earth_data <- function(file, dir = NAT_EARTH_DATA_PATH, ...) {
  #unzip(file.path(NAT_EARTH_DATA_PATH, '10m_physical.zip'), exdir = tempdir())
  if (dir == NAT_EARTH_DATA_PATH & str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(file.path(NAT_EARTH_DATA_PATH, '10m_physical'), file), verbose = FALSE, ...)
  } else if (str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(dir, file), verbose = FALSE, ...)
  } else if (str_detect(file, '.tif')) {
    data <- raster::raster(x = file.path(NAT_EARTH_DATA_PATH, file))
  }
  #file.remove(tempdir())
  return(data)
}


transform_shapefile <- function(data, transform) {
  if (transform == 'fixed') {
    return(data)
  }
  else {
    return(spTransform(data, CRSobj = CRS(transform)))
  }
}


#' @rdname acer.maps
#' @param ggobj 
#' @param projection either 'fixed' or 'robinson'
#' @param graticules either TRUE, FALSE or an integer specifying graticule spacing
#' @param topo logical, should topography be plotted?
#' @param no_bathy_clr vector of colors (length = 2): colours of land, ocean if bathy is not used
#' @param bathy logical, should bathymetry be plotted?
#' @param pnv logical, should potential natural vegetation be plotted?
#' @param lgm_ice logical, should potential ice extent for LGM (18 ka b.p.) be plotted?
#' @param black logical, black background?
#' @param cb_friendly logical, use \code{\link[viridis]{viridis}} color scheme?
#' @param bg arbitrary background as raster object
#'
#' @return gg object
shapefile_map_plot_ <- function(ggobj, 
                                projection = 'robinson', 
                                graticules = 30, 
                                grat_labs = NULL,
                                latlabs = 'outside',
                                topo = FALSE, 
                                bathy = FALSE,
                                no_bathy_clr = c('#f0e6c2', '#83b9df'),
                                pnv = FALSE, 
                                lgm_ice = FALSE, 
                                black = FALSE, 
                                cb_friendly = TRUE, 
                                bg = NULL,
                                zoom = NULL,
                                #yaxislclip = -90,
                                #xaxisrclip = 180, 
                                #xaxislclip = -180,
                                ...) {
  if (pnv & projection != 'fixed') {warning('using fixed projection for pnv = TRUE'); projection <- 'fixed'; topo <- FALSE}
  if (!exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv)) {load <- TRUE; exists <- FALSE}
  else if(exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv) & !(projection %in% names(shapefile_map_plot_data))) {load <- TRUE; exists <- TRUE}
  else {load <- FALSE}
  
  load_topo <- FALSE
  
  if (load) {
    cat(cr('caching shapefile_map_plot_data to .GlobalEnv for projection', projection, '\n'))
    map_data <- load_natural_earth_data(file = 'ne_10m_land.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    ocean_data <- load_natural_earth_data(file = 'ne_10m_ocean.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    dupl_id <- {if (is.logical(graticules)) {180/30-1} else if (90 %% graticules == 0) {180/graticules-1} else {180/graticules}}
    grat <- load_natural_earth_data(file = paste('ne_10m_graticules_', {if (is.logical(graticules)) {30}else {graticules}}, '.shp', sep = '')) %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.) %>% 
      mutate_at(vars(id), as.integer)
    grat <- grat %>% mutate_at(vars(group), as.character) %>% 
      bind_rows(grat %>% filter(id == dupl_id) %>%
                  mutate(long=-long, id = max(grat$id) + 1, group = as.character(max(grat$group %>% levels() %>% as.numeric()) + 1)))
    bbox <- load_natural_earth_data(file = 'ne_10m_wgs84_bounding_box.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    bathy_data <- lapply(sort(c(200, seq(0, 10000, 1000))), function (x, L) {
      load_natural_earth_data(file = paste0('ne_10m_bathymetry_', L[[as.character(x)]], '_', x, '.shp')) %>% 
        transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
        fortify(.)}, L = LETTERS[seq(1, 12, 1) %>% rev()] %>% setNames(., c(200, seq(0, 10000, 1000)) %>% sort())) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    lgm_ice_data <- load_natural_earth_data(dir = file.path(DIR_DATASETS_SUPPLEMENTARY, 'lgm_18k_ice_extent'), file = 'lgm.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    sfile_map_plot_data <- list(); sfile_map_plot_data[[projection]] <- list(map_data = map_data, ocean_data = ocean_data, bathy_data = bathy_data, grat = grat, bbox = bbox, lgm_ice_data = lgm_ice_data)
    if(!exists) {shapefile_map_plot_data <<- sfile_map_plot_data}else {shapefile_map_plot_data <<- merge.list(shapefile_map_plot_data, sfile_map_plot_data)}
  }
  else {
    map_data <- shapefile_map_plot_data[[projection]]$map_data; ocean_data <- shapefile_map_plot_data[[projection]]$ocean_data
    grat <- shapefile_map_plot_data[[projection]]$grat; bbox <- shapefile_map_plot_data[[projection]]$bbox
    bathy_data <- lapply(as.vector(names(shapefile_map_plot_data[[projection]]$bathy_data)), function(d, x){as_tibble(x[[d]]) %>% mutate(depth = d)},
                         x = shapefile_map_plot_data[[projection]]$bathy_data) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    if (lgm_ice) {lgm_ice_data <- shapefile_map_plot_data[[projection]]$lgm_ice_data}
  }
  
  if (load_topo) {
    stop('topo not allowed at the moment')
    cat('caching shapefile_topo to .GlobalEnv for projection', projection, '\n')
    topo_data <- raster::stack(file.path(NAT_EARTH_DATA_PATH, 'HYP_LR', 'HYP_LR.tif')) %>% as(., 'SpatialPixelsDataFrame') %>% as.data.frame(); colnames(topo_data) <- c('z', 'x', 'y')
  }
  
  if (!is.null(zoom)) {
    yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
    zoom <- c(zoom[1], 0, 0, zoom[2], zoom[3], 0, 0, zoom[4])
    zoom <- project(matrix(c(unlist(zoom)), nrow = 4, byrow = T), proj = GLOBAL_CRS[[projection]]) #;zoom <- project(matrix(c(unlist(zoom)), nrow = 2, byrow = T), proj = GLOBAL_CRS[[projection]])
    zoom <- c(zoom[1,1], zoom[2,2], zoom[3,1], zoom[4,2])
    zoom <- matrix(zoom, nrow = 2, byrow = T)
  } else {
    yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
  }
  
  gobj <- ggobj + 
    geom_polygon(data = bbox, mapping = aes(x = long, y = lat, group = group), fill = GLOBAL_GREY_DARK)
  if (topo) {
    ggobj <- ggobj + 
      geom_tile(data = topo_data, mapping = aes(x = x, y = y, fill = z))
  }
  else if (bathy) {
    gobj <- gobj + 
      geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group), fill = '#d9ebf9', na.rm = T) + 
      geom_polygon(data = bathy_data[['200']], mapping = aes(x = long, y = lat, group = group), fill = '#cae1f4', na.rm = T) + 
      geom_polygon(data = bathy_data[['1000']], mapping = aes(x = long, y = lat, group = group), fill = '#afd3ef', na.rm = T) + 
      geom_polygon(data = bathy_data[['2000']], mapping = aes(x = long, y = lat, group = group), fill = '#aacde9', na.rm = T) + 
      geom_polygon(data = bathy_data[['3000']], mapping = aes(x = long, y = lat, group = group), fill = '#96c1e3', na.rm = T) + 
      geom_polygon(data = bathy_data[['4000']], mapping = aes(x = long, y = lat, group = group), fill = '#83b9df', na.rm = T) + 
      geom_polygon(data = bathy_data[['5000']], mapping = aes(x = long, y = lat, group = group), fill = '#6fadd6', na.rm = T) + 
      geom_polygon(data = bathy_data[['6000']], mapping = aes(x = long, y = lat, group = group), fill = '#5ba2d0', na.rm = T) + 
      geom_polygon(data = bathy_data[['7000']], mapping = aes(x = long, y = lat, group = group), fill = '#589cc9', na.rm = T) + 
      geom_polygon(data = bathy_data[['8000']], mapping = aes(x = long, y = lat, group = group), fill = '#337fb2', na.rm = T) + 
      geom_polygon(data = bathy_data[['9000']], mapping = aes(x = long, y = lat, group = group), fill = '#2a77ac', na.rm = T) + 
      geom_polygon(data = bathy_data[['10000']], mapping = aes(x = long, y = lat, group = group), fill = '#2371a6', na.rm = T)
    if (!pnv) {
      gobj <- gobj + 
        geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole)) + 
        scale_fill_manual(values = c('#f0e6c2', NA), guide = FALSE)
    }
  }else if(!pnv & is.null(bg)) {
    gobj <- gobj + 
      #geom_polygon(data = bathymap, mapping = aes(x = long, y = lat, group = group, fill = hole), na.rm = T) + 
      geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group,  fill = TRUE, color = TRUE)) + 
      geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole, color = hole)) + 
      scale_fill_manual(values = no_bathy_clr, guide = FALSE) + 
      scale_color_manual(values = no_bathy_clr, guide = FALSE) + 
      new_scale_fill() + new_scale_color() # grey scale: scale_fill_manual(values = c(GLOBAL_GREY_LIGHT, GLOBAL_GREY_DARK), guide = FALSE) 
  }
  
  if (pnv & is.null(bg)) {
    pnv <- load_pnv_nc_file(NAME_POT_VEG_DISTR)
    gobj <- gobj + 
      geom_polygon(data = filter(map_data, lat < -60), mapping = aes(x = long, y = lat, group = group), fill = grey.colors(1, start = 0.94, gamma = 0.1)) +  
      geom_tile(data = pnv,
                mapping = aes(x=lat, y=lon, fill=megabiome, color=megabiome),
                show.legend = c(color = FALSE, fill = TRUE))
    if (cb_friendly) {
      gobj <- gobj + 
        scale_fill_viridis_d(begin = 0.2, guide = global_legend_horiz('megabiome', ncol = 5)) + 
        scale_color_viridis_d(begin = 0.1, guide = global_legend_horiz('megabiome', ncol = 5))
    }else {
      gobj <- gobj + 
        scale_fill_manual(values = PNV_COLORING) + 
        scale_color_manual(values = PNV_COLORING)
    }
    gobj <- gobj + 
      new_scale('fill') + new_scale('colour')
  }
  
  if (lgm_ice & is.null(bg)) {
    gobj <- gobj + 
      geom_polygon(data = lgm_ice_data, mapping = aes(x = long, y = lat, group = group), fill = NA, #fill = grey.colors(1, 0.96, gamma = 0.2),
                   alpha = 0.15, color = grey.colors(1, 0.2, gamma = 0.1), size = 0.45)
  }
  
  if (!is.null(bg)) {
    if (projection == 'fixed') {
      bg <- bg %>% as_tibble(); colnames(bg) <- famous_lat
      bg <- tibble(long = famous_lon) %>% bind_cols(bg) %>% gather(key = 'lat', value = 'slope', -long) %>% mutate_at(vars(lat), as.numeric)
      colnames(bg) <- c(colnames(bg)[1:2], 'z')
    } else if (projection != 'fixed') {
      print('#')
      bg <- raster::raster(t(bg), crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% raster::flip('y') %>% 
        raster::rasterToPolygons() %>% 
        spTransform(GLOBAL_CRS[[projection]]) %>% as('sf')
    }
    gobj <- ggplot() + 
      geom_sf(data = bg, mapping = aes(color = layer, fill = layer)) + 
      geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group), fill = NA, colour = 'black', size = 0.25)
  }
  
  if (is.logical(graticules)) {
    xgrats <- c(-90, 0, 90); ygrats <- c(-45, 0, 45); xgrats_names <- c(-90, 0, 90); ygrats_names <- c(-45, 0, 45)
  } else {
    if (!is.null(grat_labs)) {
      xgrats <- grat_labs$x
      ygrats <- grat_labs$y
    } else {
      xgrats <- seq(-180, 180, graticules) %>% ifelse(. >= xaxisrclip | . <= xaxislclip, NA, .) %>% .[!is.na(.)]
      ygrats <- seq(-90, 90, graticules)
    }
    xgrats_names <- xgrats
    ygrats_names <- ygrats
  }
  if (!projection == 'fixed') {
    graticules <- tibble(x = c(xgrats, rep(xaxislclip, length(ygrats))), y = c(rep(yaxislclip, length(xgrats)), ygrats))
    graticules <- as_tibble(project(cbind(graticules$x, graticules$y), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(x = V1, y = V2)
    xgrats <- graticules$x[1:length(xgrats)]; ygrats <- graticules$y[(length(xgrats) + 1):(length(ygrats) + length(xgrats))]
    xgratsy <- graticules$y[1:length(xgrats)]; ygratsx <- graticules$x[(length(xgrats) + 1):(length(ygrats) + length(xgrats))]# - 2*abs(earth_rad$x)
  }
  
  if (!is.logical(graticules)) {plot_grats <- TRUE} else {if (graticules) {plot_grats <- FALSE} else {plot_grats <- FALSE}}
  # add axis labels work-around
  #for (i in 1:length(ygrats)) {
  if (latlabs == 'outside') {
    latlabs_hjust <- rep(1.45, length(ygrats))
    latlabs_vjust <- rep(0.5, length(ygrats))
  } else if (latlabs == 'inside') {
    latlabs_hjust <- rep(-0.15, length(ygrats))
    latlabs_vjust <- rep(-0.25, length(ygrats))
  } else {
    stop('latlabs has to be `inside` or `outside`')
  }
  gratdat <- tibble(
    x = c(xgrats, ygratsx),
    y = c(xgratsy, ygrats),
    label = as.character(abs(c(xgrats_names, ygrats_names))),
    dir = c(case_when(xgrats_names > 0 ~ 'E',
                      xgrats_names < 0 ~ 'W',
                      xgrats_names == 0 ~ ''),
            case_when(ygrats_names < 0 ~ 'S',
                      ygrats_names > 0 ~ 'N',
                      ygrats_names == 0 ~ '')),
    hjust = c(rep(0.5, length(xgrats)),latlabs_hjust),
    vjust = c(rep(1.25, length(xgrats)),latlabs_vjust)
  ) %>% 
    rowwise() %>% 
    mutate(label = if_else(dir != '', paste0(label,'\u00B0',dir), paste0(label,'\u00B0'))) 
  # was *degree* instead \u00B0, but this way it works better with cairo_pdf -> change parse = T/F flag in annotate below accordingly
  if (latlabs == 'outside') {
    gratdat <- gratdat %>% 
      mutate(hjust = if_else(label == '90\u00B0N' | label == '90\u00B0S', hjust + sign(hjust) * 0.25, hjust)) # extra space at poles
  }
  
  gobj <- gobj + 
    {if (!plot_grats){} else {geom_path(data = grat, aes(long, lat, group=group, fill = NULL), linetype="dashed", color="grey50")}}
  
  gobj <- gobj + 
    annotate('rect', fill = 'white', color = 'white',
             xmin = min(xgrats), xmax = max(xgrats) - min(xgrats),
             ymin = min(xgratsy) - 1 * abs(zoom[1,2]), ymax = min(xgratsy))# + 0.1*abs(min(xgratsy)))
  for (i in 1:length(gratdat$x)) {
    gobj <- gobj + 
      annotate('text', x = gratdat$x[i], y = gratdat$y[i],
               label = gratdat$label[i], hjust = gratdat$hjust[i], vjust = gratdat$vjust[i], parse = F, color = 'grey50', size = GLOBAL_FONT_SIZE - 11)
  }
  
  gobj <- gobj + 
    {if (!is.null(bg)) {coord_sf(crs = GLOBAL_CRS[[projection]], xlim = zoom[,1], ylim = zoom[,2] - c(0.5 * abs(zoom[1,2]), 0))} else {coord_equal(xlim = zoom[,1], ylim = zoom[,2])}} + 
    global_title_and_axis() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom', legend.box = 'horizontal') + 
    {if (black) {theme(rect = element_rect(fill = 'black'), text = element_text(colour = 'white'),
                       panel.background = element_rect('black'), panel.border = element_rect(color = 'white'),
                       legend.background = element_rect(color = 'white'), strip.text = element_text(colour = 'white'))}}
  return(gobj)
}


shapefile_map_plot_bg_ <- function(ggobj, 
                                   projection = 'robinson', 
                                   graticules = 30, 
                                   grat_labs = NULL,
                                   latlabs = 'outside',
                                   no_bg_clr = c('#f0e6c2', '#83b9df'),
                                   bg_clr = NULL,
                                   lgm_ice = FALSE, 
                                   cb_friendly = TRUE, 
                                   bg = NULL,
                                   zoom = c(-180,-90,180,90),
                                   #yaxislclip = -90,
                                   #xaxisrclip = 180, 
                                   #xaxislclip = -180,
                                   ...) {
  if (!exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv)) {load <- TRUE; exists <- FALSE}
  else if(exists(x = 'shapefile_map_plot_data', envir = .GlobalEnv) & !(projection %in% names(shapefile_map_plot_data))) {load <- TRUE; exists <- TRUE}
  else {load <- FALSE}
  
  if (load) {
    cat(cr('caching shapefile_map_plot_data to .GlobalEnv for projection', projection, '\n'))
    map_data <- load_natural_earth_data(file = 'ne_10m_land.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    ocean_data <- load_natural_earth_data(file = 'ne_10m_ocean.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    dupl_id <- {if (is.logical(graticules)) {180/30-1} else if (90 %% graticules == 0) {180/graticules-1} else {180/graticules}}
    grat <- load_natural_earth_data(file = paste('ne_10m_graticules_', {if (is.logical(graticules)) {30}else {graticules}}, '.shp', sep = '')) %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.) %>% 
      mutate_at(vars(id), as.integer)
    grat <- grat %>% mutate_at(vars(group), as.character) %>% 
      bind_rows(grat %>% filter(id == dupl_id) %>%
                  mutate(long=-long, id = max(grat$id) + 1, group = as.character(max(grat$group %>% levels() %>% as.numeric()) + 1)))
    bbox <- load_natural_earth_data(file = 'ne_10m_wgs84_bounding_box.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    bathy_data <- lapply(sort(c(200, seq(0, 10000, 1000))), function (x, L) {
      load_natural_earth_data(file = paste0('ne_10m_bathymetry_', L[[as.character(x)]], '_', x, '.shp')) %>% 
        transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
        fortify(.)}, L = LETTERS[seq(1, 12, 1) %>% rev()] %>% setNames(., c(200, seq(0, 10000, 1000)) %>% sort())) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    lgm_ice_data <- load_natural_earth_data(dir = file.path(DIR_DATASETS_SUPPLEMENTARY, 'lgm_18k_ice_extent'), file = 'lgm.shp') %>% 
      transform_shapefile(data = ., transform = GLOBAL_CRS[[projection]]) %>% 
      fortify(.)
    sfile_map_plot_data <- list(); sfile_map_plot_data[[projection]] <- list(map_data = map_data, ocean_data = ocean_data, bathy_data = bathy_data, grat = grat, bbox = bbox, lgm_ice_data = lgm_ice_data)
    if(!exists) {shapefile_map_plot_data <<- sfile_map_plot_data}else {shapefile_map_plot_data <<- merge.list(shapefile_map_plot_data, sfile_map_plot_data)}
  }
  else {
    map_data <- shapefile_map_plot_data[[projection]]$map_data; ocean_data <- shapefile_map_plot_data[[projection]]$ocean_data
    grat <- shapefile_map_plot_data[[projection]]$grat; bbox <- shapefile_map_plot_data[[projection]]$bbox
    bathy_data <- lapply(as.vector(names(shapefile_map_plot_data[[projection]]$bathy_data)), function(d, x){as_tibble(x[[d]]) %>% mutate(depth = d)},
                         x = shapefile_map_plot_data[[projection]]$bathy_data) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
    if (lgm_ice) {lgm_ice_data <- shapefile_map_plot_data[[projection]]$lgm_ice_data}
  }
  
  if (!is.null(zoom)) {
    yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
    zoom <- c(zoom[1], 0, 0, zoom[2], zoom[3], 0, 0, zoom[4])
    zoom <- project(matrix(c(unlist(zoom)), nrow = 4, byrow = T), proj = GLOBAL_CRS[[projection]]) #;zoom <- project(matrix(c(unlist(zoom)), nrow = 2, byrow = T), proj = GLOBAL_CRS[[projection]])
    zoom <- c(zoom[1,1], zoom[2,2], zoom[3,1], zoom[4,2])
    zoom <- matrix(zoom, nrow = 2, byrow = T)
  } else {
    yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
  }
  
  gobj <- ggobj + 
    geom_polygon(data = bbox, mapping = aes(x = long, y = lat, group = group), fill = GLOBAL_GREY_DARK)
  if(is.null(bg)) {
    #map_data <- mutate_at(map_data, vars(group, id), as.numeric)
    #bathymap <- bathy_data[['0']] %>% mutate_at(vars(group, id), as.numeric) %>% 
    #  mutate(id = max(map_data$id) + 1, group = max(map_data$group) + 1) %>% 
    #  bind_rows(map_data)
    #bathy$group <- bathy$group
    gobj <- gobj + 
      #geom_polygon(data = bathymap, mapping = aes(x = long, y = lat, group = group, fill = hole), na.rm = T) + 
      geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group,  fill = TRUE, color = TRUE)) + 
      geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole, color = hole)) + 
      scale_fill_manual(values = no_bg_clr, guide = FALSE) + 
      scale_color_manual(values = no_bg_clr, guide = FALSE) + 
      new_scale_fill() + new_scale_color() # grey scale: scale_fill_manual(values = c(GLOBAL_GREY_LIGHT, GLOBAL_GREY_DARK), guide = FALSE) 
  }
  
  if (lgm_ice & is.null(bg)) {
    gobj <- gobj + 
      geom_polygon(data = lgm_ice_data, mapping = aes(x = long, y = lat, group = group), fill = NA, #fill = grey.colors(1, 0.96, gamma = 0.2),
                   alpha = 0.15, color = grey.colors(1, 0.2, gamma = 0.1), size = 0.45)
  }
  
  if (!is.null(bg)) {
    if (class(bg) != 'raster') {
      bg <- raster::raster(t(bg), crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% raster::flip('y') %>% 
        raster::rasterToPolygons() %>% 
        spTransform(GLOBAL_CRS[[projection]]) %>% as('sf')
    }
    gobj <- gobj + 
      geom_sf(data = bg, mapping = aes(color = as.character(layer), fill = as.character(layer))) + #coord_sf() + 
      scale_color_manual(values = rev(bg_clr), guide = FALSE) + 
      scale_fill_manual(values = rev(bg_clr), guide = FALSE) + 
      new_scale_color() + new_scale_fill() + 
      geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group), color = 'black',  fill = NA) 
  }
  
  if (is.logical(graticules)) {
    xgrats <- c(-90, 0, 90); ygrats <- c(-45, 0, 45); xgrats_names <- c(-90, 0, 90); ygrats_names <- c(-45, 0, 45)
  } else {
    if (!is.null(grat_labs)) {
      xgrats <- grat_labs$x
      ygrats <- grat_labs$y
    } else {
      xgrats <- seq(-180, 180, graticules) %>% ifelse(. >= xaxisrclip | . <= xaxislclip, NA, .) %>% .[!is.na(.)]
      ygrats <- seq(-90, 90, graticules)
    }
    xgrats_names <- xgrats
    ygrats_names <- ygrats
  }
  if (!projection == 'fixed') {
    graticules <- tibble(x = c(xgrats, rep(xaxislclip, length(ygrats))), y = c(rep(yaxislclip, length(xgrats)), ygrats))
    graticules <- as_tibble(project(cbind(graticules$x, graticules$y), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(x = V1, y = V2)
    xgrats <- graticules$x[1:length(xgrats)]; ygrats <- graticules$y[(length(xgrats) + 1):(length(ygrats) + length(xgrats))]
    xgratsy <- graticules$y[1:length(xgrats)]; ygratsx <- graticules$x[(length(xgrats) + 1):(length(ygrats) + length(xgrats))]# - 2*abs(earth_rad$x)
  }
  
  if (!is.logical(graticules)) {plot_grats <- TRUE} else {if (graticules) {plot_grats <- FALSE} else {plot_grats <- FALSE}}
  # add axis labels work-around
  #for (i in 1:length(ygrats)) {
  if (latlabs == 'outside') {
    latlabs_hjust <- rep(1.45, length(ygrats))
    latlabs_vjust <- rep(0.5, length(ygrats))
  } else if (latlabs == 'inside') {
    latlabs_hjust <- rep(-0.15, length(ygrats))
    latlabs_vjust <- rep(-0.25, length(ygrats))
  } else {
    stop('latlabs has to be `inside` or `outside`')
  }
  gratdat <- tibble(
    x = c(xgrats, ygratsx),
    y = c(xgratsy, ygrats),
    label = as.character(abs(c(xgrats_names, ygrats_names))),
    dir = c(case_when(xgrats_names > 0 ~ 'E',
                      xgrats_names < 0 ~ 'W',
                      xgrats_names == 0 ~ ''),
            case_when(ygrats_names < 0 ~ 'S',
                      ygrats_names > 0 ~ 'N',
                      ygrats_names == 0 ~ '')),
    hjust = c(rep(0.5, length(xgrats)),latlabs_hjust),
    vjust = c(rep(1.25, length(xgrats)),latlabs_vjust)
  ) %>% 
    rowwise() %>% 
    mutate(label = if_else(dir != '', paste0(label,'\u00B0',dir), paste0(label,'\u00B0')))
  if (latlabs == 'outside') {
    gratdat <- gratdat %>% 
      mutate(hjust = if_else(label == '90\u00B0N' | label == '90\u00B0S', hjust + sign(hjust) * 0.25, hjust)) # extra space at poles
  }
  for (i in 1:length(gratdat$x)) {
    gobj <- gobj + 
      annotate('text', x = gratdat$x[i], y = gratdat$y[i],
               label = gratdat$label[i], hjust = gratdat$hjust[i], vjust = gratdat$vjust[i], parse = F, color = 'grey50', size = GLOBAL_FONT_SIZE - 11.5)
  }
  
  gobj <- gobj + 
    {if (!plot_grats){} else {geom_path(data = grat, aes(long, lat, group=group, fill = NULL), linetype="dashed", color="grey50")}} + 
    {if (!is.null(bg)) {coord_sf(crs = GLOBAL_CRS[[projection]], xlim = zoom[,1], ylim = zoom[,2])} else {coord_equal(xlim = zoom[,1], ylim = zoom[,2])}} + 
    global_title_and_axis() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom', legend.box = 'horizontal')
  return(gobj)
}


plot_ACER_sites_on_map2 <- function(sites, sample_dating, expl_var,
                                    lgm_ice = FALSE, 
                                    names = FALSE,
                                    site_id = FALSE, 
                                    ISR = FALSE, ISR_stat = 'med', ISR_window = c(0,65000),
                                    save_plot = NULL, 
                                    legend_inside = FALSE, 
                                    #force_shapefiles = FALSE, 
                                    projection = 'fixed',
                                    bg = NULL,
                                    graticules = 30,
                                    zoom = NULL, 
                                    expl_var_tscale = 'raw',
                                    ...) {
  if (ISR) {if (!(ISR_stat %in% c('med', 'mean'))) {stop('unknown ISR_stat')}}
  if (!(projection %in% names(GLOBAL_CRS))) {
    stop('unknown projection')
  } else if (!(projection == 'fixed')) {
    sites <- arrange(sites, site_id)
    sites_transformed <- as_tibble(project(cbind(sites$long, sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
      rename(long_new = V1, lat_new = V2) %>% 
      bind_cols(select(sites, site_id))
    
    sites <- sites %>% 
      full_join(sites_transformed,
                by = 'site_id') %>% 
      select(-lat, -long) %>% 
      rename(long = long_new, lat = lat_new)
    if (!is.null(zoom)) {
      yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
    } else {
      yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
    }
    
    cat(paste('transformed locations to projection:', projection, '\n'))
  }
  
  data <- inner_join(sites, compute_age_interval_stats(samples = sample_dating,
                                                       ranges = list(list(start = ISR_window[1], stop = ISR_window[2]))) %>%
                       categorize_age_interval_stats(., ...) %>% 
                       select(-oldest_dp_in_interv, -youngest_dp_in_interv, -status_smp_per_10k, -status_ev_res_in_interv),
                     by = 'site_id') %>% 
    left_join(expl_var, by = 'site_id')
  
  #if (force_shapefiles | projection != 'fixed') {
  map_plot <- shapefile_map_plot_bg_(ggobj = ggplot(), projection = projection, grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)),
                                     graticules = graticules, no_bg_clr = c('grey55', 'grey98'), bg_clr = rev(c('#96c1e3', '#e8e8d1', '#fcfcfc')),#'#cae1f4', '#e8e8e8',# c('grey75', 'grey55', 'grey98'),
                                     lgm_ice = lgm_ice, cb_friendly = TRUE, bg = bg, zoom = zoom,
                                     yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip)
  #}else {
  #  map_plot <- bathymetric_map_plot_(ggobj = ggplot(),
  #                                    plot_bathy = bathy, resolution = 15, 
  #                                    use_contours = FALSE, plot_topo = topo)
  #}
  
  if(names){
    map_plot <- map_plot + geom_label_repel(data = data,
                                            mapping = aes(x=long, y=lat, label = paste(site_id, ' - ', site_name, sep = '')), nudge_y = 2.5, show.legend = FALSE,
                                            size = GLOBAL_FONT_SIZE - 10, fill = 'black', colour='white', label.padding = 0.15, label.r = 0,
                                            min.segment.length = 0, segment.color = 'black')
  }
  if(site_id){
    if(names){
      print('name labels overwrite site id labels')
    }
    else{
      map_plot <- map_plot + geom_label_repel(data = data,
                                              mapping = aes(x=long, y=lat, label = site_id), nudge_y = 2.5, show.legend = FALSE,
                                              size = GLOBAL_FONT_SIZE - 10, fill = 'black', colour='white', label.padding = 0.15, label.r = 0,
                                              min.segment.length = 0, segment.color = 'black')
    }
  }
  
  if (ISR) {
    if (legend_inside){
      title <- {if (ISR_stat == 'med') {paste0('Median ISR [y]\n', ISR_window[1]/1e3, '-', ISR_window[2]/1e3, ' ka BP')} else {paste0('Mean ISR [y]\n', ISR_window[1]/1e3, '-', ISR_window[2]/1e3, ' ka BP')}}
    } else {
      title <- {if (ISR_stat == 'med') {paste0('Median ISR [y] ', ISR_window[1]/1e3, '-', ISR_window[2]/1e3, ' ka BP')} else {paste0('Mean ISR [y] ', ISR_window[1]/1e3, '-', ISR_window[2]/1e3, ' ka BP')}}
    }
    map_plot <- map_plot + 
      geom_point(data = data,
                 mapping = aes(x = long, y = lat,
                               shape = {if (ISR_stat == 'med') {status_med_smp_res} else {status_mean_smp_res}}, 
                               fill = expl_var),#{if (ISR_stat == 'med') {status_med_smp_res} else {status_mean_smp_res}}), 
                 size = 5,
                 alpha = I(0.6),
                 stroke = 1) + 
      scale_shape_manual(guide = guide_legend(title = title, 
                                              order = 2,
                                              direction = {if(legend_inside){'vertical'} else {'horizontal'}},
                                              title.position = 'top'), 
                         values = c(24, 23, 21)) + 
      #scale_fill_manual(values = c('#fb1d21', '#ee9c4d', '#2e8ec9', '#8edd83'))
      scale_fill_fermenter(palette = 'Spectral', direction = -1, breaks = c(0.25,0.5,0.75,1), labels = function(x) round(x,2),
                           guide = guide_colorsteps(title = paste0('Expl. var. AP (',expl_var_tscale,')'), order = 3, show.limits = T, ticks = T, direction = 'horizontal',
                                                                                          barheight = 0.75, barwidth = 10, title.position = 'top'))
  } else {
    map_plot <- map_plot + 
      geom_point(data = data,
                 mapping = aes(x = long, y = lat, 
                               fill = expl_var), size = 3.5, alpha=I(0.6)) 
  }
  
  if(legend_inside){map_plot <- map_plot + 
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0), legend.box = 'vertical')
  } else {
    map_plot <- map_plot + 
      theme(legend.position = 'bottom', legend.box = 'horizontal')
  }
  if(!is.null(save_plot)){
    global_save_plot(save_plot = save_plot, plot = map_plot)
  }
  return(map_plot)
}
