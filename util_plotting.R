# global plotting elements (file named util_global_plotting_options.R in v1)
global_colorbar <- function(title = '', order = 1, ...){
  cb <- guide_colorbar(order = order, barwidth = 0.5, barheight = 10,
                       title = title, title.position = 'top', ...)
  return(cb)
}


global_colorbar_horiz <- function(title = '', order = 1, ...){
  cb <- guide_colorbar(order = order, barheight = 0.5, barwidth = 10, direction = 'horizontal', 
                       title = title, title.position = 'top', ...)
  return(cb)
}


global_edge_colorbar <- function(title = '', order = 1, ...){
  cb <- guide_edge_colorbar(order = order, barwidth = 0.5, barheight = 10,
                            title = title, title.position = 'top', ...)
  return(cb)
}


global_edge_colorbar_horiz <- function(title = '', order = 1, ...){
  cb <- guide_edge_colorbar(order = order, barheight = 0.5, barwidth = 10, direction = 'horizontal', 
                            title = title, title.position = 'top', ...)
  return(cb)
}


global_legend <- function(title = '', order = 1, key_color = NULL, key_fill = NULL, ...){
  if(is.null(key_color) & is.null(key_fill)){
    lg <- guide_legend(order = order, #keywidth = 0.5, keyheight = 0.5, 
                       title = title, title.position = 'top', ...)
  }
  else{
    lg <- guide_legend(order = order, #keywidth = 0.5, keyheight = 0.5, 
                       title = title, title.position = 'top', 
                       override.aes = list(colour = key_color, fill = key_fill, ...))
  }
  return(lg)
}


global_legend_horiz <- function(title = '', order = 1, key_color = NULL, key_fill = NULL, ...){
  if(is.null(key_color) & is.null(key_fill)){
    lg <- guide_legend(order = order, keywidth = 0.5, keyheight = 0.5, 
                       title = title, title.position = 'top', direction = 'horizontal', ...)
  }
  else{
    lg <- guide_legend(order = order, keywidth = 0.5, keyheight = 0.5, 
                       title = title, title.position = 'top',direction = 'horizontal', 
                       override.aes = list(colour = key_color, fill = key_fill, ...))
  }
  return(lg)
}


global_title_and_axis <- function(){
  ax <- theme_bw() + 
    theme(title = element_text(face = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE + 1, family = GLOBAL_FONT_FAMILY),
          text = element_text(face = GLOBAL_FONT_FACE_TEXT, size = GLOBAL_FONT_SIZE, family = GLOBAL_FONT_FAMILY), 
          #line = element_blank(),
          legend.direction = 'vertical', legend.box = 'vertical', legend.box.background = element_blank(), legend.background = element_rect(colour = 'black'), 
          strip.background.x = element_blank(), strip.text.x = element_text(face = GLOBAL_FONT_FACE_TITLE, hjust = 0), 
          strip.background.y = element_rect(fill = 'white'), strip.text.y = element_text(face = GLOBAL_FONT_FACE_TITLE, vjust = 0))
  return(ax)
}


signed_log_trans <- scales::trans_new("signed_log",
                                      transform = function(x) sign(x)*log(abs(x)),
                                      inverse = function(x) sign(x)*exp(abs(x)))

signed_log10_trans <- scales::trans_new("signed_log10",
                                        transform = function(x) sign(x)*log(abs(x), base = 10),
                                        inverse = function(x) sign(x)*exp(abs(x / 10)))

signed_atanh_trans <- scales::trans_new("signed_atanh",
                                        transform = function(x) sign(x)*atanh(abs(x)),
                                        inverse = function(x) sign(x)*tanh(abs(x)))


global_save_plot <- function(plot, save_plot, png_and_pdf = TRUE){
  require(RCurl)
  if (!('activate' %in% names(save_plot))) {save <- TRUE}else if (save_plot$activate) {save <- TRUE}else {save<- FALSE}
  if (save) {
    def_params <- list(filename = 'default', scale = 1) # could specify width, height
    save_plot <- merge.list(save_plot, def_params)
    
    cat('Saving image as pdf \n')
    if('width' %in% names(save_plot) & 'height' %in% names(save_plot)){
      ggsave(filename = paste(save_plot$filename, 'pdf', sep = '.'), plot = plot,# path = DIR_FIGURES,
             width = save_plot$width, height = save_plot$height, units = 'cm', dpi = 'print', limitsize = FALSE)
    }
    else {
      ggsave(filename = paste(save_plot$filename, 'pdf', sep = '.'), plot = plot,# path = DIR_FIGURES,
             scale = save_plot$scale, units = 'cm', dpi = 'print')
    }
    
    if (png_and_pdf) {
      cat('Saving image as png \n')
      if ('width' %in% names(save_plot) & 'height' %in% names(save_plot)) {
        ggsave(filename = paste(save_plot$filename, 'png', sep = '.'), plot = plot,# path = DIR_FIGURES,
               width = save_plot$width, height = save_plot$height, units = 'cm', dpi = 'print', limitsize = FALSE)
      }
      else {
        ggsave(filename = paste(save_plot$filename, 'png', sep = '.'), plot = plot,# path = DIR_FIGURES,
               scale = save_plot$scale, units = 'cm', dpi = 'print')
      }
    }
  }
}


global_save_plot_cairo <- function(plot, save_plot, png_and_pdf = TRUE){
  if (!('activate' %in% names(save_plot))) {save <- TRUE}else if (save_plot$activate) {save <- TRUE}else {save<- FALSE}
  if (save) {
    def_params <- list(filename = 'default', scale = 1) # could specify width, height
    save_plot <- merge.list(save_plot, def_params)
    
    cat('Saving image as pdf \n')
    if('width' %in% names(save_plot) & 'height' %in% names(save_plot)){
      ggsave(filename = paste(save_plot$filename, 'pdf', sep = '.'), plot = plot,# path = DIR_FIGURES,
             width = save_plot$width, height = save_plot$height, units = 'cm', dpi = 'print', limitsize = FALSE, device = cairo_pdf)
    }
    else {
      ggsave(filename = paste(save_plot$filename, 'pdf', sep = '.'), plot = plot,# path = DIR_FIGURES,
             scale = save_plot$scale, units = 'cm', dpi = 'print', device = cairo_pdf)
    }
    
    if (png_and_pdf) {
      cat('Saving image as png \n')
      if ('width' %in% names(save_plot) & 'height' %in% names(save_plot)) {
        ggsave(filename = paste(save_plot$filename, 'png', sep = '.'), plot = plot,# path = DIR_FIGURES,
               width = save_plot$width, height = save_plot$height, units = 'cm', dpi = 'print', limitsize = FALSE)
      }
      else {
        ggsave(filename = paste(save_plot$filename, 'png', sep = '.'), plot = plot,# path = DIR_FIGURES,
               scale = save_plot$scale, units = 'cm', dpi = 'print')
      }
    }
  }
}


geom_scatterpie_legend_custom <- function (radius, x, y, size = 10, face = 'plain', n = 5, labeller) 
{
  if (length(radius) > n) {
    radius <- unique(sapply(seq(min(radius), max(radius), 
                                length.out = n), signif, digits = 0))
  }
  label <- FALSE
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("labeller should be a function for converting radius")
    }
    label <- TRUE
  }
  dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x, 
                   y = y + radius - max(radius), maxr = max(radius))
  if (label) {
    dd$label <- labeller(dd$r)
  }
  else {
    dd$label <- dd$r
  }
  list(ggforce::geom_arc_bar(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, 
                                  start = ~start, end = ~end), data = dd, inherit.aes = FALSE), 
       ggplot2::geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + 
                                    r, yend = ~y + r), data = dd, inherit.aes = FALSE), 
       ggplot2::geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label), 
                          data = dd, hjust = "left", inherit.aes = FALSE, size = size, fontface = face))
}


### border for insets
element_grob.element_border <- function(element, ...)  {
  
  segmentsGrob(c(1,0),
               c(0,0),
               c(0,0),
               c(0,1), gp=gpar(lwd=1))
}
my_border <- function(...){
  structure(
    list(...), # this ... information is not used, btw
    class = c("element_border","element_blank", "element") # inheritance test workaround
  ) 
  
}


col_anomalies_ipcc_temperature <- function(n=23) {
  rev(grDevices::colorRampPalette(c(rgb(103,0,31,maxColorValue = 255),
                                    rgb(178,24,43,maxColorValue = 255),
                                    rgb(214,96,77,maxColorValue = 255),
                                    rgb(244,165,130,maxColorValue = 255),
                                    rgb(253,219,199,maxColorValue = 255),
                                    rgb(247,247,247,maxColorValue = 255),
                                    rgb(209,229,240,maxColorValue = 255),
                                    rgb(146,197,222,maxColorValue = 255),
                                    rgb(67,147,195,maxColorValue = 255),
                                    rgb(33,102,172,maxColorValue = 255),
                                    rgb(5,48,97,maxColorValue = 255)),space="rgb")(n))
}

