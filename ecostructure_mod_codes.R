# codes for running ecostructure
# these codes are modified from ecostructure package
# specifically designed to create more custom maps
# Jacob C. Cooper, April 2021

# adjusted ecos_plot_mod

ecos_plot_mod = function(omega = NULL,
                         coords = NULL,
                         bgmap_path = NULL,
                         adjust = FALSE,
                         thresh = 0.7,
                         long_lim = c(-180,180),
                         lat_lim = c(-60,90),
                         coastline_lwd = 10,
                         intensity = 1,
                         radius = 0.5,
                         sea=NULL,
                         lake=NULL,
                         country=NULL,
                         color = c("dodgerblue2","#E31A1C", "green4", "#6A3D9A","#FF7F00", "black","gold1","skyblue2","#FB9A99",
                                   "palegreen2", "#CAB2D6", "#FDBF6F",  "gray70", "khaki2", "maroon","orchid1","deeppink1",
                                   "blue1","steelblue4", "darkturquoise","green1","yellow4","yellow3", "darkorange4","brown",
                                   "red", "cornflowerblue", "cyan", "brown4", "burlywood", "darkgoldenrod1",
                                   "azure4", "green","deepskyblue","yellow", "azure1"),
                         pie_control = list(),
                         image_width = 1000,
                         image_height = 800,
                         path = "geostructure_plot.tiff"){
  
  require(sf)
  require(rnaturalearth)
  
  if(is.null(coords)){
    if(is.null(rownames(omega))){
      stop("coords not provided, omega rownames do not have latitude longitude
           information either")
    }
    latlong_chars <- rownames(omega)
    coords <- cbind.data.frame(
      as.numeric(sapply(latlong_chars, function(x) strsplit(x, "_")[[1]][1])),
      as.numeric(sapply(latlong_chars, function(x) strsplit(x, "_")[[1]][2])))
    colnames(coords) <- c("lat", "long")
  }else{
    if(dim(coords)[1] != dim(omega)[1]){
      stop("coords provided, but the number of rows in coords data does not
           match the number of rows in omega matrix")
    }
  }
  
  pie_control_default <- list(edges = 200, clockwise = TRUE, 
                              init.angle = 90, density = NULL, 
                              angle = 45, border = NULL,
                              lty = NULL, label.dist = 1.1)
  
  pie_control <- modifyList(pie_control_default, pie_control)
  
  if(is.null(sea)){
    print("Loading shapefile layers.")
    sea=ne_download(scale=110,type="ocean",category="physical")%>%
      st_as_sf()%>%st_geometry()
    lake=ne_download(scale=110,type="lakes",category="physical")%>%
      st_as_sf()%>%st_geometry()
    country=ne_download(scale=110,type="boundary_lines_land",
                        category="cultural")%>%
      st_as_sf()%>%st_geometry()
  }else{
    print("Maps preloaded.")
  }
  
  #glob <- c(xmin=long_lim[1], xmax=long_lim[2], ymin=lat_lim[1], ymax=lat_lim[2])
  #glob <- sf::st_bbox(glob)
  #glob <- structure(glob, crs = sf::st_crs(sea))
  #GlobalCoast <-suppressWarnings(suppressMessages(sf::st_intersection(GlobalCoast,
  #                                                                    sf::st_as_sfc(glob))))
  
  if(adjust){
    idx <- which(omega[,1] > thresh)
    omega <- omega[-idx,]
    coords <- coords[-idx,]
    omega <- omega[,-1]
    omega <- t(apply(omega, 1, function(x) return(x/sum(x))))
  }
  
  output_type <- strsplit(path, "[.]")[[1]][2]
  
  if(output_type == "tiff"){
    tiff(path, width = image_width, height = image_height)
  }else if(output_type == "png"){
    png(path, width = image_width, height = image_height)
  }else if(output_type == "pdf"){
    pdf(path, width = image_width, height = image_height)
  }else{
    stop("the output image may either be of  tiff, png or pdf extension")
  }
  
  plot(sea,col="#e1e1e1",axes=T,main="",lwd=coastline_lwd,
       xlim=long_lim,ylim=lat_lim)
  plot(country,add=T,lwd=0.5)
  plot(lake,col="#e1e1e1",axes=T,main="",lwd=coastline_lwd,add=T)
  
  par(lwd =.01)
  invisible(lapply(1:dim(omega)[1], function(r)
    do.call(mapplots::add.pie, append(list(
      z=as.integer(100*omega[r,]),
      x=coords[r,1], 
      y=coords[r,2], 
      labels=c("","",""),
      radius = radius,
      col=sapply(color, scales::alpha, intensity))
      , pie_control))))
  invisible(dev.off())
}

# adjusted eco_block

ecos_blocks_mod=function(omega,
                         filepath,
                         level,
                         ncluster,
                         blocker_metadata,
                         order_metadata,
                         palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                     "#0072B2", "#D55E00", "#CC79A7"),
                         structure_control = list(),
                         layout){
  
  if(!is.factor(blocker_metadata)){
    stop("the blocker_metadata must be a factor variable")
  }
  if(!is.numeric(order_metadata)){
    stop("the order_metadata must be a numeric variable")
  }
  num_levels <- length(levels(blocker_metadata))
  if(missing(layout)){
    layout <- c(floor(sqrt(num_levels)), ceiling(sqrt(num_levels)))
  }
  structure_control_default <- list(split_line=list(split_lwd = 1,
                                                    split_col = "white"),
                                    axis_tick = list(axis_ticks_length = .1,
                                                     axis_ticks_lwd_y = .1,
                                                     axis_ticks_lwd_x = .1,
                                                     axis_label_size = 6,
                                                     axis_label_face = "bold"),
                                    plot_labels = TRUE,
                                    levels_decreasing=FALSE,
                                    order_sample=TRUE,
                                    round_off=1,
                                    panel_title_size=10,
                                    panel_title_font=4,
                                    main_title="Block Structure Plot",
                                    yaxis_label = "Locality")
  
  structure_control <- modifyList(structure_control_default, structure_control)
  
  split_indices <- split(1:dim(omega)[1], as.factor(blocker_metadata))
  split_struct <- list()
  
  for(l in 1:length(split_indices)){
    
    order_split <- round(order_metadata[split_indices[[l]]], structure_control$round_off);
    omega_split <- omega[split_indices[[l]],]
    if(structure_control$levels_decreasing){
      order_split_ordered <- order_split[order(order_split, decreasing=TRUE)]
      omega_split_ordered <- omega_split[order(order_split, decreasing=TRUE),]
    }else{
      order_split_ordered <- order_split[order(order_split, decreasing=FALSE)]
      omega_split_ordered <- omega_split[order(order_split, decreasing=FALSE),]
    }
    annotation <- data.frame(
      sample_id = paste0("X", c(1:NROW(omega_split_ordered))),
      tissue_label = factor(order_split_ordered,
                            levels = unique(order_split_ordered) ) );
    
    split_struct[[l]] <- CountClust::StructureGGplot(omega = omega_split_ordered,
                                                     annotation = annotation,
                                                     figure_title = names(split_indices)[l],
                                                     palette = palette,
                                                     yaxis_label = structure_control$yaxis_label,
                                                     split_line=structure_control$split_line,
                                                     order_sample = structure_control$order_sample,
                                                     axis_tick = structure_control$axis_tick,
                                                     plot_labels=structure_control$plot_labels)
    
  }
  
  plot.struct=do.call(gridExtra::grid.arrange,
                      args = list(grobs=split_struct,
                                  ncol = layout[2],
                                  nrow = layout[1],
                                  top=grid::textGrob(structure_control$main_title,
                                                     gp=gpar(fontsize=structure_control$panel_title_size,
                                                             font=structure_control$panel_title_font))))
  
  print(plot.struct)
  ggsave(filename=paste0(filepath,level,
                         '_',ncluster,'_','blocks_plot.png'),
         plot=plot.struct,
         dpi=400)
}


# code for running over levels for multiple cluster sizes
# enables analyses within document
# level refers to taxonomic level

eco_africa=function(level,ncluster=NULL,data.x,
                    tolerance=NULL,n.trials=NULL,coords.x=NA,
                    sea=NULL,lake=NULL,country=NULL){
  if(is.null(ncluster)==T){ncluster=2}
  if(is.null(tolerance)==T){tolerance=0.1}
  if(is.null(n.trials)==T){n.trials=10}
  
  palette.x=c('#a6cee3','#1f78b4',
              '#b2df8a','#33a02c',
              '#fb9a99','#e5e5e5',
              '#e31a1c','#fdbf6f',
              '#ff7f00','#cab2d6',
              '#6a3d9a','#ffff99',
              '#b15928','#000000')
  
  x3=data.x%>%
    filter(Exclude!="Exclude")%>%
    select(-Clements,-Exclude)
  
  # case sensitive, finds column match
  xnames=x3[,which(colnames(x3)%flike%level)] %>% as.data.frame()
  
  x4=x3%>%select(-Genus,-Superspecies,-Species,
                 -Group,-Subspecies)%>%
    t()
  
  colnames(x4)=xnames[,1]
  
  u.names=unique(xnames)
  
  for(i in 1:nrow(u.names)){
    target=u.names[i,1]
    if(sum(colnames(x4)==target)<=1){
      index=which(colnames(x4)==target)
      nu.x=x4[,c(index,index)]
      if(i==1){
        x6=as.data.frame(nu.x)[,1]
      }else{
        x6=cbind(x6,nu.x[,1])
      }
    }
    if(sum(colnames(x4)==target)>1){
      xx=x4[,which(colnames(x4)==target)]
      nu.x=rowSums(xx)
      nu.x[nu.x>1]=1
      if(i==1){
        x6=as.data.frame(cbind(nu.x,nu.x))
        x6=x6[,1]
      }else{
        x6=cbind(x6,nu.x)
      }
    }
  }
  
  colnames(x6)=u.names[,1]
  
  fit=ecos_fit(x6,K=ncluster,
               tol=tolerance,num_trials=n.trials)
  
  ord.x=1:nrow(fit$omega)
  
  ecos_blocks_mod(fit$omega,blocker_metadata=as.factor('Afromontane'),
                  order_metadata = ord.x,
                  palette = palette.x,
                  filepath=filepath,
                  level=level,
                  ncluster=ncluster)
  
  # make maps
  
  features=CountClust::ExtractTopFeatures(fit$theta,
                                          top_features = 5,
                                          method="poisson",
                                          options="max")
  
  t(apply(features$indices,c(1,2),
          function(x){return(rownames(fit$theta)[x])}))
  
  # compare observed to null
  out=ecos_nullmodel(x6,K=ncluster,null.model = "richness",
                     iter_randomized = 10,option="BF")
  
  print(out)
  
  if(is.na(coords.x)==F){
    #ymin=min(coords.x$Latitude)+0.5
    #ymax=max(coords.x$Latitude)+0.5
    #xmin=min(coords.x$Longitude)+0.5
    #xmax=max(coords.x$Longitude)+0.5
    
    ecos_plot_mod(omega=fit$omega,
                  lat_lim=c(-42,14),
                  long_lim=c(-20,60),
                  coords=coords.x,
                  path=paste0(filepath,level,
                              '_',ncluster,'_','geostructure_plot.png'),
                  color = palette.x,
                  coastline_lwd = 2,
                  sea=sea,lake=lake,country=country)
  }
  
  if(is.na(coords.x)==F){
    #ymin=min(coords.x$Latitude)+0.5
    #ymax=max(coords.x$Latitude)+0.5
    #xmin=min(coords.x$Longitude)+0.5
    #xmax=max(coords.x$Longitude)+0.5
    
    ecos_plot_mod(omega=fit$omega,
                  lat_lim=c(-14,3),
                  long_lim=c(30,42),
                  coords=coords.x,
                  path=paste0(filepath,level,
                              '_',ncluster,'_','subset_geostructure_plot.png'),
                  color = palette.x,
                  coastline_lwd = 2,
                  sea=sea,lake=lake,country=country)
  }
}