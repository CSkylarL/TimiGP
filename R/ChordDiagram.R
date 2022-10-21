##' Chord diagram of Cell Interaction
##' 
##' Chord diagram reveals the functional interaction network.
##' The arrow points from favorable cell type A to unfavorable cell type B, 
##' which denotes that function of A greater than B associated with a good prognosis.
##' The width of the arrow represents -log10(Adjust.P.Value),namely,
##' the wider the arrow is, the smaller the adjusted p-value is.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param select a numeric vector of selected cell interactions according to "Index" column in resdata. 
##' Default selection is all functional interactions(adjusted p value < 0.05).
##' @param dataset a value in one of 
##' c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Xu2018"
##' ,"Newman2015", 
##' "Tirosh2016","Zheng2021","Other"). 
##' The first four options include default group and color settings. 
##' If you use other dataset or want to change group and colors, 
##' please choose "Other".
##' @param group a vector of self-defined cell groups, 
##' whose names are the all favorable and unfavorable cells types in selection. 
##' @param color a vector of self-defined cell colors, 
##' whose names are the all favorable and unfavorable cells types in selection. 
##' @return Chord diagram of Cell Interaction
##' @import circlize
##' @import RColorBrewer
##' @export 
##' @examples
##' \dontrun{
##'   data("Bindea2013c_enrich")
##'   res <- Bindea2013c_enrich
##'   TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer")
##'   TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer",select = 1:10)
##'  
##' }
##' @author Chenyang Skylar Li

TimiCellChord<-  function(resdata = NULL,
                        select = NULL,
                        dataset = NULL,
                        group = NULL,
                        color = NULL){
   # Examine required parameters ###############################################
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (is.null(dataset)){
    stop('The parameter "dataset" is required. Please choose one of c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Xu2018","Newman2015",  "Tirosh2016","Zheng2021","Other")')
  } else if(sum(dataset %in% c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Xu2018","Newman2015",  "Tirosh2016","Zheng2021","Other")) == 0){
    stop('Please choose one of c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Xu2018","Newman2015",  "Tirosh2016","Zheng2021","Other")')
    
  }
   # -log10(adjust.P.Value) ####################################################
  
  resdata$rev.p.adj <- -log10(resdata$Adjust.P.Value+1e-300)
  
  # selection ##################################################################
  if (is.null(select)) {
    se.c <- c("Favorable.Cell.Type", "Unfavorable.Cell.Type", "rev.p.adj")
    se.r <- which(resdata$Adjust.P.Value < 0.05)
    
    if (length(se.r) < 5){
      stop('There are less than 5 functional interactions.')
    } else {
      message('Using all functional interactions(Adjested P.Value < 0.05)')
      p.data <- resdata[se.r, se.c]
    }
    
  } else {
    message('Using selected cell interactions')
    se.c <- c("Favorable.Cell.Type", "Unfavorable.Cell.Type", "rev.p.adj")
    se.r <- which(resdata$Index %in%  select)
    p.data <- resdata[se.r, se.c]
  }
    
  # Default group and color setting for internal cell type annotations #########  
  selected <- unique(c(p.data$Unfavorable.Cell.Type, p.data$Favorable.Cell.Type))
  
  if (dataset == "Other"){
    if(is.null(group)){
      stop('Please define your groups of cell types')
    
    }
    if(is.null(color)){
      stop('Please define your colors of cell types')
      
    }
    cell <- selected
    
    if(sum(names(group) %in% cell) == 0) {
      stop('No Cell Type were found in given group')
    }
    
    if(sum(names(color)%in% cell) == 0) {
      stop('No Cell Type were found in given color')
    }
    
  } else{
      if(!is.null(group)){
        warning('Using default groups. If you have self-defined group, please set dataset = "Other"')
        
      }
      if(!is.null(color)){
        warning('Using default colors. If you have self-defined color, please set dataset = "Other"')
      }
    # "Bindea2013" =============================================================
      if (dataset == "Bindea2013"){
        if(length(grep(selected,pattern = "Tumor",ignore.case=TRUE)) == 1) {
          stop('You should set dataset = "Bindea2013_Cancer"')
        }
        
        cell <- c("B", 
                  "T", "CD8 T" ,"Th", "Th1", "Th2" , 
                  "Tfh","Tem", "Tcm","Tgd",
                  "Cytotoxic", 
                  "NK" , "CD56dim NK","CD56bright NK",
                  "DC" , "iDC", "aDC" ,
                  "Macrophage", 
                  "Neutrophil","Mast","Eosinophil")
        group <- c(rep("B Cell", 1),
                   rep("T Cell", 9),
                   rep("Cytotoxic Cell", 1),
                   rep("NK Cell", 3),
                   rep("DC", 3),
                   rep("Mononuclear phagocyte system", 1),
                   rep("Granulocytes", 3))
        
        names(group) <- cell
        color <- c("#66c2a4",
                   "#dadaeb","#fcc5c0",
                   "#bcbddc","#fa9fb5",
                   "#9e9ac8","#f768a1",
                   "#6a51a3","#ae017e",
                   "#54278f","#7a0177",
                   "#4eb3d3",
                   "#2b8cbe",
                   "#0868ac",
                   "#1d91c0",
                   "#225ea8",
                   "#253494",
                   "#fed976",
                   "#fc4e2a",
                   "#e31a1c",
                   "#bd0026")
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Bindea2013. Please choose "Other" and defined your groups and colors ')
        }
        
      }
      
    # "Bindea2013_Cancer" ======================================================
    
      if (dataset == "Bindea2013_Cancer"){
        
        cell <- c("B", 
                  "T", "CD8 T" ,"Th", "Th1", "Th2" , 
                  "Tfh","Tem", "Tcm","Tgd",
                  "Cytotoxic", 
                  "NK" , "CD56dim NK","CD56bright NK",
                  "DC" , "iDC", "aDC" ,
                  "Macrophage", 
                  "Neutrophil","Mast","Eosinophil",
                  "Tumor")
        group <- c(rep("B Cell", 1),
                   rep("T Cell", 9),
                   rep("Cytotoxic Cell", 1),
                   rep("NK Cell", 3),
                   rep("DC", 3),
                   rep("Mononuclear phagocyte system", 1),
                   rep("Granulocytes", 3), 
                   rep("Tumor", 1))
        
        names(group) <- cell
        color <- c("#66c2a4",
                   "#dadaeb","#fcc5c0",
                   "#bcbddc","#fa9fb5",
                   "#9e9ac8","#f768a1",
                   "#6a51a3","#ae017e",
                   "#54278f","#7a0177",
                   "#4eb3d3",
                   "#2b8cbe",
                   "#0868ac",
                   "#1d91c0",
                   "#225ea8",
                   "#253494",
                   "#fed976",
                   "#fc4e2a",
                   "#e31a1c",
                   "#bd0026",
                   "#525252")
        
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Bindea2013_Cancer. Please choose "Other" and defined your groups and colors ')
        }
      }
      # "Charoentong2017" ======================================================
      if (dataset == "Charoentong2017"){
        cell <- c("iB", "aB" , "Memory B",
                  "Th1" ,"Th2", "Th17" ,
                  "Tfh", "Treg" ,
                  "aCD4 T","aCD8 T",
                  "CD4 Tem" , "CD8 Tem",
                  "CD4 Tcm" , "CD8 Tcm",
                  "Tgd", "NKT" ,
                  "NK","CD56bright NK","CD56dim NK" ,
                  "iDC" ,"aDC" ,  "pDC",
                  "Macrophage", "Monocyte",
                  "Neutrophil","Mast", "Eosinophil" ,
                  "MDSC")
        
        group <- c(rep("B Cell", 3),
                   rep("T Cell", 13),
                   rep("NK Cell", 3),
                   rep("DC", 3),
                   rep("Mononuclear phagocyte system", 2),
                   rep("Granulocytes", 3),
                   rep("MDSC", 1))
        
        names(group) <- cell
        color <- c("#66c2a4",
                   "#2ca25f",
                   "#006d2c",
                   "#dadaeb","#fcc5c0",
                   "#bcbddc","#fa9fb5",
                   "#9e9ac8","#f768a1",
                   "#807dba","#dd3497",
                   "#6a51a3","#ae017e",
                   "#54278f","#7a0177","#3f007d",
                   "#4eb3d3",
                   "#2b8cbe",
                   "#0868ac",
                   "#1d91c0",
                   "#225ea8",
                   "#253494",
                   "#fed976",
                   "#feb24c",
                   "#fc4e2a",
                   "#e31a1c",
                   "#bd0026",
                   "#7f0000")
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Charoentong2017. Please choose "Other" and defined your groups and colors ')
        }
      }
      
    # "Xu2018" =================================================================
      if (dataset == "Xu2018"){
        cell <- c("B",
                  "T" , "CD4 T", "CD8 T",
                  "Th1" ,"Th17" , "Th2" ,"Th22" ,"Treg", 
                  "NK",
                  "DC",
                  "Macrophage" , "Monocyte",
                  "Neutrophil", "Basophil", "Eosinophil" , 
                  "MDSC")
        
        group <- c(rep("B Cell", 1),
                   rep("T Cell", 8),
                   rep("NK Cell",1),
                   rep("DC", 1),
                   rep("Mononuclear phagocyte system", 2),
                   rep("Granulocytes", 3),
                   rep("MDSC", 1))
        
        names(group) <- cell
        
        color <- c("#66c2a4",
                   "#dadaeb","#fcc5c0",
                   "#bcbddc","#fa9fb5",
                   "#9e9ac8","#f768a1",
                   "#807dba","#dd3497",
                   "#4eb3d3",
                   "#1d91c0",
                   "#fed976",
                   "#feb24c",
                   "#fc4e2a",
                   "#e31a1c",
                   "#bd0026",
                   "#7f0000")
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Xu2018. Please choose "Other" and defined your groups and colors ')
        }
        
      }
      
    # "Newman2015" ======================================================
    
    if (dataset == "Newman2015"){
      
      cell <- c("Bn", "Bm", "Plasma" ,
                "CD8 T" ,"CD4 Tn","Tfh",  "Treg", "rCD4 Tm", "aCD4 Tm" , "Tgd" ,
                "rNK", "aNK",
                "rDC", "aDC",
                "M0",   "M1",   "M2" , "Monocyte",
                "Neutrophil", "rMast" ,"aMast", "Eosinophil" )
      group <- c(rep("B Cell", 3),
                 rep("T Cell", 7),
                 rep("NK Cell", 2),
                 rep("DC", 2),
                 rep("Mononuclear phagocyte system", 4),
                 rep("Granulocytes", 4))
      
      
      names(group) <- cell
      color <- c("#66c2a4",
                 "#2ca25f",
                 "#006d2c",
                 "#dadaeb","#fcc5c0",
                 "#bcbddc","#fa9fb5",
                 "#9e9ac8","#f768a1",
                 "#807dba",
                 "#4eb3d3",
                 "#2b8cbe",
                 "#1d91c0",
                 "#225ea8",
                 "#fed976",
                 "#feb24c",
                 "#ec7014",
                 "#cc4c02",
                 "#fc4e2a",
                 "#e31a1c",
                 "#bd0026",
                 "#800026")
      
      names(color) <- cell
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Newman2015. Please choose "Other" and defined your groups and colors ')
      }
    }
    
   
    # "Tirosh2016" ===========================================================
      if (dataset == "Tirosh2016"){
        cell <- c("B"  ,"T", "Macrophage"  ,
                  "Endo","CAF",
                  "Melanoma"  )
        
        group <- c(rep("TIL", 3),
                   rep("Stroma", 2),
                   rep("Cancer", 1))
        
        names(group) <- cell
        
        color <- c("#9e9ac8",
                   "#756bb1",
                   "#54278f",
                   
                   "#fdae6b", "#fcc5c0",

                   "#969696")
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Tirosh2016. Please choose "Other" and defined your groups and colors ')
        }
        
      }
      # "Zheng2021" ============================================================
      if (dataset == "Zheng2021"){
        cell <- c("CD4+Tn",              "CD4+CXCR5+ pre-Tfh" ,
                  "CD4+ADSL+ Tn"  ,      "CD4+IL7R- Tn",      "CD4+TNF+ T",
                  "CD4+AREG+ Tm"  ,      "CD4+TIMP1+ Tm"  ,   "CD4+CREM+ Tm" ,
                  "CD4+CCL5+ Tm"  ,      "CD4+CAPG+ Tm",      "CD4+CAPG+CREM- Tm"  ,
                  "CD4+GZMK+ Tem" ,      "CD4+Temra",         "CD4+CCR6+ Th17"  ,
                  "CD4+IL26+ Th17",      "CD4+IL21+ Tfh"  ,   "CD4+IFNG+ Tfh/Th1"  ,
                  "CD4+TNFRSF9- Treg",   "CD4+S1PR1+ Treg",   "CD4+TNFRSF9+ Treg"  ,
                  "CD4+ISG+ Treg" ,      "CD4+ISG+ Th" ,      "CD4+NME1+CCR4- T",
                  "CD4+NME1+CCR4+ T" ,
                  "CD8+Tn"  ,            "CD8+IL7R+ Tm",      "CD8+ZNF683+CXCR6- Tm"  ,
                  "CD8+GZMK+ early Tem", "CD8+GZMK+ Tem"  ,   "CD8+Temra" ,
                  "CD8+KIR+EOMES+ NK-like" , "CD8+KIR+TXK+ NK-like" ,  "CD8+ZNF683+CXCR6+ Trm" ,
                  "CD8+GZMK+ Tex" , "CD8+terminal Tex"  ,     "CD8+OXPHOS- Tex" ,
                  "CD8+TCF7+ Tex" , "CD8+ISG+ CD8+ T",        "CD8+Tc17"  ,
                  "CD8+NME1+ T")
        
        group <- c(rep("CD4 T Cell", 24),
                   rep("CD8 T Cell", 16))
        
        names(group) <- cell
        
        color <- c("#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b",
                   "#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d",
                   "#e5f5e0","#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#006d2c","#00441b",
                   
                   "#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704",
                   "#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f")
        names(color) <- cell
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Zheng2021. Please choose "Other" and defined your groups and colors ')
        }
        
      }
    
    se <- which(cell %in% selected)
    cell <- cell[se]
    group <- group[se]
    color <-color[se]
      
  }
  
  

# plot
  
  circos.clear()
  circos.par(start.degree = 180, clock.wise = T, cell.padding = c(0, 0, 0, 0))
  chordDiagram(p.data, 
               order = names(group),
               direction.type = c("diffHeight", "arrows"),
               diffHeight = mm_h(4),
               directional = 1,
               link.arr.type = "big.arrow",
               preAllocateTracks = 1,
               # annotationTrack =  c("name", "grid"),
               annotationTrack =  "grid",
               annotationTrackHeight = c(0.04, mm_h(5)),
               group = factor(group, 
                              levels =  unique(group)),
               big.gap = 5,
               small.gap = 0.3,
               transparency = 0.4,
               #             scale = T,
               grid.col = color,
               #            col = arrow_col,
               link.sort = TRUE, 
               link.decreasing = TRUE,
               link.zindex = rank(p.data$rev.p.adj))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(-0.01, 0.5), 
                cex = 1.5)
  }, bg.border = NA)
  circos.clear()
}





##' Chord diagram of Gene Interaction
##' 
##' Chord diagram reveals the gene interactions of a selected cell interaction.
##' The arrow points from enriched marker F of favorable cell type A 
##' to enriched marker U of unfavorable cell type B, 
##' which denotes that 
##' the expression of F greater than that of U is associated with a good prognosis.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param select a numeric value of selected cell interactions according to "Index" column. 
##' Default selection is the top 1 cell interaction.
##' @param color a vector of self-defined marker colors, 
##' whose names are the enriched markers
##' (not pair,you need to split the genes in Shared.IMGP column) 
##' in selected cell interaction. 
##' @return Chord diagram of Gene Interaction
##' @import circlize
##' @import RColorBrewer
##' @export 
##' @examples
##' \dontrun{
##'   data("Bindea2013c_enrich")
##'   res <- Bindea2013c_enrich
##'   TimiGeneChord(resdata = res)
##'   TimiGeneChord(resdata = res,select = 2)
##'  
##' }
##' @author Chenyang Skylar Li

TimiGeneChord<-  function(resdata = NULL,
                          select = 1,
                          color = NULL){
  # Examine required parameters
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  
  # selction
  if (length(select) >1){
    stop('Please only choose 1 cell interaction at one time')
  } else if (select > max(resdata$Index)){
    stop('There are only ',max(resdata$Index),' cell interactions')
  } else {
    se.idx <- which(resdata$Index == select)
    
    xx <- unlist(strsplit(resdata$Shared.IMGP[se.idx ], split = "/"))
    xx <- unlist(strsplit(xx, split = "_"))
    favorable.gene <- xx[seq(1,length(xx),2)]
    unfavorable.gene <- xx[seq(2,length(xx),2)]
    p.data <- data.frame(favorable.gene,unfavorable.gene)
    colnames(p.data) <- resdata[se.idx,c("Favorable.Cell.Type", "Unfavorable.Cell.Type")]
  }
  
  # Color  
  if (is.null(color)){
    
    message('Using default color')
    set.seed(1234)
    red <- rand_color(floor(length(unique(p.data[,1]))/2), hue = "red", luminosity = "light")
    pink <- rand_color(ceiling(length(unique(p.data[,1]))/2), hue = "pink", luminosity = "bright")
    fav.col <- c(red,pink)
    names(fav.col) <- unique(p.data[,1])
    
    blue <- rand_color(floor(length(unique(p.data[,2]))/2), hue = "blue", luminosity = "bright")
    purple <- rand_color(ceiling(length(unique(p.data[,2]))/2), hue = "purple", luminosity = "bright")
    unfav.col <- c(blue,purple)
    names(unfav.col) <- unique(p.data[,2])
    
    color <- c(fav.col,unfav.col)
  } else if (length(color) != length(unique(c(p.data[,1],p.data[,2])))) {
    stop('The number of colors are different of unique markers enriched in the given cell interaction')
  } else {
    message('Using self-defined color')
  }
  # plot
  PV.label <- round(resdata$Adjust.P.Value,3)
  PV.label <- ifelse(PV.label < 0.001, "< 0.001",paste0("= ", PV.label))
  
  circos.clear()
  circos.par(start.degree = 180, clock.wise = T,
             cell.padding = c(0, 0, 0, 0))
  chordDiagram(p.data, 
               direction.type = c("diffHeight", "arrows"),
               diffHeight = mm_h(4),
               directional = 1,
               link.arr.type = "big.arrow",
               preAllocateTracks = 1,
               # annotationTrack =  c("name", "grid"),
               annotationTrack =  "grid",
               annotationTrackHeight = c(0.04, mm_h(5)),
               big.gap = 10,
               transparency = 0.4,
               #             scale = T,
               grid.col = color)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(-0.01, 0.5), 
                cex = 0.8)
  }, bg.border = NA)
  
  title(
    paste0(resdata$Cell.Interaction[se.idx],
           "(Adjusted P-Value ",PV.label[se.idx],")")
  )
  
  #circos.info()
  circos.clear()
}

