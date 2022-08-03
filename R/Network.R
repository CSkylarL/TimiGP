##' Generate cell network files for Cytoscape
##' 
##' It generates three files that can be used to build network in Cytoscape:
##' 1. network files: simple interaction file (network.sif);
##' 2. node attributes (node.txt);
##' 3. edge attributes (edge.txt).
##' The function also returns a list of above files that can be modified in R.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param select a numeric vector of selected cell pairs according to "Index" column in resdata. 
##' Default selection is all statistically significantly enriched cell pairs(adjusted p value < 0.05).
##' @param dataset a value in one of 
##' c("Bindea2013","Bindea2013_Cancer","Charoentong2017", 
##' "Hughes2020_Tirosh2016","Zheng2021","Xu2018","Other"). 
##' The first four options include default group and color settings. 
##' If you use other dataset or want to change group and colors, 
##' please choose "Other".
##' @param group a vector of self-defined cell groups, 
##' whose names are the all favorable and unfavorable cells types in selection. 
##' @param geneset a data.frame of cell markers, 
##' in which the 1st is cell type, 
##' the 2nd column is the marker gene 
##' and the 3rd column is the name of the dataset(optional)
##' @param export a logical value. If TRUE, it will generate network files 
##' for Cytoscape analysis. The default value is TRUE.
##' @param path a directory to export the network files when "export = TRUE".
##' @return A list of network required files(network,node,edge)
##' @import dplyr
##' @export 
##' @examples
##' \dontrun{
##'   data("Bindea2013c_enrich")
##'   res <- Bindea2013c_enrich
##'   NET <- TimiCellNetwork(resdata = res,dataset = "Bindea2013_Cancer")
##' }
##' @author Chenyang Skylar Li

TimiCellNetwork<-  function(resdata = NULL,
                            select = NULL,
                            dataset = NULL,
                            group = NULL,
                            geneset = NULL,
                            export = TRUE,
                            path = NULL){
  # Examine required parameters
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  
  if (is.null(dataset)){
    stop('The parameter "dataset" is required. Please choose one of c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Hughes2020_Tirosh2016","Zheng2021","Xu2018","Other")')
  } else if(sum(dataset %in% c("Bindea2013","Bindea2013_Cancer","Charoentong2017","Hughes2020_Tirosh2016","Zheng2021","Xu2018","Other")) == 0){
    stop('Please choose one of c("Bindea2013","Bindea2013_Cancer","Charoentong2017", "Hughes2020_Tirosh2016","Zheng2021","Xu2018","Other")')
    
  }

  
  # selection
  if (is.null(select)) {
    se.r <- which(resdata$Adjust.P.Value < 0.05)
    
    if (length(se.r) < 5){
      stop('There are less than 5 significant cell pairs.')
    } else {
      message('Using all significant cell pairs(Adjested P.Value < 0.05)')
      resdata <- resdata[se.r, ]
    }
    
  } else {
    message('Using selected cell pairs')
    se.r <- which(resdata$Index %in%  select)
    resdata <- resdata[se.r, ]
  }
  
  # Generate network file
  net <- resdata[c("Favorable.Cell.Type" , "Unfavorable.Cell.Type")]
  colnames(net) <- c("Source","Target")
  net$Interaction <- "TimiGP"
  net <- net[c("Source","Interaction", "Target")]
  
  
  # Generate edge file
  edge <- resdata
  edge$Key <- paste0(edge$Favorable.Cell.Type, " (TimiGP) ", edge$Unfavorable.Cell.Type)
  edge <- edge[c(ncol(edge),1:(ncol(edge)-1))]
 
  
  # Generate node file
  selected <- unique(c(resdata$Unfavorable.Cell.Type, resdata$Favorable.Cell.Type))
  
  if (dataset == "Other"){
    if(is.null(group)){
      stop('Please define your groups of cell types')
      
    }
    if(is.null(geneset)){
      stop('Please input the geneset of cell markers')
      
    }
    cell <- selected
    
    if(sum(names(group) %in% cell) == 0) {
      stop('No Cell Type were found in given group')
    }
    
    if(sum(unique(geneset[1]) %in% cell) == 0) {
      stop('No Cell Type were found in given geneset')
    }
    
    att <- data.frame(row.names =cell,group=group)
    node <- table(geneset[2]) %>% 
      data.frame(row.names = 1) %>% merge(att,by=0) 
    names(node) <- c("Key", "No.Markers", "Group")
    
    
    
  } else{
    if(!is.null(group)){
      warning('Using default groups. If you have self-defined group, please set dataset = "Other"')
      
    }
    
    if(!is.null(geneset)){
      warning('Using default geneset. If you have self-defined geneset, please set dataset = "Other"')
      
    }
    
    if (dataset == "Bindea2013"){
      if(length(grep(selected,pattern = "Cancer cells")) == 1) {
        stop('You should set dataset = "Bindea2013_Cancer"')
      } else {
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
        
        
        geneset <- Ann_Immune
        geneset <- geneset[which(geneset[,3] == "Bindea2013"),]
        att <- data.frame(row.names =cell,group=group)
        node <- table(geneset$CellType) %>% 
          data.frame(row.names = 1) %>% merge(att,by=0)
        names(node) <- c("Key", "No.Markers", "Group")
        
        if(sum(selected %in% cell) == 0) {
          stop('No Cell Type were found in Bindea2013. Please choose "Other" and defined your groups')
        }
        
      }
      
    }
    
    
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
      
      geneset <- Ann_Bindea2013_Cancer
      att <- data.frame(row.names =cell,group=group)
      node <- table(geneset$CellType) %>% 
        data.frame(row.names = 1) %>% 
        merge(att,by=0)
      names(node) <- c("Key", "No.Markers", "Group")
      
      
      
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Bindea2013_Cancer. Please choose "Other" and defined your groups')
      }
    }
    
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
      
      geneset <- Ann_Immune
      geneset <- geneset[which(geneset[,3] == "Charoentong2017"),]
      att <- data.frame(row.names =cell,group=group)
      node <- table(geneset$CellType) %>% 
        data.frame(row.names = 1) %>% merge(att,by=0)
      names(node) <- c("Key", "No.Markers", "Group")
      
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Charoentong2017. Please choose "Other" and defined your groups')
      }
    }
    
    
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
      
      geneset <- Ann_Immune
      geneset <- geneset[which(geneset[,3] == "Xu2018"),]
      att <- data.frame(row.names =cell,group=group)
      node <- table(geneset$CellType) %>% 
        data.frame(row.names = 1) %>% merge(att,by=0)
      names(node) <- c("Key", "No.Markers", "Group")
      
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Xu2018. Please choose "Other" and defined your groups')
      }
      
    }
    
    if (dataset == "Hughes2020_Tirosh2016"){
      cell <- c("TI B"  ,"TI T", "TI Macrophage"  ,
                "Skin B"  ,"Skin Plasma"  ,"Skin T" ,
                "Skin Myeloid"  ,"Langerhans" ,"Skin Mast" ,
                "Endo","LEC"  ,
                "Fibroblast","Hair follicle"  ,"Keratinocyte",
                "Schwann" ,"Sebocyte"  ,
                "VSMC" , "Melanocyte",
                "Melanoma"  ,"CAF")
      
      group <- c(rep("TIL", 3),
                 rep("Skin Immune", 6),
                 rep("Skin Tissue",9),
                 rep("Cancer", 2))
      
      geneset <- Ann_melanoma_TME
      att <- data.frame(row.names =cell,group=group)
      node <- table(geneset$CellType) %>% 
        data.frame(row.names = 1) %>% merge(att,by=0)
      names(node) <- c("Key", "No.Markers", "Group")
      
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Hughes2020_Tirosh2016. Please choose "Other" and defined your groups')
      }
      
    }
    
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
      
      geneset <- Ann_Tcell
      att <- data.frame(row.names =cell,group=group)
      node <- table(geneset$CellType) %>% 
        data.frame(row.names = 1) %>% merge(att,by=0)
      names(node) <- c("Key", "No.Markers", "Group")
      
      if(sum(selected %in% cell) == 0) {
        stop('No Cell Type were found in Zheng2021. Please choose "Other" and defined your groups')
      }
      
    }
    
    se <- which(node$Key %in% selected)
    node <- node[se,]

  }
  if (export == TRUE) {
    
    if (is.null(path)){
      warning('No output directory was provided. Export to working directory')
      path <- getwd()
    }
    
    write.table(node,file = paste0(path,"/cell_node.txt"),quote = F,
                row.names = F,col.names = T,sep = "\t")
    write.table(net,file = paste0(path,"/cell_network.sif"),quote = F,
                row.names = F,col.names = F,sep = "\t")
    write.table(edge,file = paste0(path,"/cell_edge.txt"),quote = F,
                row.names = F,col.names = T,sep = "\t")
  } else {
    warning('No network file has been export. If you need these file for Cytoscape, please set "export=TRUE"')
  }
  
  
  combined.list <- list(net,node,edge)
  names(combined.list ) <- c("network","node","edge")
  return(combined.list)
}



##' Generate gene network files for Cytoscape
##' 
##' It generates three files that can be used to build network in Cytoscape:
##' 1. network files: simple interaction file (network.sif);
##' 2. node attributes (node.txt);
##' 3. edge attributes (edge.txt).
##' The function also returns a list of above files that can be modified in R.
##'
##' @param resdata TimiGP cox result generated from TimiCOX
##' @param select a character vector of selected gene pairs according to rownames in resdata. 
##' Default selection is all statistically significantly enriched cell pairs(adjusted p value < 0.05).
##' @param dataset a value in one of 
##' c("Bindea2013_Cancer","Immune3", "Hughes2020_Tirosh2016","Zheng2021","Other")
##' The first four options include default group and color settings. 
##' If you use other dataset or want to change group and colors, 
##' please choose "Other".
##' @param geneset a data.frame of cell markers, 
##' in which the 1st column is cell type, 
##' the 2nd column is the marker gene 
##' and the 3rd column is the name of the dataset(optional)
##' @param export a logical value. If TRUE, it will generate network files 
##' for Cytoscape analysis. The default value is TRUE.
##' @param path a directory to export the network files when "export = TRUE".
##' @return A list of network required files(network,node,edge)
##' @import dplyr
##' @export 
##' @examples
##' \dontrun{
##'   data(Bindea2013c_COX_MP_SKCM06)
##'   cox_res <- Bindea2013c_COX_MP_SKCM06
##'   NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Bindea2013_Cancer")
##' }
##' @author Chenyang Skylar Li

TimiGeneNetwork<-  function(resdata = NULL,
                            select = NULL,
                            dataset = NULL,
                            geneset = NULL,
                            export = TRUE,
                            path = NULL){
  # Examine required parameters
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  
  if (is.null(dataset)){
    stop('The parameter "dataset" is required. Please choose one of c("Bindea2013_Cancer","Immune3", "Hughes2020_Tirosh2016","Zheng2021","Other")')
  } else if(sum(dataset %in% c("Bindea2013_Cancer","Immune3", "Hughes2020_Tirosh2016","Zheng2021","Other")) == 0){
    stop('Please choose one of c("Bindea2013_Cancer","Immune3", "Hughes2020_Tirosh2016","Zheng2021","Other")')
    
  }
  
  # selection
  if (is.null(select)) {
    se.r <- which(resdata$QV < 0.05)
    
    if (length(se.r) < 5){
      stop('There are less than 5 significant gene pairs.')
    } else {
      message('Using all significant gene pairs(Adjested P.Value < 0.05)')
      resdata <- resdata[se.r, ]
    }
    
  } else {
    message('Using selected gene pairs')
    se.r <- which(rownames(resdata) %in%  select)
    resdata <- resdata[se.r, ]
  }
  
  # Generate network file
  xx <- row.names(resdata)
  tmp <- unlist(strsplit(xx, "_"))
  nn <- length(xx)
  resdata$Favorable.Gene <- tmp[(1:nn)*2-1]
  resdata$Unfavorable.Gene <- tmp[(1:nn)*2-0]
  
  
  net <- resdata[c("Favorable.Gene" , "Unfavorable.Gene")]
  colnames(net) <- c("Source","Target")
  net$Interaction <- "TimiGP"
  net <- net[c("Source","Interaction", "Target")]

  
  # Generate edge file
  edge <- resdata
  edge$Key <- paste0(edge$Favorable.Gene, " (TimiGP) ", edge$Unfavorable.Gene)
  edge <- edge[c(ncol(edge),1:(ncol(edge)-1))]

  
  # Generate node file
  selected <- unique(c(resdata$Favorable.Gene, resdata$Unfavorable.Gene))
  
  if (dataset == "Other"){
    
    if(is.null(geneset)){
      stop('Please input the geneset of cell markers')
      
    }
    
    
    if(sum(unique(geneset[,2]) %in% selected) == 0) {
      stop('No Gene were found in given geneset')
    } else{
      se <- which(geneset[,2] %in% selected)
      node <- geneset[se,]
    }
    
    
    
  } else{
    
    if(!is.null(geneset)){
      warning('Using default geneset. If you have self-defined geneset, please set dataset = "Other"')
      
    }
    
    
    
    
    if (dataset == "Bindea2013_Cancer"){
      
      geneset <- Ann_Bindea2013_Cancer
      if(sum(selected %in% geneset[,2] ) == 0) {
        stop('No Gene were found in Bindea2013_Cancer. Please choose "Other" and use your geneset')
      } else {
        se <- which(geneset[,2] %in% selected)
        node <- geneset[se,]
      }
      
    }
    
    
    if (dataset == "Immune3"){
      
      geneset <- Ann_Immune
      if(sum(selected %in% geneset[,2] ) == 0) {
        stop('No Gene were found in Immune3(Charoentong2017_Bindea2013_Xu2018_Immune). Please choose "Other" and use your geneset')
      } else {
        se <- which(geneset[,2] %in% selected)
        node <- geneset[se,]
      }
      
      
      
    }
    
    if (dataset == "Hughes2020_Tirosh2016"){
      
      geneset <- Ann_melanoma_TME
      if(sum(selected %in% geneset[,2] ) == 0) {
        stop('No Gene were found in Hughes2020_Tirosh2016. Please choose "Other" and use your geneset')
      } else {
        se <- which(geneset[,2] %in% selected)
        node <- geneset[se,]
      }
      
    }
    
    if (dataset == "Zheng2021"){
      
      geneset <- Ann_Tcell
      if(sum(selected %in% geneset[,2] ) == 0) {
        stop('No Gene were found in Hughes2020_Tirosh2016. Please choose "Other" and use your geneset')
      } else {
        se <- which(geneset[,2] %in% selected)
        node <- geneset[se,]
      }
      
    }
    
  }
  
  if (export == TRUE) {
    
    if (is.null(path)){
      warning('No output directory was provided. Export to working directory')
      path <- getwd()
    }
    
    write.table(node,file = paste0(path,"/gene_node.txt"),quote = F,
                row.names = F,col.names = T,sep = "\t")
    write.table(net,file = paste0(path,"/gene_network.sif"),quote = F,
                row.names = F,col.names = F,sep = "\t")
    write.table(edge,file = paste0(path,"/gene_edge.txt"),quote = F,
                row.names = F,col.names = T,sep = "\t")
  } else {
    warning('No network file has been export. If you need these file for Cytoscape, please set "export=TRUE"')
  }

  
  combined.list <- list(net,node,edge)
  names(combined.list ) <- c("network","node","edge")
  return(combined.list)
}
