##' Calulate Favorability Score of cell types
##' 
##' It calculates the favorability score of each cell type 
##' based on TimiGP cell interaction network.
##' In the network, 
##' the out-degree of cell A means how many times high  A-to-other cell ratio 
##' is associated with favorable prognosis.
##' The in-degree of cell A means how many times high A-to-other cell ratio 
##' is  associated with unfavorable prognosis.
##' The favorability score  includes favorable score and unfavorable score.
##' Favorable score: out-degree of the cell/sum of out-degree of all cell*100.
##' Unfavorable score: in-degree of the cell/sum of in-degree of all cell*100.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param cutoff A cutoff of adjusted P-value used to filter cell interactions. 
##' The default cutoff is 0.05.
##' @return A matrix of Favorability Score 
##' which includes favorable score and unfavorable score.
##' Difference = Favorable.Score - Unfavorable.Score. 
##' If Difference greater than 0, the cell type is classified as F(avorable) cell. 
##' If Difference is negative, the cell type is classified as F(avorable) cell.
##' @import dplyr
##' @export 
##' @examples
##' \dontrun{
##'   data("Galon2013c_enrich")
##'   res <- Galon2013c_enrich
##'   score <- TimiFS(res)
##'   head(score)
##' }
##' @author Chenyang Skylar Li

TimiFS<-  function(resdata = NULL,
                   cutoff=0.05){
  # Examine required parameters
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  resdata <- resdata %>% filter(Adjust.P.Value < cutoff)
  # favorable cells
  total_ct <- unique(c(resdata$Favorable.Cell.Type,
                       resdata$Unfavorable.Cell.Type))

  faCell <- table(resdata$Favorable.Cell.Type) %>% data.frame() 
  se <- which(! total_ct %in% faCell$Var1)
  faCell <- faCell %>% rbind(data.frame(Var1=total_ct[se],Freq=rep(0,length(se)))) %>%
    mutate(F=Freq/sum(Freq) *100)
  
  dim(faCell)
  
  # unfavorable cell
  unfaCell <- table(resdata$Unfavorable.Cell.Type) %>% data.frame()
  se <- which(! total_ct %in% unfaCell$Var1)
  unfaCell <- unfaCell %>% rbind(data.frame(Var1=total_ct[se],Freq=rep(0,length(se)))) %>%
    mutate(U=Freq/sum(Freq) *100)
  dim(unfaCell)
  
  # all cell 
  score <- merge(faCell,unfaCell,by="Var1") %>% data.frame() %>%
    mutate(diff=F-U) %>%
    arrange(-diff) %>%
    mutate(Var1 = factor(Var1, levels=Var1)) %>%
    filter(diff != 0) %>%
    mutate( group = ifelse(diff>0, "F","U"))%>% 
    mutate(Rank = rank(-diff,ties.method = "min"))
  score$Index <- 1:nrow(score)
  score <- score[c(9,8,1:7)]
  names(score) <-c("Index", "Rank","Cell.Type","OutDegree","Favorable.Score",
                   "Indegree","Unfavorable.Score",
                   "Difference",
                   "Group")

  
  return(score)
}





##' Visualization of favorability score
##' 
##' It generate a bar plot of the favorability score of cell type 
##' to evaluate its favorable(orange,positive) or unfavorable(Blue,negative) role in prognosis.
##'
##' @param score Favorability score calculated from TimiFS
##' @param select a numeric vector of selected cell pairs according to "Index" column in resdata. 
##' Default selection is all cell type.
##' @return A barplot of Favorability Score 
##' @import dplyr
##' @export 
##' @examples
##' \dontrun{
##'   data("Galon2013c_enrich")
##'   res <- Galon2013c_enrich
##'   score <- TimiFS(res)
##'   head(score)
##'   p <- TimiFSBar(score)
##'   p
##' }
##' @author Chenyang Skylar Li

TimiFSBar<-  function(score = NULL,
                      select=NULL){
  # Examine required parameters
  if (is.null(score)){
    stop('The parameter "score" is required.')
  }
  
  if (!is.null(select)){
    se <- which(score$Index %in% select)
    if(length(se) == 0|length(se) > nrow(score)) {
      stop('There are only ',nrow(score), 'cell types')
    } else {
      score <- score[se,]
    }
    
  }
  
  
  # plot
  max <- ceiling(max(score$Favorable.Score)/10)*10
  min <- ceiling(max(score$Unfavorable.Score)/10)*10
  score$Cell.Type <- factor(score$Cell.Type,levels = score$Cell.Type)
  
  p <- ggplot() +
    geom_bar(data = score,
             mapping = aes(x = Cell.Type,y=Favorable.Score),
             width=0.8, position = position_dodge(width=0.01),
             stat="identity" ,alpha=0.5,fill="#fc4e2a") +
    geom_bar(data = score,
             mapping = aes(x = Cell.Type,y=-Unfavorable.Score ),
             width=0.8, position = position_dodge(width=0.01),
             stat="identity" ,alpha=0.5,fill="#253494") +
    theme_bw(base_size = 20, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Cell Type",y="Favorability Score",
         title=paste0("Favorability"),
         base_size = 22, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          plot.margin = unit(c(1, 1, 1, 5), "lines"),
          legend.position="right",
          legend.box = "vertical",
          legend.direction= "vertical",
          panel.grid=element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.title = element_text(face="bold", color="black",family = "serif", size=10),
          legend.text= element_text(face="bold", color="black",family = "serif", size=10),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black", size=13, 
                                     angle = 60,  vjust = 1, hjust=1),
          axis.text.y = element_text(face="bold", color="black", size=15),
          axis.title.x = element_text(face="bold", color="black", size=17),
          axis.title.y = element_text(face="bold",color="black", size=17),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))) +
    scale_y_continuous(limits=c(-min, max), breaks=seq(-min,max, 10),expand = c(0,0)) +
    geom_hline(yintercept = 0,lty=1,size=1) +
    annotate("text", x = -1.5, y = max, label = "Favorable\nScore")+
    annotate("text", x = -1.5, y = -min, label = "Unfavorable\nScore")+
    coord_cartesian(xlim = c(1,nrow(score)),
                    clip = 'off')
  return(p)
}

