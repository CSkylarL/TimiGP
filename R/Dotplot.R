##' Dotplot
##' 
##' A dotplot shows the results of selected cell pair enrichment.
##' The x-axis shows selected cell Pairs. The y-axis shows Enrichment Ratio.
##' The color of dots represents adjusted P-value. 
##' The size of dots represents the number of marker pairs 
##' shared by query pairs and annotation pairs.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param select a numeric vector of selected cell pairs according to "Index" column in resdata. 
##' Default selection is top 5 enrichment
##' @return Enrichment dotpot 
##' @import ggplot2
##' @export 
##' @examples
##' \dontrun{
##'   data("Galon2013c_enrich")
##'   res <- Galon2013c_enrich
##'   p1<-TimiDotplot(res)
##'   p1
##'   p2<-TimiDotplot(res,select = c(1:10))
##'   p2
##'  
##' }
##' @author Chenyang Skylar Li

TimiDotplot<-  function(resdata = NULL,
                     select = 1:5){
  
  # Examine required parameters
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  # plot
  resdata$Cell.Pair <- factor(resdata$Cell.Pair, levels=resdata$Cell.Pair)
  se <- which(resdata$Index %in% select)
  max <- max(ceiling(resdata$Enrichment.Ratio[se]))+1
  p<- ggplot(resdata[se,] ) +
      geom_point(
        mapping = aes(x=Cell.Pair, y=Enrichment.Ratio, color=Adjust.P.Value, 
                    size=No.Shared.Marker.Pair)) +
    scale_color_continuous(low="#9970ab", high="#5aae61", 
                           limits = c(0, 0.05),
                           name = "Adjust P-Value",
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range = c(3,8), name="No. Shared\nMarker Pair") +
    geom_hline(yintercept=0, lty=4, col="black", lwd=1)+
    
    theme_bw(base_size = 20, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Cell Pair",y="Enrichment Ratio",
         title="Cell Pair Enrichment",
         base_size = 12, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
          legend.position="right",
          legend.box = "vertical",
          legend.direction= "vertical",
          plot.margin = unit(c(1, 1, 1, 5), "lines"),
          panel.grid=element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.title = element_text(face="bold", color="black",family = "serif", size=10),
          legend.text= element_text(face="bold", color="black",family = "serif", size=10),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black", size=10, 
                                     angle = 40,  vjust = 1, hjust=1),
          axis.text.y = element_text(face="bold", color="black", size=13),
          axis.title.x = element_text(face="bold", color="black", size=15),
          axis.title.y = element_text(face="bold",color="black", size=15),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last"))) +
    scale_y_continuous(limits=c(0, max), breaks=seq(0,max, 2),expand = c(0,0)) +
    guides(shape = "none") +
    annotate("text", x = -1, y = max, label = "Good\nPrognosis")+
    coord_cartesian(xlim = c(1,nrow(resdata[se,] )),
                    clip = 'off')
  
  return(p)
}
