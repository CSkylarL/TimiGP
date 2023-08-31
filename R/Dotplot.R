##' Dotplot
##' 
##' A dotplot shows the results of selected cell interaction enrichment.
##' The x-axis shows selected cell interactions. The y-axis shows Enrichment Ratio.
##' The color of dots represents FDR. 
##' The size of dots represents the number of marker pairs 
##' shared by query pairs and annotation pairs.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param select a numeric vector of selected cell interactions according to "Index" column in resdata. 
##' Default selection is top 5 enrichment.
##' @param condition A value in one of 
##' c("P.Value","Adjust.P.Value","Permutation.FDR"),
##' which is the column name of resdata and will be represented by color for the dot plot.
##' The default value is the "Adjust.P.Value".
##' @param cutoff The maximum value for the color bar.
##' The default value is 0.05, which is used for the default condition "Adjust.P.Value". 
##' For other conditions, 0.01 is recommended for "P.Value",
##' and 0.2 is recommended for "Permutation.FDR".
##' @return Enrichment dotpot 
##' @import ggplot2
##' @export 
##' @examples
##' \dontrun{
##'   data("Bindea2013c_enrich")
##'   res <- Bindea2013c_enrich
##'   p1<-TimiDotplot(res)
##'   p1
##'   p2<-TimiDotplot(res,select = c(1:10))
##'   p2
##'  
##' }
##' @author Chenyang Skylar Li

TimiDotplot<-  function(resdata = NULL,
                select = 1:5,
                condition = "Adjust.P.Value",
                cutoff = 0.05){
  
  # Examine required parameters-------------------------------------------------
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }

  if (is.null(condition)){
    stop('The parameter "condition" is required. Please choose one of c("P.Value","Adjust.P.Value","Permutation.FDR")')
  } else if(sum(condition %in% c("P.Value","Adjust.P.Value","Permutation.FDR")) !=1){
    stop('Please choose one of c("P.Value","Adjust.P.Value","Permutation.FDR")')
  } else {
    message("The color for dotplot is ", condition)
  }

  if (is.null(cutoff)){
    stop('The parameter "cutoff" is required. Please set it like 0.05, 0.01')
  } else {
    message("The color bar limits for dotplot is ", cutoff)
  }
  # plot
  resdata$Cell.Interaction <- factor(resdata$Cell.Interaction, 
                                     levels=resdata$Cell.Interaction)
  se <- which(resdata$Index %in% select)
  max <- max(ceiling(resdata$Enrichment.Ratio[se]))+1
  p<- ggplot(resdata[se,] ) +
      geom_point(
        mapping = aes(x=Cell.Interaction, y=Enrichment.Ratio, 
                      color=get(condition), 
                    size=No.Shared.IMGP)) +
    scale_color_continuous(low="#9970ab", high="#5aae61", 
                           limits = c(0, cutoff),
                           name = "FDR",
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range = c(3,8), name="No. Shared\nMarker Pair") +
    geom_hline(yintercept=0, lty=4, col="black", lwd=1)+
    
    theme_bw(base_size = 20, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Cell Interaction",y="Enrichment Ratio",
         title="Cell Interaction Enrichment",
         base_size = 12, base_family = "serif",face="bold") +
    theme(legend.background = element_rect(linetype = 1, 
                                           linewidth = 0.5, colour = 1),
          legend.position="right",
          legend.box = "vertical",
          legend.direction= "vertical",
          plot.margin = unit(c(1, 1, 1, 5), "lines"),
          panel.grid=element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.title = element_text(face="bold", color="black",
                                      family = "serif", size=10),
          legend.text= element_text(face="bold", color="black",
                                    family = "serif", size=10),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black", size=10, 
                                     angle = 40,  vjust = 1, hjust=1),
          axis.text.y = element_text(face="bold", color="black", size=13),
          axis.title.x = element_text(face="bold", color="black", size=15),
          axis.title.y = element_text(face="bold",color="black", size=15),
          axis.line.y = element_line(arrow = grid::arrow(
            length = unit(0.3, "cm"), ends = "last"))) +
    scale_y_continuous(limits=c(0, max), breaks=seq(0,max, 2),expand = c(0,0)) +
    guides(shape = "none") +
    annotate("text", x = -1, y = max, label = "Favorable\nOutcome")+
    coord_cartesian(xlim = c(1,nrow(resdata[se,] )),
                    clip = 'off')
  
  return(p)
}
