

library(microbiome)

otus_of_interest= list()


prevalence <-
  function(ODLEPobj_cecum, ODLEPobj_ileum, ODLEPobj_feces) {
    #Ileum 
    pseq.rel_ileum <- microbiome::transform(ODLEPobj_ileum, "compositional")
    otu_relative_ileum <- as.data.frame(otu_table(pseq.rel_ileum))
    otu_relative_ileum = otu_relative_ileum[!apply(otu_relative_ileum, 1, function(x) all(x == 0)), ]
    total_samples= ncol(otu_relative_ileum) 
    total_otus= nrow(otu_relative_ileum) 
    absent=apply(otu_relative_ileum ==0, 1, sum)  
    #add a col with the percentage of samples it is present
    otu_relative_ileum$percentage_samples_present=(1-(absent/total_samples))*100 
    #add a col with OTU names
    otu_relative_ileum$OTU <- row.names(otu_relative_ileum)
    #add otus of interest to highlight in the plot
    otu_relative_ileum$of_interest <- ifelse(otu_relative_ileum$OTU %in% otus_of_interest == TRUE, "YES", "NO")
    
    #plot
    colors <- c("#075b44","#f56342")
    sizes <- c(0.05, 4)
    prev_plot_ileum=ggplot(data = otu_relative_ileum) +
      theme_classic()+
      aes(x = reorder(OTU,-percentage_samples_present,sum), y =percentage_samples_present, color=of_interest, size=of_interest) +
      geom_point() +
      scale_colour_manual(values=colors) + 
      scale_size_manual (values=sizes)+
      #geom_point(color="#343aeb", size=0.05) +
      ggtitle("Ileum") +
      #xlab("OTU") + ylab("samples where OTU is present (%)")+
      theme(legend.position = "none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    
    #Feces
    pseq.rel_feces <- microbiome::transform(ODLEPobj_feces, "compositional")
    otu_relative_feces <- as.data.frame(otu_table(pseq.rel_feces))
    otu_relative_feces = otu_relative_feces[!apply(otu_relative_feces, 1, function(x) all(x == 0)), ]
    total_samples= ncol(otu_relative_feces) 
    total_otus= nrow(otu_relative_feces) 
    absent=apply(otu_relative_feces == 0, 1, sum)  
    otu_relative_feces$percentage_samples_present=(1-(absent/total_samples))*100 
    otu_relative_feces$OTU <- row.names(otu_relative_feces)
    otu_relative_feces$of_interest <- ifelse(otu_relative_feces$OTU %in% otus_of_interest == TRUE, "YES", "NO")
    
    #plot
    colors <- c("#fcd8b6","#f56342")
    sizes <- c(0.05, 4)
    prev_plot_feces=ggplot(data = otu_relative_feces) +
      theme_classic()+
      aes(x = reorder(OTU,-percentage_samples_present,sum), y =percentage_samples_present, color=of_interest, size=of_interest) +
      geom_point() +
      scale_colour_manual(values=colors) + 
      scale_size_manual (values=sizes)+
      #geom_point(color="#16DCC9", size= 0.05) +
      ggtitle("Feces") +
      #xlab("OTU") + ylab("samples where OTU is present (%)")+
      theme(legend.position = "none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    
    #Cecum 
    pseq.rel_cecum <- microbiome::transform(ODLEPobj_cecum, "compositional")
    otu_relative_cecum <- as.data.frame(otu_table(pseq.rel_cecum))
    otu_relative_cecum = otu_relative_cecum[!apply(otu_relative_cecum, 1, function(x) all(x == 0)), ]
    total_samples= ncol(otu_relative_cecum) 
    total_otus= nrow(otu_relative_cecum) 
    absent=apply(otu_relative_cecum ==0, 1, sum) 
    otu_relative_cecum$percentage_samples_present=(1-(absent/total_samples))*100 
    otu_relative_cecum$OTU <- row.names(otu_relative_cecum)
    otu_relative_cecum$of_interest <- ifelse(otu_relative_cecum$OTU %in% otus_of_interest == TRUE, "YES", "NO")
    
    #plot
    colors <- c("#035060","#f56342")
    sizes <- c(0.05, 4)
    prev_plot_cecum=ggplot(data = otu_relative_cecum) +
      theme_classic()+
      aes(x = reorder(OTU,-percentage_samples_present,sum), y =percentage_samples_present, color=of_interest,   size=of_interest) +
      geom_point() +
      scale_colour_manual(values=colors) + 
      scale_size_manual (values=sizes)+
      #geom_point(color="#BA4FC8", size=0.05) +
      ggtitle("Cecum") +
      #xlab("OTU") + ylab("samples where OTU is present (%)")+
      theme(legend.position = "none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    
    ##  Add dashed line to each plot
    #ileum 
    otu_relative_ileum_more_50 = otu_relative_ileum[otu_relative_ileum$percentage_samples_present >= 50, c("OTU", "percentage_samples_present")]
    more50_ileum = dim(otu_relative_ileum_more_50)[1]
    prev_plot_ileum = prev_plot_ileum + geom_segment(aes(x=0, y=50, xend=more50_ileum, yend=50), linetype="dashed", color="black")
    prev_plot_ileum = prev_plot_ileum + geom_segment(aes(x=more50_ileum, y=0, xend=more50_ileum, yend=50), linetype="dashed", color="black")
    prev_plot_ileum= prev_plot_ileum + annotate("text", x = more50_ileum, y = 50, label = paste0("                        ", round(more50_ileum/dim(otu_relative_ileum)[1] * 100, 1), "% OTUs"))
    
    #cecum
    otu_relative_cecum_more_50 = otu_relative_cecum[otu_relative_cecum$percentage_samples_present >= 50, c("OTU", "percentage_samples_present")]
    more50_cecum = dim(otu_relative_cecum_more_50)[1]
    prev_plot_cecum = prev_plot_cecum + geom_segment(aes(x=0, y=50, xend=more50_cecum, yend=50), linetype="dashed", color="black")
    prev_plot_cecum = prev_plot_cecum + geom_segment(aes(x=more50_cecum, y=0, xend=more50_cecum, yend=50), linetype="dashed", color="black")
    prev_plot_cecum= prev_plot_cecum + annotate("text", x = more50_cecum, y = 50, label = paste0("                        ", round(more50_cecum/dim(otu_relative_cecum)[1] * 100, 1), "% OTUs"))
    
    #feces
    otu_relative_feces_more_50 = otu_relative_feces[otu_relative_feces$percentage_samples_present >= 50, c("OTU", "percentage_samples_present")]
    more50_feces = dim(otu_relative_feces_more_50)[1]
    prev_plot_feces = prev_plot_feces + geom_segment(aes(x=0, y=50, xend=more50_feces, yend=50), linetype="dashed", color="black")
    prev_plot_feces = prev_plot_feces + geom_segment(aes(x=more50_feces, y=0, xend=more50_feces, yend=50), linetype="dashed", color="black")
    prev_plot_feces= prev_plot_feces + annotate("text", x = more50_feces, y = 50, label = paste0("                        ", round(more50_feces/dim(otu_relative_feces)[1] * 100, 1), "% OTUs"))

    fecesper=paste0(round(more50_feces/dim(otu_relative_feces)[1] * 100, 1), "%")
    cecumper=paste0(round(more50_cecum/dim(otu_relative_cecum)[1] * 100, 1), "%")
    ileumper=paste0(round(more50_ileum/dim(otu_relative_ileum)[1] * 100, 1), "%")
    r <- list(prev_plot_cecum,prev_plot_ileum, prev_plot_feces)
    per <- list(cecumper,ileumper, fecesper)
    names(r) <- c("cecum", "ileum", "feces")
    names(per) <- c("cecum", "ileum", "feces")
    out <- list(r, per)
    return(out)
  }
    
