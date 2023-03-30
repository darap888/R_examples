#example of working functions to create volcano plots 
#from custom files

#it was intended for use in for-cycle to compare multiple selected groups of samples
#(not uploaded) - read from external output (Agilent volcano plot tables)

#not optimized and not intended for end-users



volcano_plot3 = function(volcano_tab0, strain_compounds.list) {
  #volcano_tab=`14SAvsi14SA`
  #list(ST.14,ST.14.pc)[[2]]
  vs1=str_split(volcano_tab0,"vs")[[1]][1]
  vs2=str_split(volcano_tab0,"vs")[[1]][2]
  #library(ggrepel)
  #library(tidyverse)
  #
  cols <- c("PC&CC" = "white", "rest" = "darkred", "PC (but not CC)" = "steelblue") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  
  title0=str_c("Upregulated in ",vs2," (left) vs in ",vs1," (right)")
  
  volcano_tab=get(volcano_tab0) %>% select(1,3,5,9) 
  colnames(volcano_tab)=(c("Compound0","Pcorr","Regulation0","Log2FC"))
  volcano_tab=volcano_tab %>% mutate(Regulation=if_else(Pcorr<0.05,Regulation0,"ns"),
                                     Compound=if_else(Pcorr<0.05&Compound0%in%str_subset(volcano_tab$Compound0,"^C(?=\\d)|@|^>|<",negate = TRUE),Compound0,""))
  volcano_tab=volcano_tab %>% 
    
    mutate(Detected_in=case_when(Pcorr<0.05&(Compound0%in%strain_compounds.list[[1]])~"PC&CC",Pcorr<0.05&(Compound0%in%strain_compounds.list[[2]])~"PC (but not CC)",Pcorr<0.05~"rest"))
  
  
  
  FCmax=volcano_tab %>% filter(Pcorr<0.05) %>%  pull((Log2FC)) %>% abs() %>% max()
  
  vol_plot <- volcano_tab %>%
    ggplot(aes(x = Log2FC,
               y = -log10(Pcorr), size=Regulation, alpha=Regulation, fill = Detected_in)
    ) + 
    geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
               colour = "black", alpha=0.7) + 
    theme_minimal()+
    
    
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed")  +
    
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) + # Modify point transparency
    scale_x_continuous(breaks = c(seq(-10, 10, 2)),  limits = c(-(FCmax+0.5), FCmax+0.5))+
    #scale_y_continuous(breaks = c(seq(0, 20, 2)))+
    
    geom_label_repel( # Add labels last to appear as the top layer  
      aes(label = Compound),
      force = 2,
      nudge_y = 1)+
    
    
    labs(title = title0,
         x = "log2(fold change)",
         y = "-log10(adjusted P-value)") +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = cols) 
  
  
  return(vol_plot)
}


volcano_features_stack = function(volcano_tab0,strain_compounds.list){
cols2 <- c("PC&CC" = "lightgrey", "rest" = "darkred", "PC (but not CC)" = "steelblue") 

#volcano_tab0="SCvsSA"
trt1=str_split(volcano_tab0,"vs")[[1]][1]
trt2=str_split(volcano_tab0,"vs")[[1]][2]

volcano_tab=get(volcano_tab0) %>% select(1,3,5,9) 
colnames(volcano_tab)=(c("Compound0","Pcorr","Regulation0","Log2FC"))
volcano_tab=volcano_tab %>% mutate(Regulation=if_else(Pcorr<0.05,Regulation0,"ns"),
                                   Compound=if_else(Pcorr<0.05&Compound0%in%str_subset(volcano_tab$Compound0,"^C(?=\\d)|@|^>|<",negate = TRUE),Compound0,""))
volcano_tab=volcano_tab %>% 
  
  mutate(Detected_in=case_when(Pcorr<0.05&(Compound0%in%strain_compounds.list[[1]])~"PC&CC",Pcorr<0.05&(Compound0%in%strain_compounds.list[[2]])~"PC (but not CC)",Pcorr<0.05~"rest"))

p=volcano_tab %>% filter(Pcorr<0.05&!(is.na(Regulation))) %>% mutate(log2FC.lvl=if_else(abs(Log2FC)>1,"large_FC","small_FC"),
                                              Upregulated_in=case_when(Regulation=="up"~trt1,
                                                                   Regulation=="down"~trt2,
                                                                   TRUE~"ns")) %>% 
  select(Upregulated_in, log2FC.lvl,Detected_in) %>% 
 
 table() %>% as.data.frame() %>%  #rename_at("Freq",~"Number of significantly different compounds")
  
ggplot(aes( y = Freq, fill = Detected_in, x=Upregulated_in))+
  geom_bar(stat = "identity", colour="black")+theme_bw()+facet_grid(~log2FC.lvl)+
  
  ylab("Number of differentially abundant compounds")+
 scale_fill_manual(values = cols2)

return(p)
}
