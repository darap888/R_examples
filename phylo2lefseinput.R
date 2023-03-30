#Script to prepare data
#from phyloseq object to lefse input https://huttenhower.sph.harvard.edu/galaxy/

library(dplyr)

L1=otu_table(phylo)%>%data.frame(check.names = F,stringsAsFactors=F)%>%tibble::rownames_to_column()
L2=tax_table(phylo)%>%data.frame(check.names = F,stringsAsFactors=F)%>%tibble::rownames_to_column()
L2$X16S_rRNA_Count=NULL
L3=left_join(L1,L2,by="rowname")

rank_lv=rank_names(phylo)[1:6]
#x=L3[83,]
L4=apply(L3,1,function(x){
  x1=lapply(x,rep,each=6)%>%data.frame(check.names = F,stringsAsFactors=F)
  x1$class=sapply(1:6,function(i){
    if(is.na(x1[i,rank_lv[i]])|x1[i,rank_lv[i]]==""){
      return("")
    } else {
      paste(x1[i,rank_lv][1:i],collapse="|",sep="")
    }
  })
  x2=dplyr::filter(x1,class!="")
  x2$rowname=NULL
  for(i in 1:6){x2[[rank_lv[i]]]=NULL}
  return(x2)
}) %>% bind_rows()

L5=group_by(L4,class)%>%mutate_all(funs(as.numeric))%>%summarise_all(funs(sum))
L6=rbind(c("group",as.character(get_variable(phylo,g1))),as.matrix(L5))
L7=rbind(c("days",as.character(get_variable(phylo,"Days"))),as.matrix(L6))

write.table(L7,file=sprintf("lefse0.txt"), append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F, qmethod = c("escape", "double"),fileEncoding = "")