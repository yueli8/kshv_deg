library(fgsea)
library(tidyverse)
library(data.table)
library(ggplot2)
setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/uci/deg")

#delta11_vs_mock
res<-read.table("delta11_vs_mock_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
#pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
#fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="delta11_vs_mock_kegg.txt", sep="\t", sep2=c("", " ", ""))
#fwrite(fgseaRes, file="delta11_vs_mock_all.txt", sep="\t", sep2=c("", " ", ""))


#wt_vs_mock
res<-read.table("wt_vs_mock_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
#pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
#fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="wt_vs_mock_kegg.txt", sep="\t", sep2=c("", " ", ""))
#fwrite(fgseaRes, file="delta11_vs_mock_all.txt", sep="\t", sep2=c("", " ", ""))

plotEnrichment(pathways.kegg[["KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"]],
               ranks) + labs(title="COMPLEMENT_AND_COAGULATION_CASCADES")

plotEnrichment(pathways.kegg[["KEGG_RIBOSOME"]],
               ranks) + labs(title="RIBOSOME")



#delta11_vs_wt
res<-read.table("delta11_vs_wt_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fwrite(fgseaRes_kegg, file="delta11_vs_wt_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="delta11_vs_wt_all.txt", sep="\t", sep2=c("", " ", ""))

plotEnrichment(pathways.kegg[["KEGG_HUNTINGTONS_DISEASE"]],
               ranks) + labs(title="HUNTINGTONS_DISEASE")

plotEnrichment(pathways.kegg[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) + labs(title="OXIDATIVE_PHOSPHORYLATION")

#mock_vs_wt
res<-read.table("mock_vs_wt_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fwrite(fgseaRes_kegg, file="mock_vs_wt_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="mock_vs_wt_all.txt", sep="\t", sep2=c("", " ", ""))


