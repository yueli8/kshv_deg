library(ggplot2)
setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/uci/david")
delta11_mock_bp<-read.table("delta11_vs_mock_bp01.txt", sep="\t", header=TRUE)
S1<- ggplot(delta11_mock_bp, aes(x= Enrichment_Ratio, y=Description, size=Gene_Count, color=P_Value, group=Enrichment_Ratio)) + geom_point(alpha = 25) + theme(axis.text.x = element_text(face="bold",size=20), axis.text.y = element_text(face="bold",size=20))
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
S1<-S1 + theme(axis.title.x =element_text(face="bold", size=25), axis.title.y=element_text(face="bold",size=25))
S1

delta11_mock_cc<-read.table("delta11_vs_mock_cc01.txt", sep="\t", header=TRUE)
S1<- ggplot(delta11_mock_cc, aes(x= Enrichment_Ratio, y=Description, size=Gene_Count, color=P_Value, group=Enrichment_Ratio)) + geom_point(alpha = 25) + theme(axis.text.x = element_text(face="bold",size=20), axis.text.y = element_text(face="bold",size=20))
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
S1<-S1 + theme(axis.title.x =element_text(face="bold", size=25), axis.title.y=element_text(face="bold",size=25))
S1

delta11_mock_mf<-read.table("delta11_vs_mock_mf01.txt", sep="\t", header=TRUE)
S1<- ggplot(delta11_mock_mf, aes(x= Enrichment_Ratio, y=Description, size=Gene_Count, color=P_Value, group=Enrichment_Ratio)) + geom_point(alpha = 25) + theme(axis.text.x = element_text(face="bold",size=20), axis.text.y = element_text(face="bold",size=20))
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
S1<-S1 + theme(axis.title.x =element_text(face="bold", size=25), axis.title.y=element_text(face="bold",size=25))
S1

delta11_wt_bp<-read.table("delta11_vs_wt_bp01.txt", sep="\t", header=TRUE)
S1<- ggplot(delta11_wt_bp, aes(x= Enrichment_Ratio, y=Description, size=Gene_Count, color=P_Value, group=Enrichment_Ratio)) + geom_point(alpha = 25) + theme(axis.text.x = element_text(face="bold",size=20), axis.text.y = element_text(face="bold",size=20))
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
S1<-S1 + theme(axis.title.x =element_text(face="bold", size=25), axis.title.y=element_text(face="bold",size=25))
S1
