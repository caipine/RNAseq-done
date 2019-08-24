tmp <- hmCount[hmCount$Symbol %in% namei,] 
rownames(tmp) <- paste(tmp$GeneID,tmp$Symbol, sep = ".")
tmp2 <- (rotate_df(tmp[,2:17]))
tmp3 <- cbind(sampleinfo,tmp2)

colnames(tmp3)

ggplot(tmp3, aes(x=condition, y=ENST00000572832.1.PKMYT1)) +
      geom_boxplot()  +
      geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="black")
list1 <- colnames(tmp3)[5:ncol(tmp3)]
list1

tiff(paste("boxdot_", namei,".tif", sep ='') ,width=300,height=300) 
ggplot(tmp3, aes_string("condition", list1[1])) +
      geom_boxplot()  +
      geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="black")
dev.off()


