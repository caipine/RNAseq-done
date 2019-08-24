
t1 <- resCond_vV[resCond_vV$Symbol %in% namei,] 
t1_geneid<- t1[order(t1$pvalue),][1,1]
t2 <- data.frame(normalizedCounts[rownames(normalizedCounts) %in% t1_geneid,])
colnames(t2) <- "reads"
t2$sample <- rownames(t2)
t2$gene <- namei
t2$condition <- sampleinfo$condition
t3 <- rbind(t3,t2)

