###################################################
## Main entry point for MULTI-seq classification ##
###################################################

perform_sweep <- function(bar.table, plot=F) {
    
    bar.table_sweep.list <- list()
    n <- 0
    for (q in seq(0.01, 0.99, by=0.02)) {
      n <- n + 1
      bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
      names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
    }
    threshold.results <- findThresh(call.list=bar.table_sweep.list)
    
    if (plot) {
        ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
              geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
    }
    
    return(threshold.results)
}

perform_classification <- function(bar.table) {
    
    qthresh = perform_sweep(bar.table, plot=F)
    round.calls <- classifyCells(bar.table, q=findQ(qthresh$res, qthresh$extrema))
    neg.cells <- names(round.calls)[which(round.calls == "Negative")]
        
    # remove neg cells
    bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
    all.neg.cells = neg.cells
    
    while (length(neg.cells) > 0) {
        
        qthresh = perform_sweep(bar.table, plot=F)
        round.calls <- classifyCells(bar.table, q=findQ(qthresh$res, qthresh$extrema))
        neg.cells <- names(round.calls)[which(round.calls == "Negative")]
        if (length(neg.cells) > 0) {
            bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
        }
        all.neg.cells = c(all.neg.cells, neg.cells)
        
    }
    
    final.calls <- c(round.calls, rep("Negative",length(all.neg.cells)))
    names(final.calls) <- c(names(round.calls),all.neg.cells)
    return(final.calls)
    
}

analyze_multi_sample <- function(bar.ref, cell.id.vec, R1, R2, cell.pos = c(1,16), umi.pos = c(17,26), 
                                    tag.pos = c(1,8), exp.name = NA, write=TRUE) {
    
    if (is.na(exp.name)) {
        warning("You didn't specify an experiment name, automatically setting this to the default: 'multi_analysis' ")

        exp.name = 'multi_analysis'
    }

    message("Preprocessing...")
    readTable <- MULTIseq.preProcess(R1 = R1, R2 = R2, cellIDs = cell.id.vec, cell=cell.pos, umi=umi.pos, tag=tag.pos)

    message("Aligning...")
    bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
    colnames(bar.table) = bar.ref
    rownames(bar.table) = cell.id.vec
    
    message("Classifying cells...")
    final.calls = perform_classification(bar.table)

    if (write) { 
        write.table(final.call, paste0(exp.name, '.txt'), sep='\t')
    }

    message("Producing tSNE projection")
    bar.table.n = bar.table[,1:96]
    bar.table.n = apply(bar.table.n, c(1,2), log1p)

    colnames(tsne.res) <- c("TSNE1","TSNE2")
    rownames(tsne.res) <- rownames(bar.table)
    tsne.res['label'] = final.calls[rownames(tsne.res)]

    barcoded = tsne.res[which((tsne.res$label != 'Doublet') & (tsne.res$label != 'Negative')),]

    neg = tsne.res[which(tsne.res$label == 'Negative'),]
    doub = tsne.res[which(tsne.res$label == 'Doublet'), ]

    ggplot(barcoded, aes(TSNE1, TSNE2, color=label)) + 
        geom_point(data = neg, aes(TSNE1, TSNE2, color = 'Negative'), color = 'grey') + 
        geom_point(data = doub, aes(TSNE1, TSNE2, color = 'Doublet'), color = 'black') + 
        geom_point(data = barcoded, aes(TSNE1, TSNE2, color=label)) + 
        theme_classic()
    
    ggsave(paste0(exp.name, '.tsne.pdf'))

    return(final.calls)

}