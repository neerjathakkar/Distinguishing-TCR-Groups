library(ggplot2)
library(grid)
library(gtable)

pathologies <- read.csv('datasets/mcpas_data/pathologies.csv') # TODO: subset this to what matters
small.raw <- read.csv('results/mcpas_results/small_human_results.csv')
large.raw <- read.csv('results/mcpas_results/large_human_results.csv')
# for processing, slap it all together with an extra column
all.raw <- rbind(cbind(small.raw, list(size='small')),
                 cbind(large.raw, list(size='large')))
# convert internal abbrev to label name
ps <- pathologies[pathologies$abbrev != '',]
abbrev2label <- setNames(ps$label, ps$abbrev)
all.raw$p1 <- abbrev2label[as.character(all.raw$p1)]
all.raw$p2 <- abbrev2label[as.character(all.raw$p2)]
# make symmetric
all.flipped <- all.raw; all.flipped[,'p1'] <- all.raw[,'p2']; all.flipped[,'p2'] <- all.raw[,'p1']
all.symm <- rbind(all.raw, all.flipped)

# only 0.2, 0.3, and 0.4 thresholds
all.symm <- all.symm[all.symm$threshold %in% c('0.2','0.3','0.4'),]
# order according to pathologies.csv (which groups together types of pathologies)
all.symm$p1 <- factor(all.symm$p1, levels=pathologies$label, ordered=TRUE)
all.symm$p2 <- factor(all.symm$p2, levels=pathologies$label, ordered=TRUE)

# separate three copies (unidentified, correct, incorrect), then merge back with $prediction indicating source
unid <- all.symm; unid[,'pct'] <- 100*unid[,'unidentified_frac']; unid[,'prediction'] <- 'unidentified'
cor <- all.symm; cor[,'pct'] <- 100*cor['correct_frac']; cor[,'prediction'] <- 'correct'
inc <- all.symm; inc[,'pct'] <- 100*inc['incorrect_frac']; inc[,'prediction'] <- 'incorrect'
combined <- rbind(unid,cor,inc)[,c('p1','p2','size','threshold','prediction','pct')]

# my best stab at colors from other figs
col.unid <- '#0F7E12'; col.cor <- '#0A20E4'; col.inc <- '#AE19AE'

# function to do the plotting for large or small
# (uses global variables from above)
plot_triangle <- function(size, width) {
  data <- combined[combined$size==size,]
  labels <- pathologies$label[pathologies$group==size]
  # subset to focus on lower triangle
  data <- data[data$p1 %in% labels[2:length(labels)] & data$p2 %in% labels[1:(length(labels)-1)],]
  gg <- ggplot(data, aes(x=threshold, weight=pct, fill=prediction)) +
    geom_bar() +
    scale_fill_manual(values=c('unidentified'=col.unid, 'correct'=col.cor, 'incorrect'=col.inc)) +
    facet_grid(rows=vars(p1), cols=vars(p2)) +
    theme_classic() + theme(strip.background=element_rect(colour="white", fill="white")) +
    labs(x='threshold', y='% predictions')
  # remove upper triangle, thanks to https://stackoverflow.com/questions/49521848/remove-unused-facet-combinations-in-2-way-facet-grid
  # (which also suggests maybe should've done this a different way to start with, but oh well)
  grob <- ggplotGrob(gg)
  for (row in 1:(length(labels)-2)) {
    for (col in (row+1):(length(labels)-1)) {
      idx <- which(grob$layout$name==paste0('panel-',row,'-',col))
      grob$grobs[[idx]] <- nullGrob()
    }
  }
  grid.newpage()
  grid.draw(grob)
  ggsave(plot=grob, filename=paste0('results/mcpas_results/compare-',size,'.pdf'), width=width, units='in')
}

plot_triangle('small',  8)
plot_triangle('large', 7)
