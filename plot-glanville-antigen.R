library(ggplot2)
library(grid)
library(gtable)

# TODO: copied & modified from plot-mcpas -- generalize

raw <- read.csv('results/glanville_results/curated_separate_pp65_results_hla.csv') # edited to include HLA with each antigen
# make symmetric
flipped <- raw; flipped[,'p1'] <- raw[,'p2']; flipped[,'p2'] <- raw[,'p1']
symm <- rbind(raw, flipped)
# keep order from their table, which puts same HLA adjacent
ags <- c('pp50_A1','NP44_A1','pp65_A2','BMLF1_A2','M1_A2','pp65_B7','NP177_B7')
symm$p1 <- factor(symm$p1, levels=ags, ordered=TRUE)
symm$p2 <- factor(symm$p2, levels=ags, ordered=TRUE)

# only 0.2, 0.3, and 0.4 thresholds
symm <- symm[symm$threshold %in% c('0.2','0.3','0.4'),]

# separate three copies (unidentified, correct, incorrect), then merge back with $prediction indicating source
unid <- symm; unid[,'pct'] <- 100*unid[,'unidentified_frac']; unid[,'prediction'] <- 'unidentified'
cor <- symm; cor[,'pct'] <- 100*cor['correct_frac']; cor[,'prediction'] <- 'correct'
inc <- symm; inc[,'pct'] <- 100*inc['incorrect_frac']; inc[,'prediction'] <- 'incorrect'
combined <- rbind(unid,cor,inc)[,c('p1','p2','threshold','prediction','pct')]

# my best stab at colors from other figs
col.unid <- '#0F7E12'; col.cor <- '#0A20E4'; col.inc <- '#AE19AE'

# subset to focus on lower triangle
combined <- combined[combined$p1 %in% ags[2:length(ags)] & combined$p2 %in% ags[1:(length(ags)-1)],]
gg <- ggplot(combined, aes(x=threshold, weight=pct, fill=prediction)) +
  geom_bar() +
  scale_fill_manual(values=c('unidentified'=col.unid, 'correct'=col.cor, 'incorrect'=col.inc)) +
  facet_grid(rows=vars(p1), cols=vars(p2)) +
  theme_classic() + theme(strip.background=element_rect(colour="white", fill="white")) +
  labs(x='threshold', y='% predictions')
# remove upper triangle, thanks to https://stackoverflow.com/questions/49521848/remove-unused-facet-combinations-in-2-way-facet-grid
# (which also suggests maybe should've done this a different way to start with, but oh well)
grob <- ggplotGrob(gg)
for (row in 1:(length(ags)-2)) {
  for (col in (row+1):(length(ags)-1)) {
    idx <- which(grob$layout$name==paste0('panel-',row,'-',col))
    grob$grobs[[idx]] <- nullGrob()
  }
}
grid.newpage()
grid.draw(grob)
ggsave(plot=grob, filename=paste0('results/glanville_results/curated-antigen_hla-pairs.pdf'), width=8, height=7, units='in')
