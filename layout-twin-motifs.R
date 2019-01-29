require(png)

layout <- read.csv('results/twin_results/motif_figure/twin_motifs.csv')
groups <- c('A1-specific', 'twin-specific', 'twin-distinct', 'common') # TODO: good order? reverse?
individuals <- c('A1','A2','C1','C2','D1','D2') # force this order for clarity

# leaving a blank row between groups
xpad <- 0.4 # estimate of width of logo, to pad x axis
plot(c(0,xpad+max(layout$A1.distance)), c(0.5, length(groups)*(1+length(individuals))), type='n', xlab='distance to A1', ylab='individual', yaxt='n', cex.lab=1.5, cex.axis=1.25)
yticks <- NULL
for (g in 0:(length(groups)-1)) {
  yticks <- c(yticks, seq(1 + g*(1+length(individuals)), (g+1)*(1+length(individuals))-1))
}
axis(2, at=yticks, labels=rep(individuals, length(groups)), las=2)
y <- 1
for (group in groups) { 
  mtext(group, side=4, at=y+3)
  for (individual in individuals) { 
    row <- layout[layout$group==group & layout$individual==individual,] # inefficient, but oh well
    x <- row$A1.distance
    img <- readPNG(paste0('results/twin_results/motif_figure/',row$file,'.png'))
    w <- dim(img)[2]; h <- dim(img)[1]
    rasterImage(img, x,y, x+h/w, y+1)
    y <- y+1
  }
  y <- y+1 # blank between groups
}
dev.print(pdf, width=6, height=10, file='results/twin_results/motif_figure/twin-motifs.pdf')
