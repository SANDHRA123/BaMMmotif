################################################################################
#
# Documentations
# 
################################################################################

# Plot higher-order BaMM logo
# using files with BaMM probabilities as obtained from BaMM!motif

# This routine is built on code from the Bioconductor package seqLogo
# bioconductor.org/packages/release/bioc/html/seqLogo.html

################################################################################
#
# R Libraries
# 
################################################################################

library( grid )

################################################################################
#
# R Sources
# 
################################################################################

# main directory
mainDir <- file.path( "", "home", "user", "workspace", "BaMMmotif", "R" )

# function(s) to plot BaMM logo
source( file.path( mainDir, "plotBaMM.R" ) )
# function(s) to read in BaMM file
source( file.path( mainDir, "readBaMM.R" ) )

################################################################################
#
# Parameter settings
# 
################################################################################

# directory with BaMM files
bammDir <- file.path( mainDir, "BaMMs" )

# directory for resulting logo
logoDir <- file.path( mainDir, "Logos" )

# file with BaMM probabilities as obtained from BaMM!motif
# omit filename extensions
filename <- "wgEncodeSydhTfbsHepg2Mafkab50322IggrabAlnRep0_summits102_restr5000_k2-K2-a1-20-60"

# BaMM logo order (maximum order = 5)
order <- 0

# use RNA alphabet (ACGU) instead of DNA (ACGT) alphabet
RNASeqLogo <- FALSE

# use nucleotide frequencies from <filename>.freqs file to calculate 0'th-order
# background frequencies (default: uniformly distributed nucleotide frequencies)
useFreqs <- FALSE

# scale columns according to information content (instead of probability)
icColumnScale <- TRUE # icLetterScale and not icColumnScale is not implemented

# scale nucleotides according to information content (instead of probability)
icLetterScale <- TRUE # icLetterScale and not icColumnScale is not implemented

# alpha channel for transparency
alpha <- .3 # icLetterScale only

# plot x-axis
plot.xaxis <- TRUE

# plot x-axis title: Model position [nt]
plot.xlab <- TRUE

# plot labels of first and last BaMM positions only
plot.border.labels <- TRUE

# BaMM position that defines an anchor point (on the x-axis)
# e.g. transcription start site (TSS)
xanchor <- NULL
# xanchor <- 51

# numeric or character vector of length 1 or 3 giving anchor point label(s)
xanchor.labels <- FALSE
# xanchor.labels <- "TSS"

# plot x-axis anchor labels
plot.xanchor <- FALSE # numeric xanchor.labels only

# plot y-axis
plot.yaxis <- TRUE

# plot y-axis title: Information content
plot.ylab <- TRUE

# numeric vector of length 2 giving the y coordinate range.
ylim <- c( -.5, 2 )
# ylim <- c( -.3, .5 )

# the points at which tick-marks are to be drawn.
yat <- c( -.5, 0, 1, 2 )
# yat <- c( -.3, 0, .5 )

# a numeric value specifying the vertical justification of the y-axis title
yaxis.vjust <- .5

# a numerical value giving the amount by which plotting text and symbols are
# magnified relative to the default
cex <- 2.7

# line width
lwd <- 4

# the width of the svg device in inches
width <- 7

# the height of the svg device in inches
height <- 7

# the default pointsize of plotted text (in big points)
pointsize <- 12

################################################################################
#
# Code to plot BaMM logo
# 
################################################################################

if( !( file.exists( logoDir ) ) ) dir.create( logoDir )

svg( file.path( logoDir, paste( filename, "_logo-k", order, ".svg", sep="" ) ),
     width=width, height=height, pointsize=pointsize )
opar <- par( no.readonly=TRUE )

hoSeqLogo( filename=file.path( bammDir, filename ), order=order, rna=RNASeqLogo,
           useFreqs=useFreqs, icColumnScale=icColumnScale,
           icLetterScale=icLetterScale, alpha=alpha, plot.xaxis=plot.xaxis,
           plot.xlab=plot.xlab, plot.border.labels=plot.border.labels,
           xanchor=xanchor, xanchor.label=xanchor.labels,
           plot.xanchor=plot.xanchor, plot.yaxis=plot.yaxis,
           plot.ylab=plot.ylab, ylim=ylim, yat=yat, yaxis.vjust=yaxis.vjust,
           cex=cex, lwd=lwd )

par( opar )
dev.off()
