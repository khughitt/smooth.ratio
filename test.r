#!/usr/bin/env Rscript
################################################################################
# Box Smooth Test
# Keith Hughitt <khughitt@umd.edu>
# Chengxi Ye <cxy@cs.umd.edu>
# 2013/04/28
#
# This file is a translation of Matlab code originally written by Chengxi Ye.
#
################################################################################
library(ggplot2)

# load data
source('smooth.r')
load(file='data/sample.rda')
dat = fast.smooth(x$cpg0x2Dsites1, x$methylation1, x$coverage1)

# plot coarse scale
dat = stack(as.data.frame(smoothed_signal))
dat$x <- rep(seq_len(nrow(smoothed_signal)), ncol(smoothed_signal))
ggplot(dat, aes(x, values)) + geom_line(aes(colour=ind))

