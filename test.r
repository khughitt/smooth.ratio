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

# load data
source('smooth.r')
load(file='data/sample.rda')
result = fast.smooth(cpgsites, methylation, coverage)

# plot coarse scale
library(ggplot2)
dat = stack(as.data.frame(result$small))
dat$x <- rep(seq_len(nrow(result$small)), ncol(result$small))
ggplot(dat, aes(x, values)) + geom_line(aes(colour=ind))

