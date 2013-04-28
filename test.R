#!/usr/bin/env Rscript
#################################################
# Box Smooth Test
# Keith Hughitt <khughitt@umd.edu>
# 2013/04/19
#
# This file is a translation of test.m written by
# Chengxi ye.
#################################################
library(R.matlab)

# load data
data = readMat('test.mat')
attach(data)

# smoothing parameters
bins = 350
radius = 350

# create matrix to store result
smoothed_percentage = matrix(0, nrow(coverage), ncol(coverage))
num_rows = nrow(smoothed_percentage)

# smoothing
for (i in 1:ncol(coverage)) {
    c_cum = cumsum(coverage[,i])
    m_cum = cumsum(methylation[,i])

    smoothed_percentage[1:radius,i] = m_cum[(1+radius):(2*radius)] / 
                                      c_cum[(1+radius):(2*radius)]

    smoothed_percentage[(num_rows - radius + 1):num_rows, i] = (matrix(rep(m_cum[num_rows], radius)) - m_cum[(num_rows-2*radius+1):(num_rows-radius)]) / 
                                                               (matrix(rep(m_cum[num_rows], radius)) - c_cum[(num_rows-2*radius+1):(num_rows-radius)])

    smoothed_percentage[(radius+1):(num_rows-radius),i] = (m_cum[(1+2*radius):num_rows] - m_cum[1:(num_rows-2*radius)]) /
                                                          (c_cum[(1+2*radius):num_rows] - c_cum[1:(num_rows-2*radius)])
}



