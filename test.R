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
library(R.matlab)
library(signal)
library(ggplot2)

# load data
data = readMat('test.mat')
attach(data)

# Domain sigma
sigma_d    = max(cpg0x2Dsites1) / 100
sampling_d = sigma_d
derived_sigma = sigma_d / sampling_d

xi    = round(cpg0x2Dsites1 / sampling_d) + 1
max_x = max(xi)
numerator   = matrix(0, max_x, ncol(methylation1))
denominator = matrix(0, max_x, ncol(methylation1))

# Kernel
kernel_width = 2 * derived_sigma + 1

kernel = 0:(kernel_width - 1)
kernel = kernel - floor(kernel_width / 2)
kernel = kernel^2 / (derived_sigma * derived_sigma)
kernel = exp( -0.5 * kernel ) # gaussian

for (i in 1:max_x) {
    mask = (xi == i)
    numerator[i,]   = apply(methylation1[mask,], 2, sum)
    denominator[i,] = apply(coverage1[mask,], 2, sum)
}

#
# Centered convolution
#
# Performs a convolution and returns the the center part that is the same size
# as the original input. This is similar to using the 'same' option for the
# conv function in Octave.
#
# See: http://www.inside-r.org/packages/cran/signal/docs/conv
#
conv_same = function(a, b) {
    x = conv(a, b)
    offset = (length(x) - length(a)) / 2
    return(x[(1 + ceiling(offset)):(length(x) - floor(offset))])
}

# Convolve kernel
yi              = matrix(0, nrow(methylation1), ncol(methylation1))
smoothed_signal = matrix(0, max_x, ncol(methylation1))

for (col in 1:ncol(methylation1)) {
    numerator[,col]   = conv_same(numerator[,col], kernel)
    denominator[,col] = conv_same(denominator[,col], kernel)
    
    # Work-around 2013/04/29
    # conv currently has a bug relating to precision when convolving over a
    # range of zeros. As a temporary work-around, any values very close to 0
    # will be set to 0.
    numerator[,col][numerator[,col] <= 1E-10] = 0
    denominator[,col][denominator[,col] <= 1E-10] = 0
        
    mask = (denominator[,col] == 0)
    
    numerator[mask,col]   = 0
    denominator[mask,col] = 1
    
    smoothed_signal[,col] = numerator[,col] / denominator[,col]

    yi[,col] = interp1(1:max_x, smoothed_signal[,col], 
                       (cpg0x2Dsites1 / sampling_d) +1)
}

# plot coarse scale
dat = stack(as.data.frame(smoothed_signal))
dat$x <- rep(seq_len(nrow(smoothed_signal)), ncol(smoothed_signal))
ggplot(dat, aes(x, values)) + geom_line(aes(colour=ind))

# plot interpolated result
# plot(yi)

