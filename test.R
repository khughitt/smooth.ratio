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

# load data
data = readMat('test.mat')
attach(data)

# Domain sigma
sigma_d = max(cpg0x2Dsites1) / 100
sampling_d = sigma_d

derived_sigma = sigma_d / sampling_d

xi = round(cpg0x2Dsites1 / sampling_d) + 1
max_x = max(xi)
numerator = matrix(0, 1, max_x)
denominator = matrix(0, 1, max(xi))

# Kernel
kernel_width = 2 * derived_sigma + 1

halfKernelWidth = floor( kernel_width / 2 )
kernel = 0:(kernel_width - 1)
kernel = kernel - halfKernelWidth;
kernel = kernel^2 / (derived_sigma * derived_sigma )

kernel = exp( -0.5 * kernel ) # gaussian

#############################################################
# NOT FINISHED PORTING -- Code below is incomplete...
# Keith 2013/04/29
#############################################################
for (i in 1:max_x) {
    mask = (xi == i)
    numerator[i] = apply(methylation1[mask,], 2, sum)
    denominator[i] = apply(coverage1[mask,], 2, sum)
}

# convolve kernel
numerator = conv(numerator, kernel)
denominator = conv(denominator, kernel)
numerator[denominator == 0] = 0
denominator[denominator == 0] = 1

smoothed_signal = numerator / denominator

# interpolate and plot
yi = interp1(1:max_x, smoothed_signal, (cpg0x2Dsites1 / sampling_d) +1)

# coarse scale
plot(smoothed_signal)

# interpolated result
plot(yi)
