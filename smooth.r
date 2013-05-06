#' Fast smoothing functions for R
#' 
#' This file contains smoothing methods for R.
#' 
#' @author Keith Hughitt, based on a Matlab version written by Chengxi Ye.
library(signal)
library(ggplot2)

#' Approximate bilateral smoothing
#' 
#' @param spatial    mx1 Spatial component of matrix to be smoothed
#' @param intensity  mxn Intensity component of matrix to be smoothed
#' @param confidence mxn Confidence measure for vector to be smooted
#' @param sigma_d    Standard deviation for Gaussian kernel approximation
#' @param sampling_d Factor to divide spatial range by to determine number
#'                   of bins to use for down-sampling.
#' @return list containing smoothed matrix along with intermediate downsampled
#'         version of the matrix.
fast.smooth = function(spatial, intensity, confidence, 
                       sigma_d=(max(spatial) - min(spatial)) / 10000,
                       sampling_d=sigma_d) {

    # Convert any dataframe input to matrices
    spatial    = as.matrix(spatial)
    intensity  = as.matrix(intensity)
    confidence = as.matrix(confidence)
    
    # Down-sampling parameter
    derived_sigma = sigma_d / sampling_d
    
    # Bin data for down-sampling
    xi    = round((spatial - min(spatial)) / sampling_d) + 1
    max_x = max(xi)

    # Generate kernel (Gaussian approximation)
    kernel_width = 2 * derived_sigma + 1

    kernel = 0:(kernel_width - 1)
    kernel = kernel - floor(kernel_width / 2)
    kernel = kernel^2 / (derived_sigma * derived_sigma)
    kernel = exp(-0.5 * kernel)

    # Down-sample data
    numerator   = matrix(0, max_x, ncol(intensity))
    denominator = matrix(0, max_x, ncol(confidence))

    # @TODO: Would it be possible to collapse repeats in xi to pairs of 
    # coefficients to use in below calculation?
    for (i in 1:nrow(intensity)) {
        numerator[xi[i],]   = numerator[xi[i],] + intensity[i,]
        denominator[xi[i],] = denominator[xi[i],] + confidence[i,]
    }

    #numerator = aggregate(intensity, list(xi[[1]]), "sum")
    #denominator = aggregate(confidence, list(xi[[1]]), "sum")

    # Instantiate matrices to hold smooth down-sampled and interpolated values
    # smoothed_signal will contain the smooth down-sampled matrix while yi
    # will be used to store the final version which has been interpolated back
    # up to its original size.
    smoothed_signal = matrix(0, max_x, ncol(intensity))
    yi = matrix(0, nrow(intensity), ncol(intensity))

    # Iterate through samples and convolve down-sampled data
    for (col in 1:ncol(intensity)) {
        # Perform convolution
        numerator[,col]   = conv_same(numerator[,col], kernel)
        denominator[,col] = conv_same(denominator[,col], kernel)

        # Work-around 2013/04/29
        # conv currently has a bug relating to precision when convolving over a
        # range of zeros. As a temporary work-around, any values very close to 0
        # will be set to 0.
        numerator[,col][numerator[,col] <= 1E-10] = 0
        denominator[,col][denominator[,col] <= 1E-10] = 1

        mask = which(denominator[,col] == 0)
        
        numerator[mask,col]   = 0
        denominator[mask,col] = 1
        
        # Scale the smoothed result by the smoothed confidence vector
        smoothed_signal[,col] = numerator[,col] / denominator[,col]

        # interpolate back up to the original vector length
        yi[,col] = interp1(1:max_x, smoothed_signal[,col], 
                           (spatial / sampling_d) +1)
    }

    return(list(y=yi, small=smoothed_signal))
}

#' Centered convolution
#'
#' Performs a convolution and returns the the center part that is the same size
#' as the original input. This is similar to using the 'same' option for the
#' conv function in Octave.
#' 
#' @param a First numeric sequence to be convolved
#' @param b Second numeric sequence to be convolved
#' @return Vector resulting from convolution of a and b of the same length
#'         as input vectors.
#' @seealso: \url{http://www.inside-r.org/packages/cran/signal/docs/conv}
conv_same = function(a, b) {
    x = conv(a, b)
    offset = (length(x) - length(a)) / 2
    return(x[(1 + ceiling(offset)):(length(x) - floor(offset))])
}
