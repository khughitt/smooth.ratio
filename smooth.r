#' Fast smoothing functions for R
#' 
#' This file contains smoothing methods for R.
#' 
#' @author Keith Hughitt, based on a Matlab version written by Chengxi Ye.
library(signal)
library(matrixcalc)
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
                       sigma_d=(max(spatial) - min(spatial)) / 1000,
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

    # Sum data in bins
    for (i in 1:nrow(xi)) {
        numerator[xi[i],]   = numerator[xi[i],] + intensity[i,]
        denominator[xi[i],] = denominator[xi[i],] + confidence[i,]
    }
    
    # Instantiate matrices to hold smooth down-sampled and interpolated values
    # smoothed_signal will contain the smooth down-sampled matrix while yi
    # will be used to store the final version which has been interpolated back
    # up to its original size.
    smoothed_signal = matrix(0, max_x, ncol(intensity))
    yi = matrix(0, nrow(intensity), ncol(intensity))

    numerator = conv_fast(numerator, kernel)
    denominator = conv_fast(denominator, kernel)

    mask = denominator == 0
    numerator[mask]   = 0
    denominator[mask] = 1

    smoothed_signal = numerator/denominator
    
    for (col in 1:ncol(yi)) {
        # Interpolate back up to the original vector length
        yi[,col] = interp1(1:max_x, smoothed_signal[,col], 
                           ((spatial - min(spatial)) / sampling_d) +1)
    }

    return(new("SmoothedData", spatial=spatial, intensity=intensity, 
               confidence=confidence, smoothed=yi))
}

# faster convolution
# a = mxn matrix
# b = 1x3 column vector
# LIMITATION: sampling_d must equal sigma_d, otherwise kernel_width != 3
conv_fast = function(a, b) {
    x1 = b[2] * a
    x2 = b[1] * shift.up(a, 1)
    x3 = b[3] * shift.down(a, 1)
    return(x1 + x2 + x3)
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

#'
#' SmoothedData class definition
#'
setClass('SmoothedData', representation(spatial='matrix', intensity='matrix', 
                                        confidence='matrix', smoothed='matrix'))
setMethod('plot', 'SmoothedData', function(object, x, y, col) {
    #col=1
    indices = which((object@spatial > 830000) & (object@spatial <= 850000))
    
    # create dataframes for ggplot
    df1 = data.frame(x=object@spatial[indices], 
                     y=object@intensity[indices,col] / object@confidence[indices,col])
    
    # first column
    df2 = data.frame(x=object@spatial[indices], y=object@smoothed[indices,col])
    
    ggplot(df1, aes(x=x, y=y)) + geom_point(color="#5A5A5A") +
        geom_line(data=df2, aes(x=x, y=y), color='red') +
        xlab("CpG site (nt)") +
        ylab("% Methylation") +
        ggtitle("Smooted data fit")
    #ggsave('output_sample1.png', width=12, height=9, dpi=96)
})
