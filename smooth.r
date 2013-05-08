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
                       sigma_d=(max(spatial) - min(spatial)) / 1000,
                       sampling_d=sigma_d) {

    t0 = Sys.time()
    
    # Convert any dataframe input to matrices
    spatial    = as.matrix(spatial)
    intensity  = as.matrix(intensity)
    confidence = as.matrix(confidence)
    
    # Down-sampling parameter
    derived_sigma = sigma_d / sampling_d
    
    t1 = Sys.time()
    print(sprintf("t1: %f", t1-t0))
    
    # Bin data for down-sampling
    xi    = round((spatial - min(spatial)) / sampling_d) + 1
    max_x = max(xi)
    
    t2 = Sys.time()
    print(sprintf("t2: %f", t2-t1))

    # Generate kernel (Gaussian approximation)
    kernel_width = 2 * derived_sigma + 1

    kernel = 0:(kernel_width - 1)
    kernel = kernel - floor(kernel_width / 2)
    kernel = kernel^2 / (derived_sigma * derived_sigma)
    kernel = exp(-0.5 * kernel)

    # Down-sample data
    numerator   = matrix(0, max_x, ncol(intensity))
    denominator = matrix(0, max_x, ncol(confidence))
    
    t3 = Sys.time()
    print(sprintf("t3: %f", t3-t2))

    # @TODO: Would it be possible to collapse repeats in xi to pairs of 
    # coefficients to use in below calculation?
    for (i in 1:nrow(xi)) {
        numerator[xi[i],]   = numerator[xi[i],] + intensity[i,]
        denominator[xi[i],] = denominator[xi[i],] + confidence[i,]
    }
    
    t4 = Sys.time()
    print(sprintf("t4: %f", t4-t3))

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
                           ((spatial - min(spatial)) / sampling_d) +1)
    }
    
    t5 = Sys.time()
    print(sprintf("t5: %f", t5-t4))

    #return(list(y=yi, small=smoothed_signal))
    return(new("SmoothedData", spatial=spatial, intensity=intensity, confidence=confidence, smoothed=yi))
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
