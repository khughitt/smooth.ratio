#' Fast smoothing functions for R
#'
#' This file contains smoothing methods for R.
#'
#' @TODO: fix trailing zero issue for subset input
#' @TODO: have SmoothedData return smoothed data vector when cast to vector.
#'

#' Approximate bilateral smoothing
#'
#' @author Keith Hughitt, based on a Matlab version written by Chengxi Ye.
#'
#' @param x An mx1 vector of predictors.
#' @param y An mxn Matrix of response variables. If more than one column is
#'          specified, each additional column is treated as a separate sample.
#' @param weights    An mxn weight matrix for the above response matrix.
#' @param sigma_d    Standard deviation for Gaussian kernel approximation
#' @param sampling_d Factor to divide x range by to determine number
#'                   of bins to use for down-sampling.
#' @return SmoothedData instance
#'
smooth.ratio = function(x, y, weights,
                        sigma_d=max(round((max(x) - min(x)) / 1e5), 100)) {
    require(signal)

    # Sampling window
    # 2013/05/16: hard-coding until convolution can be generalized
    # to alternative sizes.
    sampling_d=sigma_d

    # Convert any dataframe input to matrices
    x = as.matrix(x)
    y = as.matrix(y)
    weights = as.matrix(weights)

    # Down-sampling parameter
    derived_sigma = sigma_d / sampling_d

    # Assign each row in the input data to a bin
    xi    = round((x - min(x)) / sampling_d) + 1
    max_x = max(xi)

    # Generate kernel (Gaussian approximation)
    kernel_width = 2 * derived_sigma + 1

    kernel = 0:(kernel_width - 1)
    kernel = kernel - floor(kernel_width / 2)
    kernel = kernel^2 / (derived_sigma * derived_sigma)
    kernel = exp(-0.5 * kernel)

    # Down-sample data
    numerator   = matrix(0, max_x, ncol(y))
    denominator = matrix(0, max_x, ncol(weights))

    # Sum data in bins
    index_start = 1

    for (index_end in 1:nrow(xi)) {
        v1 = xi[index_start]
        v2 = xi[index_end]

        if (v1 != v2) {
            indices = index_start:(index_end - 1)

            if (length(indices) == 1) {
                numerator[v1,]   = y[indices,]
                denominator[v1,] = weights[indices,]
            } else if (ncol(y) == 1) {
                numerator[v1,]   = sum(y[indices,])
                denominator[v1,] = sum(weights[indices,])
            } else {
                numerator[v1,]   = colSums(y[indices,])
                denominator[v1,] = colSums(weights[indices,])
            }
            index_start = index_end
        }
    }

    # last row
    indices = index_start:index_end

    if (length(indices) == 1) {
        numerator[v1,]   = y[indices,]
        denominator[v1,] = weights[indices,]
    } else if (ncol(y) == 1) {
        numerator[v1,]   = sum(y[indices,])
        denominator[v1,] = sum(weights[indices,])
    } else {
        numerator[v1,]   = colSums(y[indices,])
        denominator[v1,] = colSums(weights[indices,])
    }

    # Cleaner but slower way of binning the data...
    #     for (i in 1:nrow(xi)) {
    #         numerator[xi[i],]   = numerator[xi[i],] + y[i,]
    #         denominator[xi[i],] = denominator[xi[i],] + weights[i,]
    #     }

    # Instantiate matrices to hold smooth down-sampled and interpolated values
    # smoothed_signal will contain the smooth down-sampled matrix while yi
    # will be used to store the final version which has been interpolated back
    # up to its original size.
    smoothed = matrix(0, max_x, ncol(y))
    yi = matrix(0, nrow(y), ncol(y))

    numerator   = conv_fast(numerator, kernel)
    denominator = conv_fast(denominator, kernel)

    mask = (denominator == 0)
    numerator[mask]   = 0
    denominator[mask] = 1

    smoothed = numerator/denominator

    for (col in 1:ncol(yi)) {
        # Interpolate back up to the original vector length
        yi[,col] = interp1(1:max_x, smoothed[,col],
                          ((x - min(x)) / sampling_d) +1)
    }

    return(new("SmoothedData", x=x, y=y, weights=weights, fitted=yi))
}

#' <TITLE>
#'
#' <DESCRIPTION>
#'
#' @author Mahfuza Sharmin, based on a C++ version written by Khoa Trinh.
#'
#' @param x An mx1 vector of predictors.
#' @param y An mxn Matrix of response variables. If more than one column is
#'          specified, each additional column is treated as a separate sample.
#' @param weights    An mxn weight matrix for the above response matrix.
#' @param window     ...
#' @param a          ...
#' @param b          ...
#'
#' @return SmoothedData instance
#'
smooth.ratio2 = function(x, y, weights, window=70, a=0.5, b=0.5) {
    require(Rcpp)
    require(ggplot2)

    # Load C++ code
    sourceCpp("../src/smooth.cpp")

    # Convert any dataframe input to matrices
    x       = as.matrix(x)
    y       = as.matrix(y)
    weights = as.matrix(weights)

    d_weights = weights
    d_weights[which(d_weights == 0)] = 1

    # Pre-procesing
    M = apply(y, 2, cumsum)
    C = apply(weights, 2, cumsum)
    S = apply(y / d_weights, 2, cumsum)

    max_weights = apply(weights, 2, max)

    smoothed_data = smoothing(window, a, b, d_weights, max_weights,
                              x, y, weights, M, C, S)

    return(new("SmoothedData", x=x, y=y, weights=weights, fitted=smoothed_data))
}

# faster convolution
# a = mxn matrix
# b = 1x3 column vector
# LIMITATION: sampling_d must equal sigma_d, otherwise kernel_width != 3
conv_fast = function(a, b) {
    require(matrixcalc)

    x1 = b[2] * a
    if (ncol(a) != 1) {
        x2 = b[1] * shift.up(a, 1)
        x3 = b[3] * shift.down(a, 1)
    } else {
        x2 = b[1] * matrix(c(a[2:length(a)], 0))
        x3 = b[3] * matrix(c(0, a[1:(length(a) -1)]))
    }
    return(x1 + x2 + x3)
}

#'
#' SmoothedData class definition
#'
setClass('SmoothedData', representation(x='matrix', y='matrix',
                                        weights='matrix', fitted='matrix'))
setMethod('plot', 'SmoothedData', function(object, x, y,
                                           columns=1:ncol(object@fitted),
                                           locfit=FALSE,
                                           title='Smoothed data fit') {
    require(reshape2)

    # raw data
    df1 = data.frame(x=object@x,
                     y=object@y[,1] / object@weights[,1],
                     weights=object@weights[,1])

    # smoothed curves
    stacked = cbind(data.frame(x=object@x), y=object@fitted[,columns])
    df2 = melt(stacked, id='x')

    # optional comparison with locfit
    if (locfit == TRUE) {
        require(locfit)
        locplt = stat_smooth(method='locfit', data=df1, formula=y~lp(x, nn=0.08))
    } else {
        locplt = NULL
    }

    # construct plot
    ggplot(df1, aes(x=x, y=y)) + geom_point(color="#5A5A5A", aes(size=weights)) +
        scale_size_continuous(range=c(1,7)) +
        geom_line(data=df2, aes(x=x, y=value, color=variable)) +
        locplt +
        ylim(0, 1) +
        xlab("CpG site (nt)") +
        ylab("% Methylation") +
        ggtitle(title)
})
