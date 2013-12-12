#
# Generate plots for smooth.ratio paper
# Keith Hughitt <khughitt@umd.edu>
# 2013/05/15
#
# @NOTE: Eventually this should be combined with rest of manuscript using
# knitr + LaTeX.
#
library(bsseq)
library(locfit)
library(reshape2)
library(ggplot2)
source('R/smooth.r')

# Raw data
load('data/sample.rda')

indices = which((cpgsites > 830000) & (cpgsites < 850000))

cpgsites = cpgsites[indices]
methylation = methylation[indices,]
coverage = coverage[indices,]

#################################################
# PLOT1: smooth.ratio: method=downsampling
#################################################
result1 = fast.smooth(cpgsites, methylation, coverage, 1000)
plot(result1, title='smooth.ratio (Method #1: down-sampling)')
ggsave("figure1.eps")

#################################################
# PLOT2: smooth.ratio: method=downsampling
#################################################

# Load Khoa's output (using output directly from C++ version of code 
# until a bug in the R port can be worked out...)
# awk '{$1=$1}1' smooth1.txt > smooth1_stripped.txt
result2 = read.csv(file='~/smooth1_70_stripped.txt', sep=' ', header=FALSE)
result2 = result2[indices,]

# raw data
df1 = data.frame(x=cpgsites, 
                 y=methylation[,1] / coverage[,1],
                 weights=coverage[,1])

stacked = cbind(data.frame(x=cpgsites), y=result2)
df2 = melt(stacked, id='x')

# construct plot
ggplot(df1, aes(x=x, y=y)) + geom_point(color="#5A5A5A", aes(size=weights)) +
    scale_size_continuous(range=c(1,7)) +
    geom_line(data=df2, aes(x=x, y=value, color=variable)) +
    ylim(0, 1) +
    xlab("CpG site (nt)") +
    ylab("% Methylation") +
    ggtitle('smooth.ratio (Method #2: convex combination)')
ggsave("figure2.eps")


#################################################
# PLOT3: smooth.ratio vs. BSmooth vs. locfit
#################################################
bsmooth_data = BSseq(M=methylation, Cov=coverage,
                     chr=c('chr1'), pos=1:nrow(methylation))

bsmooth_result = BSmooth(bsmooth_data, mc.cores=4, ns=70, h=1000, verbose=TRUE)
result3 = bsmooth_result@trans(bsmooth_result@coef)

# locfit
#result4 = fitted.values(locfit(y~lp(x, nn=0.08), data=df1))
result4 = fitted.values(locfit(y~lp(x, nn=0.7), data=df1))

# Timing
#for (time_amt in c(1, 0.5)) {
percentage = methylation / coverage
for (time_amt in c(10, 5, 1, 0.7, 0.5, 0.3, 0.1, 0.05)) {
    t1 = Sys.time()
    for (i in 1:ncol(methylation)) {
        x=fitted.values(locfit(percentage[,i] ~ lp(cpgsites, nn=time_amt)))
    }
    t2 = Sys.time()
    print(sprintf("Time elapsed (nn=%f): %0.2f", time_amt, t2 - t1))
}


# @WORKAROUND: output currently not the same size as input!
#size_diff = nrow(df1) - length(result4)
#pad_left = rep(0.5, ceiling(size_diff / 2))
#pad_right = rep(0.5, floor(size_diff / 2))
pad_left = rep(0.5, 3)
pad_right = rep(0.5, 12)
result4 = c(pad_left, result4, pad_right)

# Combine smoothed output for first sample of all three approaches
combined_fits = cbind(result1@fitted[,1], result2[,1], result3[,1], result4)

# smoothed curves
stacked = cbind(data.frame(x=cpgsites), y=combined_fits)
names(stacked) = c('x', 'smooth.ratio #1', 'smooth.ratio #2', 'BSmooth', 'locfit')
df2 = melt(stacked, id='x')

# construct plot
ggplot(df1, aes(x=x, y=y)) + geom_point(color="#5A5A5A", aes(size=weights)) +
    scale_size_continuous(range=c(1,7)) +
    geom_line(data=df2, aes(x=x, y=value, color=variable)) +
    ylim(0, 1) +
    xlab("CpG site (nt)") +
    ylab("% Methylation") +
    ggtitle('Multiple smoothing method comparison')
ggsave("figure3.eps")
