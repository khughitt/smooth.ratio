#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector smoothing(int window, double a, double b, NumericMatrix d_coverage, NumericVector maxCoverage,
                            NumericVector cpgsites, NumericMatrix methylation, NumericMatrix coverage,
                            NumericMatrix M, NumericMatrix C, NumericMatrix S) {
    int n = cpgsites.size();
    int num_samples = methylation.ncol();
    NumericMatrix smoothed_data(n, num_samples);
    for(int j=0;j<num_samples; j++){
      double L_denominator = 0, R_denominator = 0;
      double L_numerator = 0, R_numerator = 0;
      for (int i=1; i<=window; i++){
        R_denominator += cpgsites[i] - cpgsites[0];
        R_numerator += 1.0 * (cpgsites[i] - cpgsites[0]) * (1.0*methylation(i,j) / d_coverage(i,j));
      }
          
      for (int i=0; i<n; i++){
        int left = 0;// = max(0, i - window);
        int right = n-1;// = min(n-1, i + window);
        if((i-window)>0){
          left = i - window;
        }
        if((i+window)<(n-1)){
          right = i + window;
        }

          
        double D = cpgsites[right] - cpgsites[left];
        double S_c = C(right,j) - C(left,j) + coverage(left,j);
          
        double denominator = a*S_c/maxCoverage[j] + (b/D)*((right-left+1)*D - (L_denominator+R_denominator));
        
        double numerator = (a/maxCoverage[j]) * (M(right,j) - M(left,j) + methylation(left,j))
                            + b * (S(right,j) - S(left,j) + 1.0*methylation(left,j)/d_coverage(left,j))
                            - (b / D) * (L_numerator + R_numerator);
          
        smoothed_data(i,j) = numerator / denominator;
          
          
        // update
        if (i+1<n){
          double delta = cpgsites[i+1] - cpgsites[i];
          
          if (right+1 < n){
            R_denominator = R_denominator + (cpgsites[right+1] - cpgsites[i+1]);
          }
          R_denominator = R_denominator - 1.0*(right - i)*delta;
          
          if (i - left == window){
            L_denominator = L_denominator - (cpgsites[i] - cpgsites[left]);
          }
          else{
            L_denominator += delta;
          }
          L_denominator = L_denominator + 1.0*(i - left)*delta;
          
          if (right+1 < n){
            R_numerator = R_numerator + (cpgsites[right+1] - cpgsites[i+1])
                            *(1.0*methylation(right+1,j) / d_coverage(right+1,j));
          }
          R_numerator = R_numerator - 1.0 * delta * (S(right,j) - S(i,j));
            
          if (i - left == window){
            L_numerator = L_numerator - (cpgsites[i] - cpgsites[left])
                          *(1.0 * methylation(left,j) / d_coverage(left,j));
          }
          L_numerator += delta * (1.0 * methylation(i,j) / d_coverage(i,j));
          if (i){
            L_numerator = L_numerator + 1.0 * delta * ( S(i-1,j) - ( i<window ? 0:S(left,j)) );
          }
        }
      }
    }
    return smoothed_data;
}
