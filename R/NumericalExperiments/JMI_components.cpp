#include <Rcpp.h> 
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericVector JMI_components (IntegerVector X, IntegerVector Y, IntegerMatrix Z) {
  
  int n = X.size(), s = Z.ncol();     // n - No of observations, s - no of conditioning variables
  int bins = 2;                       // bins - No of distinct values of the predictors
  double prob_x[bins][bins*s], prob_y[2][bins*s], prob_xy[2*bins][bins*s], prob_z[bins*s];
                                      // vectors/tables of probabilities of X|Z, Y|Z, (X, Y)|Z, Z
  
  NumericVector result(n);            // sum_i log p(x, y, z_i)*p(z_i)/p(x, z_i)*p(y, z_i), i - row (observation) index
  
  for(int i = 0; i < bins*s; i++){
    for(int j = 0; j < bins; j++){
      prob_x[j][i] = 0;
      prob_xy[j][i] = 0;
      prob_xy[bins + j][i] = 0;
    }
    prob_z[i] = 0;
    prob_y[0][i] = 0;
    prob_y[1][i] = 0;
  }
  
  // counts in tables, for example:
  // prob_x - rows:    P(X=0 or 1, Z_i=.)
  // prob_x - columns: bins (2) first are for Z_1=0, Z_1=1, ..., Z_1=bins-1, 
  //                   bins (2) next are for Z_2=0, Z_2=1, ..., Z_2=bins-1 and so on (Z_s is the last one).
  
  for(int i = 0; i < n; i++){
    int z = 0;
    for(int j = 0; j < s; j++){
      z = bins*j + Z(i, j);
      prob_x[X[i]][z] += 1;
      prob_y[Y[i]][z] += 1;
      prob_xy[X[i] + bins*Y[i]][z] += 1;
      prob_z[z] += 1;
    }
  }
  
  double value_tmp = 0;
  for(int i = 0; i < n; i++){
    value_tmp = 0;
    for(int j = 0; j < s; j++){
      value_tmp += log((prob_xy[X[i] + bins*Y[i]][j*bins + Z(i,j)] * prob_z[j*bins + Z(i,j)]) / (prob_x[X[i]][j*bins + Z(i,j)] * prob_y[Y[i]][j*bins + Z(i,j)]));
    }
    result[i] = value_tmp;
  }

  return result;
}
