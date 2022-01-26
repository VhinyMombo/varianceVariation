
#include <Rcpp.h>
using namespace Rcpp;


template<typename T> 
struct square 
{
  T operator()(const T& Left, const T& Right) const
  {
    return (Left + Right*Right);
  }
};


NumericVector subset_range(NumericVector x,
                                 int start , int end ) {

  return x[Rcpp::Range(start, end)];
}



int which_max_cpp(NumericVector v) {
  double m = Rcpp::max(v);
  std::vector< int > res;
  
  int i;
  for( i = 0; i < v.size(); ++i) {
    if( v[i] == m ) {
      res.push_back( i );
    }
  }
  Rcpp::IntegerVector iv( res.begin(), res.end() );
  return iv[0];
}

// [[Rcpp::export]]
NumericVector css_statistic(NumericVector &y){
  NumericVector stat;
  int n = y.size();
  
  for(int t=1; t<y.size()+1; t++){
    float sum1 = 0.0;
    float sum2 = 0.0;
    float c = static_cast< float >(t) / static_cast< float >(n);
    
    sum1 = std::accumulate(y.begin(), y.begin()+t, 0.0,square<float>());
    sum2 = std::accumulate(y.begin(), y.end(), 0.0,square<float>());
    
    float res = sqrt(n/2)*std::abs((sum1/sum2)-c);
    stat.push_back(res);
  }
  return stat;
}

// [[Rcpp::export]]
Rcpp::List BS(int s, int e, NumericVector &y, float penality){
  if(e-s ==1){
    Rcpp::List res;
    return res;
  }
  else{
    NumericVector sub_vec = subset_range(y, s, e-1);
    NumericVector css_stats = css_statistic(sub_vec);
    int cp = which_max_cpp(css_stats)+s;
    
    float m = max(css_stats);

    if(m>penality){
      NumericVector v; 
      v.push_back(cp);
      return Rcpp::List::create(BS(s, cp, y, penality), v, BS(cp-1, e, y, penality));
    }
    else{
      Rcpp::List res;
      return res;
    }
  }
}

