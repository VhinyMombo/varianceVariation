#include <iostream>
#include <vector>
#include <math.h>

#include <Rcpp.h> //to use the NumericVector object
using namespace Rcpp; //to use the NumericVector object
using namespace std;



// [[Rcpp::export]]

float penalty_fun(int n, int params){
    double beta;
    beta = (params+1) * log(n);
    return beta;
}



// [[Rcpp::export]]
float variance( std::vector<float> &v, float &mu){
    int n;
    n = v.size();
    double result;
    result = 0;

    for (int i = 0; i < n; i++)
    {
        result += (v[i]*v[i]-mu)/n;
    }

    return result;
}

// [[Rcpp::export]]
float cost_function(std::vector<float> &v){
    int n = v.size();
    float mu = 0;
    float C = 0;
    float var;
    if (n==0){
        C = 0;
    }
    else{
        var = variance(v, mu);
        C = n * (log(2*M_PI) + log(var) + 1);
    }
    return(C);
}

// Template class to slice a std::vector
// from range X to Y
// [[Rcpp::export]]
std::vector<float> slicing(std::vector<float> const& v,
                  int X, int Y)
{

    // Begin and End iterator
    auto first = v.begin() + X;
    auto last = v.begin() + Y + 1;

    // Copy the element
    std::vector<float> vector(first, last);

    // Return the results
    return vector;
}

std::vector<int> slicing(std::vector<int> const& v,
                  int X, int Y)
{

    // Begin and End iterator
    auto first = v.begin() + X;
    auto last = v.begin() + Y +1;

    // Copy the element
    std::vector<int> vector(first, last);

    // Return the results
    return vector;
}
// [[Rcpp::export]]
int argmin(std::vector<float> &v){
 int res = std::min_element(v.begin(), v.end()) - v.begin();
 return (res);
}
/*
template <typename T>
void printResult(std::vector<T> const& v)
{

    // Traverse the vector v
    for (auto i : v) {
        cout << i << ' ';
    }
    cout << '\n';
}
 */

void print(std::vector<float> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        std::cout << input.at(i) << ' ';
    }
    std::cout << '\n';
}


// [[Rcpp::export]]
std::vector<int> OP_cpp(std::vector<float> &data, int params = 1){
    int n = data.size();
    float beta;
    float F_temp;
    int i_temp = 0;
    float F_tau;
    int v;
    std::vector<int> P;
    std::vector<float> y_0t;
    std::vector<float> y_st;
    std::vector<int> cp(n,0);
    std::vector<float> F_comp(n,0);
    std::vector<int>::iterator ip;


    beta = penalty_fun(n,params);
    F_comp[0] = -beta ;
    for (int tau = 1; tau < n ; tau++)
    {
        y_0t = slicing(data,0,tau);
        F_temp = cost_function(y_0t);
        i_temp = 0;
        for (int i = 1; i < tau; i++)
        {
            y_st = slicing(data,i,tau-1);
            F_tau = F_comp[i-1] + cost_function(y_st) + beta;

            if (F_tau < F_temp)
            {
                F_temp = F_tau;
                i_temp = i-1;
            }

        }
        F_comp[tau] = F_temp;
        cp[tau] = i_temp;

    }

    v  = cp[n-1];
    P.push_back(cp[n-1]);

    while (v > 0)
    {
        P.push_back(cp[v]);
        v = cp[v];
    }
    std::reverse(P.begin(), P.end());
    //* ip = std::unique(cp.begin(),cp.begin() + n);
    //* cp.resize
    P = slicing(P,1,P.size()-1);
    return (P);
}

// [[Rcpp::export]]
std::vector<int> PELT_cpp(std::vector<float> &data, int params = 1, int K = 0){
    int n = data.size();
    float beta;
    float F_temp;
    int i_temp = 0;
    float F_tau;
    int v;
    std::vector<int> P;
    std::vector<float> y_0t;
    std::vector<float> y_st;
    std::vector<int> cp(n,0);
    std::vector<float> F_comp(n,0);
    std::vector<int>::iterator ip;
    std::vector<int> R;
    std::vector<float> F_tau_K;
    std::vector<int> R_new;
    beta = penalty_fun(n,params);
    F_comp[0] = -beta ;
    R = {0};
    for (int tau = 1; tau < n ; tau++)
    { 
        y_0t = slicing(data,0,tau);
        F_temp = cost_function(y_0t);
        i_temp = 0;

        for(int k = 0; k < R.size(); k++) 
        {
            int i = R[k];
            y_st = slicing(data,i+1,tau);
            F_tau = F_comp[i] + cost_function(y_st) + beta;
            F_tau_K.push_back(F_comp[i] + cost_function(y_st) + K);

            if (F_tau < F_temp)
            {
                F_temp = F_tau;
                i_temp = i;
            }

        }
        F_comp[tau] = F_temp;
        cp[tau] = i_temp;
        for (int j = 0; j < F_tau_K.size(); j++)
        {
            if(F_tau_K[j] < F_temp){
                R_new.push_back(R[j]);
            }
        }
        F_tau_K.clear();
        R = R_new;
        R.push_back(tau);
        R_new.clear();

    }
    v  = cp[n-1];
    P.push_back(cp[n-1]);

    while (v > 1)
    {
        P.push_back(cp[v]);
        v = cp[v];
    }
    std::reverse(P.begin(), P.end());
    P = slicing(P,0,P.size()-1);
    return (P);
}



/*
int main(){

    int params = 2;
    int n;
    double beta;
    float var;
    float mu = 0.0;
    float C;
    int s = 0;
    int t = 6;
    float res_OP;

    std::vector<float> v = {1,2,3,4,6,45.6,4};
    std::vector<float> v1;

    n = OP(v);
    v1 = slicing(v,s,t);
    beta = penalty_fun(n,params);
    var = variance(v, mu);
    C = cost_function(v);
    res_OP = Optimal_Part(v);



    cout << "size is :" << n;
    cout << "beta is :" << beta;
    printf("%f\n",beta);
    printf(" variance  is %f\n",var);
    printf("cost function is %f\n",C );
    // printResult(v1);
    printf("OP is %f\n",res_OP);

    return 0;

}
*/
