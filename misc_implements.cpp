#include "misc_implements.h"

#include<math.h>
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>

Date::Date(){ day = 0; month = 0; year = 0; }

Date::Date(std::string string_form){
    // Assume that the string has the format mm.dd.yyyy
    month = std::stod(string_form.substr(0, 2));
    day = std::stod(string_form.substr(3, 2));
    year = std::stod(string_form.substr(6, 4));
}

void Date::print_date(){
    std::cout << month << "/" << day << "/" << year << "\n";
}

/*Overloading and improving functionality of Date class*/
bool operator==(const Date& lhs, const Date& rhs){
    if ((rhs.day == lhs.day) && (rhs.month == lhs.month) && (rhs.year == lhs.year)){
        return true;
    }
    else
        return false;
}

// Function that returns the difference in time (in years) between two dates
// Later date must go as the second argument. (End - Beginning) / 365 = time difference in years
double get_time_diff(const Date& begin, const Date& end){
    return ((end.year * 365 + end.month * 30 + end.day) - (begin.year * 365 + begin.month * 30 + begin.day)) / 365.00;
}

std::ostream& operator<<(std::ostream& os, Date t){
    os << t.month << "/" << t.day << "/" << t.year << " ";
    return os;
}

/* VECTOR RELATED FUNCTIONS BELOW */
// Since all the classes and member functions use vectors for computational and storage purposes, their functionalities are improved by the following implementations
// The function names are self explanatory, and will mostly rely on lambda functions and generic algorithms for implementations

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
    std::for_each(v.begin(), v.end(), [](const T& a){os << a << " "; });
    return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v){
    std::for_each(v.begin(), v.end(), [&](const double& a)->void{os << a << " "; });
    return os;
}

template<typename T>
// print the first n elements of v, or last n elements of v (if n is negative)
void print_vec(const std::vector<T>& v, int n){
    if (n >= 0){ // From the front
        for (size_t i = 0; i < n; i++){
            std::cout << v[i] << " ";
        }
    }
    if (n < 0){ // From the back
        for (size_t i = v.size() + n; i < v.size(); i++){
            std::cout << v[i] << " ";
        }
    }
}

double find_arithmatic_average(const std::vector<double>& v){
    double total = 0;
    std::for_each(v.begin(), v.end(), [&](double a){total += a; });
    return (total / v.size());
}

double find_geometric_average(const std::vector<double>& v){
    double product = 1;
    std::for_each(v.begin(), v.end(), [&](double a){product *= a; });
    return pow(product, 1.0 / v.size());
}

double find_max(const std::vector<double>& v){
    double max = v[0];
    std::for_each(v.begin(), v.end(), [&](double a){(a > max) ? max = a : max = max; });
    return max;
}

double find_min(const std::vector<double>& v){
    double min = v[0];
    std::for_each(v.begin(), v.end(), [&](double a){(a < min) ? min = a : min = min; });
    return min;
}

template<typename T>
bool search_vector(const std::vector<T>& v, T element){
    for (const auto& i : v){
        if (i == element){
            return true;
        }
    }
    return false;
}

/* SUPPLEMENTAL FUNCTIONS IMPLEMENTED BELOW*/

// Recursive factorial function
inline int factorial(int n){
    return (n == 0) ? 1 : n*(factorial(n - 1));
}

/* IMPLEMENTATION OF SOME PROBABILITY DISTRIBUTIONS AND TOOLS */
// Implementation of binomial probability mass function with parameters n, p, k
double bin(const int& n, const double& p, const int& k){
    return (factorial(n) / (factorial(n - k)*factorial(k)))*pow(p, k)*pow(1 - p, n - k);
}

// Implementation of the standard Gaussian CDF
double normalCDF(double value){
    return 0.5 *(1 + erf(value / sqrt(2)));
}

// Numerical approximation of the inverse CDF of the standard Gaussian Distribution
/* From Stegun and Abramowitz's Handbook of Mathematical Functions: The inverse CCDF can be approximated for values of p < 0.5 by the solving
a polynomial, but due to the symmetry of the standard Gaussian, flipping the sign of the result produces the inv CDF. To handle cases where p > 0.5
symmetry allows us to take the complement of p, (that is, 1-p) and apply the same technique, but this time do not flip the sign.*/
// First implement the approximation method for the inverse CCDF
double inv_CCDF(const double& p){
    double t = sqrt(-2.0*log(p));
    const double c0 = 2.511516;
    const double c1 = .802853;
    const double c2 = .010328;
    const double d1 = 1.432788;
    const double d2 = .189269;
    const double d3 = .001308;
    return t - ((c0 + c1*t + c2*t*t) / (1 + d1*t + d2*t*t + d3*t*t*t));
}

double norm_inv(const double& p){
    if (p <= 0.5){
        return -inv_CCDF(p);
    }
    else
        return inv_CCDF(1 - p);
}
