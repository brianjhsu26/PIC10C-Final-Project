// Contains all necessary supplementary classes and operator overloads (for vectors)

#ifndef MISC_H
#define MISC_H
#include<math.h>
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>

class Date{
public:
    // should have find time difference function and output a time in date
    Date();
    Date(std::string string_form);
    void print_date();
    friend bool operator==(const Date& lhs, const Date& rhs);
    friend std::ostream& operator<<(std::ostream& os, Date t);
    friend double get_time_diff(const Date& begin, const Date& end);
private:
    // Private variables represent three times, this allows you to find the difference between dates
    size_t day;
    size_t month;
    size_t year;
};

bool operator==(const Date& lhs, const Date& rhs);
double get_time_diff(const Date& begin, const Date& end);
std::ostream& operator<<(std::ostream& os, Date t);

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);
template<>
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);

template<typename T>
void print_vec(const std::vector<T>& v, int n);

double find_arithmatic_average(const std::vector<double>& v);
double find_geometric_average(const std::vector<double>& v);
double find_max(const std::vector<double>& v);
double find_min(const std::vector<double>& v);

template<typename T>
bool search_vector(const std::vector<T>& v, T element);

inline int factorial(int n);

double bin(const int& n, const double& p, const int& k);
double normalCDF(double value);
double inv_CCDF(const double& p);
double norm_inv(const double& p);

#endif

