#ifndef MATRIX_H
#define MATRIX_H

#include<vector>
#include<fstream>
#include<regex>
#include<string>
#include<iostream>
#include<math.h>
#include<algorithm>
#include<functional>
#include"misc_implements.h"

static const double E = 2.7182818284590452353602874713527;

class Custom_Matrix{
public:
    Custom_Matrix(){ rows = 0; columns = 0; }
    void write_to_header(const std::vector<std::string>& header);
    inline void write_to_data(const double& data);
    inline void assign_dimensions(const size_t& n, const size_t& m);
    // Overloading () operator to act as an accessor and get the element in the (a,b) position of the matrix
    // IMPORTANT: THIS ACCESOR SHOULD NOT BE USED WITH B = 1, IT WILL RETURN THE DIFFERENCE IN TIME BETWEEN DATE IN VECTOR AND ROOT_DATE
    inline double operator()(const size_t& a, const size_t& b);
    // Implementing getter functions
    // Outputs dimensional information of the matrix
    void dim();
    std::vector<Date> dates;
    std::vector<std::string> column_labels; // Note that the header row does NOT count as the first row
    std::vector<double> matrix_data;
    size_t rows;
    size_t columns;
};

inline void Custom_Matrix::assign_dimensions(const size_t& n, const size_t& m){
    rows = n; columns = m;
}

// For some reason, this function overload cannot be implemented outside
inline double Custom_Matrix::operator()(const size_t& a, const size_t& b){
    if ((b > columns) || (b < 1)){
        std::cout << "Error: Invalid Matrix Index" << "\n";
        return 0;
    }
    if ((a == 1) && (b != 1)){
        return matrix_data[b - 2];
    }
    else if ((a != 1) && (b != 1)){
        return matrix_data[(a - 1)*(columns - 1) + b - 2]; // Note dates are in the FIRST column, so subtract one from the column#
    }
    else if (b == 1){
        //std::cout << dates[a - 1] << "\n";
        Date root_date("01.01.1960");
        return get_time_diff(root_date, dates[a - 1]);
    }
}

class Option_Matrix : public Custom_Matrix{
public:
    Option_Matrix(){
        stock_name = "";
        current_date = Date("11.04.2016");
        current_price = 120.75;
        Custom_Matrix calls();
        Custom_Matrix puts();
    }
    Option_Matrix(std::fstream& call_data, std::fstream& put_data){
        // Read in each data set seperately, beginning with calls
        std::regex date("\\d{2}\.\\d{2}\.\\d{4}");
        std::string x;
        int row_counter = 0;
        while (!call_data.eof()){
            call_data >> x;
            if (x == "Call"){ // Once the option type is described, we assume the next line provides the headers for the option
                row_counter = 0;
                call_data >> x;
                calls.column_labels.push_back(x);
                while ((call_data.peek() != '\n') && (call_data >> x)){ // For the next line
                    calls.column_labels.push_back(x); // Read all the column titles into this vector
                }
                call_data >> x; // Shift the fstream one forward
            }
            //std::cout << "check: " << x << "\n";
            // Now continue to parse the text file, look for instances where the string matches a date, if it does, convert it and store it in the exp_dates vector
            if (std::regex_match(x, date)){
                calls.dates.push_back(Date(x));
                call_data >> x;
            }
            calls.matrix_data.push_back(std::stod(x)); // If it's not a date, then push it into the matrix data
            if (call_data.peek() == '\n'){
                ++row_counter;
            }
        }
        calls.assign_dimensions(row_counter, calls.column_labels.size());
        while (!put_data.eof()){
            put_data >> x;
            if (x == "Put"){ // Repeat the process for put options
                row_counter = 0;
                put_data >> x;
                puts.column_labels.push_back(x);
                while ((put_data.peek() != '\n') && (put_data >> x)){ // For the next line
                    puts.column_labels.push_back(x); // Read all the column titles into this vector
                }
                put_data >> x;
            }
            if (std::regex_match(x, date)){
                puts.dates.push_back(Date(x));
                put_data >> x; // Shift the fstream forward by one or else it's going to read the part after the decimal as one word
            }
            puts.matrix_data.push_back(std::stod(x));
            if (put_data.peek() == '\n'){
                ++row_counter;
            }
        }
        puts.assign_dimensions(row_counter, puts.column_labels.size());
        // Optional feature: do checking on the formatting of the data
    }

    void set_characteristics(const std::string& name, const std::string& date, double price){
        stock_name = name;
        current_date = Date(date);
        current_price = price;
    }

    /* // TEST CODE
    void print_matrix_characteristics(){
        std::cout << "Stock Name: " << stock_name << "\n";
        std::cout << "Call option headers: " << calls.column_labels << std::endl;
        std::cout << "Rows: " << calls.rows << "\n";
        std::cout << "Columns: " << calls.columns << "\n" << std::endl;
        std::cout << "Put option headers: " << puts.column_labels << std::endl;
        std::cout << "Rows: " << puts.rows << "\n";
        std::cout << "Columns: " << puts.columns << "\n";
    } */

    // Prints a list of strike prices in the vector
    std::vector<double> get_strike_price_vector(){
        std::vector<double> v;
        for (size_t i = 1; i < calls.rows-1; i++){
            v.push_back(calls(i, 2));
        }
        v.push_back(calls(calls.rows, 2));
        return v;
    }

    // By put call parity, we know r = (1/-T) * ln[(C-P-S)/-K]
    // Pass in a function pointer so that this function can return either an average (arithmatic), max, or min
    double find_r(std::function<double(std::vector<double>)> func){
        std::vector<double> rfrates;
        for (size_t i = 1; i != calls.rows; i++){
            double time_until_exp = get_time_diff(current_date, calls.dates[i]); // <-- 1 is arbitrary
            rfrates.push_back((1 / -time_until_exp) * log((calls(i, 3) - puts(i, 3) - current_price) / (-calls(i, 2))));
        }
        return func(rfrates);
    }

    // We first need a function that finds the derivative of the call price with respect to volatility (this derivative is called vega)
    // The closed form of vega (as derived from the Black-Scholes equation) is v = S*sqrt(T)*N(d1)
    // This is so that we can continuously plug in values of volatility for Newton-Raphson until it stabilizes, hence it is only used in implied volatility function
    // Note that the vega for a call option is equal to the vega for a put option, so it should be the same whether the option is call or put

    double find_vega(const double& r, double strike, double vol){
        size_t option_index = 1;
        while (calls(option_index, 2) != strike){
            ++option_index;
        }
        double d1 = (log(current_price / calls(option_index, 2)) +
            (r - 0.5*vol*vol)*get_time_diff(current_date, calls.dates[option_index])) /
            (vol*sqrt(get_time_diff(current_date, calls.dates[option_index])));
        return current_price*sqrt(get_time_diff(current_date, calls.dates[option_index]))*normalCDF(d1);
    }

    // Computes the option price with strike K given a volatility (which can be either historical, or an estimate of implied volatility)
    double black_scholes_price(const double& vol, const double& r, const int& strike, std::string option_type){
        // First find the row with the corresponding strike price
        if (option_type == "Call"){
            size_t option_index = 1;
            while (calls(option_index, 2) != strike){
                ++option_index;
            }
            double T = get_time_diff(current_date, calls.dates[option_index]);
            double d1 = (log(current_price / strike) + (r + 0.5*vol*vol)*T) / (vol*sqrt(T));
            double d2 = d1 - vol*sqrt(T);
            /* For error checking...
            std::cout << "risk free rate: " << r << "\n";
            std::cout << "Time until expiration: " << T << "\n";
            std::cout << "d1 = " << d1 << "\n";
            std::cout << "d2 = " << d2 << "\n"; */
            return current_price*normalCDF(d1) - strike*pow(E, -r*T)*normalCDF(d2);
        }
    }

    // Computes the option price based on the binomial tree method
    std::vector<double> binomial_tree_price(const double& vol, const double& r, const int& strike, std::string option_type){ // Input for option type should be "Call" or "Put"
        // Use the Jarrow-Rudd model, where u = e^[(r-0.5vol*vol)*h + vol*sqrt(h)] , d = e^[(r-0.5vol*vol)*h - vol*sqrt(h)]
        // h is determined by the length of each timestep. Get the time difference between expiry date and current day, then h = time_to_exp/timesteps
        size_t timesteps = 1; // timesteps determines the number of nodes we have
        double time_to_exp = get_time_diff(current_date, calls.dates[0]); // time until expiration
        double h = time_to_exp / timesteps;
        std::cout << "exp time: " << time_to_exp << "\n";
        // First compute the probability of an upward movement, commonly denoted as p = (e^(-rt) - D) / (U - D)
        std::vector<double> option_price;
        option_price.push_back(0); // Initialize the vector to zero to allow the while loop to activate
        double U = exp((r)*h + vol*sqrt(h));
        double D = exp((r)*h - vol*sqrt(h));
        double p = (exp(r*h) - D) / (U - D);
        size_t iter = 1;
        double payoff = 0;
        if (option_type == "Call"){
            double call_price = 0;
            call_price = exp(-r*h)*(std::max(U*current_price - strike, 0.00)*p + std::max(D*current_price - strike, 0.00)*(1 - p));
            std::cout << "Call option value for " << timesteps << " timesteps is: " << call_price << "\n";
            option_price.push_back(exp(-r*h)*(std::max(U*current_price - strike, 0.00)*p + std::max(D*current_price - strike, 0.00)*(1 - p))); // Push back the one period binomial tree price
            while ((timesteps < 12) && (abs(option_price[iter] - option_price[iter - 1]) > .001)){ // Repeat until convergence or 1000 timesteps
                ++timesteps;
                ++iter;
                //h = time_to_exp / timesteps; // Update the size of the timestep
                h = time_to_exp / timesteps;
                U = exp((r)*h + vol*sqrt(h)); // Update U to fit the new timestep
                D = exp((r)*h - vol*sqrt(h)); // Update D to fit the new timestep
                p = (exp(r*h) - D) / (U - D); // Update p to fit the new timestep
                payoff = 0; // Reset and calculate the risk-neutral payoffs of each end branch of the binomial tree
                for (size_t k = 0; k <= timesteps; k++){
                    payoff += bin(timesteps, p, k)*std::max(current_price*pow(U, k)*pow(D, timesteps - k) - strike, 0.00);
                }

                call_price = (exp(-r*time_to_exp))*payoff; // Discount the payoffs to present time for present value
                std::cout << "Call option value for " << timesteps << " timesteps is: " << call_price << "\n";
                option_price.push_back(call_price);
            }
            return option_price;
        }
        if (option_type == "Put"){
            return option_price;
        }
    }

    /* To estimate the implied volatility (Ivol), we use Newton-Raphson starting with the initial approximation of volatility of
    v = sqrt(2*pi/T)*(C/S)(Brenner, 1988) with the at the money option. Loop until the error converges or hits a maximum of 100 iterations
    INPUTS: strike price (strike), risk free interest rate (r), accepted error (err) */
    double estimate_Ivolatility(const int& strike, const double& r, const double& err = .0001){
        size_t option_index = 1;
        while (calls(option_index, 2) != strike){
            ++option_index;
        }
        size_t iter = 1;
        std::vector<double> computed_vols;
        computed_vols.push_back(0);
        computed_vols.push_back(sqrt(2 * 3.1415 / get_time_diff(current_date, calls.dates[option_index]))*(calls(option_index, 3) / current_price)); // volatility at time t minus 1
        // Now that there's an initial guess for volatility, we use the Newton-Raphson method to pinpoint it
        while ((iter != 100) && (abs(computed_vols[iter] - computed_vols[iter - 1]) > err)){
            double vega = this->find_vega(r, strike, computed_vols[iter]);
            computed_vols.push_back(computed_vols[iter] - (black_scholes_price(computed_vols[iter], r, strike, "Call") - calls(option_index, 3)) / vega);
            ++iter;
        }
        std::cout << "Estimated implied volatility: " << computed_vols[iter] << "\n";
        std::cout << "Steps taken to convergence: " << iter << "\n";
        std::cout << "Computed volatilities: " << computed_vols << "\n";
        return computed_vols[iter];
    }

    // Return the market price of the call option at a specified strike price
    double get_actual_call_price(double strike){
            size_t i = 1;
            while (calls(i, 2) != strike){
                i++;
            }
            return calls(i, 3);
        }
private:
    std::string stock_name;
    Date current_date = Date("11.04.2016");
    double current_price = 120.75;
    Custom_Matrix calls;
    Custom_Matrix puts;
};



class Stock_Price_Matrix : public Custom_Matrix{
public:
    Stock_Price_Matrix(){
        earliest_date = Date("1.1.1990");
        latest_date = Date("1.1.1990");
    }
    Stock_Price_Matrix(std::fstream& price_data){
        size_t row_counter = 0;
        std::regex date("\\d{2}\.\\d{2}\.\\d{4}");
        std::string x;
        //Assuming the first line of the data file has the headers
        while ((price_data.peek() != '\n') && (price_data >> x)){
            this->column_labels.push_back(x);
        }

        price_data >> x; //Push the fstream by one word
        while (!price_data.eof()){
            if (std::regex_match(x, date)){ //check if the item is a date
                if (dates.size() == 0){
                    earliest_date = Date(x);
                }
                dates.push_back(Date(x));
                price_data >> x;
            }
            matrix_data.push_back(std::stod(x)); // If not a date, then push into the matrix
            price_data >> x;
            if (price_data.peek() == '\n'){ // Keep track of the number of rows we read through
                row_counter++;
            }
        }
        assign_dimensions(row_counter, column_labels.size());
        latest_date = this->dates[dates.size() - 1];
    }

    /* Estimate the historical volatility of the stock by computing the standard deviation of the % change in stock price (daily) and multiplying by
    15.937, which is the square root of the number of trading days in a year
    UPDATE: This function will take in two dates and figure out the historical volatility based on that time frame
            The idea is to be able to get different historical volatilities based on how far we look (Relative to present time) */
    double estimate_Hvolatility(const Date begin, const Date end, const std::string timeframe = "null"){
        Date begin_time = begin;
        Date end_time = end;
        double time_length = get_time_diff(begin_time, end_time);
        /*if (( get_time_diff(earliest_date, begin_time) < 0) || !(search_vector(dates, begin_time))){ // If the supplied date is not a trading day or is earlier than the earliest date
        std::cout << (get_time_diff(earliest_date, begin_time) < 0) << "\n";
        std::cout << !(search_vector(dates, begin_time)) << "\n";
        std::cout << "Error: Invalid date supplied as Argument";
        return 0;
        } */
        size_t date_index_b = 1; // initialize time index for start date
        for (const auto& i : dates){
            if (i == begin_time){
                break;
            }
            ++date_index_b;
        }
        size_t date_index_e = date_index_b; // initialize time index for end date
        for (size_t i = date_index_b; i < dates.size(); i++){
            if (dates[i] == end_time){
                break;
            }
            ++date_index_e;
        }
        std::cout << "date index (begin): " << date_index_b << " \n";
        std::cout << "date index (end): " << date_index_e << " \n";
        std::vector<double> daily_changes;
        size_t trading_days = 1;
        double average_change = 0;
        double change_variance = 0;
        double change_SD = 0;
        // Find the column index that contains the closing prices "Close"
        size_t close_index = 0;
        /* for (const auto& i : column_labels){
        if (i != "Close"){
        close_index++;
        }
        } */
        std::for_each(column_labels.begin(), column_labels.end(), [&](std::string a){(a != "Close") ? ++close_index : close_index = close_index; });
        // iteration should start at a value >= 2
        for (size_t i = date_index_b + 1; i < date_index_e; i++){
            // For each closing prices, starting at the inputted date, compute the percentage change and store it in the vector
            double single_day_change = 100 * ((*this)(i, close_index) - (*this)(i - 1, close_index)) / (*this)(i - 1, close_index);
            daily_changes.push_back(single_day_change);
            ++trading_days;
        }

        // Find the standard deviation of the daily_changes vector
        average_change = find_arithmatic_average(daily_changes);
        for (size_t i = 0; i < daily_changes.size(); i++){
            change_variance += (daily_changes[i] - average_change)*(daily_changes[i] - average_change);
        }

        change_variance /= daily_changes.size();
        // TEST CODE
        std::cout << "time difference: " << date_index_e - date_index_b << "\n";
        return sqrt(change_variance)*sqrt((date_index_e)-(date_index_b)) / 100.0; // end-begin index describes the number of trading days
    }
    /* Wrap the estimate_Hvolatility in another function such that the user can input a time frame (inputs are "Week, Month, Quarter, Year, All Time")
    */
    double Stock_Price_Matrix::estimate_Hvolatility_TF(std::string timeframe){
        if (timeframe == "Week"){ // Look back 5 days
            Date start = dates[dates.size() - 5];
            Date end = dates[dates.size() - 1];
            return estimate_Hvolatility(start, end)*sqrt(52);
        }
        else if (timeframe == "Month"){ // Look back 20 days
            Date start = dates[dates.size() - 20];
            Date end = dates[dates.size() - 1];
            return estimate_Hvolatility(start, end)*sqrt(12);
        }
        else if (timeframe == "Quarter"){ // Look back 63 days
            Date start = dates[dates.size() - 63];
            Date end = dates[dates.size() - 1];
            return estimate_Hvolatility(start, end)*sqrt(4);
        }
        else if (timeframe == "Year"){ // Start from the beginning
            Date start = dates[0];
            Date end = dates[dates.size() - 1];
            return estimate_Hvolatility(start, end);
        }
    }
    void dim();

private:
    Date earliest_date;
    Date latest_date;
};


#endif
