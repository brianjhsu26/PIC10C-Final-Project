#include "matrix_classes.h"

//Custom_Matrix::Custom_Matrix(){ rows = 0; columns = 0; }

void Custom_Matrix::write_to_header(const std::vector<std::string>& header){
    {
        for (const auto& i : header){
            column_labels.push_back(i);
        }
    }
}

inline void Custom_Matrix::write_to_data(const double& data){
    matrix_data.push_back(data);
}



void Custom_Matrix::dim(){
    std::cout << "Column Headers: ";
    for (const auto& i : column_labels){
        std::cout << i << " ";
    }
    std::cout << "\n";
    std::cout << "Rows: " << rows << "\n";
    std::cout << "Columns: " << columns << "\n";
}



// OPTION MATRIX IMPLEMENTATIONS



std::vector<double> Option_Matrix::list_strike_prices(){
    std::vector<double> strike_prices;
    std::cout << "Call option strike prices: ";
    for (size_t i = 1; i < calls.rows - 1; i++){
        strike_prices.push_back(calls(i, 2));
    }
    strike_prices.push_back(calls(calls.rows, 2));
    return strike_prices;
}






// STOCK MATRIX FUNCTION IMPLEMENTATIONS




void Stock_Price_Matrix::dim(){
    std::cout << "Stock Historical Price Matrix" << "\n";
    std::cout << "Earliest Date: " << earliest_date << "\n";
    this->Custom_Matrix::dim();
}

