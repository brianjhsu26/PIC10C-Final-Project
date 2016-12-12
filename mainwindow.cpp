#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDesktopServices>
#include <QTextStream>
#include <QFileInfo>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    // Select Call option text files
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("Text files (*.txt)"));
    QFile file(fileName);
    QFileInfo fileinfo = QFileInfo(fileName);
    QDir::setCurrent(fileinfo.absolutePath()); // Changes current directory to one containing put option .txt file
    this->calls = std::fstream(fileinfo.fileName().toStdString());
}

void MainWindow::on_pushButton_2_clicked()
{
    // Select Put option text files
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("Text files (*.txt)"));
    QFile file(fileName);
    QFileInfo fileinfo = QFileInfo(fileName);
    QDir::setCurrent(fileinfo.absolutePath()); // Changes current directory to one containing put option .txt file
    this->puts = std::fstream(fileinfo.fileName().toStdString());
}


void MainWindow::on_pushButton_3_clicked()
{
    // Select stock history text file
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr("Text files (*.txt)"));
    QFile file(fileName);
    QFileInfo fileinfo = QFileInfo(fileName);
    QDir::setCurrent(fileinfo.absolutePath()); // Changes current directory to one containing put option .txt file
    this->stocks = std::fstream(fileinfo.fileName().toStdString());

    // After putting this info in, we can assemble the matrixes and create the dropdown menu for the strike price
}

void MainWindow::on_pushButton_4_clicked()
{
    // Assemble the option matrix and stock price matrix
    OM = Option_Matrix(calls, puts);
    SPM = Stock_Price_Matrix(stocks);
    std::vector<double> v = OM.get_strike_price_vector();

    // Label the option matrix and stock
    OM.set_characteristics(this->ui->lineEdit->text().toStdString(),
                           this->ui->lineEdit_2->text().toStdString(),
                           this->ui->lineEdit_3->text().toDouble());

    // Fill the dropdown menu with a list of strike prices as determined from the option file
    for(size_t i = 0; i < v.size(); i++){
        double element = v[i];
        QString valueAsString = QString::number(element);
        this->ui->comboBox->addItem(valueAsString);
    }

    this->ui->comboBox_2->addItem(QString("Week"));
    this->ui->comboBox_2->addItem(QString("Month"));
    this->ui->comboBox_2->addItem(QString("Quarter"));
    this->ui->comboBox_2->addItem(QString("Year"));
}


void MainWindow::on_pushButton_5_clicked()
{
    // Clear the screen
    this->ui->textEdit->setText("");

    // Find all Necessary constants
    double rf_rate = OM.find_r(find_max);
    std::string lookback_period = this->ui->comboBox_2->currentText().toStdString();
    double hist_vol = SPM.estimate_Hvolatility_TF(lookback_period);
    double imp_vol = OM.estimate_Ivolatility(this->ui->comboBox->currentText().toDouble(), rf_rate);
    double strikeprice_ = this->ui->comboBox->currentText().toDouble();
    double market_price = OM.get_actual_call_price(strikeprice_);
    double bsprice_impvol = OM.black_scholes_price(imp_vol, rf_rate, strikeprice_, "Call");
    double bsprice_histvol = OM.black_scholes_price(hist_vol, rf_rate, strikeprice_, "Call");
    double btreeprice = OM.binomial_tree_price(imp_vol, rf_rate, strikeprice_, "Call")[OM.binomial_tree_price(imp_vol, rf_rate, strikeprice_, "Call").size()-1];

    // Make them into std::string objects
    std::string rfrate_info = "Risk Free Rate: " + std::to_string(rf_rate) + "\n";
    std::string histvol_info = "Historical volatility: " + std::to_string(hist_vol) + "\n";
    std::string impvol_info = "Implied volatility: " + std::to_string(imp_vol) + "\n";
    std::string marketprice_info = "Actual market price: " + std::to_string(market_price) + "\n";
    std::string bs_info_imp = "Option price as generated by the Black Scholes Equation using the implied volatility: " + std::to_string(bsprice_impvol) + "\n";
    std::string bs_info_hist = "Option price as generated by the Black Scholes Equation using the historical volatility: " + std::to_string(bsprice_histvol) + "\n";
    std::string btreeprice_info = "Option price as generated by implied volatility Binomial Tree Convergence: " + std::to_string(btreeprice) + "\n";


    // Convert and display the created strings as qstrings in the text box
    this->ui->textEdit->append(QString::fromStdString(rfrate_info));
    this->ui->textEdit->append(QString::fromStdString(histvol_info));
    this->ui->textEdit->append(QString::fromStdString(impvol_info));
    this->ui->textEdit->append(QString::fromStdString(marketprice_info));
    this->ui->textEdit->append(QString::fromStdString(bs_info_imp));
    this->ui->textEdit->append(QString::fromStdString(bs_info_hist));
    this->ui->textEdit->append(QString::fromStdString(btreeprice_info));


}