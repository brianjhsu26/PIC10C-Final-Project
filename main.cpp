#include "mainwindow.h"
#include "matrix_classes.h"
#include <QApplication>
#include <QDir>

int main(int argc, char *argv[])
{
    // Googl: $760.54 at 11.20.2016
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
    return 0;
}
