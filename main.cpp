#include "widget.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MyWidget w;
    w.InitData();
    w.BuildSpline();
    w.show();

    return a.exec();
}
