#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

class MyWidget : public QWidget
{
    Q_OBJECT
public:
    QVector <QPointF> SplinePoints;
    QVector <QPointF> StartPoints;
    QVector <qreal> Bk;

    int SplinePointsIterator=0;
    int N;
    bool change_point_flag=false;
    int change_point_index=0;



    QPointF cursor;

    qreal Hk (int k);
    qreal Fxk (int k);
    qreal Alfa(int k);
    qreal Beta(int k);
    qreal Gamma(int k);
    qreal Delta(int k);

    qreal A(int k);
    void Create_Bk();
    qreal B(int k);
    qreal C(int k);
    qreal D(int k);

    qreal P (qreal x, int k);
    void BuildSpline(); //заполняет массив SplinePoints достаточным для отрисовки количеством точек
    void InitData(); //инициализация данных
    void PrintInfo();

    MyWidget(QWidget* p=0) {}
    virtual ~MyWidget() {}
protected:
     void paintEvent(QPaintEvent*e);
     void mouseReleaseEvent(QMouseEvent* pe);
     void mousePressEvent(QMouseEvent *pe);
     void mouseMoveEvent(QMouseEvent *pe);
};


#endif // WIDGET_H
