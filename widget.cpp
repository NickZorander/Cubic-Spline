#include "widget.h"
#include <QApplication>
#include<QPainter>
#include<QPaintEvent>
#include<QtWidgets>
#include<cstdlib>
#include <QtDebug>



void MyWidget::paintEvent(QPaintEvent* e){

    QPainter painter(this);
    int i;

    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(Qt::black));

    QVector<QPointF>::iterator it=SplinePoints.begin();
    painter.drawPolyline(it, SplinePoints.size());

    painter.setPen(QPen(Qt::yellow, 6));
    painter.drawPoint(cursor);

    for (int i=0; i<=N-1; ++i)
    {
        painter.drawPoint(StartPoints[i]);
    }
}





void MyWidget::mouseReleaseEvent(QMouseEvent* pe)
{
    QPointF Q(pe->x(),pe->y());

    if (change_point_flag)
        {
            StartPoints[change_point_index]=Q;
            BuildSpline();
            repaint();
        }
    change_point_flag=false;
}



qreal distance (const QPointF& p1, const QPointF& p2)
{
    return qSqrt(  (p2.x()-p1.x())*(p2.x()-p1.x()) + (p2.y()-p1.y())*(p2.y()-p1.y()) );
}





void MyWidget::mousePressEvent(QMouseEvent *pe)
{
    change_point_flag=false;

    QPointF Q(pe->x(),pe->y());


    for (int i=0; i<=N-1; ++i)
        {
            if (distance(Q, StartPoints[i]) <=7)
                {
                    change_point_index=i;
                    change_point_flag=true;
                    break;
                }
        }
}


void MyWidget::mouseMoveEvent(QMouseEvent *pe)
{
    cursor=QPointF(pe->x(), pe->y());
    repaint();
}





void MyWidget::InitData()
{
    QFile inputfile("D:\\QT\\Bezier_Curve\\points.txt");

    if (!inputfile.open(QFile::ReadOnly | QFile::Text)){
        QMessageBox::information(this,"Error", "path not found!");
        return;
    } //проверили на открытие

    QTextStream stream(&inputfile);

    QString A =stream.readAll(); //сделали содержимое файла строчкой

    QStringList list = A.split(' ', QString::SkipEmptyParts); //распилили строку в список где каждый элемент - одна координата точки
    bool ok;

    N=list[0].toInt(&ok); //первым числом передается число точек

    StartPoints.resize(N);
    SplinePoints.resize((N-1)*100);


    for (int i=0; i<=N-1; ++i)
        {
            QPointF p(list[i*2+1].toInt(&ok), list[i*2+2].toInt(&ok));
            StartPoints[i]=p;
        } //создали массив точек


    change_point_flag=false;


}


//______________________________________________________________________________________________________________________________________________________________________________


bool operator<(const QPointF& p1, const QPointF& p2 )
{
      return (p1.x()<p2.x());
}


qreal MyWidget::Hk(int k)
{
    return StartPoints[k].x()-StartPoints[k-1].x();
}

qreal MyWidget::Alfa(int k)
{
    return 1/Hk(k);
}

qreal MyWidget::Fxk (int k) //k~x+1
{
    qDebug() << "Fxk [" << k-1 << "," << k <<"]";
    return (StartPoints[k-1].y()-StartPoints[k-2].y())/(StartPoints[k-1].x()-StartPoints[k-2].x());
}

qreal MyWidget::Beta (int k)
{
    qDebug() << "Beta k:" << k;
    return 2*( (1/Hk(k)) + (1/Hk(k+1)) );
}

qreal MyWidget::Gamma (int k)
{
    return 1/Hk(k+1);
}

qreal MyWidget::Delta(int k)
{
    qDebug() << "Delta k:" << k;
    return 3*( Fxk(k+2)/Hk(k+1) + Fxk(k+1)/Hk(k) );
}



/**
     * n - число уравнений (строк матрицы)
     * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
     * c - главная диагональ матрицы A (нумеруется: [0;n-1])
     * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
     * f - правая часть (столбец)
     * x - решение, массив x будет содержать ответ
     */
void solveMatrix (int n, qreal *a, qreal *c, qreal *b, qreal *f, qreal *x)
{
    qreal m;
    for (int i = 1; i < n; i++)
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*b[i-1];
        f[i] = f[i] - m*f[i-1];
    }

    x[n-1] = f[n-1]/c[n-1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}


void MyWidget::Create_Bk()
{
    QVector <qreal> alfa;
    QVector <qreal> beta;
    QVector <qreal> gamma;
    QVector <qreal> delta;


    alfa.resize(N-2);
    alfa[0]=0;
    for (int i=1; i<=N-3;++i)
        {
            alfa[i]=Alfa(i+1);
        }
    qDebug()<<"check1";

    beta.resize(N-2);
    for (int i=0; i<=N-3;++i)
        {
            qDebug()<<"here?";
            beta[i]=Beta(i+1);
            qDebug()<<"hmmm";
        }
    qDebug()<<"check2";

    delta.resize(N-2);
    for (int i=0; i<=N-3;++i)
        {
            qDebug()<<"i:="<<i;
            delta[i]=Delta(i+1);
        }
    qDebug()<<"check3";

    gamma.resize(N-2);
    gamma[N-3]=0;
    for (int i=0; i<=N-4;++i)
        {
            gamma[i]=Gamma(i+1);
        }
    qDebug()<<"check4";

     QVector <qreal>::iterator alfa_it=alfa.begin();
     QVector <qreal>::iterator beta_it=beta.begin();
     QVector <qreal>::iterator gamma_it=gamma.begin();
     QVector <qreal>::iterator delta_it=delta.begin();

     Bk.resize(N-2);
     QVector <qreal>::iterator Bk_it=Bk.begin();

     qDebug()<<"check5";
     solveMatrix(N-2, alfa_it, beta_it, gamma_it, delta_it, Bk_it);
     qDebug()<<"check6";
}



qreal MyWidget::A(int k)
{
    return StartPoints[k-1].y();
}

qreal MyWidget::B(int k)
{
    if (k==1 || k==N )
        return 0;
    else
        {
            return Bk[k-2];
        }
}

qreal MyWidget::C(int k)
{
    return ( 3*Fxk(k+1) - B(k+1) -2*B(k) )/Hk(k);
}

qreal MyWidget::D(int k)
{
    return ( B(k) + B(k+1) -2*Fxk(k+1) )/ ( Hk(k) * Hk(k) );
}





qreal MyWidget::P (qreal x, int k)
{
    return (A(k) + B(k)*(x-StartPoints[k-1].x()) + C(k)*(x-StartPoints[k-1].x())*(x-StartPoints[k-1].x()) + D(k)*(x-StartPoints[k-1].x())*(x-StartPoints[k-1].x())*(x-StartPoints[k-1].x()));
}




void MyWidget::BuildSpline()
{
    SplinePointsIterator=0;
    qSort(StartPoints);
    qDebug()<<"N:"<< N;
    Create_Bk();
    qDebug()<<"!";
    for (int i=1; i<=N-1; i++)
        {
            qreal x_start=StartPoints[i-1].x();
            qreal x_end=StartPoints[i].x();

            for (qreal t=0; t<=1; t+=0.01)
                {
                    qDebug()<<"?";
                    SplinePoints[SplinePointsIterator]=QPointF(x_start+t*(x_end-x_start), P(x_start+t*(x_end-x_start), i));
                     ++SplinePointsIterator;
                }
        }
}


































