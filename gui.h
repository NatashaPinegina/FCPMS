//
// Created by natasha on 29.12.22.
//

#ifndef UNTITLED5_GUI_H
#define UNTITLED5_GUI_H

#include <QWidget>
#include "ui_gui.h"
#include "qcustomplot/qcustomplot.h"
#include <vector>
//#include <QMainWindow>

namespace ui_gui {
    QT_BEGIN_NAMESPACE
    namespace Ui { class gui; }
    QT_END_NAMESPACE

    struct ParamSignal
    {
        double f0;
        double fdisk;
        double snr;
        double t;
        double duration;
        double bitrate;
        double outputT;
        double p;
        double koef_dur;
        int typeShum;
    };

    struct InfoList
    {
        std::vector<std::vector<double>> MassOtrisovka;
        std::vector<std::vector<double>> MassOtshetX;
    };

    class gui : public QWidget {
        Q_OBJECT

    public:
        explicit gui(QWidget *parent = nullptr);
        QVector<QVector<double>> XY;
        QVector<double> x,y;
        int ind;
        int size;
        double bufmaxY;

        ~gui() override;

    private:
        Ui::gui *ui;

    private slots:
        void getParam(ParamSignal& param);
        void on_getShortSigma_clicked();
        void on_getLongSigma_clicked();
        void on_calculate_clicked();
        void on_listWidget_graf_itemClicked(QListWidgetItem *item);
        void on_whiteGauss_clicked(bool checked);
        void on_FlikerShum_clicked(bool checked);
        void on_raceNoise_clicked(bool checked);
        void find_max_min(QVector<double>& x, QVector<double>y,double& maxX, double& minX, double& maxY, double& minY);
    };
} // ui_gui


/*namespace Ui{
    class MainWindow;
}
class Ui::MainWindow: public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
private:
    Ui::MainWindow * ui;
};
*/
#endif //UNTITLED5_GUI_H
