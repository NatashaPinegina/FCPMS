//
// Created by natasha on 29.12.22.
//

// You may need to build the project (run Qt uic code generator) to get "ui_gui.h" resolved

#include "gui.h"
#include "GenerationSignal.h"


namespace ui_gui {

    GenerationSignal n;
    ParamSignal param;

    gui::gui(QWidget *parent) : QWidget(parent), ui(new Ui::gui) {
        ui->setupUi(this);
        ui->f0->setText(QString::number(25000));
        ui->fdisk->setText(QString::number(250000));
        ui->SNR->setText(QString::number(20));
        ui->bitrate->setText(QString::number(4800));
        ui->duration->setText(QString::number(0.001));
        ui->outputT->setText(QString::number(30));
        ui->p->setText(QString::number(7));
        ui->koef_dur->setText(QString::number(0.5));
        ui->t->setText(QString::number(0));
    }

    gui::~gui() {
        delete ui;
    }

    int kol = 0;

    void gui::find_max_min(QVector<double>& x, QVector<double>y,double& maxX, double& minX, double& maxY, double& minY)
    {
        for(int i=0;i<y.size();i++)
        {
            if(maxY<y[i]) maxY=y[i];
        }
        minY=-maxY/50;
        minX=x[0];
        maxX=x[x.size()-1];
    }

    void gui::getParam(ParamSignal& param)
    {
        param.f0 = ui->f0->text().toDouble();
        param.fdisk = ui->fdisk->text().toDouble();
        param.snr = ui->SNR->text().toDouble();
        param.bitrate = ui->bitrate->text().toDouble();
        param.duration = ui->duration->text().toDouble();
        param.outputT = ui->outputT->text().toDouble();
        param.p = ui->p->text().toDouble();
        param.koef_dur = ui->koef_dur->text().toDouble();
        param.t = ui->t->text().toDouble();
    }

    void gui::on_getShortSigma_clicked(){

        n.info.MassOtrisovka.clear();
        n.info.MassOtshetX.clear();
        n.listt.clear();
        ui->listWidget_graf->clear();
        //n.~GenerationSignal();

        getParam(param);
        XY = n.GenerateShortSignal(param);
        size = XY[0].size();
        double maxX=0, minX=0, maxY=0, minY=0;

        find_max_min(XY[0], XY[1],maxX, minX,  maxY,  minY);

        // создаем график и добавляем данные:
        ui->shortSigma->addGraph();
        ui->shortSigma->graph(0)->setData(XY[0], XY[1]);
        // задаем имена осей координат
        ui->shortSigma->xAxis->setLabel("duration");
        ui->shortSigma->yAxis->setLabel("ampl");
        // задаем размеры осей
        ui->shortSigma->xAxis->setRange(minX, maxX);
        ui->shortSigma->yAxis->setRange(minY, maxY);
        ui->shortSigma->replot();
    }

    void gui::on_getLongSigma_clicked(){
        XY.clear();
        XY = n.GenerateLongSignal(param, size);
        double maxX=0, minX=0, maxY=0, minY=0;

        find_max_min(XY[0], XY[1],maxX, minX,  maxY,  minY);

        // создаем график и добавляем данные:
        ui->longSigma->addGraph();
        ui->longSigma->graph(0)->setData(XY[0], XY[1]);
        // задаем имена осей координат
        ui->longSigma->xAxis->setLabel("duration");
        ui->longSigma->yAxis->setLabel("ampl");
        // задаем размеры осей
        ui->longSigma->xAxis->setRange(minX, maxX);
        ui->longSigma->yAxis->setRange(minY, maxY);
        ui->longSigma->replot();

    }

    void gui::on_calculate_clicked() {
        InfoList info = n.Calculate(param);
        for(int i=0;i<n.listt.size();i++)
        {
            ui->listWidget_graf->addItem(n.listt[i]);
        }
    }
    void gui::on_listWidget_graf_itemClicked(QListWidgetItem *item)
    {
        ind=ui->listWidget_graf->row(item);
        y.resize(n.info.MassOtrisovka[ind].size());
        x.resize(n.info.MassOtshetX[ind].size());
        double min=0;
        for(int i=0;i<n.info.MassOtrisovka[ind].size();i++)
        {
            y[i] = n.info.MassOtrisovka[ind][i];
            x[i] = n.info.MassOtshetX[ind][i];
            if (min>y[i]) min=y[i];
        }
        double maxX=0, minX=0, maxY=0, minY=0;

        find_max_min(x, y,maxX, minX,  maxY,  minY);

        // создаем график и добавляем данные:
        ui->vse->addGraph();
        ui->vse->graph(0)->setData(x, y);
        // задаем имена осей координат
        ui->vse->xAxis->setLabel("duration");


        ui->vse->yAxis->setLabel("ampl");
        // задаем размеры осей
        ui->vse->xAxis->setRange(minX, maxX);
        ui->vse->yAxis->setRange(min, maxY);
        ui->vse->replot();
    }

    void gui::on_whiteGauss_clicked(bool checked)
    {
        if(checked) param.typeShum=1;
    }
    void gui::on_FlikerShum_clicked(bool checked)
    {
        if(checked) param.typeShum=2;
    }
    void gui::on_raceNoise_clicked(bool checked)
    {
        if(checked) param.typeShum=3;
    }

} // ui_gui
