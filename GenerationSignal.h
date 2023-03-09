//
// Created by natasha on 31.01.23.
//

#ifndef UNTITLED5_GENERATIONSIGNAL_H
#define UNTITLED5_GENERATIONSIGNAL_H

#include <vector>
#include <complex>
#include <fstream>
#include "gui.h"

using namespace std;
using namespace ui_gui;

struct Signal{
    vector<double> signal;
};

class GenerationSignal
{
public:

    /**Variables*/
    double zader=0.3;
    int kol_bit;
    double sizeSPM;
    int KolSat = 3;
    int KolSour = 2;


    /**Function*/
    double GenerateSin();
    bool RandomBit(double low_chance);
    void Rxx(vector<vector<double>>& R, vector<double>& signal1, vector<double>& signal2);
    void Newcorrelate(vector<double>& signal1, vector<double>& signal2, vector<double>& correlation);
    void correlate(vector<double>& base_signal, vector<double>& correlation);
    double DoubleRand(double _max, double _min);
    void addNoise(vector<double>& buf, double SNR, ParamSignal& param);
    void GetSigma(vector<double>& InputSignal, vector<double>rxx, vector<double>& Sigma, ParamSignal& param);
    void GetLongSignal(vector<double>& InputSignal);
    void correlate(vector<double>& base_signal, vector<double>& analyzed_signal, vector<double>& correlation, double& sredSigma, double&sredLongSigma);
    QVector<QVector<double>> GenerateShortSignal(ParamSignal& param);
    void newFFT(vector<complex<double>>& in, int direction);
    void GetK(vector<complex<double>>& K, vector<double>& SPNSigma, vector<double>& SPNShum, vector<complex<double>>& HSopr, vector<double>& HMod);
    void GetSPM(vector<complex<double>>& Pw, vector<complex<double>>& PShum, vector<double>& PSh, vector<double>& PSg, ParamSignal& param);
    void newcorrelate(vector<complex<double>>& base_signal, vector<complex<double>>& analyzed_signal, vector<double>& PSg, vector<double>& PSh, vector<complex<double>>& correlation);
    QVector<QVector<double>> GenerateLongSignal(ParamSignal& param, int& size);
    InfoList Calculate(ParamSignal& param);
    void coherentSumm(vector<vector<double>> &sgs, vector<double> &result);
    void transformSignal(vector<double>& base_signal, double delay, double duration, double fshift, double scale, double SNR, vector<double>& ret_sig, ParamSignal& param);

    ~GenerationSignal();

    /**Vector*/
    vector<double>Sigma;
    vector<vector<double>>LongSigma;
    vector<vector<double>> InvR;
    vector<vector<double>> InputSignal;
    vector<double> LongInputSignal;
    vector<double> ClearSignal;
    vector<double> SHUM;
    vector<vector<double>> bitPosl;
    QVector<double> x,y;
    InfoList info;
    QVector<QString> listt;
    vector<vector<double>> SummSignal;
    vector<double> sign;

    vector<double> sdvigTime;
    vector<double> sdvigFreq;
};


#endif //UNTITLED5_GENERATIONSIGNAL_H
