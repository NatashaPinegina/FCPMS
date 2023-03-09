//
// Created by natasha on 31.01.23.
//
#include "GenerationSignal.h"
#include "gui.h"

using namespace ui_gui;

bool Ymin = false;
double dur=0;

bool GenerationSignal::RandomBit(double low_chance)
{
    return (rand() < RAND_MAX * low_chance);
}

double GenerationSignal::DoubleRand(double _max, double _min)
{
    return _min + double(rand()) / RAND_MAX * (_max - _min);
}

void GenerationSignal::addNoise(vector<double>& buf, double SNR, ParamSignal& param)
{
    //srand(0);
    double energy = std::accumulate(buf.begin(), buf.end(), 0.0,
                                    [&](double a, double b)
                                    {
                                        return a + b * b;
                                    });

    vector <double> noise(buf.size(), 0.0);
    vector <double> ampl(buf.size(), 0.0);
    std::transform(ampl.begin(), ampl.end(), ampl.begin(),
                   [&](double a)
                   {
                       int count_rep = 20;
                       a = 0.0;
                       for (int rep = 0; rep < count_rep; rep++)
                       {
                           a += ((double)rand() / RAND_MAX) * 2.0 - 1.0;
                       }
                       a /= count_rep;
                       return a;
                   });
    std::transform(noise.begin(), noise.end(), ampl.begin(), noise.begin(),
                   [&](double a, double ampl)
                   {
                       double fi = (double)rand() / RAND_MAX * 2.0 * 3.1415926535897;
                       return ampl * cos(fi);
                   });

    if (param.typeShum==1) ;
    else if (param.typeShum==2)
    {
        /*double time =0;
        double timeshag=1/param.fdisk;
        vector<double> H(noise.size(),0);
        for (int i = 0; i < noise.size(); i++)
        {
            double Re = 0;
            if (i > 0) Re = 1/sqrt(time);
            H[i] = Re;
            time += timeshag;
        }

        vector<double> convolution;
        for (int n = 0; n < (int)(noise.size()); n++)
        {
            double counter = 0;
            for (int m = 0; m < H.size(); m++)
            {
                if ((-m + n) >= 0)
                    counter += noise[-m + n] * H[m];
                else
                    continue;
            }
            convolution.push_back(counter);
        }
        noise=convolution;
        convolution.clear();
        H.clear();*/

        vector<complex<double>> n;
        for(int i=0;i< noise.size();i++)
        {
            complex<double> s (noise[i],0);
            n.push_back(s);
        }
        newFFT(n, -1);

        double freq =0;
        double freqshag = param.fdisk / n.size();
        vector<complex<double>> H(n.size(),0);
        for (int i = 0; i < noise.size(); i++)
        {
            double Re = 0;
            if (i > 0) Re = 1/sqrt(freq);
            complex<double> s (Re, 0);
            H[i] = s;
            freq += freqshag;
        }

        vector<double>convolution;

        for(int i = 0; i < n.size(); i++)
        {
            vector<complex<double>> Peremn;
            for (int j = 0; j < H.size(); j++)
            {
                //int sdvig=i-j;
                //if (sdvig >= 0)
                    Peremn.push_back(H[j] * conj(n[j]));
               // else
                //    Peremn.push_back(0);
            }
            newFFT(Peremn, 1);
            convolution.push_back(Peremn[0].real());
        }

        noise=convolution;
        convolution.clear();
        H.clear();

    }
    else if (param.typeShum==3)
    {
        double K=12;
        double koeff= pow(10,-K/20);

        vector<double> buf_;
        buf_.resize(buf.size());

        for(int i=0;i<buf.size();i++)
            buf_[i]=koeff*buf[i];


        double t1=DoubleRand( dur/3,0);
        double t2=DoubleRand( 2*dur/3, dur/3);
        double t3=DoubleRand(dur, 2*dur/3);

        vector<double> zader_buf1;
        zader_buf1.resize(buf_.size());
        vector<double> zader_buf2;
        zader_buf2.resize(buf_.size());
        vector<double> zader_buf3;
        zader_buf3.resize(buf_.size());

        double time=0;
        double time_shag=1/param.fdisk;
        for(int i=0;i<buf_.size();i++)
        {
            if(time>t1) zader_buf1[i]=buf_[i];
            else zader_buf1[i]=0;
            if(time>t2) zader_buf2[i]=buf_[i];
            else zader_buf2[i]=0;
            if(time>t3) zader_buf3[i]=buf_[i];
            else zader_buf3[i]=0;

            double koeff_zamir = zader_buf1[i]+zader_buf2[i]+zader_buf3[i];
            buf_[i]*=koeff_zamir;
            time+=time_shag;
        }

        for(int i=0;i<buf.size();i++)
        {
            buf[i]+=buf_[i];
        }

        energy = std::accumulate(buf.begin(), buf.end(), 0.0,
                                 [&](double a, double b)
                                 {
                                     return a + b * b;
                                 });
    }

    double noise_energy = std::accumulate(noise.begin(), noise.end(), 0.0,
                                          [&](double a, double b)
                                          {
                                              return a + b* b;
                                          });
    double norm_coef = energy * pow(10.0, -SNR / 10.0) / noise_energy;


    std::transform(buf.begin(), buf.end(), noise.begin(), buf.begin(),
                   [&](double a, double b)
                   {
                       return (a + sqrt(norm_coef) * b);
                   });
    if (SNR == 10)
        SHUM = noise;
}

double** matrix(int n, int m)
{
    double** matr = new double* [n];
    for (int i = 0; i < n; i++)
        matr[i] = new double[m];
    return matr;
}

int svd_hestenes(int m_m, int n_n, double* a, double* u, double* v, double* sigma)
{
    double thr = 1.E-4f, nul = 1.E-16f;
    int n, m, i, j, l, k, lort, iter, in, ll, kk;
    double alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
    n = n_n;
    m = m_m;
    for (i = 0; i < n; i++)
    {
        in = i * n;
        for (j = 0; j < n; j++)
            if (i == j) v[in + j] = 1.;
            else v[in + j] = 0.;
    }
    for (i = 0; i < m; i++)
    {
        in = i * n;
        for (j = 0; j < n; j++)
        {
            u[in + j] = a[in + j];
        }
    }

    iter = 0;
    while (1)
    {
        lort = 0;
        iter++;
        for (l = 0; l < n - 1; l++)
            for (k = l + 1; k < n; k++)
            {
                alfa = 0.; betta = 0.; hamma = 0.;
                for (i = 0; i < m; i++)
                {
                    in = i * n;
                    ll = in + l;
                    kk = in + k;
                    alfa += u[ll] * u[ll];
                    betta += u[kk] * u[kk];
                    hamma += u[ll] * u[kk];
                }

                if (sqrt(alfa * betta) < nul)	continue;
                if (fabs(hamma) / sqrt(alfa * betta) < thr) continue;

                lort = 1;
                eta = (betta - alfa) / (2.f * hamma);
                t = (double)((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta * eta)));
                cos0 = (double)(1. / sqrt(1. + t * t));
                sin0 = t * cos0;

                for (i = 0; i < m; i++)
                {
                    in = i * n;
                    buf = u[in + l] * cos0 - u[in + k] * sin0;
                    u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
                    u[in + l] = buf;

                    if (i >= n) continue;
                    buf = v[in + l] * cos0 - v[in + k] * sin0;
                    v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
                    v[in + l] = buf;
                }
            }

        if (!lort) break;
    }

    for (i = 0; i < n; i++)
    {
        s = 0.;
        for (j = 0; j < m; j++)	s += u[j * n + i] * u[j * n + i];
        s = (double)sqrt(s);
        sigma[i] = s;
        if (s < nul)	continue;
        for (j = 0; j < m; j++)	u[j * n + i] /= s;
    }
    //======= Sortirovka ==============
    for (i = 0; i < n - 1; i++)
        for (j = i; j < n; j++)
            if (sigma[i] < sigma[j])
            {
                s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
                for (k = 0; k < m; k++)
                {
                    s = u[i + k * n]; u[i + k * n] = u[j + k * n]; u[j + k * n] = s;
                }
                for (k = 0; k < n; k++)
                {
                    s = v[i + k * n]; v[i + k * n] = v[j + k * n]; v[j + k * n] = s;
                }
            }

    return iter;
}

double** Composition_Matrix_Two(double** Inital_Matrix, double** Inverse_Matrix, int n)
{
    double** Composition = matrix(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Composition[i][j] = 0;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                Composition[i][j] += Inital_Matrix[i][k] * Inverse_Matrix[k][j];
    return Composition;

}

vector<double> Composition_Matrix_Stroka(vector<vector<double>> Inital_Matrix, vector<double> Inverse_Matrix, int n, int m)
{
    vector<double> Composition;
    Composition.resize(m);
    for (int i = 0; i < m; i++)
    {
        Composition[i] = 0.;
    }

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            Composition[i] += Inital_Matrix[i][j] * Inverse_Matrix[j];
    return Composition;
}

void GenerationSignal::GetSigma(vector<double>& InputSignal, vector<double>rx, vector<double>& Sigma, ParamSignal& param)
{
    Ymin = false;
    //частота дискретизации
    double samplingFrequency = param.fdisk;
    //строим автокорреляционную матрицу
    int P = param.p;
    double* rxx = new double[P];
    memset(rxx, 0, (P) * sizeof(double));
    double summ;
    int iter = 0;
    vector<vector<double>> cor;
    cor.clear();
    cor.resize(1);
    for (int m = 0; m < P; m++)
    {
        /*summ = 0;
        for (int k = 0; k < sign.size(); k++)
        {
            if(k+m>=0 && (m + k) < sign.size())
            summ += sign[k] * sign[k + m];
        }
        rxx[iter] = summ / (sign.size());*/
        rxx[iter] = rx[iter];
        cor[0].push_back(rx[iter]);
        iter++;
    }

    double** Rij = matrix(P, P);
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            Rij[i][j] = 0;
        }
    }
    /*int t = p-1;
    iter = 1;
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            Rij[i][j] = rxx[t-j];
        }
        t = t + iter;
    }*/

    int t = 0;
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            if (i - j < 0)
            {
                t = abs(i - j);
                Rij[i][j] = rxx[t];
            }
            else
                Rij[i][j] = rxx[i - j];
        }
    }


    double* MassInStroka = new double[P * P];
    memset(MassInStroka, 0, (P * P) * sizeof(double));
    int k = 0;
    for (int j = 0; j < P; j++)
    {
        for (int i = 0; i < P; i++)
        {
            MassInStroka[k] = Rij[i][j];
            k++;
        }
    }

    double* U = new double[P * P];
    double* V = new double[P * P];
    double* G = new double[P];

    memset(U, 0, P * P * sizeof(double));
    memset(V, 0, P * P * sizeof(double));
    memset(G, 0, P * sizeof(double));

    svd_hestenes(P, P, MassInStroka, U, V, G);

    double* GKrest = new double[P];
    memset(GKrest, 0, P * sizeof(double));
    double max = 0;
    for (int i = 0; i < P; i++)
        if (max < G[i]) max = G[i];
    double koef = 0;
    if (param.snr == 0) koef = 0.1 / 100 * max;
    else koef = abs(param.snr) / 100 * max;
    for (int i = 0; i < P; i++)
    {
        if (G[i] != 0 && G[i] > koef)
            GKrest[i] = 1 / G[i];
        else
            GKrest[i] = G[i];
    }

    //псевдообратная матрица
    double** UKrest = matrix(P, P);
    double** VKrest = matrix(P, P);
    /*memset(UKrest, 0, P* P * sizeof(double));
    memset(VKrest, 0, P* P * sizeof(double));*/
    int s = 0;
    int x = 0;
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            UKrest[i][j] = U[s];
            VKrest[i][j] = V[x];
            s++;
            x++;
        }
    }

    double** UKrest_transpose = matrix(P, P);
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            UKrest_transpose[i][j] = UKrest[j][i];
        }
    }

    double** Akrest = matrix(P, P);

    double** GGkrest = matrix(P, P);
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            if (i != j) GGkrest[i][j] = 0;
            else GGkrest[i][j] = GKrest[i];
        }
    }

    double** Composition_Two = Composition_Matrix_Two(GGkrest, UKrest_transpose, P);
    Akrest = Composition_Matrix_Two(VKrest, Composition_Two, P);

    double** buf1 = matrix(P, P);

    double** A = matrix(P, P);

    iter = 0;
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++)
        {
            A[i][j] = MassInStroka[iter];
            iter++;
        }
    iter = 0;

    buf1 = Composition_Matrix_Two(Akrest, A, P);

    double** buf2 = matrix(P, P);
    buf2 = Composition_Matrix_Two(A, buf1, P);

    double dbl_index = 0;

    iter = 0;
    Sigma.clear();
    double M = param.outputT;

    InvR.clear();
    InvR.resize(param.p);
    for (int i = 0; i < param.p; i++)
    {
        InvR[i].resize(param.p);
        for (int j = 0; j < param.p; j++)
        {
            if (i == j)InvR[i][j] = Akrest[i][j]/* + DoubleRand(0.001, 0.05)*/;
            else
                InvR[i][j] = Akrest[i][j];
        }
    }
    double sampling_period = 1. / samplingFrequency;// период дискретизации

    int ti = 0;
    double startTimestamp = 0;
    //продолжительность
    double Duration = param.duration;
    while (dbl_index < startTimestamp + param.koef_dur*Duration)
    {
        vector<double>r;
        r.resize(param.p);

        int kol = 0;
        int Pi = 0;
        for (int i = ti; i < param.p + ti; i++)
        {
            double counter = 0;
            for (int j = 0; j < InputSignal.size(); j++)
            {
                if (j < M && j + i < InputSignal.size())
                    counter += InputSignal[j] * InputSignal[j + i];
                else continue;
            }
            r[Pi] = counter / M;
            Pi++;
        }

        vector<double> exp2 = Composition_Matrix_Stroka(InvR, r, param.p, param.p);

        double buf = 0;
        for (int i = 0; i < exp2.size(); i++)
        {
            buf += r[i] * exp2[i];
        }

        Sigma.push_back(buf);
        dbl_index += sampling_period;
        ti++;
    }
}

QVector<QVector<double>> GenerationSignal::GenerateShortSignal(ParamSignal& param)
{
    Ymin = false;
    //частота дискретизации
    double samplingFrequency = param.fdisk;
    //метка времени начала
    double startTimestamp = 0;
    //продолжительность
    double Duration = param.duration;
    dur=param.duration;;
    //начальная фаза
    double startPhase = 0;
    double nSamples = 0;
    //скорость передачи данных
    double Bitrate = param.bitrate;
    //дополнительный параметр
    double additionalParameter = 0;

    int samples_in_bit = (int)(samplingFrequency / Bitrate);// кол-во отсчетов на 1 бит
    double sampling_period = 1. / samplingFrequency;// период дискретизации

    vector<double> sign;
    vector<double> Bit;

    double dbl_index = startTimestamp;
    int iter = 0;
    double phase;
    double Ampl = 1.;
    param.koef_dur = 1;
    bitPosl.clear();
    InputSignal.clear();
    InputSignal.resize(KolSour);
    bitPosl.resize(KolSour);
    for(int i=0;i<KolSour;i++) {
        while (dbl_index < startTimestamp + param.koef_dur * Duration) {
            double bit = RandomBit(0.5);

            for (int n_sample = 0; n_sample < samples_in_bit; n_sample++) {
                if (bit == 0) phase = 0.;
                    //deltaF = 100000;
                else phase = 3.1415926535897;
                //deltaF=0;
                InputSignal[i].push_back(Ampl * sin(2 * 3.1415926535897 * (param.f0) * dbl_index + phase));
                sign.push_back(Ampl * cos(2 * 3.1415926535897 * param.f0 * dbl_index));
                dbl_index += sampling_period;
                Bit.push_back((int) bit);
                bitPosl[i].push_back(bit);
            }
        }
        dbl_index = 0;
        //ClearSignal = InputSignal;
        //addNoise(sign, 10, param);
        //addNoise(InputSignal, 10, param);
    }

    QVector<QVector<double>> XY;
    XY.resize(2);
    y.resize(InputSignal[0].size());
    x.resize(InputSignal[0].size());
    dbl_index=0;
    for(int i=0;i<InputSignal[0].size();i++)
    {
        y[i]=InputSignal[0][i];
        x[i]=dbl_index += sampling_period;
    }
    XY[0]=x;
    XY[1]=y;
    return XY;
}

void GenerationSignal::coherentSumm(vector<vector<double>> &sgs, vector<double> &result)
{
    if (sgs.size() == 0)
    {
        return;
    }
    else if (sgs.size() == 1)
    {
        result = sgs[0];
        return;
    }
    result.clear();

    result = sgs[0];
    for (unsigned int i = 1; i < sgs.size(); i++)
    {
        for (unsigned int j = 0; j < sgs[i].size(); j++)
        {
            result[j] += sgs[i][j];
        }
    }
}

void GenerationSignal::transformSignal(vector<double>& base_signal, double delay, double duration, double fshift, double scale, double SNR, vector<double>& ret_sig, ParamSignal& param)
{
    double time_shift = delay * param.fdisk;
    int size_of_sample = (int)(duration * param.fdisk);
    if (time_shift < base_signal.size() && (time_shift + size_of_sample) < base_signal.size())
    {
        for (int i = 0; i < (size_of_sample); i++)
        {
            ret_sig.push_back(
                    base_signal[(unsigned int)((double)i / scale + time_shift)] *
                   cos( 2. * M_PI * fshift * (1. / param.fdisk) * (double)i));
        }
    }
    addNoise(ret_sig, SNR, param);
}

QVector<QVector<double>> GenerationSignal::GenerateLongSignal(ParamSignal& param, int& size)
{
    vector<vector<vector<double>>> LongSignal;
    LongSignal.clear();
    Ymin = false;
    //частота дискретизации
    double samplingFrequency = param.fdisk;
    //метка времени начала
    double startTimestamp = 0;
    //продолжительность
    double Duration = param.duration;
    //начальная фаза
    double startPhase = 0;
    double nSamples = 0;
    //скорость передачи данных
    double Bitrate = param.bitrate;
    //дополнительный параметр
    double additionalParameter = 0;

    int samples_in_bit = (int)(samplingFrequency / Bitrate);// кол-во отсчетов на 1 бит
    double sampling_period = 1. / samplingFrequency;// период дискретизации
    vector<double> sign;
    sign.clear();
    //vector<double> Bit;

    double dbl_index = startTimestamp;
    int iter = 0;
    double phase;
    double deltaF = 0.;
    double Ampl = 1.;
    dur=param.koef_dur *Duration;
    param.koef_dur = 5;

    for(int i=0;i<KolSour;i++)
    {
        for(int j=0;j<KolSat;j++)
        {
            sdvigTime.push_back(1. + 0.5 * ((double)rand() / RAND_MAX - 0.5));
        }
    }

    for(int i=0;i<KolSour;i++)
    {
        for(int j=0;j<KolSat;j++)
        {
            sdvigFreq.push_back(DoubleRand(1000,-1000));
        }
    }

    SummSignal.resize(KolSat);

    kol_bit = 0;

    LongSignal.resize(KolSat);
    LongSigma.resize(KolSat);
    double d = Duration / 100;
    for(int i=0;i<KolSat;i++)
    {
        sign.clear();
        LongSignal[i].resize(KolSour);
        for(int j=0;j<KolSour;j++)
        {
            transformSignal(InputSignal[j], sdvigTime[i + KolSat * j] * d, d, sdvigFreq[i + KolSat * j], 1, param.snr, LongSignal[i][j], param);
        }

        /*vector<vector<double>> vc;
        vc.resize(KolSour);
        for(int j=0;j<KolSour;j++)
        {
            GetSigma(LongSignal[i][j], sign, vc[j], param);
        }

        if (KolSour == 1)
        {
            SummSignal[i] = vc[0];
        }

        for(int j=1;j<KolSour;j++)
        {
            if (j == 1) {
                vector<vector<double>> sgs;
                sgs.resize(2);
                sgs[0] = vc[0];
                sgs[1] = vc[j];
                coherentSumm(sgs, SummSignal[i]);
            }
            else {
                vector<vector<double>> sgs;
                sgs.resize(2);
                sgs[0] = SummSignal[i];
                sgs[1] = vc[j];
                coherentSumm(sgs, SummSignal[i]);
            }
        }*/

       // GetSigma(SummSignal[i], sign, LongSigma[i], param);
    }
    bitPosl.clear();

    Ymin = true;

    QVector<QVector<double>> XY;
    XY.resize(2);
    y.resize(LongSignal[0][0].size());
    x.resize(LongSignal[0][0].size());
    dbl_index=0;
    for(int i=0;i<LongSignal[0][0].size();i++)
    {
        y[i]=LongSignal[0][0][i];
        x[i]=dbl_index += sampling_period;
    }
    XY[0]=x;
    XY[1]=y;


    return XY;
}

void GenerationSignal::newFFT(vector<complex<double>>& in, int direction)
{
    //out = in;
    unsigned int pts = 2;
    while (in.size() > pts)
    {
        pts *= 2;
    }

    int pts_to_add = pts - in.size();

    for (int ttt = 0; ttt < pts_to_add; ttt++)
    {
        in.push_back(complex<double>(0, 0));
    }
    int n = in.size();

    int i, j, istep;
    int m, mmax;
    double r, r1, theta, w_r, w_i, temp_r, temp_i;

    r = M_PI * direction;
    j = 0;

    for (i = 0; i < n; i++)
    {
        if (i < j)
        {
            temp_r = in[j].real();
            temp_i = in[j].imag();
            in[j] = in[i];
            in[i] = complex<double>(temp_r, temp_i);
        }
        m = n >> 1;
        while (j >= m)
        {
            j -= m;
            m = (m + 1) / 2;
        }
        j += m;
    }
    mmax = 1;
    while (mmax < n)
    {
        istep = mmax << 1;
        r1 = r / (double)mmax;
        for (m = 0; m < mmax; m++)
        {
            theta = r1 * m;
            w_r = (double)cos((double)theta);
            w_i = (double)sin((double)theta);
            for (i = m; i < n; i += istep)
            {
                j = i + mmax;
                temp_r = w_r * in[j].real() - w_i * in[j].imag();
                temp_i = w_r * in[j].imag() + w_i * in[j].real();
                in[j] = complex<double>((in[i].real() - temp_r), (in[i].imag() - temp_i));
                in[i] += complex<double>(temp_r, temp_i);
            }
        }
        mmax = istep;
    }
    if (direction > 0)
    {
        double sqn = n;
        for (i = 0; i < n; i++)
        {
            in[i] /= sqn;
        }
    }
}

void GenerationSignal::GetSPM(vector<complex<double>>& Pw, vector<complex<double>>& PShum, vector<double>& PSh, vector<double>& PSg, ParamSignal& param)
{
    Ymin = false;
    //частота дискретизации
    double samplingFrequency = param.fdisk;
    //метка времени начала
    double startTimestamp = 0;
    //продолжительность
    double Duration = param.duration;
    //начальная фаза
    double startPhase = 0;
    double nSamples = 0;
    //скорость передачи данных
    double Bitrate = param.bitrate;
    //дополнительный параметр
    double additionalParameter = 0;

    int samples_in_bit = (int)(samplingFrequency / Bitrate);// кол-во отсчетов на 1 бит
    double sampling_period = 1. / samplingFrequency;// период дискретизации

    Pw.resize(sizeSPM);
    PShum.resize(sizeSPM);

    vector<double> Re(sizeSPM);
    vector<double> Im(sizeSPM);

    vector<double> ReS(sizeSPM);
    vector<double> ImS(sizeSPM);

    PSg.resize(sizeSPM);
    PSh.resize(sizeSPM);

    double kol_Usr=10;
    for (int i = 0; i < kol_Usr; i++)
    {
        vector<double> IsslSignal;
        vector<double> IssSigma;
        vector<double> sign;
        vector<double> Bit;

        double dbl_index = startTimestamp;
        int iter = 0;
        double phase;
        double Ampl = 1.;
        param.koef_dur = 1;

        while (dbl_index < startTimestamp + param.koef_dur * Duration)
        {
            double bit = RandomBit(0.5);
            for (int n_sample = 0; n_sample < samples_in_bit; n_sample++)
            {
                if (bit == 0) phase = 0.;
                    //deltaF = 100000;
                else phase = 3.1415926535897;
                //deltaF=0;
                IsslSignal.push_back(Ampl * sin(2 * 3.1415926535897 * (param.f0)*dbl_index + phase));
                sign.push_back(Ampl * cos(2 * 3.1415926535897 * param.f0 * dbl_index));
                dbl_index += sampling_period;
                Bit.push_back((int)bit);
            }
        }
        vector<double> buf = IsslSignal;
        vector<double> bufsign = sign;
        addNoise(bufsign, 10, param);
        addNoise(buf, 10, param);

        GetSigma(IsslSignal, sign, IssSigma, param);

        vector<complex<double>> IsslH;

        for (int i = 0; i < IssSigma.size(); i++)
        {
            double Re = IssSigma[i];
            complex<double> comp(Re, 0);
            IsslH.push_back(comp);
        }
        newFFT(IsslH, -1);

        vector<complex<double>> IsslHShum;

        for (int i = 0; i < IssSigma.size(); i++)
        {
            double Re = SHUM[i];
            complex<double> comp(Re, 0);
            IsslHShum.push_back(comp);
        }
        newFFT(IsslHShum, -1);

        vector<vector<double>> cor;
        cor.resize(1);
        for(int a=0;a< IsslHShum.size();a++)
            cor[0].push_back(abs(IsslHShum[a]));

        /**GraphPen.clear();
        GraphPen.push_back(new CPen(PS_SOLID, 3, RGB(178, 102, 255)));

        GraphType type = GraphType::Graphic;
        DrawGraph2(cor, startTimestamp, IsslHShum.size(), GraphPen, PicDc_LongSignal, Pic_LongSignal, type);*/



        for (int k = 0; k < IsslH.size(); k++)
        {
            Re[k] += IsslH[k].real();
            Im[k] += IsslH[k].imag();
            PSg[k] += (abs(IsslH[k]));

            ReS[k] += IsslHShum[k].real();
            ImS[k] += IsslHShum[k].imag();
            PSh[k] += (abs(IsslHShum[k]));
        }
    }

    for (int i = 0; i < Re.size(); i++)
    {
        double re = Re[i]/kol_Usr;
        double im = Im[i]/kol_Usr;
        complex<double> comp(re, im);
        Pw[i]=comp;

        double res = ReS[i]/kol_Usr;
        double ims = ImS[i]/kol_Usr;
        complex<double> compS(res, ims);
        PShum[i] = compS;
    }

    for (int i = 0; i < Pw.size(); i++)
    {
        Pw[i] /= kol_Usr;
        PShum[i] /= kol_Usr;
    }
}

void GenerationSignal::GetK(vector<complex<double>>& K, vector<double>& SPNSigma, vector<double>& SPNShum, vector<complex<double>>& HSopr, vector<double>& HMod)
{
    double alpha = 0.5;
    K.resize(HSopr.size());
    double shum = 0;
    for (int i = 0; i < SPNShum.size(); i++)
    {
        shum += SPNShum[i];
    }
    shum/= SPNShum.size();

    for (int i = 0; i < K.size(); i++)
    {
        double koef = HMod[i]/**HMod[i]*/ + alpha * shum/ SPNSigma[i];
       // if (i == 0) K[i] = 0;
        //else
        K[i] =HSopr[i]/koef;
    }
}

void GenerationSignal::newcorrelate(vector<complex<double>>& base_signal, vector<complex<double>>& analyzed_signal, vector<double>& PSg, vector<double>& PSh,vector<complex<double>>& correlation)
{

    newFFT(base_signal, -1);

    vector<complex<double>> K;

    vector<complex<double>> HSopr;
    vector<double> ModH;
    for (int i = 0; i < base_signal.size(); i++)
    {
        //if (i > 30 && i < base_signal.size() - 29)HSopr.push_back(0.);
        //else
        HSopr.push_back(conj(base_signal[i]));
        ModH.push_back(abs(base_signal[i]));
        //double mod = base_signal[i].real() * base_signal[i].real() + base_signal[i].imag() + base_signal[i].imag();
        //ModH.push_back(mod);
    }

    GetK(K, PSg, PSh, HSopr, ModH);


    vector<double> cor;
    vector<double> x;
    for (int i = 0; i < K.size(); i++) {
        cor.push_back((K[i].real()));
        x.push_back(i);
    }

    info.MassOtrisovka.push_back(cor);
    info.MassOtshetX.push_back(x);
    listt.push_back("К");

    for (unsigned int n = 0; n < (analyzed_signal.size() - K.size()); n++)
    {
        vector<complex<double>> Hu;
        for (int i = n; i < K.size() + n; i++)
            Hu.push_back(analyzed_signal[i]);
        newFFT(Hu, -1);
        vector<complex<double>> Peremn;
        for (int i = 0; i < K.size(); i++)
        {
            Peremn.push_back(K[i] * Hu[i]);
        }
        newFFT(Peremn,1);
        correlation.push_back(Peremn[0]);
    }
}

void GenerationSignal::correlate(vector<double>& base_signal, vector<double>& analyzed_signal, vector<double>& correlation, double& sredSigma, double& sredLongSigma)
{
   for (int n = -(int)(analyzed_signal.size() - 1); n <= (int)(base_signal.size() - 1); n++)
    {
        double counter = 0;
        for (unsigned int m = 0; m < analyzed_signal.size(); m++)
        {
            if ((m + n) >= 0 && (m + n) < base_signal.size())
            {
                counter += (analyzed_signal[m]) * (base_signal[m + n] );
            }
            else continue;
        }
        counter /= (double)analyzed_signal.size();//cnt;//
        correlation.push_back(counter);
    }
}

InfoList GenerationSignal::Calculate(ParamSignal& param)
{
    for(int k=0;k<KolSat;k++) {

        vector<double> correlation;

        if(k+1==KolSat)
        {
            double sredLongSigma1 = 0, sredLongSigma2 = 0;
            for (int i = 0; i < SummSignal.size(); i++) {
                sredLongSigma1 += SummSignal[k][i];
                sredLongSigma2 += SummSignal[0][i];
            }
            sredLongSigma1 /= SummSignal[k].size();
            sredLongSigma2 /= SummSignal[0].size();
            correlate(SummSignal[k],SummSignal[0],correlation, sredLongSigma1, sredLongSigma2);
        }
        else
        {
            double sredLongSigma1 = 0, sredLongSigma2 = 0;
            for (int i = 0; i < SummSignal.size(); i++) {
                sredLongSigma1 += SummSignal[k][i];
                sredLongSigma2 += SummSignal[k+1][i];
            }
            sredLongSigma1 /= SummSignal[k].size();
            sredLongSigma2 /= SummSignal[k+1].size();
            correlate(SummSignal[k],SummSignal[k+1],correlation, sredLongSigma1, sredLongSigma2);
        }

        info.MassOtrisovka.push_back(correlation);
        vector<double> x;
        double sampling_period = 1. / param.fdisk;// период дискретизации
        double t = 0;
        for (int i = 0; i < correlation.size(); i++) {
            x.push_back(t);
            t += sampling_period;
        }
        info.MassOtshetX.push_back(x);
        listt.push_back("Корреляция");
    }

    return info;
}

GenerationSignal::~GenerationSignal()
{
    Sigma.clear();
    LongSigma.clear();
    InvR.clear();
    InputSignal.clear();
    LongInputSignal.clear();
    ClearSignal.clear();
    SHUM.clear();
    bitPosl.clear();
    x.clear();
    y.clear();
    info.MassOtshetX.clear();
    info.MassOtrisovka.clear();
    listt.clear();
}