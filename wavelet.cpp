/************************************************************************/
/* wavelet.cpp
* Author: Collin
* Date: 2008/12/01
*/
/************************************************************************/
#include <vector>
#include <string>
#include <iostream>
#include "wavelet.h"
using namespace std;
using namespace Wavelet;
C_L Wavelet::WaveDec(const vector<double>& signal,
                    const int nMaxLevel,                 //小波分解层数
                    const char* strWaveName
                    )
{
    const WaveFilter& filters = WFilters(strWaveName, 'd');      //获取小波滤波器组
    int len = signal.size();                 //信号的长度
    C_L cl;
    cl.L.push_back(len);
    WaveCoeff waveCoeff;
    waveCoeff.app = signal;
    vector<double>::iterator itC;
    vector<int>::iterator itL;
    for (int i = 0; i < nMaxLevel; ++i){
        waveCoeff = DWT(waveCoeff.app, filters.Low, filters.High);  //进行单层小波变化
        itC = cl.C.begin();
        cl.C.insert(itC, waveCoeff.det.begin(), waveCoeff.det.end());
        itL = cl.L.begin();
        cl.L.insert(itL, waveCoeff.det.size());
    }
    itC = cl.C.begin();
    cl.C.insert(itC, waveCoeff.app.begin(), waveCoeff.app.end());
    itL = cl.L.begin();
    cl.L.insert(itL, waveCoeff.app.size());
    return cl;
}

vector<double> Wavelet::WRCoef(const char a_or_d,
                               const vector<double>& C, 
                               const vector<int>& L, 
                               const char* strWaveName,
                               const int nLevel
                               )
{
    vector<double> Coef;
    const WaveFilter& filter = WFilters(strWaveName, 'r');
    int nMax = L.size() - 2;
    int nMin;
    char type = tolower(a_or_d);
    if ('a' == type)
        nMin = 0;
    else if ('d' == type)
        nMin = 1;
    else {
        cerr << "bad parameter: a_or_d: "<< a_or_d << "\n";
        exit(1);
    }
    if (nLevel < nMin || nLevel > nMax){
        cerr << "bad parameter for level\n";
        exit(1);
    }
    vector<double> F1;
    switch (type){
        case 'a':
            Coef = AppCoef(C, L, strWaveName, nLevel);
            if (0 == nLevel)
                return Coef;
            F1 = filter.Low;
            break;
        case 'd':
            Coef = DetCoef(C, L, nLevel);
            F1 = filter.High;
            break;
        default:
            ;
    }
    int iMin = L.size() - nLevel;
    Coef = UpsConv1(Coef, F1, L[iMin], "db");
    for (int k = 1; k < nLevel; ++k){
        Coef = UpsConv1(Coef, filter.Low, L[iMin + k], "db");
    }
    return Coef;
}
vector<double> Wavelet::UpsConv1(const vector<double>& signal,
                                const vector<double>& filter, 
                                const int nLen, 
                                const char* strMode
                                )
{
    //implement dyadup(y,0)
    vector<double> y(2 * signal.size() - 1);
    y[0] = signal[0];
    for (int i = 1; i < signal.size(); ++i){
        y[2*i - 1] = 0;
        y[2*i] = signal[i];
    }
    y = Conv(y, filter);

    //extract the central portion
    vector<double>::iterator it = y.begin();
    return vector<double>(it + (y.size() - nLen) / 2, it + (y.size() + nLen) / 2);
}
vector<double> Wavelet::Conv(const vector<double>& vecSignal, const vector<double>& vecFilter){
    vector<double> signal(vecSignal);
    vector<double> filter(vecFilter);
    if (signal.size() < filter.size())
        signal.swap(filter);    
    int lenSignal = signal.size();
    int lenFilter = filter.size();
    vector<double> result(lenSignal + lenFilter - 1);
    for (int i = 0; i < lenFilter; i++){
        for (int j = 0; j <= i; j++)
            result[i] += signal[j] * filter[i - j];
    }
    for (int i = lenFilter; i < lenSignal; i++){
        for (int j = 0; j <lenFilter; j++)
            result[i] += signal[i - j] * filter[j];
    }
    for (int i = lenSignal; i < lenSignal + lenFilter - 1; i++){
        for (int j = i - lenSignal + 1; j < lenFilter; j++)
            result[i] += signal[i - j] * filter[j];
    }
    return result;    
}
vector<double> Wavelet::DetCoef(const vector<double>& C,
                       const vector<int>& L,
                       const int nLevel
                       )
{
    if (nLevel < 1 || nLevel > L.size() - 2){
        cerr << "bad level parameter\n";
        exit(1);
    }

    int nlast = 0, nfirst = 0;
    vector<int>::const_reverse_iterator it = L.rbegin();
    ++it;
    for (int i = 1; i < nLevel; ++i){
        nlast += *it;
        ++it;
    }
    nfirst = nlast + *it;
    return vector<double>(C.end() - nfirst, C.end() - nlast);
}

/***********一维单层离散小波分解************/
WaveCoeff Wavelet::DWT(const vector<double>& signal,
            const vector<double>& Lo_D,
            const vector<double>& Hi_D
            )
{
    int nLenExt = Lo_D.size() - 1;    //信号边界延拓的长度为滤波器的长度减1
    vector<double> y;
    y = WExtend(signal, nLenExt, "db");  //信号边界延拓，matlab的边界延拓采用对称延拓方式
    vector<double> z;
    z = WConv1(y, Lo_D, "valid");         //一维卷积计算近似系数
    WaveCoeff coeff;
	/****二抽法获取近似系数***/
    for (int i = 1; i < z.size(); i += 2){
        coeff.app.push_back(z[i]);
    }
    z = WConv1(y, Hi_D, "valid");        //一维卷积计算细节系数
	/****二抽法获取细节系数***/
    for (int i = 1; i < z.size(); i += 2){
        coeff.det.push_back(z[i]);
    }
    return coeff;
}

/**************************构建小波变换滤波器组*************************************/
const WaveFilter& Wavelet::WFilters(const char* strWaveName,
                           const char d_or_r
                           )
{
    char type = tolower(d_or_r);           //tolower为库函数，功 能: 把字符转换成小写字母,非字母字符不做出处理   头文件：在VC6.0可以是ctype.h或者stdlib.h，常用ctype.h
    if (!strcmp(strWaveName, "db24")){
        switch(type){
        case 'd':
            return Wavelet::db24_d;
            break;
        case 'r':
            return Wavelet::db24_r;
            break;
        default:
            cerr << "bad parameter for d_or_r\n";
            exit(1);
        }
    }
    else {
        cerr << "not implement \n";
        exit(1);
    }
}

vector<double> Wavelet::AppCoef(
                                const vector<double>& C,
                                const vector<int>& L,
                                const char* strWaveName,
                                const int nLevel
                                )
{
    int nMaxLevel = L.size() - 2;
    if (nLevel < 0 || nLevel > nMaxLevel){
        cerr << "bad parameter for level\n";
        exit(1);
    }
    const WaveFilter& filters = WFilters(strWaveName, 'r');
    vector<double> app(C.begin(), C.begin() + L[0]); //app for the last level 
    vector<double> det;
    for (int i = 0; i < nMaxLevel - nLevel; ++i){
        det = DetCoef(C, L, nMaxLevel - i);
        app = IDWT(app, det, filters.Low, filters.High, L[i + 2]);
    }
    return app;
}
vector<double> Wavelet::IDWT(
                    const vector<double>& app, 
                    const vector<double>& det, 
                    const vector<double>& Lo_R, 
                    const vector<double>& Hi_R, 
                    const int nLenCentral
                    )
{
    vector<double> app1, app2;
    app1 = UpsConv1(app, Lo_R, nLenCentral, "sym");
    app2 = UpsConv1(det, Hi_R, nLenCentral, "sym");
    for (int i = 0; i < nLenCentral; ++i){
        app1[i] += app2[i];
    }
    return app1;
}

/***************************信号的边界延拓*********************************/
vector<double> Wavelet::WExtend(
                       const vector<double>& signal,
                       const int nLenExt,
                       const char* mode
                       )
{
    int signalLen = signal.size();                         //信号的长度
    vector<double> result(signalLen + 2 * nLenExt);        //延拓之后的长度为 信号长度+2*（滤波器长度-1）
    for (int i = 0, idx = nLenExt; idx < signalLen + nLenExt; ++i, ++idx){
        result[idx] = signal[i];                           //原信号放在 nlenExt~signalLen + nLenExt-1之间，共signalLen个
    }
    for (int idx = nLenExt - 1, bFlag = 1, signalIdx = 0; idx >= 0; --idx){
        result[idx] = signal[signalIdx];                  //左延拓，将原信号对称到左边
		/*********************
		*当滤波器的长度大于信号的长度时，左延拓的元素个数将大于原信号个数，这时在对称一次之后，继续左对称过去
		*|....a(n-2),a(n-1)|a(n-1)....a(2),a(1),a(0)|a(0),a(1),a(2)......a(n-1)|.....|
		**********************/
        if (bFlag && ++signalIdx == signalLen){
            bFlag = 0;
            signalIdx = signalLen - 1;
        }
        else if (!bFlag && --signalIdx == -1) {
            bFlag = 1;
            signalIdx = 0;
        }
    }
    for (int idx = nLenExt + signalLen, bFlag = 0, signalIdx = signalLen - 1; idx < 2 * nLenExt + signalLen; ++idx){
        result[idx] = signal[signalIdx];                   //右延拓，将原信号对称到右边
		/*********************
		*当滤波器的长度大于信号的长度时，右延拓的元素个数将大于原信号个数，这时在对称一次之后，继续右对称过去
		*|......|a(0),a(1),a(2)......a(n-1)|a(n-1)，a(n-2)..........a(2),a(1),a(0)|a(0),a(1),a(2)...|
		**********************/
        if (bFlag && ++signalIdx == signalLen){
            bFlag = 0;
            signalIdx = signalLen - 1;
        }
        else if (!bFlag && --signalIdx == -1) {
            bFlag = 1;
            signalIdx = 0;
        }
    }
    return result;
}
vector<double> Wavelet::WConv1(
                    const vector<double>& signal,
                    const vector<double>& filter,
                    const char* shape
                    )
{
    vector<double> y;
    y = Conv(signal, filter);
    int nLenExt = filter.size() - 1;
    return vector<double>(y.begin() + nLenExt, y.end() - nLenExt);
}
