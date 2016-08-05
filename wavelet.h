#ifndef WAVELET_H
#define WAVELET_H
#include <vector>
#include "db.h"
namespace Wavelet{
    using std::vector;
    struct C_L 
    {
        vector<double> C;
        vector<int> L;
    };
	/**************************
	* The structure is organized as:
    * C      = [app. coef.(N)|det. coef.(N)|... |det. coef.(1)]
    * L(1)   = length of app. coef.(N)
    * L(i)   = length of det. coef.(N-i+2) for i = 2,...,N+1
    * L(N+2) = length(X).
	**************************/
	
    struct WaveFilter{
        vector<double> Low;        //��ͨ�˲���
        vector<double> High;       //��ͨ�˲���
    };
    struct WaveCoeff{
        vector<double> app;   //����ϵ��
        vector<double> det;   //ϸ��ϵ��
    };
	
    const static WaveFilter sym4_d = {vector<double>(sym4_Lo_D, sym4_Lo_D + 8), vector<double>(sym4_Hi_D, sym4_Hi_D + 8)};
    const static WaveFilter sym4_r = {vector<double>(sym4_Lo_R, sym4_Lo_R + 8), vector<double>(sym4_Hi_R, sym4_Hi_R + 8)};
	
	const static WaveFilter db24_d = {vector<double>(db24_Lo_D, db24_Lo_D + 48), vector<double>(db24_Hi_D, db24_Hi_D + 48)};
    const static WaveFilter db24_r = {vector<double>(db24_Lo_R, db24_Lo_R + 48), vector<double>(db24_Hi_R, db24_Hi_R + 48)};

    const WaveFilter& WFilters(const char* strWaveName,const char d_or_r); 
    C_L WaveDec(const vector<double>& signal, const int nMaxLevel, const char* strWaveName); //һά���С���ֽ�
    WaveCoeff DWT(const vector<double>& signal,const vector<double>& Lo_D,const vector<double>& Hi_D);//һά����С���ֽ�
    vector<double> WRCoef(const char a_or_d,const vector<double>& C,   const vector<int>& L,  const char* strWaveName, const int nLevel);//С���ع�ϵ��
    vector<double> AppCoef(const vector<double>& C, const vector<int>& L, const char* strWaveName,  const int nLevel );//����ϵ��
    vector<double> DetCoef( const vector<double>& C,  const vector<int>& L,  const int nLevel );//ϸ��ϵ��
    
	//upsample and convolution
    vector<double> UpsConv1(const vector<double>& signal,const vector<double>& filter,const int nLen, const char* strMode = "db" );  //�ϲ����;�������
    vector<double> Conv( const vector<double>& vecSignal,  const vector<double>& vecFilter); //��������
    vector<double> IDWT( const vector<double>& app,   const vector<double>& det,  const vector<double>& Lo_R,  const vector<double>& Hi_R,   const int nLenCentral ); //һά����С���ع�
    vector<double> WExtend(const vector<double>& signal,  const int nLenExt, const char* mode = "db" );//С������
    vector<double> WConv1( const vector<double>& signal, const vector<double>& filter, const char* shape = "valid" );
}
#endif