//
// Created by Ghost Gloomy on 2020/9/14.
//

#ifndef CMAKE_TN_FNDS_MAIN_H
#define CMAKE_TN_FNDS_MAIN_H

#include <cmath>
#include "world/fft.h"
#include "world/matlabfunctions.h"
#include <iostream>

using std::min;
using std::max;

typedef uint32_t DWORD;

#define PI 3.1415926535897932384
#define DEFAULT_F0 500.0//150.0   tn_fnds v0.0.3 大きくしてみたら無声子音がきれいになった
#define FLOOR_F0 90.0

// 71は，fs: 44100においてFFT長を2048にできる下限．
// 70 Hzにすると4096点必要になる．
// DEFAULT_F0は，0.0.4での新機能．調整の余地はあるが，暫定的に決定する．

double* wavread(char* filename, int* fs, int* Nbit, int* waveLength, int* offset, int* endbr);

// F0推定法 DIO : Distributed Inline-filter Operation
void dio(double *x, int xLen, int fs, double framePeriod, double *timeAxis, double *f0);

int getSamplesForDIO(int fs, int xLen, double framePeriod);

// スペクトル包絡推定法 STAR : Synchronous Technique and Adroit Restoration
int getFFTLengthForStar(int fs);

//tn_fnds v0.0.3 にて追加
int pt101(double *x, int xLen, int fs, double *timeAxis, double *f0, double ***residualSpecgram, int **residualSpecgramLength, int *residualSpecgramIndex);

//tn_fnds v0.0.4 にて追加
void PulseResidualWindow(double **residualSpecgram, int *residualSpecgramLength, int pCount);

//tn_fnds v0.0.3 にて追加
void synthesisPt101(double fixedDefault_f0, double *f0, int tLen, double **aperiodicity, int *ResidualSpecgramLength, int *fixedResidualSpecgramIndex, double *volume, int fftl, double framePeriod, int fs, double *synthesisOut, int xLen);

#endif //CMAKE_TN_FNDS_MAIN_H
