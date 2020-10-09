//tn_fnds v0.0.5   2012/3/17
//追加されているコメントには誤りがあるかもしれません。
#include "main.h"

#include <stdio.h>  // for debug
#include <stdlib.h>

void synthesisPt101(double fixedDefault_f0, double* f0, int tLen, double** aperiodicity, int* ResidualSpecgramLength, int* fixedResidualSpecgramIndex, double* volume, int fftl, double framePeriod, int fs, double* synthesisOut, int xLen) {
    int i, j;
    double currentTime = 0;
    int currentPosition = 0;  //currentTime / framePeriod;
    int currentFrame = 0;

    for (i = 0;; i++) {
        for (j = 0; j < ResidualSpecgramLength[fixedResidualSpecgramIndex[currentFrame]]; j++) {
            if (j + currentPosition >= xLen) {
                break;
            }
            synthesisOut[max(0, j + currentPosition)] += aperiodicity[fixedResidualSpecgramIndex[currentFrame]][j] * volume[currentFrame];
        }

        // 更新
        currentTime += 1.0 / (f0[currentFrame] == 0.0 ? fixedDefault_f0 : f0[currentFrame]);  //時刻を1周期分進める
        currentFrame = ( int )(currentTime / (framePeriod / 1000.0) + 0.5);                   //次に継ぎ足すデータ位置は次の時刻に最も近いフレーム
        currentPosition = ( int )(currentTime * ( double )fs);                                //周期単位で継ぎ足していく
        if (j + currentPosition >= xLen || currentFrame >= tLen) {
            break;
        }
    }
    return;
}
