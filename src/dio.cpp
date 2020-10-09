//tn_fnds v0.0.5   2012/3/17
//追加されているコメントには誤りがあるかもしれません。
#include "main.h"

// 内部関数(ユーザは触らないほうが良い)
//tn_fnds v0.0.4
void rawEventByDio(double boundaryF0, double fs, fft_complex* xSpec, int xLength, int fftl, double shiftTime, double f0Floor, double f0Ceil, double* timeAxis, int tLen,
                   double* f0Deviations, double* interpolatedF0,
                   fft_plan& forwardFFT, fft_plan& inverseFFT, double* equivalentFIR, fft_complex* eSpec);
void zeroCrossingEngine(double* x, int xLen, double fs,
                        double* eLocations, double* iLocations, double* intervals, int* iLen);
void filterForDecimate(double* x, int xLen, double* y, int r);
void nuttallWindow(int yLen, double* y);
//tn_fnds v0.0.4
void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double** f0Map,
                    double** stabilityMap, double* bestF0, double* f0);
void histc(double* x, int xLen, double* y, int yLen, int* index);
void filterForDecimate(double* x, int xLen, double* y, int r) {
    double w[3], wt;
    w[0] = w[1] = w[2] = 0.0;
    double a[3], b[2];  // �t�B���^�W�� (r�ˑ�)

    switch (r) {
    case 11:  // fs : 44100
        a[0] = 2.450743295230728;
        a[1] = -2.06794904601978;
        a[2] = 0.59574774438332101;
        b[0] = 0.0026822508007163792;
        b[1] = 0.0080467524021491377;
        break;
    case 12:  // fs : 48000
        a[0] = 2.4981398605924205;
        a[1] = -2.1368928194784025;
        a[2] = 0.62187513816221485;
        b[0] = 0.0021097275904709001;
        b[1] = 0.0063291827714127002;
        break;
    case 8:  // fs : 32000
        a[0] = 2.2357462340187593;
        a[1] = -1.7780899984041358;
        a[2] = 0.49152555365968692;
        b[0] = 0.0063522763407111993;
        b[1] = 0.019056829022133598;
        break;
    case 6:  // fs : 24000 and 22050
        a[0] = 1.9715352749512141;
        a[1] = -1.4686795689225347;
        a[2] = 0.3893908434965701;
        b[0] = 0.013469181309343825;
        b[1] = 0.040407543928031475;
        break;
    case 4:  // fs : 16000
        a[0] = 1.4499664446880227;
        a[1] = -0.98943497080950582;
        a[2] = 0.24578252340690215;
        b[0] = 0.036710750339322612;
        b[1] = 0.11013225101796784;
        break;
    case 2:  // fs : 8000
        a[0] = 0.041156734567757189;
        a[1] = -0.42599112459189636;
        a[2] = 0.041037215479961225;
        b[0] = 0.16797464681802227;
        b[1] = 0.50392394045406674;
    }

    for (int i = 0; i < xLen; i++) {
        wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];

        y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];

        w[2] = w[1];
        w[1] = w[0];
        w[0] = wt;
    }
}

// F0軌跡の要素数を得る（事前にユーザがメモリ確保できるように）
// framePeriod の単位はmsec
int getSamplesForDIO(int fs, int xLen, double framePeriod) {
    return ( int )(( double )xLen / ( double )fs / (framePeriod / 1000.0)) + 1;
}

long decimateForF0(double* x, int xLen, double* y, int r) {
    //	int r = 11;
    int nfact = 9;  // ��������͌Œ��OK
    double *tmp1, *tmp2;
    tmp1 = ( double* )malloc(sizeof(double) * (xLen + nfact * 2));
    tmp2 = ( double* )malloc(sizeof(double) * (xLen + nfact * 2));

    int i;
    for (i = 0; i < nfact; i++)
        tmp1[i] = 2 * x[0] - x[nfact - i];
    for (i = nfact; i < nfact + xLen; i++)
        tmp1[i] = x[i - nfact];
    for (i = nfact + xLen; i < 2 * nfact + xLen; i++)
        tmp1[i] = 2 * x[xLen - 1] - x[xLen - 2 - (i - (nfact + xLen))];

    filterForDecimate(tmp1, 2 * nfact + xLen, tmp2, r);
    for (i = 0; i < 2 * nfact + xLen; i++)
        tmp1[i] = tmp2[2 * nfact + xLen - i - 1];
    filterForDecimate(tmp1, 2 * nfact + xLen, tmp2, r);
    for (i = 0; i < 2 * nfact + xLen; i++)
        tmp1[i] = tmp2[2 * nfact + xLen - i - 1];

    int nout = ( int )(xLen / r) + 1;
    int nbeg = r - (r * nout - xLen);
    int count;

    for (i = nbeg, count = 0; i < xLen + nfact; i += r, count++)
        y[count] = tmp1[i + nfact - 1];

    free(tmp1);
    free(tmp2);
    return count;
}

// DIO (Distributed Inline filter Operation) によるF0推定
// x	: 入力信号
// xLen : 信号長 [sample].
// f0	: 推定結果
void dio(double* x, int xLen, int fs, double framePeriod,
         double* timeAxis, double* f0) {
    int i, j;

    // 初期条件 (改良したい人はここから頑張って)
    double f0Floor = 80;
    double f0Ceil = 640;
    double channelsInOctave = 2;
    double targetFs = 4000;

    // 基礎パラメタの計算
    int decimationRatio = ( int )(fs / targetFs);
    double fss = ( double )fs / ( double )decimationRatio;
    int nBands = ( int )(log(( double )f0Ceil / ( double )f0Floor) / log(2.0) * channelsInOctave);

    // ここも基礎パラメタ
    double* boundaryF0List = ( double* )malloc(sizeof(double) * (nBands + 1));
    for (i = 0; i <= nBands; i++)
        boundaryF0List[i] = f0Floor * pow(2.0, i / channelsInOctave);

    // fft Lengthの計算
    int yLen = (1 + ( int )(xLen / decimationRatio));
    int fftl = ( int )pow(2.0, 1.0 + ( int )(log(( double )yLen +
                                                 ( double )(4 * ( int )(1.0 + ( double )fs / boundaryF0List[0] / 2.0))) /
                                             log(2.0)));
    double* y = ( double* )malloc(sizeof(double) * fftl);

    // ダウンサンプリング
    decimateForF0(x, xLen, y, decimationRatio);

    // 直流成分の除去 y = y - mean(y)
    double meanY = 0.0;
    for (i = 0; i < yLen; i++)
        meanY += y[i];
    meanY /= ( double )yLen;
    for (i = 0; i < yLen; i++)
        y[i] -= meanY;
    for (i = yLen; i < fftl; i++)
        y[i] = 0.0;

    // 中間データの保存用
    int tLen;                                        // F0軌跡のサンプル数
    tLen = getSamplesForDIO(fs, xLen, framePeriod);  // debug
    int lengthInMs = 1 + ( int )(( double )xLen / ( double )fs * 1000.0);
    double **stabilityMap, **f0Map;  // f0mapに候補が全て入るので，結果に納得できない場合は，f0Mapを直接操作する．
    stabilityMap = ( double** )malloc(sizeof(double*) * (nBands + 1));
    f0Map = ( double** )malloc(sizeof(double*) * (nBands + 1));
    for (i = 0; i <= nBands; i++) {
        stabilityMap[i] = ( double* )malloc(sizeof(double) * tLen);
        f0Map[i] = ( double* )malloc(sizeof(double) * tLen);
    }

    // 波形のスペクトルを事前に計算（ここは高速化の余地有り）
    fft_plan forwardFFT;  // FFTセット
    fft_complex* ySpec;   // スペクトル
    ySpec = ( fft_complex* )malloc(sizeof(fft_complex) * fftl);
    forwardFFT = fft_plan_dft_r2c_1d(fftl, y, ySpec, FFT_ESTIMATE);
    fft_execute(forwardFFT);  // FFTの実行

    // temporary values
    double* interpolatedF0;
    double* f0Deviations;
    interpolatedF0 = ( double* )malloc(sizeof(double) * lengthInMs);
    f0Deviations = ( double* )malloc(sizeof(double) * lengthInMs);

    for (i = 0; i < tLen; i++)
        timeAxis[i] = ( double )i * framePeriod / 1000.0;

    //tn_fnds v0.0.4 FFTWのプランを再利用するため、DIO側でメモリ確保、プランの作成を行う
    fft_destroy_plan(forwardFFT);
    double* equivalentFIR;
    equivalentFIR = ( double* )malloc(sizeof(double) * fftl);
    //	fft_plan	forwardFFT;				// FFTセット
    fft_plan inverseFFT;
    fft_complex* eSpec;  // スペクトル
    eSpec = ( fft_complex* )malloc(sizeof(fft_complex) * fftl);
    forwardFFT = fft_plan_dft_r2c_1d(fftl, equivalentFIR, eSpec, FFT_ESTIMATE);
    inverseFFT = fft_plan_dft_c2r_1d(fftl, eSpec, equivalentFIR, FFT_ESTIMATE);

    // イベントの計算 (4つのゼロ交差．詳しくは論文にて)
    for (i = 0; i <= nBands; i++) {
        //		rawEventByDio(boundaryF0List[i], fss, ySpec, yLen, fftl, framePeriod/1000.0, f0Floor, f0Ceil, timeAxis, tLen,
        //			f0Deviations, interpolatedF0);
        //tn_fnds v0.0.4 FFTWのプランを再利用する
        rawEventByDio(boundaryF0List[i], fss, ySpec, yLen, fftl, framePeriod / 1000.0, f0Floor, f0Ceil, timeAxis, tLen,
                      f0Deviations, interpolatedF0, forwardFFT, inverseFFT, equivalentFIR, eSpec);
        for (j = 0; j < tLen; j++) {
            stabilityMap[i][j] = f0Deviations[j] / (interpolatedF0[j] + 0.00000001);
            f0Map[i][j] = interpolatedF0[j];
        }
    }

    free(equivalentFIR);
    free(eSpec);

    // ベスト候補の選定 (基本波らしさを使い一意に決める)
    double* bestF0;
    //	bestF0 = (double *)malloc(sizeof(double) * (int)((double)xLen / (double)fs / (framePeriod/1000.0) ) + 1);
    bestF0 = ( double* )malloc(sizeof(double) * tLen);  // 2010/6/14 修正 (死にたい)
    // 後処理 (第一候補と候補マップから最適なパスを探す)
    //	postprocessing(framePeriod/1000.0, f0Floor, nBands+1, xLen, fs, f0Map, bestF0, f0);
    //tn_fnds v0.0.4 F0補正にstabilityMapを使用するようにした
    postprocessing(framePeriod / 1000.0, f0Floor, nBands + 1, xLen, fs, f0Map, stabilityMap, bestF0, f0);

    // お片づけ(メモリの開放)
    free(bestF0);
    free(interpolatedF0);
    free(f0Deviations);
    fft_destroy_plan(forwardFFT);
    fft_destroy_plan(inverseFFT);
    free(ySpec);
    for (i = 0; i <= nBands; i++) {
        free(stabilityMap[i]);
        free(f0Map[i]);
    }
    free(stabilityMap);
    free(f0Map);
    free(boundaryF0List);
    free(y);
}

// イベント数があったか判定
// longの範囲を超えてしまったので苦肉の策
int checkEvent(int x) {
    if (x > 0)
        return 1;
    return 0;
}

// 後処理（4ステップ）
//void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double **f0Map,
//					double *bestF0, double *f0)
void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double** f0Map,
                    double** stabilityMap, double* bestF0, double* f0) {
    int i, j, k;
    int voiceRangeMinimum = ( int )(0.5 + 1.0 / framePeriod / f0Floor);
    int f0Len = ( int )(( double )xLen / ( double )fs / framePeriod) + 1;
    //	double allowedRange = 0.1; // これは5 msecの基準なのでframePeriodに併せて調整する．
    double allowedRange = 0.1 * framePeriod / 0.005;  // これは5 msecの基準なのでframePeriodに併せて調整する．

    //tn_fnds v0.0.4 ベストなF0を抽出（DIO本体からこちらに移動）
    double tmp;
    double* bestF0Stab;
    bestF0Stab = ( double* )malloc(sizeof(double) * f0Len);
    for (i = 0; i < f0Len; i++) {
        tmp = stabilityMap[0][i];
        bestF0[i] = 0.0;
        bestF0Stab[i] = stabilityMap[0][i];
        for (j = 1; j < candidates; j++) {
            if (tmp > stabilityMap[j][i]) {
                tmp = stabilityMap[j][i];
                bestF0[i] = f0Map[j][i];
                bestF0Stab[i] = stabilityMap[j][i];
            }
        }
    }

    //tn_fnds v0.0.4 安定性の低いF0を除外
    int addCount = 0;
    double addValue = 0.0;
    for (i = 0; i < f0Len; i++) {
        if (bestF0Stab[i] < 0.05) {
            addCount++;
            addValue += bestF0Stab[i];
        }
    }
    addValue = addValue * 2.0 / addCount;
    for (i = 0; i < f0Len; i++)
        if (bestF0Stab[i] > addValue)
            bestF0[i] = 0.0;

    // メモリ節約はできるけど，どうせ少量なのでデバッグのしやすさを優先
    double* f0Base;
    f0Base = ( double* )malloc(sizeof(double) * f0Len);
    double* f0Step1;
    f0Step1 = ( double* )malloc(sizeof(double) * f0Len);
    double* f0Step2;
    f0Step2 = ( double* )malloc(sizeof(double) * f0Len);
    double* f0Step3;
    f0Step3 = ( double* )malloc(sizeof(double) * f0Len);
    double* f0Step4;
    f0Step4 = ( double* )malloc(sizeof(double) * f0Len);

    // まずは初期化
    for (i = 0; i < voiceRangeMinimum; i++)
        f0Base[i] = 0;
    for (; i < f0Len - voiceRangeMinimum; i++)
        f0Base[i] = bestF0[i];
    for (; i < f0Len; i++)
        f0Base[i] = 0;
    for (i = 0; i < f0Len; i++)
        f0Step1[i] = 0.0;

    // 第一のステップ (F0の跳躍防止)
    for (i = voiceRangeMinimum; i < f0Len; i++)
        if (fabs((f0Base[i] - f0Base[i - 1]) / (0.00001 + f0Base[i])) < allowedRange)
            f0Step1[i] = f0Base[i];

    // 第二のステップ (無声区間の切り離し)
    for (i = 0; i < f0Len; i++)
        f0Step2[i] = f0Step1[i];
    for (i = voiceRangeMinimum; i < f0Len; i++) {
        for (j = 0; j < voiceRangeMinimum; j++) {
            if (f0Step1[i - j] == 0) {
                f0Step2[i] = 0.0;
                break;
            }
        }
    }

    // tn_fnds v0.0.4  第二のステップ (無声区間の切り離し)逆方向
    for (i = f0Len - 1 - voiceRangeMinimum; i >= 0; i--) {
        for (j = 0; j < voiceRangeMinimum; j++) {
            if (f0Step1[i + j] == 0) {
                f0Step2[i] = 0.0;
                break;
            }
        }
    }

    // 島数の検出
    int *positiveIndex, *negativeIndex;
    positiveIndex = ( int* )malloc(sizeof(int) * f0Len);
    negativeIndex = ( int* )malloc(sizeof(int) * f0Len);
    int positiveCount, negativeCount;
    positiveCount = negativeCount = 0;
    for (i = 1; i < f0Len; i++) {
        if (f0Step2[i] == 0 && f0Step2[i - 1] != 0)
            negativeIndex[negativeCount++] = i - 1;
        else if (f0Step2[i - 1] == 0 && f0Step2[i] != 0)
            positiveIndex[positiveCount++] = i;
    }

    // ステップ3（前向き補正）
    double refValue1, refValue2, bestError, errorValue;
    for (i = 0; i < f0Len; i++)
        f0Step3[i] = f0Step2[i];
    for (i = 0; i < negativeCount; i++) {
        for (j = negativeIndex[i]; j < f0Len - 1; j++) {
            if (f0Step3[j + 1] != 0)
                break;
            refValue1 = f0Step3[j] * 2 - f0Step3[j - 1];
            refValue2 = f0Step3[j];
            //			bestError = fabs(refValue - f0Map[0][j+1]);
            bestError = min(fabs(refValue1 - f0Map[0][j + 1]), fabs(refValue2 - f0Map[0][j + 1]));
            for (k = 1; k < candidates; k++) {
                //				errorValue = fabs(refValue - f0Map[k][j+1]);
                errorValue = min(fabs(refValue1 - f0Map[k][j + 1]), fabs(refValue2 - f0Map[k][j + 1]));
                //				if(errorValue < bestError)
                if (errorValue < bestError && stabilityMap[k][j + 1] < 0.1)  //tn_fnds v0.0.4 安定性の低いF0は使用しない
                {
                    bestError = errorValue;
                    f0Step3[j + 1] = f0Map[k][j + 1];
                }
            }
            //			if(bestError / (refValue+0.0001) > allowedRange)
            if (min(bestError / (refValue1 + 0.0001), bestError / (refValue2 + 0.0001)) > allowedRange) {
                f0Step3[j + 1] = 0.0;
                break;
            }
            if (i != negativeCount && j == positiveIndex[i + 1] - 1) {
                negativeIndex[j] = j;
                break;
            }
        }
    }

    // ステップ4（後向き補正）
    for (i = 0; i < f0Len; i++)
        f0Step4[i] = f0Step3[i];
    for (i = positiveCount - 1; i >= 0; i--) {
        for (j = positiveIndex[i] /*+1*/; j > 1; j--)  //tn_fnds v0.0.4
        {
            if (f0Step4[j - 1] != 0)
                break;
            refValue1 = f0Step4[j] * 2 - f0Step4[j - 1];
            refValue2 = f0Step4[j];
            //			refValue = f0Step4[j]*2 - f0Step4[j+1];
            bestError = min(fabs(refValue1 - f0Map[0][j + 1]), fabs(refValue2 - f0Map[0][j + 1]));
            //			bestError = fabs(refValue - f0Map[0][j-1]);
            for (k = 1; k < candidates; k++) {
                errorValue = min(fabs(refValue1 - f0Map[k][j - 1]), fabs(refValue2 - f0Map[k][j - 1]));
                //				errorValue = fabs(refValue - f0Map[k][j-1]);
                //				if(min(bestError / (refValue1+0.0001), bestError / (refValue2+0.0001)) > allowedRange)
                if (min(bestError / (refValue1 + 0.0001), bestError / (refValue2 + 0.0001)) > allowedRange && stabilityMap[k][j - 1] < 0.1)  //tn_fnds v0.0.4 安定性の低いF0は使用しない
                                                                                                                                             //				if(errorValue < bestError)
                {
                    bestError = errorValue;
                    f0Step4[j - 1] = f0Map[k][j - 1];
                }
            }
            if (min(bestError / (refValue1 + 0.0001), bestError / (refValue2 + 0.0001)) > allowedRange)
            //			if(bestError / (refValue+0.0001) > allowedRange)
            {
                f0Step4[j - 1] = 0.0;
                break;
            }
            if (i != 0 && j == negativeIndex[i - 1] + 1)
                break;
        }
    }

    // コピー
    for (i = 0; i < f0Len; i++)
        f0[i] = f0Step4[i];

    // メモリの開放
    free(bestF0Stab);
    free(f0Base);
    free(f0Step1);
    free(f0Step2);
    free(f0Step3);
    free(f0Step4);
}

// イベントを計算する内部関数 (内部変数なので引数・戻り値に手加減なし)
void rawEventByDio(double boundaryF0, double fs, fft_complex* xSpec, int xLength, int fftl, double framePeriod, double f0Floor, double f0Ceil, double* timeAxis, int tLen,
                   double* f0Deviations, double* interpolatedF0,
                   fft_plan& forwardFFT, fft_plan& inverseFFT, double* equivalentFIR, fft_complex* eSpec) {
    int i;
    int halfAverageLength = ( int )(fs / boundaryF0 / 2 + 0.5);
    int indexBias = halfAverageLength * 2;
    //	double *equivalentFIR;
    //	equivalentFIR = (double *)malloc(sizeof(double) * fftl);
    for (i = halfAverageLength * 2; i < fftl; i++)
        equivalentFIR[i] = 0.0;
    nuttallWindow(halfAverageLength * 4, equivalentFIR);

    //tn_fnds v0.0.4 FFTWプランはDIO本体で作成することにした
    //	fft_plan			forwardFFT;				// FFTセット
    //	fft_complex		*eSpec;	// スペクトル
    //	eSpec = (fft_complex *)malloc(sizeof(fft_complex) * fftl);
    //	forwardFFT = fft_plan_dft_r2c_1d(fftl, equivalentFIR, eSpec, FFT_ESTIMATE);
    fft_execute(forwardFFT);  // FFTの実行

    // 複素数の掛け算
    double tmp;
    //	for(i = 0;i <= fftl-1;i++)　
    for (i = 0; i <= fftl >> 1; i++)  //tn_fnds v0.0.4 FFT長の半分でいいはず
    {
        tmp = xSpec[i][0] * eSpec[i][0] - xSpec[i][1] * eSpec[i][1];
        eSpec[i][1] = xSpec[i][0] * eSpec[i][1] + xSpec[i][1] * eSpec[i][0];
        eSpec[i][0] = tmp;
    }

    // 低域通過フィルタリング
    //tn_fnds v0.0.4 FFTWプランはDIO本体で作成することにした
    //	fft_plan	 inverseFFT;
    //	inverseFFT = fft_plan_dft_c2r_1d(fftl, eSpec, equivalentFIR, FFT_ESTIMATE);
    fft_execute(inverseFFT);
    // バイアス（低域通過フィルタによる遅延）の除去
    for (i = 0; i < xLength; i++)
        equivalentFIR[i] = equivalentFIR[i + indexBias];

    // ４つのゼロ交差(構造体のほうがいいね) e:event, i:interval
    double *nELocations, *pELocations, *dnELocations, *dpELocations;
    double *nILocations, *pILocations, *dnILocations, *dpILocations;
    double *nIntervals, *pIntervals, *dnIntervals, *dpIntervals;
    int nLen, pLen, dnLen, dpLen;
    nELocations = ( double* )malloc(sizeof(double) * xLength);  // xLengthはかなりの保険
    pELocations = ( double* )malloc(sizeof(double) * xLength);
    dnELocations = ( double* )malloc(sizeof(double) * xLength);
    dpELocations = ( double* )malloc(sizeof(double) * xLength);
    nILocations = ( double* )malloc(sizeof(double) * xLength);
    pILocations = ( double* )malloc(sizeof(double) * xLength);
    dnILocations = ( double* )malloc(sizeof(double) * xLength);
    dpILocations = ( double* )malloc(sizeof(double) * xLength);
    nIntervals = ( double* )malloc(sizeof(double) * xLength);
    pIntervals = ( double* )malloc(sizeof(double) * xLength);
    dnIntervals = ( double* )malloc(sizeof(double) * xLength);
    dpIntervals = ( double* )malloc(sizeof(double) * xLength);

    zeroCrossingEngine(equivalentFIR, xLength, fs,
                       nELocations, nILocations, nIntervals, &nLen);

    for (i = 0; i < xLength; i++)
        equivalentFIR[i] = -equivalentFIR[i];
    zeroCrossingEngine(equivalentFIR, xLength, fs,
                       pELocations, pILocations, pIntervals, &pLen);

    for (i = 0; i < xLength - 1; i++)
        equivalentFIR[i] = equivalentFIR[i] - equivalentFIR[i + 1];
    zeroCrossingEngine(equivalentFIR, xLength - 1, fs,
                       dnELocations, dnILocations, dnIntervals, &dnLen);

    for (i = 0; i < xLength - 1; i++)
        equivalentFIR[i] = -equivalentFIR[i];
    zeroCrossingEngine(equivalentFIR, xLength - 1, fs,
                       dpELocations, dpILocations, dpIntervals, &dpLen);

    int usableChannel;
    usableChannel = checkEvent(nLen - 2) * checkEvent(pLen - 2) *
                    checkEvent(dnLen - 2) * checkEvent(dpLen - 2);

    double* interpolatedF0Set[4];
    if (usableChannel <= 0) {  // ノー候補でフィニッシュです
        for (i = 0; i < tLen; i++) {
            f0Deviations[i] = 100000.0;
            interpolatedF0[i] = 0.0;
        }
    } else {
        for (i = 0; i < 4; i++)
            interpolatedF0Set[i] = ( double* )malloc(sizeof(double) * tLen);
        // 4つのゼロ交差
        interp1(nILocations, nIntervals, nLen, timeAxis, tLen, interpolatedF0Set[0]);
        interp1(pILocations, pIntervals, pLen, timeAxis, tLen, interpolatedF0Set[1]);
        interp1(dnILocations, dnIntervals, dnLen, timeAxis, tLen, interpolatedF0Set[2]);
        interp1(dpILocations, dpIntervals, dpLen, timeAxis, tLen, interpolatedF0Set[3]);

        for (i = 0; i < tLen; i++) {
            interpolatedF0[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] +
                                 interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) /
                                4.0;

            f0Deviations[i] = sqrt(((interpolatedF0Set[0][i] - interpolatedF0[i]) * (interpolatedF0Set[0][i] - interpolatedF0[i]) + (interpolatedF0Set[1][i] - interpolatedF0[i]) * (interpolatedF0Set[1][i] - interpolatedF0[i]) + (interpolatedF0Set[2][i] - interpolatedF0[i]) * (interpolatedF0Set[2][i] - interpolatedF0[i]) + (interpolatedF0Set[3][i] - interpolatedF0[i]) * (interpolatedF0Set[3][i] - interpolatedF0[i])) / 3.0);

            if (interpolatedF0[i] > boundaryF0 || interpolatedF0[i] < boundaryF0 / 2.0 || interpolatedF0[i] > f0Ceil || interpolatedF0[i] < FLOOR_F0)  // 70 Hz以下はF0としない．
            {
                interpolatedF0[i] = 0.0;
                f0Deviations[i] = 100000.0;
            }
        }

        for (i = 0; i < 4; i++)
            free(interpolatedF0Set[i]);
    }

    // メモリの開放
    free(nELocations);
    free(pELocations);
    free(dnELocations);
    free(dpELocations);
    free(nILocations);
    free(pILocations);
    free(dnILocations);
    free(dpILocations);
    free(nIntervals);
    free(pIntervals);
    free(dnIntervals);
    free(dpIntervals);
    //	fft_destroy_plan(inverseFFT);
    //	fft_destroy_plan(forwardFFT);
    //	free(eSpec);
    //	free(equivalentFIR);
}

// ゼロ交差を計算
void zeroCrossingEngine(double* x, int xLen, double fs,
                        double* eLocations, double* iLocations, double* intervals, int* iLen) {
    int i;
    int* negativeGoingPoints;
    negativeGoingPoints = ( int* )malloc(sizeof(int) * xLen);

    int tmp1, tmp2;
    for (i = 0; i < xLen - 1; i++)  // 毎回余りを計算するのは無駄
    {
        tmp1 = x[i] * x[i + 1] < 0 ? 1 : 0;
        tmp2 = x[i + 1] < x[i] ? 1 : 0;
        negativeGoingPoints[i] = (i + 1) * tmp1 * tmp2;
    }
    negativeGoingPoints[xLen - 1] = 0;

    // 有効イベントの検出
    int* edges;
    edges = ( int* )malloc(sizeof(int) * xLen);
    int count = 0;
    for (i = 0; i < xLen; i++) {
        if (negativeGoingPoints[i] > 0)
            edges[count++] = negativeGoingPoints[i];
    }
    // 最終戻り値の計算準備
    double* fineEdges;
    fineEdges = ( double* )malloc(sizeof(double) * count);
    for (i = 0; i < count; i++) {
        fineEdges[i] = ( double )edges[i] - x[edges[i] - 1] / (x[edges[i]] - x[edges[i] - 1]);
    }

    *iLen = count - 1;
    for (i = 0; i < *iLen; i++) {
        intervals[i] = fs / (fineEdges[i + 1] - fineEdges[i]);
        iLocations[i] = (fineEdges[i] + fineEdges[i + 1]) / 2.0 / fs;
        eLocations[i] = fineEdges[i] / fs;
    }
    if (count != 0)
        eLocations[count - 1] = fineEdges[count - 1] / fs;  //0の場合を考慮

    free(fineEdges);
    free(edges);
    free(negativeGoingPoints);
}

// ナットール窓．マジックナンバーのように見えるけどこれが正解．
void nuttallWindow(int yLen, double* y) {
    int i;
    double tmp;
    for (i = 0; i < yLen; i++) {
        tmp = (( double )(i + 1) - ( double )(yLen + 1) / 2.0) / ( double )(yLen + 1);
        y[i] = 0.355768 + 0.487396 * cos(2 * PI * tmp) + 0.144232 * cos(4 * PI * tmp) + 0.012604 * cos(6 * PI * tmp);
    }
}
