#include "main.h"

// �����֐�(���[�U�͐G��Ȃ��ق����ǂ�)
//void rawEventByDio(double boundaryF0, double fs, fft_complex *xSpec, int xLength, int fftl, double shiftTime, double f0Floor, double f0Ceil, double *timeAxis, int tLen, 
//				   double *f0Deviations, double *interpolatedF0);
//tn_fnds v0.0.4
void rawEventByDio(double boundaryF0, double fs, fft_complex *xSpec, int xLength, int fftl, double shiftTime, double f0Floor, double f0Ceil, double *timeAxis, int tLen, double *f0Deviations, double *interpolatedF0, fft_plan &forwardFFT, fft_plan &inverseFFT, double *equivalentFIR, fft_complex *eSpec);

void zeroCrossingEngine(double *x, int xLen, double fs, double *eLocations, double *iLocations, double *intervals, int *iLen);

void nuttallWindow(int yLen, double *y);

//void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs,
//					double **f0Map, double *bestF0, double *f0);
//tn_fnds v0.0.4
void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double **f0Map, double **stabilityMap, double *bestF0, double *f0);

void histc(double *x, int xLen, double *y, int yLen, int *index);

// F0�O�Ղ̗v�f���𓾂�i���O�Ƀ��[�U���������m�ۂł���悤�Ɂj
// framePeriod �̒P�ʂ�msec
int getSamplesForDIO(int fs, int xLen, double framePeriod) {
    return (int) ((double) xLen / (double) fs / (framePeriod / 1000.0)) + 1;
}

// DIO (Distributed Inline filter Operation) �ɂ��F0����
// x	: ���͐M��
// xLen : �M���� [sample].
// f0	: ���茋��
void dio(double *x, int xLen, int fs, double framePeriod,
         double *timeAxis, double *f0) {
    int i, j;

    // �������� (���ǂ������l�͂�������撣����)
    double f0Floor = 80;
    double f0Ceil = 640;
    double channelsInOctave = 2;
    double targetFs = 4000;

    // ��b�p�����^�̌v�Z
    int decimationRatio = (int) (fs / targetFs);
    double fss = (double) fs / (double) decimationRatio;
    int nBands = (int) (log((double) f0Ceil / (double) f0Floor) / log(2.0) * channelsInOctave);

    // ��������b�p�����^
    double *boundaryF0List = (double *) malloc(sizeof(double) * (nBands + 1));
    for (i = 0; i <= nBands; i++)
        boundaryF0List[i] = f0Floor * pow(2.0, i / channelsInOctave);

    // fft Length�̌v�Z
    int yLen = (1 + (int) (xLen / decimationRatio));
    int fftl = (int) pow(2.0, 1.0 + (int) (log((double) yLen +
                                               (double) (4 * (int) (1.0 + (double) fs / boundaryF0List[0] / 2.0))) / log(2.0)));
    auto *y = (double *) malloc(sizeof(double) * fftl);

    // ���������̏��� y = y - mean(y)
    double meanY = 0.0;
    for (i = 0; i < yLen; i++) meanY += y[i];
    meanY /= (double) yLen;
    for (i = 0; i < yLen; i++) y[i] -= meanY;
    for (i = yLen; i < fftl; i++) y[i] = 0.0;

    // ���ԃf�[�^�̕ۑ��p
    int tLen; // F0�O�Ղ̃T���v����
    tLen = getSamplesForDIO(fs, xLen, framePeriod); // debug
    int lengthInMs = 1 + (int) ((double) xLen / (double) fs * 1000.0);
    double **stabilityMap, **f0Map; // f0map�Ɍ�₪�S�ē���̂ŁC���ʂɔ[���ł��Ȃ��ꍇ�́Cf0Map�𒼐ڑ��삷��D
    stabilityMap = (double **) malloc(sizeof(double *) * (nBands + 1));
    f0Map = (double **) malloc(sizeof(double *) * (nBands + 1));
    for (i = 0; i <= nBands; i++) {
        stabilityMap[i] = (double *) malloc(sizeof(double) * tLen);
        f0Map[i] = (double *) malloc(sizeof(double) * tLen);
    }

    // �g�`�̃X�y�N�g�������O�Ɍv�Z�i�����͍������̗]�n�L��j
    fft_plan forwardFFT;                // FFT�Z�b�g
    fft_complex *ySpec;    // �X�y�N�g��
    ySpec = (fft_complex *) malloc(sizeof(fft_complex) * fftl);
    forwardFFT = fft_plan_dft_r2c_1d(fftl, y, ySpec, FFT_ESTIMATE);
    fft_execute(forwardFFT); // FFT�̎��s

    // temporary values
    double *interpolatedF0;
    double *f0Deviations;
    interpolatedF0 = (double *) malloc(sizeof(double) * lengthInMs);
    f0Deviations = (double *) malloc(sizeof(double) * lengthInMs);

    for (i = 0; i < tLen; i++)
        timeAxis[i] = (double) i * framePeriod / 1000.0;

    //tn_fnds v0.0.4 FFTW�̃v�������ė��p���邽�߁ADIO���Ń������m�ہA�v�����̍쐬���s��
    fft_destroy_plan(forwardFFT);
    double *equivalentFIR;
    equivalentFIR = (double *) malloc(sizeof(double) * fftl);
//	fft_plan	forwardFFT;				// FFT�Z�b�g
    fft_plan inverseFFT;
    fft_complex *eSpec;    // �X�y�N�g��
    eSpec = (fft_complex *) malloc(sizeof(fft_complex) * fftl);
    forwardFFT = fft_plan_dft_r2c_1d(fftl, equivalentFIR, eSpec, FFT_ESTIMATE);
    inverseFFT = fft_plan_dft_c2r_1d(fftl, eSpec, equivalentFIR, FFT_ESTIMATE);

    // �C�x���g�̌v�Z (4�̃[�������D�ڂ����͘_���ɂ�)
    for (i = 0; i <= nBands; i++) {
        rawEventByDio(boundaryF0List[i], fss, ySpec, yLen, fftl, framePeriod / 1000.0, f0Floor, f0Ceil, timeAxis, tLen,f0Deviations, interpolatedF0, forwardFFT, inverseFFT, equivalentFIR, eSpec);
        for (j = 0; j < tLen; j++) {
            stabilityMap[i][j] = f0Deviations[j] / (interpolatedF0[j] + 0.00000001);
            f0Map[i][j] = interpolatedF0[j];
        }
    }

    free(equivalentFIR);
    free(eSpec);

    // �x�X�g���̑I�� (��{�g�炵�����g����ӂɌ��߂�)
    double *bestF0;
//	bestF0 = (double *)malloc(sizeof(double) * (int)((double)xLen / (double)fs / (framePeriod/1000.0) ) + 1);
    bestF0 = (double *) malloc(sizeof(double) * tLen); // 2010/6/14 �C�� (���ɂ���)
/* tn_fnds v0.0.4 �x�X�g���̑I��͑S��postprocessing�ōs���悤�ɂ���
	double tmp;
	for(i = 0;i < tLen;i++)
	{
		tmp = stabilityMap[0][i];
		bestF0[i] = (stabilityMap[0][i] < 0.002) ? f0Map[0][i] : 0.0;
		for(j = 1;j <= nBands;j++)
		{
			if(tmp > stabilityMap[j][i] && stabilityMap[j][i] < 0.002)
			{
				tmp = stabilityMap[j][i];
				bestF0[i] = f0Map[j][i];
			}
		}
	}
*/
    // �㏈�� (�����ƌ��}�b�v����œK�ȃp�X��T��)
//	postprocessing(framePeriod/1000.0, f0Floor, nBands+1, xLen, fs, f0Map, bestF0, f0);
    //tn_fnds v0.0.4 F0�␳��stabilityMap���g�p����悤�ɂ���
    postprocessing(framePeriod / 1000.0, f0Floor, nBands + 1, xLen, fs, f0Map, stabilityMap, bestF0, f0);

    // ���ЂÂ�(�������̊J��)
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

// �C�x���g����������������
// long�͈̔͂𒴂��Ă��܂����̂ŋ���̍�
int checkEvent(int x) {
    if (x > 0) return 1;
    return 0;
}

// �㏈���i4�X�e�b�v�j
//void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double **f0Map,
//					double *bestF0, double *f0)
void postprocessing(double framePeriod, double f0Floor, int candidates, int xLen, int fs, double **f0Map,
                    double **stabilityMap, double *bestF0, double *f0) {
    int i, j, k;
    int voiceRangeMinimum = (int) (0.5 + 1.0 / framePeriod / f0Floor);
    int f0Len = (int) ((double) xLen / (double) fs / framePeriod) + 1;
//	double allowedRange = 0.1; // �����5 msec�̊�Ȃ̂�framePeriod�ɕ����Ē�������D
    double allowedRange = 0.1 * framePeriod / 0.005; // �����5 msec�̊�Ȃ̂�framePeriod�ɕ����Ē�������D

    //tn_fnds v0.0.4 �x�X�g��F0�𒊏o�iDIO�{�̂��炱����Ɉړ��j
    double tmp;
    double *bestF0Stab;
    bestF0Stab = (double *) malloc(sizeof(double) * f0Len);
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

    //tn_fnds v0.0.4 ���萫�̒ႢF0�����O
    int addCount = 0;
    double addValue = 0.0;
    for (i = 0; i < f0Len; i++) {
        if (bestF0Stab[i] < 0.05) {
            addCount++;
            addValue += bestF0Stab[i];
        }
    }
    addValue = addValue * 2.0 / addCount;
    for (i = 0; i < f0Len; i++) if (bestF0Stab[i] > addValue) bestF0[i] = 0.0;

    // �������ߖ�͂ł��邯�ǁC�ǂ������ʂȂ̂Ńf�o�b�O�̂��₷����D��
    double *f0Base;
    f0Base = (double *) malloc(sizeof(double) * f0Len);
    double *f0Step1;
    f0Step1 = (double *) malloc(sizeof(double) * f0Len);
    double *f0Step2;
    f0Step2 = (double *) malloc(sizeof(double) * f0Len);
    double *f0Step3;
    f0Step3 = (double *) malloc(sizeof(double) * f0Len);
    double *f0Step4;
    f0Step4 = (double *) malloc(sizeof(double) * f0Len);

    // �܂��͏�����
    for (i = 0; i < voiceRangeMinimum; i++) f0Base[i] = 0;
    for (; i < f0Len - voiceRangeMinimum; i++) f0Base[i] = bestF0[i];
    for (; i < f0Len; i++) f0Base[i] = 0;
    for (i = 0; i < f0Len; i++) f0Step1[i] = 0.0;

    // ���̃X�e�b�v (F0�̒����h�~)
    for (i = voiceRangeMinimum; i < f0Len; i++)
        if (fabs((f0Base[i] - f0Base[i - 1]) / (0.00001 + f0Base[i])) < allowedRange)
            f0Step1[i] = f0Base[i];

    // ���̃X�e�b�v (������Ԃ̐؂藣��)
    for (i = 0; i < f0Len; i++) f0Step2[i] = f0Step1[i];
    for (i = voiceRangeMinimum; i < f0Len; i++) {
        for (j = 0; j < voiceRangeMinimum; j++) {
            if (f0Step1[i - j] == 0) {
                f0Step2[i] = 0.0;
                break;
            }
        }
    }

    // tn_fnds v0.0.4  ���̃X�e�b�v (������Ԃ̐؂藣��)�t����
    for (i = f0Len - 1 - voiceRangeMinimum; i >= 0; i--) {
        for (j = 0; j < voiceRangeMinimum; j++) {
            if (f0Step1[i + j] == 0) {
                f0Step2[i] = 0.0;
                break;
            }
        }
    }

    // �����̌��o
    int *positiveIndex, *negativeIndex;
    positiveIndex = (int *) malloc(sizeof(int) * f0Len);
    negativeIndex = (int *) malloc(sizeof(int) * f0Len);
    int positiveCount, negativeCount;
    positiveCount = negativeCount = 0;
    for (i = 1; i < f0Len; i++) {
        if (f0Step2[i] == 0 && f0Step2[i - 1] != 0)
            negativeIndex[negativeCount++] = i - 1;
        else if (f0Step2[i - 1] == 0 && f0Step2[i] != 0)
            positiveIndex[positiveCount++] = i;
    }

    // �X�e�b�v3�i�O�����␳�j
    double refValue1, refValue2, bestError, errorValue;
    for (i = 0; i < f0Len; i++) f0Step3[i] = f0Step2[i];
    for (i = 0; i < negativeCount; i++) {
        for (j = negativeIndex[i]; j < f0Len - 1; j++) {
            if (f0Step3[j + 1] != 0) break;
            refValue1 = f0Step3[j] * 2 - f0Step3[j - 1];
            refValue2 = f0Step3[j];
//			bestError = fabs(refValue - f0Map[0][j+1]);
            bestError = min(fabs(refValue1 - f0Map[0][j + 1]), fabs(refValue2 - f0Map[0][j + 1]));
            for (k = 1; k < candidates; k++) {
//				errorValue = fabs(refValue - f0Map[k][j+1]);
                errorValue = min(fabs(refValue1 - f0Map[k][j + 1]), fabs(refValue2 - f0Map[k][j + 1]));
//				if(errorValue < bestError)
                if (errorValue < bestError && stabilityMap[k][j + 1] < 0.1) //tn_fnds v0.0.4 ���萫�̒ႢF0�͎g�p���Ȃ�
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

    // �X�e�b�v4�i������␳�j
    for (i = 0; i < f0Len; i++) f0Step4[i] = f0Step3[i];
    for (i = positiveCount - 1; i >= 0; i--) {
        for (j = positiveIndex[i]/*+1*/; j > 1; j--) //tn_fnds v0.0.4
        {
            if (f0Step4[j - 1] != 0) break;
            refValue1 = f0Step4[j] * 2 - f0Step4[j - 1];
            refValue2 = f0Step4[j];
//			refValue = f0Step4[j]*2 - f0Step4[j+1];
            bestError = min(fabs(refValue1 - f0Map[0][j + 1]), fabs(refValue2 - f0Map[0][j + 1]));
//			bestError = fabs(refValue - f0Map[0][j-1]);
            for (k = 1; k < candidates; k++) {
                errorValue = min(fabs(refValue1 - f0Map[k][j - 1]), fabs(refValue2 - f0Map[k][j - 1]));
//				errorValue = fabs(refValue - f0Map[k][j-1]);
//				if(min(bestError / (refValue1+0.0001), bestError / (refValue2+0.0001)) > allowedRange)
                if (min(bestError / (refValue1 + 0.0001), bestError / (refValue2 + 0.0001)) > allowedRange && stabilityMap[k][j - 1] < 0.1) //tn_fnds v0.0.4 ���萫�̒ႢF0�͎g�p���Ȃ�
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
            if (i != 0 && j == negativeIndex[i - 1] + 1) break;
        }
    }

    // �R�s�[
    for (i = 0; i < f0Len; i++) f0[i] = f0Step4[i];
/* �X�e�b�v5�́C���\���オ��Ȃ��̂ňꎞ�I�ɍ폜
	// �X�e�b�v5�i�Ǘ����̐؂藣�� 2��ځj
	int voiceRangeMinimum2 = 2+(int)(voiceRangeMinimum/2);
	for(i = 0;i < f0Len;i++) f0[i] = f0Step4[i];
	for(i = voiceRangeMinimum2; i < f0Len-voiceRangeMinimum2;i++)
	{
		for(j = 0;j < voiceRangeMinimum2;j++)
		{
			if(f0Step4[i-j] == 0)
				break;
		}
		for(k = 0;k < voiceRangeMinimum2;k++)
		{
			if(f0Step4[i+k] == 0)
				break;
		}
		f0[i] = j != voiceRangeMinimum2 && k != voiceRangeMinimum2 ? 
			0 : f0Step4[i];
	}
*/

//���_�ƕ␳���ʂ̃t�@�C���o�́@tn_fnds�f�o�b�O�p
/*
	FILE *file;

	file = fopen("/tmp/f0map.txt", "w");

	for(k = 1;k < candidates;k++) fprintf(file,",map%d",k);
	for(k = 1;k < candidates;k++) fprintf(file,",stab%d",k);
	fprintf(file,",BEST,base,step1,step2,step3,step4\n");
	for(i = 0; i < f0Len; i++)
	{
		fprintf(file,"%d",i);
		for(k = 1;k < candidates;k++) fprintf(file,",%f",f0Map[k][i]);
		for(k = 1;k < candidates;k++) fprintf(file,",%f",stabilityMap[k][i]);
		fprintf(file,",%f", bestF0[i]);
		fprintf(file,",%f", f0Base[i]);
		fprintf(file,",%f", f0Step1[i]);
		fprintf(file,",%f", f0Step2[i]);
		fprintf(file,",%f", f0Step3[i]);
		fprintf(file,",%f\n",f0Step4[i]);
	}
	fclose(file);
*/
    // �������̊J��
    free(bestF0Stab);
    free(f0Base);
    free(f0Step1);
    free(f0Step2);
    free(f0Step3);
    free(f0Step4);
}

// �C�x���g���v�Z��������֐� (�����ϐ��Ȃ̂ň����E�߂�l�Ɏ�����Ȃ�)
void rawEventByDio(double boundaryF0, double fs, fft_complex *xSpec, int xLength, int fftl, double framePeriod, double f0Floor, double f0Ceil, double *timeAxis, int tLen,
                   double *f0Deviations, double *interpolatedF0,
                   fft_plan &forwardFFT, fft_plan &inverseFFT, double *equivalentFIR, fft_complex *eSpec) {
    int i;
    int halfAverageLength = (int) (fs / boundaryF0 / 2 + 0.5);
    int indexBias = halfAverageLength * 2;
//	double *equivalentFIR;
//	equivalentFIR = (double *)malloc(sizeof(double) * fftl);
    for (i = halfAverageLength * 2; i < fftl; i++) equivalentFIR[i] = 0.0;
    nuttallWindow(halfAverageLength * 4, equivalentFIR);

//tn_fnds v0.0.4 FFTW�v������DIO�{�̂ō쐬���邱�Ƃɂ���
//	fft_plan			forwardFFT;				// FFT�Z�b�g
//	fft_complex		*eSpec;	// �X�y�N�g��
//	eSpec = (fft_complex *)malloc(sizeof(fft_complex) * fftl);
//	forwardFFT = fft_plan_dft_r2c_1d(fftl, equivalentFIR, eSpec, FFT_ESTIMATE);
    fft_execute(forwardFFT); // FFT�̎��s

    // ���f���̊|���Z
    double tmp;
//	for(i = 0;i <= fftl-1;i++)�@
    for (i = 0; i <= fftl >> 1; i++)    //tn_fnds v0.0.4 FFT���̔����ł����͂�
    {
        tmp = xSpec[i][0] * eSpec[i][0] - xSpec[i][1] * eSpec[i][1];
        eSpec[i][1] = xSpec[i][0] * eSpec[i][1] + xSpec[i][1] * eSpec[i][0];
        eSpec[i][0] = tmp;
    }

    // ���ʉ߃t�B���^�����O
//tn_fnds v0.0.4 FFTW�v������DIO�{�̂ō쐬���邱�Ƃɂ���
//	fft_plan	 inverseFFT;
//	inverseFFT = fft_plan_dft_c2r_1d(fftl, eSpec, equivalentFIR, FFT_ESTIMATE);
    fft_execute(inverseFFT);
    // �o�C�A�X�i���ʉ߃t�B���^�ɂ��x���j�̏���
    for (i = 0; i < xLength; i++) equivalentFIR[i] = equivalentFIR[i + indexBias];

    // �S�̃[������(�\���̂̂ق���������) e:event, i:interval
    double *nELocations, *pELocations, *dnELocations, *dpELocations;
    double *nILocations, *pILocations, *dnILocations, *dpILocations;
    double *nIntervals, *pIntervals, *dnIntervals, *dpIntervals;
    int nLen, pLen, dnLen, dpLen;
    nELocations = (double *) malloc(sizeof(double) * xLength); // xLength�͂��Ȃ�̕ی�
    pELocations = (double *) malloc(sizeof(double) * xLength);
    dnELocations = (double *) malloc(sizeof(double) * xLength);
    dpELocations = (double *) malloc(sizeof(double) * xLength);
    nILocations = (double *) malloc(sizeof(double) * xLength);
    pILocations = (double *) malloc(sizeof(double) * xLength);
    dnILocations = (double *) malloc(sizeof(double) * xLength);
    dpILocations = (double *) malloc(sizeof(double) * xLength);
    nIntervals = (double *) malloc(sizeof(double) * xLength);
    pIntervals = (double *) malloc(sizeof(double) * xLength);
    dnIntervals = (double *) malloc(sizeof(double) * xLength);
    dpIntervals = (double *) malloc(sizeof(double) * xLength);

    zeroCrossingEngine(equivalentFIR, xLength, fs,
                       nELocations, nILocations, nIntervals, &nLen);

    for (i = 0; i < xLength; i++) equivalentFIR[i] = -equivalentFIR[i];
    zeroCrossingEngine(equivalentFIR, xLength, fs,
                       pELocations, pILocations, pIntervals, &pLen);

    for (i = 0; i < xLength - 1; i++) equivalentFIR[i] = equivalentFIR[i] - equivalentFIR[i + 1];
    zeroCrossingEngine(equivalentFIR, xLength - 1, fs,
                       dnELocations, dnILocations, dnIntervals, &dnLen);

    for (i = 0; i < xLength - 1; i++) equivalentFIR[i] = -equivalentFIR[i];
    zeroCrossingEngine(equivalentFIR, xLength - 1, fs,
                       dpELocations, dpILocations, dpIntervals, &dpLen);


    int usableChannel;
    usableChannel = checkEvent(nLen - 2) * checkEvent(pLen - 2) *
                    checkEvent(dnLen - 2) * checkEvent(dpLen - 2);

    double *interpolatedF0Set[4];
    if (usableChannel <= 0) { // �m�[���Ńt�B�j�b�V���ł�
        for (i = 0; i < tLen; i++) {
            f0Deviations[i] = 100000.0;
            interpolatedF0[i] = 0.0;
        }
    } else {
        for (i = 0; i < 4; i++)
            interpolatedF0Set[i] = (double *) malloc(sizeof(double) * tLen);
        // 4�̃[������
        interp1(nILocations, nIntervals, nLen, timeAxis, tLen, interpolatedF0Set[0]);
        interp1(pILocations, pIntervals, pLen, timeAxis, tLen, interpolatedF0Set[1]);
        interp1(dnILocations, dnIntervals, dnLen, timeAxis, tLen, interpolatedF0Set[2]);
        interp1(dpILocations, dpIntervals, dpLen, timeAxis, tLen, interpolatedF0Set[3]);

        for (i = 0; i < tLen; i++) {
            interpolatedF0[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] +
                                 interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) / 4.0;

            f0Deviations[i] = sqrt(((interpolatedF0Set[0][i] - interpolatedF0[i]) * (interpolatedF0Set[0][i] - interpolatedF0[i])
                                    + (interpolatedF0Set[1][i] - interpolatedF0[i]) * (interpolatedF0Set[1][i] - interpolatedF0[i])
                                    + (interpolatedF0Set[2][i] - interpolatedF0[i]) * (interpolatedF0Set[2][i] - interpolatedF0[i])
                                    + (interpolatedF0Set[3][i] - interpolatedF0[i]) * (interpolatedF0Set[3][i] - interpolatedF0[i])) / 3.0);

            if (interpolatedF0[i] > boundaryF0 || interpolatedF0[i] < boundaryF0 / 2.0
                || interpolatedF0[i] > f0Ceil || interpolatedF0[i] < FLOOR_F0) // 70 Hz�ȉ���F0�Ƃ��Ȃ��D
            {
                interpolatedF0[i] = 0.0;
                f0Deviations[i] = 100000.0;
            }
        }

        for (i = 0; i < 4; i++) free(interpolatedF0Set[i]);
    }


    // �������̊J��
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

// �[���������v�Z
void zeroCrossingEngine(double *x, int xLen, double fs,
                        double *eLocations, double *iLocations, double *intervals, int *iLen) {
    int i;
    int *negativeGoingPoints;
    negativeGoingPoints = (int *) malloc(sizeof(int) * xLen);

    int tmp1, tmp2;
    for (i = 0; i < xLen - 1; i++) // ����]����v�Z����͖̂���
    {
        tmp1 = x[i] * x[i + 1] < 0 ? 1 : 0;
        tmp2 = x[i + 1] < x[i] ? 1 : 0;
        negativeGoingPoints[i] = (i + 1) * tmp1 * tmp2;
    }
    negativeGoingPoints[xLen - 1] = 0;

    // �L���C�x���g�̌��o
    int *edges;
    edges = (int *) malloc(sizeof(int) * xLen);
    int count = 0;
    for (i = 0; i < xLen; i++) {
        if (negativeGoingPoints[i] > 0) edges[count++] = negativeGoingPoints[i];
    }
    // �ŏI�߂�l�̌v�Z����
    double *fineEdges;
    fineEdges = (double *) malloc(sizeof(double) * count);
    for (i = 0; i < count; i++) {
        fineEdges[i] = (double) edges[i] - x[edges[i] - 1] / (x[edges[i]] - x[edges[i] - 1]);
    }

    *iLen = count - 1;
    for (i = 0; i < *iLen; i++) {
        intervals[i] = fs / (fineEdges[i + 1] - fineEdges[i]);
        iLocations[i] = (fineEdges[i] + fineEdges[i + 1]) / 2.0 / fs;
        eLocations[i] = fineEdges[i] / fs;
    }
    if (count != 0) eLocations[count - 1] = fineEdges[count - 1] / fs;  //0�̏ꍇ���l��

    free(fineEdges);
    free(edges);
    free(negativeGoingPoints);
}

// �i�b�g�[�����D�}�W�b�N�i���o�[�̂悤�Ɍ����邯�ǂ��ꂪ�����D
void nuttallWindow(int yLen, double *y) {
    int i;
    double tmp;
    for (i = 0; i < yLen; i++) {
        tmp = ((double) (i + 1) - (double) (yLen + 1) / 2.0) / (double) (yLen + 1);
        y[i] = 0.355768 + 0.487396 * cos(2 * PI * tmp) + 0.144232 * cos(4 * PI * tmp) + 0.012604 * cos(6 * PI * tmp);
    }
}




