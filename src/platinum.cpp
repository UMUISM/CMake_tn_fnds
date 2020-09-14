//tn_fnds v0.0.6   2012/3/31
//�ǉ�����Ă���R�����g�ɂ͌�肪���邩������܂���B

#include "main.h"

void getOneFrameResidualSignal(double *x, int xLen, int fs, int positionIndex, double framePeriod, double f0, int fftl, double *pulseLocations, int pCount,
                               double *residualSpec) {
    int i;
    double T0;
    int index, tmpIndex, wLen;

    double tmp;
    double tmpValue = 100000.0; // safeGuard
    for (i = 0; i < pCount; i++) {
        tmp = fabs(pulseLocations[i] - (double) positionIndex * framePeriod);//�t���[���ɍł��߂��p���X�H��T��
        if (tmp < tmpValue) {
            tmpValue = tmp;
            tmpIndex = i;
        }
        index = 1 + (int) (0.5 + pulseLocations[tmpIndex] * fs);//�ł��߂��p���X�H�̃T���v���ʒu?
    }

    T0 = (double) fs / f0;//�P�����̃T���v�����i�����j
    wLen = (int) (0.5 + T0 * 2.0);//2�����̃T���v����

    if (wLen + index - (int) (0.5 + T0) >= xLen) {
        for (i = 0; i < fftl; i++) residualSpec[i] = 0.0;
        return;
    }

    for (i = 0; i < wLen; i++) {
        tmpIndex = i + index - (int) (0.5 + T0);//�ł��߂��p���X�H�̑O��1�������̃C���f�b�N�X
        residualSpec[i] = x[min(xLen - 1, max(0, tmpIndex))] *
                          (0.5 - 0.5 * cos(2.0 * PI * (double) (i + 1) / ((double) (wLen + 1))));//�����|����
    }
    for (; i < fftl / 2; i++) {
        residualSpec[i] = 0.0;
    }
}

void getOnePulseResidualSignal(double *x, int xLen, int fs, double framePeriod, double *f0, int fftl, double *pulseLocations, int pCount,
                               double **residualSpecgram, int *residualSpecgramLength) {
    int i, j;
    double T0;
    int index, tmpIndex, wLen;
    double ff0;
    double f0fi;
    int f0si, f0ei;

    residualSpecgramLength[pCount - 1] = 0;
    for (j = 0; j < pCount - 1; j++) {
        //�T���v���ʒu���v�Z
        index = 1 + (int) (0.5 + pulseLocations[j] * fs);

        //f0���v�Z
        f0fi = pulseLocations[j] / (double) framePeriod;
        f0si = (int) f0fi;
        f0ei = f0si + 1;

        ff0 = (f0[f0si] == 0.0 || f0[f0ei] == 0.0) ? DEFAULT_F0 :
              (f0[f0si] == 0) ? f0[f0ei] :
              (f0[f0ei] == 0) ? f0[f0si] :
              f0[f0si] + (f0[f0ei] - f0[f0si]) * (double) (f0fi - f0si);//�߂��̃t���[���̑O�ォ�璼���⊮



        T0 = (double) fs / ff0;//�P�����̃T���v�����i�����j
        wLen = (int) (0.5 + T0 * 2.0);//2�����̃T���v����

        if (wLen + index - (int) (0.5 + T0) >= xLen) {
            residualSpecgramLength[j] = 0;
            continue;
        }
        residualSpecgramLength[j] = min(fftl - 1, wLen);

        for (i = 0; i < residualSpecgramLength[j]; i++) {
            tmpIndex = i + index - (int) (0.5 + T0);//�ł��߂��p���X�H�̑O��1�������̃C���f�b�N�X
            residualSpecgram[j][i] = x[min(xLen - 1, max(0, tmpIndex))];
//			residualSpecgram[j][i] = x[min( xLen-1, max(0, tmpIndex))] * 
//				(0.5 - 0.5*cos(2.0*PI*(double)(i+1)/((double)(wLen+1))));//�����|����
        }
        for (i = residualSpecgramLength[j]; i < fftl; i++) residualSpecgram[j][i] = 0.0;
    }
}

void PulseResidualWindow(double **residualSpecgram, int *residualSpecgramLength, int pCount) {
    int i, j;
    for (i = 0; i < pCount - 1; i++) {
        for (j = 0; j < residualSpecgramLength[i]; j++) {
            residualSpecgram[i][j] = residualSpecgram[i][j] *
                                     (0.5 - 0.5 * cos(2.0 * PI * (double) (j + 1) / ((double) (residualSpecgramLength[i] + 1))));//�����|����
        }
    }
}

void getFrameResidualIndex(int tLen, int pCount, double framePeriod, double *pulseLocations, int *residualSpecgramIndex) {
    int i, j;
    double tmp;
    double tmpValue;
    int tmpIndex;

    for (j = 0; j < tLen; j++) {
        tmpValue = 100000.0; // safeGuard
        tmpIndex = pCount - 1;
        for (i = 0; i < pCount - 1; i++) {
            tmp = fabs(pulseLocations[i] + -(double) j * framePeriod);//�t���[���ɍł��߂��p���X�H��T��
            if (tmp < tmpValue) {
                tmpValue = tmp;
                tmpIndex = i;
            }
        }
        residualSpecgramIndex[j] = tmpIndex;//�ł��߂��p���X�H�̃C���f�b�N�X
    }
}

int getPulseLocations(double *x, int xLen, double *totalPhase, int vuvNum, int *stList, int *edList, int fs, double framePeriod, int *wedgeList, double *pulseLocations) {
    int i, j;
    int stIndex, edIndex;

    int pCount = 0;
    int numberOfLocation;
    double *tmpPulseLocations, *basePhase;
    tmpPulseLocations = (double *) malloc(sizeof(double) * xLen);
    basePhase = (double *) malloc(sizeof(double) * xLen);


    double tmp;
    for (i = 0; i < vuvNum; i++) {
        stIndex = max(0, (int) ((double) fs * (stList[i]) * framePeriod / 1000.0));  //���̐擪�̃T���v���ʒu
        edIndex = min(xLen - 1, (int) ((double) fs * (edList[i] + 1) * framePeriod / 1000.0 + 0.5) - 1);//���̖����̃T���v���ʒu

        tmp = totalPhase[wedgeList[i]];

        for (j = stIndex; j < edIndex; j++) {
//			basePhase[j] = fmod(totalPhase[j+1]-tmp, 2*PI) - fmod(totalPhase[j]-tmp, 2*PI);//�t���[���Ԃ̈ʑ����H
            basePhase[j] = fmod(totalPhase[j] - tmp + PI * 0.5, 2 * PI);  //tn_fnds �e�T���v���̈ʑ���␳
        }

        numberOfLocation = 0;
        for (j = stIndex; j < edIndex - 1; j++) {
            if (basePhase[j + 1] < basePhase[j])  //�ʑ���2*PI�𒴂���0�ɖ߂���
            {
                tmpPulseLocations[numberOfLocation++] = (double) j / (double) fs;//�[���N���X�ʒu�̎���
            }
        }
        for (j = 0; j < numberOfLocation; j++) pulseLocations[pCount++] = tmpPulseLocations[j];//
    }

    free(basePhase);
    free(tmpPulseLocations);
    return pCount;

}

void getWedgeList(double *x, int xLen, int vuvNum, int *stList, int *edList, int fs, double framePeriod, double *f0, int *wedgeList) {
    int i, j;
    double LowestF0 = 40.0;
    int center, T0;
    double peak;
    int peakIndex = 0;
    double *tmpWav;
    double currentF0;
    tmpWav = (double *) malloc(sizeof(double) * (int) (fs * 2 / LowestF0));

    for (i = 0; i < vuvNum; i++) {
        center = (int) ((stList[i] + edList[i] + 1) / 2);           //���̒����̃t���[���ʒu
        currentF0 = f0[center] == 0.0 ? DEFAULT_F0 : f0[center];//���̒�����F0 �m�C�Y�̈�̏ꍇ�̓f�t�H���g
        T0 = (int) ((fs / currentF0) + 0.5);                //���̒�����F0�̂P�����̃T���v����
//		peakIndex = (int)(((1+center)*framePeriod*fs/1000.0)+0.5);//���̒����̃T���v���ʒu
        peakIndex = (int) (((center) * framePeriod * fs / 1000.0) + 0.5);//���̒����̃T���v���ʒu
//		for(j = 0;j < T0*2;j++)
        for (j = 0; j < T0 * 2 + 1; j++) {
//			tmpWav[j] = x[peakIndex-T0+j-1];
            tmpWav[j] = x[max(0, min(xLen - 1, peakIndex - T0 + j - 1))];//���̒����̂Q�������̃f�[�^
        }
        peak = 0.0;
        peakIndex = 0;
        for (j = 0; j < T0 * 2 + 1; j++)//�g�`�̃s�[�N�̃T���v���ʒu�����o
        {
            if (fabs(tmpWav[j]) > peak) {
                peak = tmpWav[j];
                peakIndex = j;
            }
        }
//		wedgeList[i] = max(0, min(xLen-1, (int)(0.5 + ((center+1)*framePeriod*fs/1000.0)-T0+peakIndex+1.0) - 1));//���̒����̃t���[���̃s�[�N�̃T���v���ʒu
        wedgeList[i] = max(0, min(xLen - 1, (int) (0.5 + ((center) * framePeriod * fs / 1000.0) - T0 + peakIndex + 1.0) - 1));//���̒����̃t���[���̃s�[�N�̃T���v���ʒu
    }
    free(tmpWav);
}

// PLATINUM Version 0.0.4. ���炭���̎d�l�Ŋm��ł��D
// Aperiodicity estimation based on PLATINUM

void pt100(double *x, int xLen, int fs, double *timeAxis, double *f0,
           double **residualSpecgram) {
    int i, j, index;
    double framePeriod = (timeAxis[1] - timeAxis[0]) * 1000.0;

    int fftl = (int) pow(2.0, 1.0 + (int) (log(3.0 * fs / FLOOR_F0 + 1) / log(2.0)));
    int tLen = getSamplesForDIO(fs, xLen, framePeriod);

    int vuvNum;
    vuvNum = 0;
    for (i = 1; i < tLen; i++) {
        if (f0[i] != 0.0 && f0[i - 1] == 0.0) vuvNum++;
    }
    vuvNum += vuvNum - 1; // �����̒��� (�L�����Ɩ�����)
    if (f0[0] == 0) vuvNum++;
    if (f0[tLen - 1] == 0) vuvNum++;

    int stCount, edCount;
    int *stList, *edList;
    stList = (int *) malloc(sizeof(int) * vuvNum);
    edList = (int *) malloc(sizeof(int) * vuvNum);
    edCount = 0;

    stList[0] = 0;
    stCount = 1;
    index = 1;
    if (f0[0] != 0) {
        for (i = 1; i < tLen; i++) {
            if (f0[i] == 0 && f0[i - 1] != 0) {
                edList[0] = i - 1;
                edCount++;
                stList[1] = i;
                stCount++;
                index = i;
            }
        }
    }

    edList[vuvNum - 1] = tLen - 1;
    for (i = index; i < tLen; i++) {
        if (f0[i] != 0.0 && f0[i - 1] == 0.0) {
            edList[edCount++] = i - 1;
            stList[stCount++] = i;
        }
        if (f0[i] == 0.0 && f0[i - 1] != 0.0) {
            edList[edCount++] = i - 1;
            stList[stCount++] = i;
        }
    }

    int *wedgeList;
    wedgeList = (int *) malloc(sizeof(int) * vuvNum);
    getWedgeList(x, xLen, vuvNum, stList, edList, fs, framePeriod, f0, wedgeList);//�������̃s�[�N�ʒu���擾

    double *signalTime, *f0interpolatedRaw, *totalPhase;
    double *fixedF0;
    fixedF0 = (double *) malloc(sizeof(double) * tLen);
    signalTime = (double *) malloc(sizeof(double) * xLen);
    f0interpolatedRaw = (double *) malloc(sizeof(double) * xLen);
    totalPhase = (double *) malloc(sizeof(double) * xLen);

    for (i = 0; i < tLen; i++) fixedF0[i] = f0[i] == 0 ? DEFAULT_F0 : f0[i]; //F0��0�Ȃ�f�t�H���g�ɕ␳
    for (i = 0; i < xLen; i++) signalTime[i] = (double) i / (double) fs;       //�T���v���ʒu�̎���
    interp1(timeAxis, fixedF0, tLen, signalTime, xLen, f0interpolatedRaw);//�e�T���v����F0
    totalPhase[0] = f0interpolatedRaw[0] * 2 * PI / (double) fs;                 //�e�T���v���̈ʑ�
    for (i = 1; i < xLen; i++) totalPhase[i] = totalPhase[i - 1] + f0interpolatedRaw[i] * 2 * PI / (double) fs;

    double *pulseLocations;
    pulseLocations = (double *) malloc(sizeof(double) * xLen);
    int pCount;
    pCount = getPulseLocations(x, xLen, totalPhase, vuvNum, stList, edList, fs, framePeriod, wedgeList, pulseLocations);//�ʑ����傫�������T���v���ʒu�̎���

    double *tmpResidualSpec;
    tmpResidualSpec = (double *) malloc(sizeof(double) * fftl);
    double currentF0;

    for (j = 0; j < fftl / 2; j++) residualSpecgram[0][j] = 0.0;
    for (i = 1; i < tLen; i++) {
        currentF0 = f0[i] <= FLOOR_F0 ? DEFAULT_F0 : f0[i];  //�t���[����F0 �����Ȃ�f�t�H���g�ɕ␳
        getOneFrameResidualSignal(x, xLen, fs, i, framePeriod / 1000.0, currentF0, fftl, pulseLocations, pCount, //�ł��߂��p���X�H�O��1�������̔g�`�ɑ������������̂��擾
                                  tmpResidualSpec);
        for (j = 0; j < fftl / 2; j++) residualSpecgram[i][j] = tmpResidualSpec[j];
    }

    free(fixedF0);
    free(tmpResidualSpec);
    free(pulseLocations);
    free(totalPhase);
    free(f0interpolatedRaw);
    free(signalTime);
    free(wedgeList);
    free(edList);
    free(stList);
    return;
}

//residualSpecgram �g�`�̃R�s�[������B���������m�ۂ����ɓn���B���̊֐��ɂ��K�v�ȃ��������m�ۂ����B
//residualSpecgramLength �g�`�̒���������B���������m�ۂ����ɓn���B���̊֐��ɂ��K�v�ȃ��������m�ۂ����B
//residualSpecgramIndex�@�e�t���[���̔g�`���w�肷��C���f�b�N�X������B
//�߂�l�g�`�̐�(pCount)
int pt101(double *x, int xLen, int fs, double *timeAxis, double *f0,
          double ***residualSpecgram, int **residualSpecgramLength, int *residualSpecgramIndex) {
    int i, index;
    double framePeriod = (timeAxis[1] - timeAxis[0]) * 1000.0;

    int fftl = (int) pow(2.0, 1.0 + (int) (log(3.0 * fs / FLOOR_F0 + 1) / log(2.0)));
    int tLen = getSamplesForDIO(fs, xLen, framePeriod);

    int vuvNum;
//	vuvNum = 0;
    vuvNum = 1;    //tn_fuds
    for (i = 1; i < tLen; i++) {
        if (f0[i] != 0.0 && f0[i - 1] == 0.0) vuvNum++;    //�������L��
        if (f0[i] == 0.0 && f0[i - 1] != 0.0) vuvNum++;    //�L��������  tn_fnds
    }
//	vuvNum+=vuvNum-1; // �����̒��� (�L�����Ɩ�����)  tn_fnds �R�����g�A�E�g
//	if(f0[0] == 0) vuvNum++;  tn_fnds �R�����g�A�E�g
//	if(f0[tLen-1] == 0) vuvNum++;  tn_fnds �R�����g�A�E�g

    int stCount, edCount;
    int *stList, *edList;
    stList = (int *) malloc(sizeof(int) * vuvNum);
    edList = (int *) malloc(sizeof(int) * vuvNum);
    edCount = 0;

    stList[0] = 0;
    stCount = 1;
    index = 1;
    if (f0[0] != 0)    //�L������n�܂�ꍇ
    {
        for (i = 1; i < tLen; i++) {
            if (f0[i] == 0 && f0[i - 1] != 0)    //�L��������
            {
                edList[0] = i - 1;
                edCount++;
                stList[1] = i;
                stCount++;
                index = i;

                break;    //tn_fnds
            }
        }
    }

    edList[vuvNum - 1] = tLen - 1;
    for (i = index; i < tLen; i++) {
        if (f0[i] != 0.0 && f0[i - 1] == 0.0) //�������L��
        {
            edList[edCount++] = i - 1;
            stList[stCount++] = i;
        }
        if (f0[i] == 0.0 && f0[i - 1] != 0.0) //�L��������
        {
            edList[edCount++] = i - 1;
            stList[stCount++] = i;
        }
    }

    int *wedgeList;
    wedgeList = (int *) malloc(sizeof(int) * vuvNum);
    getWedgeList(x, xLen, vuvNum, stList, edList, fs, framePeriod, f0, wedgeList);//�������̃s�[�N�ʒu���擾

    double *signalTime, *f0interpolatedRaw, *totalPhase;
    double *fixedF0;
    fixedF0 = (double *) malloc(sizeof(double) * tLen);
    signalTime = (double *) malloc(sizeof(double) * xLen);
    f0interpolatedRaw = (double *) malloc(sizeof(double) * xLen);
    totalPhase = (double *) malloc(sizeof(double) * xLen);

    for (i = 0; i < tLen; i++) fixedF0[i] = f0[i] == 0 ? DEFAULT_F0 : f0[i]; //F0��0�Ȃ�f�t�H���g�ɕ␳
    for (i = 0; i < xLen; i++) signalTime[i] = (double) i / (double) fs;       //�T���v���ʒu�̎���
    interp1(timeAxis, fixedF0, tLen, signalTime, xLen, f0interpolatedRaw);//�e�T���v����F0
    totalPhase[0] = f0interpolatedRaw[0] * 2 * PI / (double) fs;                 //�e�T���v���̈ʑ�
    for (i = 1; i < xLen; i++) totalPhase[i] = totalPhase[i - 1] + f0interpolatedRaw[i] * 2 * PI / (double) fs;

    double *pulseLocations;
    pulseLocations = (double *) malloc(sizeof(double) * xLen);
    int pCount;
    pCount = getPulseLocations(x, xLen, totalPhase, vuvNum, stList,
                               edList, fs, framePeriod, wedgeList, pulseLocations);//�ʑ����傫�������T���v���ʒu�̎���

//tn_fnds�f�o�b�O
//zeroXToFile(x, xLen, f0interpolatedRaw, totalPhase, pCount, pulseLocations, fs, vuvNum, wedgeList);

    pCount++;//�����O�̃_�~�[�p���X��ǉ�

    *residualSpecgram = (double **) malloc(sizeof(double *) * pCount);
    for (i = 0; i < pCount; i++) (*residualSpecgram)[i] = (double *) malloc(sizeof(double) * fftl);
    *residualSpecgramLength = (int *) malloc(sizeof(int) * pCount);

    getOnePulseResidualSignal(x, xLen, fs, framePeriod / 1000.0, f0, fftl, pulseLocations, pCount,
                              *residualSpecgram, *residualSpecgramLength);

    getFrameResidualIndex(tLen, pCount, framePeriod / 1000, pulseLocations, residualSpecgramIndex);


//	pulseToFile(pCount, pulseLocations, *residualSpecgramLength);

    free(fixedF0);
    free(pulseLocations);
    free(totalPhase);
    free(f0interpolatedRaw);
    free(signalTime);
    free(wedgeList);
    free(edList);
    free(stList);
    return pCount;
}
