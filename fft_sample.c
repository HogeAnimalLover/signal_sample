//math.h�̒���M_PI���`���邽�߂̃p�����[�^
//�R���p�C������ŕK�v��������s�v�������肷��B
//�������~�����萔��Ǝ��ɒ�`���Ďg���Ă����͂Ȃ��B
//���łɎO�p�֐��͎��O�Ƀe�[�u�������Ƃ�葬���Ȃ�B
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include "fft.h"


#define N 1000

int main()
{
	double sample_r[N], sample_i[N];
	double result_r[N], result_i[N];
	double re_r[N], re_i[N];

	//�e�X�g�p����
	for (int i = 0; i < N; i++) {
		sample_r[i] = cos(2 * M_PI * i / N) + 10 * cos(2 * M_PI * i / N * 15);
		sample_i[i] = 0;
	}


	fft1(result_r, result_i, sample_r, sample_i, N, 1);
	ifft1(re_r, re_i, result_r, result_i, N, 1);

	printf("�ԍ�\t���͎���\t���͋���\t�ϊ������\t�ϊ��㋕��\t�t�ϊ������\t�t�ϊ��㋕��\n");

	for (int i = 0; i < N; i++){
		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, sample_r[i], sample_i[i], result_r[i],result_i[i], re_r[i], re_i[i]);
	}

	//������ł̓E�C���h���������邽��
	fgetc(stdin);
}