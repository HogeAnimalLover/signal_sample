//math.hの中でM_PIを定義するためのパラメータ
//コンパイラ次第で必要だったり不要だったりする。
//もちろん円周率定数を独自に定義して使っても問題はない。
//ついでに三角関数は事前にテーブルを作るとより速くなる。
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

	//テスト用入力
	for (int i = 0; i < N; i++) {
		sample_r[i] = cos(2 * M_PI * i / N) + 10 * cos(2 * M_PI * i / N * 15);
		sample_i[i] = 0;
	}


	fft1(result_r, result_i, sample_r, sample_i, N, 1);
	ifft1(re_r, re_i, result_r, result_i, N, 1);

	printf("番号\t入力実部\t入力虚部\t変換後実部\t変換後虚部\t逆変換後実部\t逆変換後虚部\n");

	for (int i = 0; i < N; i++){
		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, sample_r[i], sample_i[i], result_r[i],result_i[i], re_r[i], re_i[i]);
	}

	//環境次第ではウインドがすぐ閉じるため
	fgetc(stdin);
}