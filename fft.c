//math.hの中でM_PIを定義するためのパラメータ
//コンパイラ次第で必要だったり不要だったりする。
//もちろん円周率定数を独自に定義して使っても問題はない。
//ついでに三角関数は事前にテーブルを作るとより速くなる。
#define _USE_MATH_DEFINES
#include <math.h>


//離散フーリエ変換
//引数：出力実部、出力虚部、入力実部、入力虚部、配列長
void dft(double dst_real[], double dst_image[], double src_real[], double src_image[], int length)
{
	for (int i = 0; i < length; i++) {
		dst_real[i] = dst_image[i] = 0.0f;
		for (int j = 0; j < length; j++) {
			double temp_cos = cos(2 * M_PI * i * j / length), temp_sin = sin(2 * M_PI * i * j / length);

			dst_real[i] += src_real[j] * temp_cos + src_image[j] * temp_sin;
			dst_image[i] += src_image[j] * temp_cos - src_real[j] * temp_sin;
		}
	}
}

//逆離散フーリエ変換（スケールしないもの）
//引数：出力実部、出力虚部、入力実部、入力虚部、配列長
void idft(double dst_real[], double dst_image[], double src_real[], double src_image[], int length)
{
	for (int i = 0; i < length; i++) {
		dst_real[i] = dst_image[i] = 0.0f;
		for (int j = 0; j < length; j++) {
			double temp_cos = cos(2 * M_PI * i * j / length), temp_sin = sin(2 * M_PI * i * j / length);

			dst_real[i] += src_real[j] * temp_cos - src_image[j] * temp_sin;
			dst_image[i] += src_image[j] * temp_cos + src_real[j] * temp_sin;
		}
	}
}

//高速フーリエ変換
//引数：出力実部、出力虚部、入力実部、入力虚部、配列長、ステップ幅
//ステップ幅は再帰用パラメータ、外部から呼び出す際には1を指定すること
void fft1(double dst_real[], double dst_image[], double src_real[], double src_image[], int length, int step)
{
	if (length <= step || length % (2 * step) ) {
		//２分割できない場合は諦めて離散フーリエ変換（高速でないもの）
		for (int i = 0; i < length / step; i++) {
			dst_real[i] = dst_image[i] = 0.0f;
			for (int j = 0; j < length; j += step) {
				double temp_cos = cos(2 * M_PI * i * j / length), temp_sin = sin(2 * M_PI * i * j / length);

				dst_real[i] += src_real[j] * temp_cos + src_image[j] * temp_sin;
				dst_image[i] += src_image[j] * temp_cos - src_real[j] * temp_sin;
			}
		}
	}
	else {
		//２分割できる場合は偶数番目と奇数番目をそれぞれ再帰処理する。
		fft1(&dst_real[                0], &dst_image[                0], &src_real[   0], &src_image[   0], length, 2 * step);
		fft1(&dst_real[length / 2 / step], &dst_image[length / 2 / step], &src_real[step], &src_image[step], length, 2 * step);

		//分割して処理したものをバタフライ演算で結合する。
		for (int i = 0; i < length / 2 / step; i++) {
			double temp_cos = cos(2 * M_PI * i  * step / length), temp_sin = sin(2 * M_PI * i * step / length);
			double temp_real = dst_real[i + length / 2 / step] * temp_cos + dst_image[i + length / 2 / step] * temp_sin;
			double temp_image = dst_image[i + length / 2 / step] * temp_cos - dst_real[i + length / 2 / step] * temp_sin;

			dst_real[i + length / 2 / step] = dst_real[i] - temp_real;
			dst_image[i + length / 2 / step] = dst_image[i] - temp_image;
			dst_real[i] += temp_real;
			dst_image[i] += temp_image;
		}
	}
}



//高速逆フーリエ変換（スケールしないもの）
//引数：出力実部、出力虚部、入力実部、入力虚部、配列長、ステップ幅
//ステップ幅は再帰用パラメータ、外部から呼び出す際には1を指定すること
void ifft1(double dst_real[], double dst_image[], double src_real[], double src_image[], int length, int step)
{
	if (length <= step || length % (2 * step)) {
		//２分割できない場合は諦めて離散フーリエ変換（高速でないもの）
		for (int i = 0; i < length / step; i++) {
			dst_real[i] = dst_image[i] = 0.0f;
			for (int j = 0; j < length; j += step) {
				double temp_cos = cos(2 * M_PI * i * j / length), temp_sin = sin(2 * M_PI * i * j / length);

				dst_real[i] += src_real[j] * temp_cos - src_image[j] * temp_sin;
				dst_image[i] += src_image[j] * temp_cos + src_real[j] * temp_sin;
			}
		}
	}
	else {
		//２分割できる場合は偶数番目と奇数番目をそれぞれ再帰処理する。
		ifft1(&dst_real[0], &dst_image[0], &src_real[0], &src_image[0], length, 2 * step);
		ifft1(&dst_real[length / 2 / step], &dst_image[length / 2 / step], &src_real[step], &src_image[step], length, 2 * step);

		//分割して処理したものをバタフライ演算で結合する。
		for (int i = 0; i < length / 2 / step; i++) {
			double temp_cos = cos(2 * M_PI * i  * step / length), temp_sin = sin(2 * M_PI * i * step / length);
			double temp_real = dst_real[i + length / 2 / step] * temp_cos - dst_image[i + length / 2 / step] * temp_sin;
			double temp_image = dst_image[i + length / 2 / step] * temp_cos + dst_real[i + length / 2 / step] * temp_sin;

			dst_real[i + length / 2 / step] = dst_real[i] - temp_real;
			dst_image[i + length / 2 / step] = dst_image[i] - temp_image;
			dst_real[i] += temp_real;
			dst_image[i] += temp_image;
		}
	}
}

