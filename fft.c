//math.h�̒���M_PI���`���邽�߂̃p�����[�^
//�R���p�C������ŕK�v��������s�v�������肷��B
//�������~�����萔��Ǝ��ɒ�`���Ďg���Ă����͂Ȃ��B
//���łɎO�p�֐��͎��O�Ƀe�[�u�������Ƃ�葬���Ȃ�B
#define _USE_MATH_DEFINES
#include <math.h>


//���U�t�[���G�ϊ�
//�����F�o�͎����A�o�͋����A���͎����A���͋����A�z��
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

//�t���U�t�[���G�ϊ��i�X�P�[�����Ȃ����́j
//�����F�o�͎����A�o�͋����A���͎����A���͋����A�z��
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

//�����t�[���G�ϊ�
//�����F�o�͎����A�o�͋����A���͎����A���͋����A�z�񒷁A�X�e�b�v��
//�X�e�b�v���͍ċA�p�p�����[�^�A�O������Ăяo���ۂɂ�1���w�肷�邱��
void fft1(double dst_real[], double dst_image[], double src_real[], double src_image[], int length, int step)
{
	if (length <= step || length % (2 * step) ) {
		//�Q�����ł��Ȃ��ꍇ�͒��߂ė��U�t�[���G�ϊ��i�����łȂ����́j
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
		//�Q�����ł���ꍇ�͋����ԖڂƊ�Ԗڂ����ꂼ��ċA��������B
		fft1(&dst_real[                0], &dst_image[                0], &src_real[   0], &src_image[   0], length, 2 * step);
		fft1(&dst_real[length / 2 / step], &dst_image[length / 2 / step], &src_real[step], &src_image[step], length, 2 * step);

		//�������ď����������̂��o�^�t���C���Z�Ō�������B
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



//�����t�t�[���G�ϊ��i�X�P�[�����Ȃ����́j
//�����F�o�͎����A�o�͋����A���͎����A���͋����A�z�񒷁A�X�e�b�v��
//�X�e�b�v���͍ċA�p�p�����[�^�A�O������Ăяo���ۂɂ�1���w�肷�邱��
void ifft1(double dst_real[], double dst_image[], double src_real[], double src_image[], int length, int step)
{
	if (length <= step || length % (2 * step)) {
		//�Q�����ł��Ȃ��ꍇ�͒��߂ė��U�t�[���G�ϊ��i�����łȂ����́j
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
		//�Q�����ł���ꍇ�͋����ԖڂƊ�Ԗڂ����ꂼ��ċA��������B
		ifft1(&dst_real[0], &dst_image[0], &src_real[0], &src_image[0], length, 2 * step);
		ifft1(&dst_real[length / 2 / step], &dst_image[length / 2 / step], &src_real[step], &src_image[step], length, 2 * step);

		//�������ď����������̂��o�^�t���C���Z�Ō�������B
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

