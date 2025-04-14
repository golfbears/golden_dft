#include <math.h>
typedef struct
{
    double x;  //REAL
    double y;  //IMG
}Complex;
static Complex dft_matrix[FFT_LENGTH][FFT_LENGTH];
static void Add(Complex z1, Complex z2, Complex* answer)  
{
    answer->x = z1.x + z2.x;
    answer->y = z1.y + z2.y;
}
static void Sub(Complex z1, Complex z2, Complex* answer)  
{
    answer->x = z1.x - z2.x;
    answer->y = z1.y - z2.y;
}
static void Mul(Complex z1, Complex z2, Complex* answer)  
{
    answer->x = z1.x * z2.x - z1.y * z2.y;
    answer->y = z1.x * z2.y + z1.y * z2.x;
}

void kdft(double* pr, double* pi, int n, int k, double* fr, double* fi, int l, int il) {

    double p = 6.283185306 / (1.0 * FFT_LENGTH);
    double p_k;
    Complex tmp;
    Complex dft_k;

    if (l == 0) {
        for (int i = 0; i < FFT_LENGTH; i++) {
            p_k = p * i;
            for (int j = 0; j < FFT_LENGTH; j++) {
                dft_matrix[i][j].x = cos(p_k * j);
                dft_matrix[i][j].y = -sin(p_k * j);
            }
        }
    }
    else {
        for (int i = 0; i < FFT_LENGTH; i++) {
            p_k = p * i;
            for (int j = 0; j < FFT_LENGTH; j++) {
                dft_matrix[i][j].x = cos(p_k * j) / FFT_LENGTH;
                dft_matrix[i][j].y = sin(p_k * j) / FFT_LENGTH;
            }
        }
    }
    for (int i = 0; i < FFT_LENGTH; i++) {
        dft_k.x = 0.0f; dft_k.y = 0.0f;
        for (int j = 0; j < FFT_LENGTH; j++) {
            tmp.x = pr[j];
            tmp.y = pi[j];
            Mul(dft_matrix[i][j], tmp, &tmp);
            Add(tmp, dft_k, &dft_k);
        }
        fr[i] = dft_k.x;
        fi[i] = dft_k.y;
    }
    if (il != 0) {
        for (int i = 0; i <= n - 1; i++)
        {
            pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
            if (fabs(fr[i]) < 0.000001 * fabs(fi[i]))
            {
                if ((fi[i] * fr[i]) > 0) pi[i] = 90.0;
                else pi[i] = -90.0;
            }
            else
                pi[i] = atan(fi[i] / fr[i]) * 360.0 / 6.283185306;
        }
    }

    return;
}

#define FFT_FORWARD 0
#define FFT_BACKWARD 1
#define N0_NEED_MODULE_ANGLE 0
#define NEED_MODULE_ANGLE 1
#define FRAME_LENGTH 2000
void check_fft(short* in, short* out, int fft_len) {
    double dft_fre_real_frame[FRAME_LENGTH];
    double dft_fre_img_frame[FRAME_LENGTH];
    double ref_window_frame[FRAME_LENGTH];
    double in_window_frame[FRAME_LENGTH];
    double out_window_frame[FRAME_LENGTH];
    double dft_ref_window_frame[FRAME_LENGTH];
    double dft_out_window_frame[FRAME_LENGTH];
    for (int i = 0;i < FRAME_LENGTH;i++) {
        in_window_frame[i] = (double)in[i] / 32768.0;
        ref_window_frame[i] = 0.0f;
        dft_fre_real_frame[i] = 0.0f;
        dft_fre_img_frame[i] = 0.0f;
    }

    kdft(in_window_frame, ref_window_frame, FFT_LENGTH, FFT_PWR_INDEX, dft_fre_real_frame, dft_fre_img_frame, FFT_FORWARD, NEED_MODULE_ANGLE);
    kdft(dft_fre_real_frame, dft_fre_img_frame, FFT_LENGTH, FFT_PWR_INDEX, dft_out_window_frame, dft_ref_window_frame, FFT_BACKWARD, NEED_MODULE_ANGLE);
    return;
}
