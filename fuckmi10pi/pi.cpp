/* Please check the following macros before compiling */

#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#ifndef FFT_ERROR_MARGIN
#define FFT_ERROR_MARGIN 0.4  /* must be < 0.5 */
#endif


#ifndef dgt_int
#ifdef USE_DGT_LONG_INT
#define dgt_int long long int /* 64 bit int */
#define DGT_INT_MAX LLONG_MAX /* 64 bit int max */
#else
#ifdef USE_DGT_NORMAL_INT
#define dgt_int int           /* 32 bit int */
#define DGT_INT_MAX INT_MAX   /* 32 bit int max */
#else
#define dgt_int short int     /* 16 bit int */
#define DGT_INT_MAX SHRT_MAX  /* 16 bit int max */
#endif
#endif
#endif

#ifndef fft_int
#ifdef USE_FFT_LONG_INT
#define fft_int long long int /* 64 bit int */
#define FFT_INT_MAX LLONG_MAX /* 64 bit int max */
#else
#define fft_int int           /* 32 bit int */
#define FFT_INT_MAX INT_MAX   /* 32 bit int max */
#endif
#endif

#ifndef fft_float
#ifdef USE_FFT_LONG_DOUBLE
#define fft_float long double
#define FFT_FLOAT_EPSILON LDBL_EPSILON
#define xlint long int
#define XLINT_MIN LONG_MIN
#define XLINT_MAX LONG_MAX
#define xlfloat long double
#define XLFLOAT_EPSILON LDBL_EPSILON
#define XLFLOAT_NAN (-LDBL_MAX)
fft_float sqrtxl(fft_float);
fft_float atanxl(fft_float);
fft_float cosxl(fft_float);
fft_float sinxl(fft_float);
#define fft_sqrt(x) sqrtxl(x)
#define fft_atan(x) atanxl(x)
#define fft_cos(x) cosxl(x)
#define fft_sin(x) sinxl(x)
#define FC_HALF ((fft_float) 0.5L)
#define FC_PI_2 ((fft_float) 1.570796326794896619231321691639751442098584699687552910487472296153908203143104499314L)
/* -------- WRN=cos(N*pi/20000), WIN=sin(N*pi/20000) -------- */
#define WR5000 ((fft_float) 0.707106781186547524400844362104849039284835937688474036588339868995366239231053519425L)
#define WR2500 ((fft_float) 0.923879532511286756128183189396788286822416625863642486115097731280535007501102358714L)
#define WI2500 ((fft_float) 0.382683432365089771728459984030398866761344562485627041433800635627546033960089692237L)
#define WR1250 ((fft_float) 0.980785280403230449126182236134239036973933730893336095002916088545306513549605063915L)
#define WI1250 ((fft_float) 0.195090322016128267848284868477022240927691617751954807754502089494763318785924580225L)
#define WR3750 ((fft_float) 0.831469612302545237078788377617905756738560811987249963446124590227637920144642343981L)
#define WI3750 ((fft_float) 0.555570233019602224742830813948532874374937190754804045924153528202949247577480068383L)
#else
#define fft_float double
#define FFT_FLOAT_EPSILON DBL_EPSILON
#define fft_sqrt(x) sqrt(x)
#define fft_atan(x) atan(x)
#define fft_cos(x) cos(x)
#define fft_sin(x) sin(x)
#define FC_HALF ((fft_float) 0.5)
#define FC_PI_2 ((fft_float) 1.570796326794896619231321691639751442098584)
#define WR5000 ((fft_float) 0.707106781186547524400844362104849039284835)
#define WR2500 ((fft_float) 0.923879532511286756128183189396788286822416)
#define WI2500 ((fft_float) 0.382683432365089771728459984030398866761344)
#define WR1250 ((fft_float) 0.980785280403230449126182236134239036973933)
#define WI1250 ((fft_float) 0.195090322016128267848284868477022240927691)
#define WR3750 ((fft_float) 0.831469612302545237078788377617905756738560)
#define WI3750 ((fft_float) 0.555570233019602224742830813948532874374937)
#endif
#endif


#ifdef USE_CDFT_PTHREADS
#ifndef USE_CDFT_THREADS
#define USE_CDFT_THREADS
#ifndef CDFT_THREADS_BEGIN_N
#define CDFT_THREADS_BEGIN_N 8192
#endif
#ifndef CDFT_4THREADS_BEGIN_N
#define CDFT_4THREADS_BEGIN_N 65536
#endif
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define cdft_thread_t pthread_t
#define cdft_thread_create(thp,func,argp) { \
    if (pthread_create(thp, NULL, func, (void *) argp) != 0) { \
        fprintf(stderr, "cdft thread error\n"); \
        exit(1); \
    } \
}
#define cdft_thread_wait(th) { \
    if (pthread_join(th, NULL) != 0) { \
        fprintf(stderr, "cdft thread error\n"); \
        exit(1); \
    } \
}
#endif
#endif /* USE_CDFT_PTHREADS */


#ifdef USE_CDFT_WINTHREADS
#ifndef USE_CDFT_THREADS
#define USE_CDFT_THREADS
#ifndef CDFT_THREADS_BEGIN_N
#define CDFT_THREADS_BEGIN_N 32768
#endif
#ifndef CDFT_4THREADS_BEGIN_N
#define CDFT_4THREADS_BEGIN_N 524288
#endif
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#define cdft_thread_t HANDLE
#define cdft_thread_create(thp,func,argp) { \
    DWORD thid; \
    *(thp) = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) func, (LPVOID) argp, 0, &thid); \
    if (*(thp) == 0) { \
        fprintf(stderr, "cdft thread error\n"); \
        exit(1); \
    } \
}
#define cdft_thread_wait(th) { \
    WaitForSingleObject(th, INFINITE); \
    CloseHandle(th); \
}
#endif
#endif /* USE_CDFT_WINTHREADS */


#ifndef CDFT_LOOP_DIV  /* control of the CDFT's speed & tolerance */
#define CDFT_LOOP_DIV 32
#endif

#ifndef RDFT_LOOP_DIV  /* control of the RDFT's speed & tolerance */
#define RDFT_LOOP_DIV 64
#endif

#ifndef DCST_LOOP_DIV  /* control of the DCT,DST's speed & tolerance */
#define DCST_LOOP_DIV 64
#endif


#ifndef F_CMUL_BUF
#define F_CMUL_BUF 512
#endif

#ifndef SWAP_FILE_FFT1
#define SWAP_FILE_FFT1 "tmp1fft.swp"
#endif
#ifndef SWAP_FILE_VAR1
#define SWAP_FILE_VAR1 "tmp1var.swp"
#endif
#ifndef SWAP_FILE_VAR2
#define SWAP_FILE_VAR2 "tmp2var.swp"
#endif
#ifndef SWAP_FILE_TMPA
#define SWAP_FILE_TMPA "tmpavar.swp"
#endif
#ifndef SWAP_FILE_TMPC
#define SWAP_FILE_TMPC "tmpcvar.swp"
#endif


#define PI_FFTC_VER "ver. LG1.1.2-MP1.5.2af"


void mp_load_0(fft_int n, fft_int radix, fft_int out[]);
void mp_load_1(fft_int n, fft_int radix, fft_int out[]);
void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[]);
fft_int mp_cmp(fft_int n, fft_int radix, fft_int in1[], fft_int in2[]);
void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
void mp_imul(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[]);
fft_int mp_idiv(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[]);
void mp_idiv_2(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
fft_float mp_mul_radix_test(fft_int n, fft_int radix, fft_int nfft,
	fft_float tmpfft[], fft_int ip[], fft_float w[]);
void mp_mul(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int tmp[], fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
	fft_float tmp3fft[], fft_int ip[], fft_float w[]);
void mp_squ(fft_int n, fft_int radix, fft_int in[], fft_int out[], fft_int tmp[],
	fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
	fft_int ip[], fft_float w[]);
void mp_mulhf(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int tmp[], fft_int nfft, fft_float in1fft[], fft_float tmpfft[],
	fft_int ip[], fft_float w[]);
void mp_mulhf_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[], fft_int in2[],
	fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
	fft_int ip[], fft_float w[]);
void mp_squhf_use_infft(fft_int n, fft_int radix, fft_float infft[], fft_int in[],
	fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
	fft_int ip[], fft_float w[]);
void mp_mulh(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int nfft, fft_float in1fft[], fft_float outfft[],
	fft_int ip[], fft_float w[]);
void mp_squh(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int nfft, fft_float outfft[], fft_int ip[], fft_float w[]);
fft_int mp_inv(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
fft_int mp_sqrt(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
fft_int mp_invisqrt(fft_int n, fft_int radix, fft_int in, fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
void mp_sprintf(fft_int n, fft_int log10_radix, fft_int in[], char out[]);
void mp_sscanf(fft_int n, fft_int log10_radix, const char in[], fft_int out[]);
void mp_fprintf(fft_int n, fft_int log10_radix, fft_int in[], FILE* fout);
fft_int mp_chksum(fft_int n, fft_int in[]);


int main()
{
	fft_int nfft, log2_nfft, radix, log10_radix, n, npow, nprc, sum;
	fft_float err;
	fft_int* a, * b, * c, * e, * i1, * i2, * ip;
	fft_float* d1, * d2, * d3, * w;
	double x, memsize;
	time_t t_1, t_2;
	FILE* f_out;
#ifdef PI_OUT_LOGFILE
	FILE* f_log;
	f_log = fopen("pi.log", "w");
#endif

	printf("Calculation of PI using FFT and AGM, %s\n", PI_FFTC_VER);
#ifdef PI_OUT_LOGFILE
	fprintf(f_log, "Calculation of PI using FFT and AGM, %s\n", PI_FFTC_VER);
#endif
	printf("length of FFT =?\n");
	scanf("%lg", &x);
	nfft = (fft_int)x;

	printf("initializing...\n");
	for (log2_nfft = 1; (1 << log2_nfft) < nfft; log2_nfft++);
	nfft = 1 << log2_nfft;
	n = nfft + 2;
	ip = (fft_int*)malloc((3 + (fft_int)sqrt(0.5 * nfft)) * sizeof(fft_int));
	memsize = (3 + (fft_int)sqrt(0.5 * nfft)) * sizeof(fft_int);
	w = (fft_float*)malloc(nfft / 2 * sizeof(fft_float));
	memsize += nfft / 2 * sizeof(fft_float);
	a = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	b = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	c = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	e = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	i1 = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	i2 = (fft_int*)malloc((n + 2) * sizeof(fft_int));
	memsize += 6 * ((n + 2) * sizeof(fft_int));
	d1 = (fft_float*)malloc((nfft + 2) * sizeof(fft_float));
	d2 = (fft_float*)malloc((nfft + 2) * sizeof(fft_float));
	d3 = (fft_float*)malloc((nfft + 2) * sizeof(fft_float));
	memsize += 3 * ((nfft + 2) * sizeof(fft_float));
	if (d3 == NULL) {
		printf("Allocation Failure!\n");
		exit(1);
	}
	ip[0] = 0;
	/* ---- radix test ---- */
	log10_radix = 1;
	radix = 10;
	err = mp_mul_radix_test(n, radix, nfft, d1, ip, w);
	err += FFT_FLOAT_EPSILON * (n * radix * radix / 4);
	while (100 * err < FFT_ERROR_MARGIN && radix <= FFT_INT_MAX / 20) {
		err *= 100;
		log10_radix++;
		radix *= 10;
	}
	printf("nfft= %.0f\nradix= %.0f\nerror_margin= %g\nmem_alloc_size= %.0f\n", (double)nfft, (double)radix, (double)err, memsize);
	printf("calculating %.0f digits of PI...\n", (double)(log10_radix * (n - 2)));
#ifdef PI_OUT_LOGFILE
	fprintf(f_log, "nfft= %.0f\nradix= %.0f\nerror_margin= %g\nmem_alloc_size= %.0f\n", (double)nfft, (double)radix, (double)err, memsize);
	fprintf(f_log, "calculating %.0f digits of PI...\n", (double)(log10_radix * (n - 2)));
#endif
	/* ---- time check ---- */
	time(&t_1);
	/*
	 * ---- a formula based on the AGM (Arithmetic-Geometric Mean) ----
	 *   c = sqrt(0.125);
	 *   a = 1 + 3 * c;
	 *   b = sqrt(a);
	 *   e = b - 0.625;
	 *   b = 2 * b;
	 *   c = e - c;
	 *   a = a + e;
	 *   npow = 4;
	 *   do {
	 *       npow = 2 * npow;
	 *       e = (a + b) / 2;
	 *       b = sqrt(a * b);
	 *       e = e - b;
	 *       b = 2 * b;
	 *       c = c - e;
	 *       a = e + b;
	 *   } while (e > SQRT_SQRT_EPSILON);
	 *   e = e * e / 4;
	 *   a = a + b;
	 *   pi = (a * a - e - e / 2) / (a * c - e) / npow;
	 * ---- modification ----
	 *   This is a modified version of Gauss-Legendre formula
	 *   (by T.Ooura). It is faster than original version.
	 * ---- reference ----
	 *   1. E.Salamin,
	 *      Computation of PI Using Arithmetic-Geometric Mean,
	 *      Mathematics of Computation, Vol.30 1976.
	 *   2. R.P.Brent,
	 *      Fast Multiple-Precision Evaluation of Elementary Functions,
	 *      J. ACM 23 1976.
	 *   3. D.Takahasi, Y.Kanada,
	 *      Calculation of PI to 51.5 Billion Decimal Digits on
	 *      Distributed Memoriy Parallel Processors,
	 *      Transactions of Information Processing Society of Japan,
	 *      Vol.39 No.7 1998.
	 *   4. T.Ooura,
	 *      Improvement of the PI Calculation Algorithm and
	 *      Implementation of Fast Multiple-Precision Computation,
	 *      Information Processing Society of Japan SIG Notes,
	 *      98-HPC-74, 1998.
	 */
	 /* ---- c = 1 / sqrt(8) ---- */
	mp_invisqrt(n, radix, 8, c, i1, i2, nfft, d1, d2, ip, w);
	/* ---- a = 1 + 3 * c ---- */
	mp_imul(n, radix, c, 3, e);
	mp_sscanf(n, log10_radix, "1", a);
	mp_add(n, radix, a, e, a);
	/* ---- b = sqrt(a) ---- */
	mp_sqrt(n, radix, a, b, i1, i2, nfft, d1, d2, ip, w);
	/* ---- e = b - 0.625 ---- */
	mp_sscanf(n, log10_radix, "0.625", e);
	mp_sub(n, radix, b, e, e);
	/* ---- b = 2 * b ---- */
	mp_add(n, radix, b, b, b);
	/* ---- c = e - c ---- */
	mp_sub(n, radix, e, c, c);
	/* ---- a = a + e ---- */
	mp_add(n, radix, a, e, a);
	/* ---- time check ---- */
	time(&t_2);
	sum = mp_chksum(n, c);
	printf("AGM iteration,\ttime= %.0f,\tchksum= %x\n", difftime(t_2, t_1), (int)sum);
#ifdef PI_OUT_LOGFILE
	fprintf(f_log, "AGM iteration,\ttime= %.0f,\tchksum= %x\n", difftime(t_2, t_1), (int)sum);
	fflush(f_log);
#endif
	npow = 4;
	do {
		npow *= 2;
		/* ---- e = (a + b) / 2 ---- */
		mp_add(n, radix, a, b, e);
		mp_idiv_2(n, radix, e, e);
		/* ---- b = sqrt(a * b) ---- */
		mp_mul(n, radix, a, b, a, i1, nfft, d1, d2, d3, ip, w);
		mp_sqrt(n, radix, a, b, i1, i2, nfft, d1, d2, ip, w);
		/* ---- e = e - b ---- */
		mp_sub(n, radix, e, b, e);
		/* ---- b = 2 * b ---- */
		mp_add(n, radix, b, b, b);
		/* ---- c = c - e ---- */
		mp_sub(n, radix, c, e, c);
		/* ---- a = e + b ---- */
		mp_add(n, radix, e, b, a);
		/* ---- convergence check ---- */
		nprc = -e[1];
		if (e[0] == 0) {
			nprc = n;
		}
		/* ---- time check ---- */
		time(&t_2);
		sum = mp_chksum(n, c);
		printf("precision= %.0f,\ttime= %.0f,\tchksum= %x\n", (double)(4 * nprc * log10_radix), difftime(t_2, t_1), (int)sum);
#ifdef PI_OUT_LOGFILE
		fprintf(f_log, "precision= %.0f,\ttime= %.0f,\tchksum= %x\n", (double)(4 * nprc * log10_radix), difftime(t_2, t_1), (int)sum);
		fflush(f_log);
#endif
	} while (4 * nprc <= n);
	/* ---- e = e * e / 4 (half precision) ---- */
	mp_idiv_2(n, radix, e, e);
	mp_squh(n, radix, e, e, nfft, d1, ip, w);
	/* ---- a = a + b ---- */
	mp_add(n, radix, a, b, a);
	/* ---- a = (a * a - e - e / 2) / (a * c - e) / npow ---- */
	mp_mulhf(n, radix, a, c, c, i1, nfft, d1, d2, ip, w);
	mp_sub(n, radix, c, e, c);
	mp_inv(n, radix, c, b, i1, i2, nfft, d2, d3, ip, w);
	mp_squhf_use_infft(n, radix, d1, a, a, i1, nfft, d2, ip, w);
	mp_sub(n, radix, a, e, a);
	mp_idiv_2(n, radix, e, e);
	mp_sub(n, radix, a, e, a);
	mp_mul(n, radix, a, b, a, i1, nfft, d1, d2, d3, ip, w);
	mp_idiv(n, radix, a, npow, a);
	/* ---- time check ---- */
	time(&t_2);
	sum = mp_chksum(n, a);
	/* ---- output ---- */
	f_out = fopen("pi.dat", "w");
	printf("writing pi.dat...\n");
	mp_fprintf(n - 1, log10_radix, a, f_out);
	fclose(f_out);
	free(d3);
	free(d2);
	free(d1);
	free(i2);
	free(i1);
	free(e);
	free(c);
	free(b);
	free(a);
	free(w);
	free(ip);
	/* ---- difftime ---- */
	printf("Total %.0f sec. (real time),\tchksum= %x\n", difftime(t_2, t_1), (int)sum);
#ifdef PI_OUT_LOGFILE
	fprintf(f_log, "Total %.0f sec. (real time),\tchksum= %x\n", difftime(t_2, t_1), (int)sum);
	fclose(f_log);
#endif
	return 0;
}


/* -------- multiple precision routines -------- */


#include <math.h>
#include <float.h>
#include <stdio.h>

/* ---- floating point format ----
	data := data[0] * pow(radix, data[1]) *
			(data[2] + data[3]/radix + data[4]/radix/radix + ...),
	data[0]       : sign (1;data>0, -1;data<0, 0;data==0)
	data[1]       : exponent (0;data==0)
	data[2...n+1] : digits
   ---- function prototypes ----
	void mp_load_0(fft_int n, fft_int radix, fft_int out[]);
	void mp_load_1(fft_int n, fft_int radix, fft_int out[]);
	void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[]);
	fft_int mp_cmp(fft_int n, fft_int radix, fft_int in1[], fft_int in2[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_imul(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[]);
	fft_int mp_idiv(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[]);
	void mp_idiv_2(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	fft_float mp_mul_radix_test(fft_int n, fft_int radix, fft_int nfft,
			fft_float tmpfft[], fft_int ip[], fft_float w[]);
	void mp_mul(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
			fft_int tmp[], fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
			fft_float tmp3fft[], fft_int ip[], fft_float w[]);
	void mp_squ(fft_int n, fft_int radix, fft_int in[], fft_int out[], fft_int tmp[],
			fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
			fft_int ip[], fft_float w[]);
	void mp_mulhf(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
			fft_int tmp[], fft_int nfft, fft_float in1fft[], fft_float tmpfft[],
			fft_int ip[], fft_float w[]);
	void mp_mulhf_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[], fft_int in2[],
			fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
			fft_int ip[], fft_float w[]);
	void mp_squhf_use_infft(fft_int n, fft_int radix, fft_float infft[], fft_int in[],
			fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
			fft_int ip[], fft_float w[]);
	void mp_mulh(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
			fft_int nfft, fft_float in1fft[], fft_float outfft[],
			fft_int ip[], fft_float w[]);
	void mp_squh(fft_int n, fft_int radix, fft_int in[], fft_int out[],
			fft_int nfft, fft_float outfft[], fft_int ip[], fft_float w[]);
	fft_int mp_inv(fft_int n, fft_int radix, fft_int in[], fft_int out[],
			fft_int tmp1[], fft_int tmp2[], fft_int nfft,
			fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
	fft_int mp_sqrt(fft_int n, fft_int radix, fft_int in[], fft_int out[],
			fft_int tmp1[], fft_int tmp2[], fft_int nfft,
			fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
	fft_int mp_invisqrt(fft_int n, fft_int radix, fft_int in, fft_int out[],
			fft_int tmp1[], fft_int tmp2[], fft_int nfft,
			fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[]);
	void mp_sprintf(fft_int n, fft_int log10_radix, fft_int in[], char out[]);
	void mp_sscanf(fft_int n, fft_int log10_radix, char in[], fft_int out[]);
	void mp_fprintf(fft_int n, fft_int log10_radix, fft_int in[], FILE *fout);
   ----
*/


/* -------- mp_load routines -------- */


void mp_load_0(fft_int n, fft_int radix, fft_int out[])
{
	fft_int j;

	for (j = 0; j <= n + 1; j++) {
		out[j] = 0;
	}
}


void mp_load_1(fft_int n, fft_int radix, fft_int out[])
{
	fft_int j;

	out[0] = 1;
	out[1] = 0;
	out[2] = 1;
	for (j = 3; j <= n + 1; j++) {
		out[j] = 0;
	}
}


void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[])
{
	fft_int j;

	for (j = 0; j <= n + 1; j++) {
		out[j] = in[j];
	}
}


void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[])
{
	fft_int j, x;

	if (m < n) {
		for (j = n + 1; j > m + 2; j--) {
			inout[j] = 0;
		}
		x = 2 * inout[m + 2];
		inout[m + 2] = 0;
		if (x >= radix) {
			for (j = m + 1; j >= 2; j--) {
				x = inout[j] + 1;
				if (x < radix) {
					inout[j] = x;
					break;
				}
				inout[j] = 0;
			}
			if (x >= radix) {
				inout[2] = 1;
				inout[1]++;
			}
		}
	}
}


/* -------- mp_add routines -------- */


fft_int mp_cmp(fft_int n, fft_int radix, fft_int in1[], fft_int in2[])
{
	fft_int mp_unsgn_cmp(fft_int n, fft_int in1[], fft_int in2[]);

	if (in1[0] > in2[0]) {
		return 1;
	}
	else if (in1[0] < in2[0]) {
		return -1;
	}
	return in1[0] * mp_unsgn_cmp(n, &in1[1], &in2[1]);
}


void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[])
{
	fft_int mp_unsgn_cmp(fft_int n, fft_int in1[], fft_int in2[]);
	fft_int mp_unexp_add(fft_int n, fft_int radix, fft_int expdif,
		fft_int in1[], fft_int in2[], fft_int out[]);
	fft_int mp_unexp_sub(fft_int n, fft_int radix, fft_int expdif,
		fft_int in1[], fft_int in2[], fft_int out[]);
	fft_int outsgn, outexp, expdif;

	expdif = in1[1] - in2[1];
	outexp = in1[1];
	if (expdif < 0) {
		outexp = in2[1];
	}
	outsgn = in1[0] * in2[0];
	if (outsgn >= 0) {
		if (outsgn > 0) {
			outsgn = in1[0];
		}
		else {
			outsgn = in1[0] + in2[0];
			outexp = in1[1] + in2[1];
			expdif = 0;
		}
		if (expdif >= 0) {
			outexp += mp_unexp_add(n, radix, expdif,
				&in1[2], &in2[2], &out[2]);
		}
		else {
			outexp += mp_unexp_add(n, radix, -expdif,
				&in2[2], &in1[2], &out[2]);
		}
	}
	else {
		outsgn = mp_unsgn_cmp(n, &in1[1], &in2[1]);
		if (outsgn >= 0) {
			expdif = mp_unexp_sub(n, radix, expdif,
				&in1[2], &in2[2], &out[2]);
		}
		else {
			expdif = mp_unexp_sub(n, radix, -expdif,
				&in2[2], &in1[2], &out[2]);
		}
		outexp -= expdif;
		outsgn *= in1[0];
		if (expdif == n) {
			outsgn = 0;
		}
	}
	if (outsgn == 0) {
		outexp = 0;
	}
	out[0] = outsgn;
	out[1] = outexp;
}


void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[])
{
	fft_int mp_unsgn_cmp(fft_int n, fft_int in1[], fft_int in2[]);
	fft_int mp_unexp_add(fft_int n, fft_int radix, fft_int expdif,
		fft_int in1[], fft_int in2[], fft_int out[]);
	fft_int mp_unexp_sub(fft_int n, fft_int radix, fft_int expdif,
		fft_int in1[], fft_int in2[], fft_int out[]);
	fft_int outsgn, outexp, expdif;

	expdif = in1[1] - in2[1];
	outexp = in1[1];
	if (expdif < 0) {
		outexp = in2[1];
	}
	outsgn = in1[0] * in2[0];
	if (outsgn <= 0) {
		if (outsgn < 0) {
			outsgn = in1[0];
		}
		else {
			outsgn = in1[0] - in2[0];
			outexp = in1[1] + in2[1];
			expdif = 0;
		}
		if (expdif >= 0) {
			outexp += mp_unexp_add(n, radix, expdif,
				&in1[2], &in2[2], &out[2]);
		}
		else {
			outexp += mp_unexp_add(n, radix, -expdif,
				&in2[2], &in1[2], &out[2]);
		}
	}
	else {
		outsgn = mp_unsgn_cmp(n, &in1[1], &in2[1]);
		if (outsgn >= 0) {
			expdif = mp_unexp_sub(n, radix, expdif,
				&in1[2], &in2[2], &out[2]);
		}
		else {
			expdif = mp_unexp_sub(n, radix, -expdif,
				&in2[2], &in1[2], &out[2]);
		}
		outexp -= expdif;
		outsgn *= in1[0];
		if (expdif == n) {
			outsgn = 0;
		}
	}
	if (outsgn == 0) {
		outexp = 0;
	}
	out[0] = outsgn;
	out[1] = outexp;
}


/* -------- mp_add child routines -------- */


fft_int mp_unsgn_cmp(fft_int n, fft_int in1[], fft_int in2[])
{
	fft_int j, cmp;

	cmp = 0;
	for (j = 0; j <= n && cmp == 0; j++) {
		cmp = in1[j] - in2[j];
	}
	if (cmp > 0) {
		cmp = 1;
	}
	else if (cmp < 0) {
		cmp = -1;
	}
	return cmp;
}


fft_int mp_unexp_add(fft_int n, fft_int radix, fft_int expdif,
	fft_int in1[], fft_int in2[], fft_int out[])
{
	fft_int j, x, carry;

	carry = 0;
	if (expdif == 0 && in1[0] + in2[0] >= radix) {
		x = in1[n - 1] + in2[n - 1];
		carry = x >= radix ? -1 : 0;
		for (j = n - 1; j > 0; j--) {
			x = in1[j - 1] + in2[j - 1] - carry;
			carry = x >= radix ? -1 : 0;
			out[j] = x - (radix & carry);
		}
		out[0] = -carry;
	}
	else {
		if (expdif > n) {
			expdif = n;
		}
		for (j = n - 1; j >= expdif; j--) {
			x = in1[j] + in2[j - expdif] - carry;
			carry = x >= radix ? -1 : 0;
			out[j] = x - (radix & carry);
		}
		for (j = expdif - 1; j >= 0; j--) {
			x = in1[j] - carry;
			carry = x >= radix ? -1 : 0;
			out[j] = x - (radix & carry);
		}
		if (carry != 0) {
			for (j = n - 1; j > 0; j--) {
				out[j] = out[j - 1];
			}
			out[0] = -carry;
		}
	}
	return -carry;
}


fft_int mp_unexp_sub(fft_int n, fft_int radix, fft_int expdif,
	fft_int in1[], fft_int in2[], fft_int out[])
{
	fft_int j, x, borrow, ncancel;

	if (expdif > n) {
		expdif = n;
	}
	borrow = 0;
	for (j = n - 1; j >= expdif; j--) {
		x = in1[j] - in2[j - expdif] + borrow;
		borrow = x < 0 ? -1 : 0;
		out[j] = x + (radix & borrow);
	}
	for (j = expdif - 1; j >= 0; j--) {
		x = in1[j] + borrow;
		borrow = x < 0 ? -1 : 0;
		out[j] = x + (radix & borrow);
	}
	ncancel = 0;
	for (j = 0; j < n && out[j] == 0; j++) {
		ncancel = j + 1;
	}
	if (ncancel > 0 && ncancel < n) {
		for (j = 0; j < n - ncancel; j++) {
			out[j] = out[j + ncancel];
		}
		for (j = n - ncancel; j < n; j++) {
			out[j] = 0;
		}
	}
	return ncancel;
}


/* -------- mp_imul routines -------- */


void mp_imul(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[])
{
	void mp_unsgn_imul(fft_int n, fft_float dradix, fft_int in1[], fft_float din2,
		fft_int out[]);

	if (in2 > 0) {
		out[0] = in1[0];
	}
	else if (in2 < 0) {
		out[0] = -in1[0];
		in2 = -in2;
	}
	else {
		out[0] = 0;
	}
	mp_unsgn_imul(n, radix, &in1[1], in2, &out[1]);
	if (out[0] == 0) {
		out[1] = 0;
	}
}


fft_int mp_idiv(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[])
{
	void mp_load_0(fft_int n, fft_int radix, fft_int out[]);
	void mp_unsgn_idiv(fft_int n, fft_float dradix, fft_int in1[], fft_float din2,
		fft_int out[]);

	if (in2 == 0) {
		return -1;
	}
	if (in2 > 0) {
		out[0] = in1[0];
	}
	else {
		out[0] = -in1[0];
		in2 = -in2;
	}
	if (in1[0] == 0) {
		mp_load_0(n, radix, out);
		return 0;
	}
	mp_unsgn_idiv(n, radix, &in1[1], in2, &out[1]);
	return 0;
}


void mp_idiv_2(fft_int n, fft_int radix, fft_int in[], fft_int out[])
{
	fft_int j, ix, carry, shift;

	out[0] = in[0];
	shift = 0;
	if (in[2] == 1) {
		shift = 1;
	}
	out[1] = in[1] - shift;
	carry = -shift;
	for (j = 2; j <= n + 1 - shift; j++) {
		ix = in[j + shift] + (radix & carry);
		carry = -(ix & 1);
		out[j] = ix >> 1;
	}
	if (shift > 0) {
		out[n + 1] = (radix & carry) >> 1;
	}
}


/* -------- mp_imul child routines -------- */


void mp_unsgn_imul(fft_int n, fft_float dradix, fft_int in1[], fft_float din2,
	fft_int out[])
{
	fft_int j, carry, shift;
	fft_float x, d1_radix;

	d1_radix = ((fft_float)1) / dradix;
	carry = 0;
	for (j = n; j >= 1; j--) {
		x = din2 * in1[j] + carry + FC_HALF;
		carry = (fft_int)(d1_radix * x);
		out[j] = (fft_int)(x - dradix * carry);
	}
	shift = 0;
	x = carry + FC_HALF;
	while (x > 1) {
		x *= d1_radix;
		shift++;
	}
	out[0] = in1[0] + shift;
	if (shift > 0) {
		while (shift > n) {
			carry = (fft_int)(d1_radix * carry + FC_HALF);
			shift--;
		}
		for (j = n; j >= shift + 1; j--) {
			out[j] = out[j - shift];
		}
		for (j = shift; j >= 1; j--) {
			x = carry + FC_HALF;
			carry = (fft_int)(d1_radix * x);
			out[j] = (fft_int)(x - dradix * carry);
		}
	}
}


void mp_unsgn_idiv(fft_int n, fft_float dradix, fft_int in1[], fft_float din2,
	fft_int out[])
{
	fft_int j, ix, carry, shift;
	fft_float x, d1_in2;

	d1_in2 = ((fft_float)1) / din2;
	shift = 0;
	x = 0;
	do {
		shift++;
		x *= dradix;
		if (shift <= n) {
			x += in1[shift];
		}
	} while (x < din2 - FC_HALF);
	x += FC_HALF;
	ix = (fft_int)(d1_in2 * x);
	carry = (fft_int)(x - din2 * ix);
	out[1] = ix;
	shift--;
	out[0] = in1[0] - shift;
	if (shift >= n) {
		shift = n - 1;
	}
	for (j = 2; j <= n - shift; j++) {
		x = in1[j + shift] + dradix * carry + FC_HALF;
		ix = (fft_int)(d1_in2 * x);
		carry = (fft_int)(x - din2 * ix);
		out[j] = ix;
	}
	for (j = n - shift + 1; j <= n; j++) {
		x = dradix * carry + FC_HALF;
		ix = (fft_int)(d1_in2 * x);
		carry = (fft_int)(x - din2 * ix);
		out[j] = ix;
	}
}


/* -------- mp_mul routines -------- */


fft_float mp_mul_radix_test(fft_int n, fft_int radix, fft_int nfft,
	fft_float tmpfft[], fft_int ip[], fft_float w[])
{
	void mp_mul_csqu(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[]);
	fft_float mp_mul_d2i_test(fft_int radix, fft_int nfft, fft_float din[]);
	fft_int j, ndata, radix_2;

	ndata = (nfft >> 1) + 1;
	if (ndata > n) {
		ndata = n;
	}
	tmpfft[nfft + 1] = radix - 1;
	for (j = nfft; j > ndata; j--) {
		tmpfft[j] = 0;
	}
	radix_2 = (radix + 1) / 2;
	for (j = ndata; j > 2; j--) {
		tmpfft[j] = radix_2;
	}
	tmpfft[2] = radix;
	tmpfft[1] = radix - 1;
	tmpfft[0] = 0;
	mp_mul_csqu(nfft, tmpfft, ip, w);
	return 2 * mp_mul_d2i_test(radix, nfft, tmpfft);
}


void mp_mul(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int tmp[], fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
	fft_float tmp3fft[], fft_int ip[], fft_float w[])
{
	void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul_nt_out(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_cmul_nt_d2(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_cmul_nt_d1_add(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_float d3[], fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h, shift;

	shift = (nfft >> 1) + 1;
	while (n > shift) {
		if (in1[shift + 2] + in2[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp3fft = (upper) in1 * (lower) in2 ---- */
	mp_mul_i2d(n, radix, nfft, 0, in1, tmp1fft);
	mp_mul_i2d(n, radix, nfft, shift, in2, tmp3fft);
	mp_mul_cmul_nt_out(nfft, tmp1fft, tmp3fft, ip, w);
	/* ---- tmp = (upper) in1 * (upper) in2 ---- */
	mp_mul_i2d(n, radix, nfft, 0, in2, tmp2fft);
	mp_mul_cmul_nt_d2(nfft, tmp2fft, tmp1fft, ip, w);
	mp_mul_d2i(n, radix, nfft, tmp1fft, tmp);
	/* ---- tmp3fft += (upper) in2 * (lower) in1 ---- */
	mp_mul_i2d(n, radix, nfft, shift, in1, tmp1fft);
	mp_mul_cmul_nt_d1_add(nfft, tmp2fft, tmp1fft, tmp3fft, ip, w);
	/* ---- out = tmp + tmp3fft ---- */
	mp_mul_d2i(n_h, radix, nfft, tmp3fft, out);
	if (out[0] != 0) {
		mp_add(n, radix, out, tmp, out);
	}
	else {
		mp_copy(n, radix, tmp, out);
	}
}


void mp_squ(fft_int n, fft_int radix, fft_int in[], fft_int out[], fft_int tmp[],
	fft_int nfft, fft_float tmp1fft[], fft_float tmp2fft[],
	fft_int ip[], fft_float w[])
{
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_csqu_nt_d1(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h, shift;

	shift = (nfft >> 1) + 1;
	while (n > shift) {
		if (in[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp = (upper) in * (lower) in ---- */
	mp_mul_i2d(n, radix, nfft, 0, in, tmp1fft);
	mp_mul_i2d(n, radix, nfft, shift, in, tmp2fft);
	mp_mul_cmul(nfft, tmp1fft, tmp2fft, ip, w);
	mp_mul_d2i(n_h, radix, nfft, tmp2fft, tmp);
	/* ---- out = 2 * tmp + ((upper) in)^2 ---- */
	mp_mul_csqu_nt_d1(nfft, tmp1fft, ip, w);
	mp_mul_d2i(n, radix, nfft, tmp1fft, out);
	if (tmp[0] != 0) {
		mp_add(n_h, radix, tmp, tmp, tmp);
		mp_add(n, radix, out, tmp, out);
	}
}


void mp_mulhf(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int tmp[], fft_int nfft, fft_float in1fft[], fft_float tmpfft[],
	fft_int ip[], fft_float w[])
{
	void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_cmul_nt_d1(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h, shift;

	shift = (nfft >> 1) + 1;
	while (n > shift) {
		if (in2[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp = (upper) in1 * (upper) in2 ---- */
	mp_mul_i2d(n, radix, nfft, 0, in1, in1fft);
	mp_mul_i2d(n, radix, nfft, 0, in2, tmpfft);
	mp_mul_cmul(nfft, in1fft, tmpfft, ip, w);
	mp_mul_d2i(n, radix, nfft, tmpfft, tmp);
	/* ---- out = tmp + (upper) in1 * (lower) in2 ---- */
	mp_mul_i2d(n, radix, nfft, shift, in2, tmpfft);
	mp_mul_cmul_nt_d1(nfft, in1fft, tmpfft, ip, w);
	mp_mul_d2i(n_h, radix, nfft, tmpfft, out);
	if (out[0] != 0) {
		mp_add(n, radix, out, tmp, out);
	}
	else {
		mp_copy(n, radix, tmp, out);
	}
}


void mp_mulhf_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[], fft_int in2[],
	fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
	fft_int ip[], fft_float w[])
{
	void mp_copy(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul_nt_d1(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h, shift;

	shift = (nfft >> 1) + 1;
	while (n > shift) {
		if (in2[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp = (upper) in1fft * (upper) in2 ---- */
	mp_mul_i2d(n, radix, nfft, 0, in2, tmpfft);
	mp_mul_cmul_nt_d1(nfft, in1fft, tmpfft, ip, w);
	mp_mul_d2i(n, radix, nfft, tmpfft, tmp);
	/* ---- out = tmp + (upper) in1 * (lower) in2 ---- */
	mp_mul_i2d(n, radix, nfft, shift, in2, tmpfft);
	mp_mul_cmul_nt_d1(nfft, in1fft, tmpfft, ip, w);
	mp_mul_d2i(n_h, radix, nfft, tmpfft, out);
	if (out[0] != 0) {
		mp_add(n, radix, out, tmp, out);
	}
	else {
		mp_copy(n, radix, tmp, out);
	}
}


void mp_squhf_use_infft(fft_int n, fft_int radix, fft_float infft[], fft_int in[],
	fft_int out[], fft_int tmp[], fft_int nfft, fft_float tmpfft[],
	fft_int ip[], fft_float w[])
{
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul_nt_d1(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_csqu_nt_d1(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h, shift;

	shift = (nfft >> 1) + 1;
	while (n > shift) {
		if (in[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp = (upper) infft * (lower) in ---- */
	mp_mul_i2d(n, radix, nfft, shift, in, tmpfft);
	mp_mul_cmul_nt_d1(nfft, infft, tmpfft, ip, w);
	mp_mul_d2i(n_h, radix, nfft, tmpfft, tmp);
	/* ---- out = tmp + ((upper) infft)^2 ---- */
	mp_mul_csqu_nt_d1(nfft, infft, ip, w);
	mp_mul_d2i(n, radix, nfft, infft, out);
	if (tmp[0] != 0) {
		mp_add(n, radix, out, tmp, out);
	}
}


void mp_mulh(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
	fft_int nfft, fft_float in1fft[], fft_float outfft[], fft_int ip[], fft_float w[])
{
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);

	mp_mul_i2d(n, radix, nfft, 0, in1, in1fft);
	mp_mul_i2d(n, radix, nfft, 0, in2, outfft);
	mp_mul_cmul(nfft, in1fft, outfft, ip, w);
	mp_mul_d2i(n, radix, nfft, outfft, out);
}


void mp_mulh_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[],
	fft_int shift, fft_int in2[], fft_int out[], fft_int nfft, fft_float outfft[],
	fft_int ip[], fft_float w[])
{
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_cmul_nt_d1(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);
	fft_int n_h;

	while (n > shift) {
		if (in2[shift + 2] != 0) {
			break;
		}
		shift++;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	mp_mul_i2d(n, radix, nfft, shift, in2, outfft);
	mp_mul_cmul_nt_d1(nfft, in1fft, outfft, ip, w);
	mp_mul_d2i(n_h, radix, nfft, outfft, out);
}


void mp_squh(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int nfft, fft_float outfft[], fft_int ip[], fft_float w[])
{
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_csqu(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);

	mp_mul_i2d(n, radix, nfft, 0, in, outfft);
	mp_mul_csqu(nfft, outfft, ip, w);
	mp_mul_d2i(n, radix, nfft, outfft, out);
}


void mp_squh_save_infft(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int nfft, fft_float infft[], fft_float outfft[],
	fft_int ip[], fft_float w[])
{
	void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
		fft_int in[], fft_float dout[]);
	void mp_mul_csqu_save_d1(fft_int nfft, fft_float d1[], fft_float d2[],
		fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);

	mp_mul_i2d(n, radix, nfft, 0, in, infft);
	mp_mul_csqu_save_d1(nfft, infft, outfft, ip, w);
	mp_mul_d2i(n, radix, nfft, outfft, out);
}


void mp_squh_use_in1fft(fft_int n, fft_int radix, fft_float inoutfft[], fft_int out[],
	fft_int nfft, fft_int ip[], fft_float w[])
{
	void mp_mul_csqu_nt_d1(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[]);
	void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[]);

	mp_mul_csqu_nt_d1(nfft, inoutfft, ip, w);
	mp_mul_d2i(n, radix, nfft, inoutfft, out);
}


/* -------- mp_mul child routines -------- */


void mp_mul_i2d(fft_int n, fft_int radix, fft_int nfft, fft_int shift,
	fft_int in[], fft_float dout[])
{
	fft_int j, x, carry, ndata, radix_2, topdgt;

	ndata = 0;
	topdgt = 0;
	if (n > shift) {
		topdgt = in[shift + 2];
		ndata = (nfft >> 1) + 1;
		if (ndata > n - shift) {
			ndata = n - shift;
		}
	}
	dout[nfft + 1] = in[0] * topdgt;
	for (j = nfft; j > ndata; j--) {
		dout[j] = 0;
	}
	/* ---- abs(dout[j]) <= radix/2 (to keep FFT precision) ---- */
	if (ndata > 1) {
		radix_2 = radix / 2;
		carry = 0;
		for (j = ndata + 1; j > 3; j--) {
			x = in[j + shift] - carry;
			carry = x >= radix_2 ? -1 : 0;
			dout[j - 1] = x - (radix & carry);
		}
		dout[2] = in[shift + 3] - carry;
	}
	dout[1] = topdgt;
	dout[0] = in[1] - shift;
}


void mp_mul_cmul(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_int ip[], fft_float w[])
{
	void makect(fft_int nc, fft_int * ip, fft_float * c); /* in fft*g.c */
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcmul(fft_int n, fft_float * a, fft_float * b, fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d1[1], ip, w);
	cdft(nfft, 1, &d2[1], ip, w);
	if (nfft > (ip[1] << 2)) {
		makect(nfft >> 2, ip, &w[ip[0]]);
	}
	d2[0] += d1[0];
	xr = d1[1] * d2[1] + d1[2] * d2[2];
	xi = d1[1] * d2[2] + d1[2] * d2[1];
	d2[1] = xr;
	d2[2] = xi;
	if (nfft > 2) {
		mp_mul_rcmul(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
	}
	d2[nfft + 1] *= d1[nfft + 1];
	cdft(nfft, -1, &d2[1], ip, w);
}


void mp_mul_cmul_nt_d1(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_int ip[], fft_float w[])
{
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcmul_nt_in1(fft_int n, fft_float * a, fft_float * b, fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d2[1], ip, w);
	d2[0] += d1[0];
	xr = d1[1] * d2[1] + d1[2] * d2[2];
	xi = d1[1] * d2[2] + d1[2] * d2[1];
	d2[1] = xr;
	d2[2] = xi;
	if (nfft > 2) {
		mp_mul_rcmul_nt_in1(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
	}
	d2[nfft + 1] *= d1[nfft + 1];
	cdft(nfft, -1, &d2[1], ip, w);
}


void mp_mul_cmul_nt_d2(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_int ip[], fft_float w[])
{
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcmul_nt_in2(fft_int n, fft_float * a, fft_float * b, fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d1[1], ip, w);
	d2[0] += d1[0];
	xr = d1[1] * d2[1] + d1[2] * d2[2];
	xi = d1[1] * d2[2] + d1[2] * d2[1];
	d2[1] = xr;
	d2[2] = xi;
	if (nfft > 2) {
		mp_mul_rcmul_nt_in2(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
	}
	d2[nfft + 1] *= d1[nfft + 1];
	cdft(nfft, -1, &d2[1], ip, w);
}


void mp_mul_cmul_nt_out(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_int ip[], fft_float w[])
{
	void makect(fft_int nc, fft_int * ip, fft_float * c); /* in fft*g.c */
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcmul_nt_out(fft_int n, fft_float * a, fft_float * b,
		fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d1[1], ip, w);
	cdft(nfft, 1, &d2[1], ip, w);
	if (nfft > (ip[1] << 2)) {
		makect(nfft >> 2, ip, &w[ip[0]]);
	}
	d2[0] += d1[0];
	xr = d1[1] * d2[1] + d1[2] * d2[2];
	xi = d1[1] * d2[2] + d1[2] * d2[1];
	d2[1] = xr;
	d2[2] = xi;
	if (nfft > 2) {
		mp_mul_rcmul_nt_out(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
	}
	d2[nfft + 1] *= d1[nfft + 1];
}


void mp_mul_cmul_nt_d1_add(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_float d3[], fft_int ip[], fft_float w[])
{
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcmul_nt_in1_add(fft_int n, fft_float * a, fft_float * b, fft_float * badd,
		fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d2[1], ip, w);
	xr = d1[1] * d2[1] + d1[2] * d2[2];
	xi = d1[1] * d2[2] + d1[2] * d2[1];
	d3[1] += xr;
	d3[2] += xi;
	if (nfft > 2) {
		mp_mul_rcmul_nt_in1_add(nfft, &d1[1], &d2[1], &d3[1],
			ip[1], &w[ip[0]]);
	}
	d3[nfft + 1] += d1[nfft + 1] * d2[nfft + 1];
	cdft(nfft, -1, &d3[1], ip, w);
}


void mp_mul_csqu(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[])
{
	void makect(fft_int nc, fft_int * ip, fft_float * c); /* in fft*g.c */
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcsqu(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d1[1], ip, w);
	if (nfft > (ip[1] << 2)) {
		makect(nfft >> 2, ip, &w[ip[0]]);
	}
	d1[0] *= 2;
	xr = d1[1] * d1[1] + d1[2] * d1[2];
	xi = 2 * d1[1] * d1[2];
	d1[1] = xr;
	d1[2] = xi;
	if (nfft > 2) {
		mp_mul_rcsqu(nfft, &d1[1], ip[1], &w[ip[0]]);
	}
	d1[nfft + 1] *= d1[nfft + 1];
	cdft(nfft, -1, &d1[1], ip, w);
}


void mp_mul_csqu_save_d1(fft_int nfft, fft_float d1[], fft_float d2[],
	fft_int ip[], fft_float w[])
{
	void makect(fft_int nc, fft_int * ip, fft_float * c); /* in fft*g.c */
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcsqu_save(fft_int n, fft_float * a, fft_float * b, fft_int nc, fft_float * c);
	fft_float xr, xi;

	cdft(nfft, 1, &d1[1], ip, w);
	if (nfft > (ip[1] << 2)) {
		makect(nfft >> 2, ip, &w[ip[0]]);
	}
	d2[0] = 2 * d1[0];
	xr = d1[1] * d1[1] + d1[2] * d1[2];
	xi = 2 * d1[1] * d1[2];
	d2[1] = xr;
	d2[2] = xi;
	if (nfft > 2) {
		mp_mul_rcsqu_save(nfft, &d1[1], &d2[1], ip[1], &w[ip[0]]);
	}
	d2[nfft + 1] = d1[nfft + 1] * d1[nfft + 1];
	cdft(nfft, -1, &d2[1], ip, w);
}


void mp_mul_csqu_nt_d1(fft_int nfft, fft_float d1[], fft_int ip[], fft_float w[])
{
	void cdft(fft_int n, fft_int isgn, fft_float * a, fft_int * ip, fft_float * w);
	void mp_mul_rcsqu_nt_in(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_float xr, xi;

	d1[0] *= 2;
	xr = d1[1] * d1[1] + d1[2] * d1[2];
	xi = 2 * d1[1] * d1[2];
	d1[1] = xr;
	d1[2] = xi;
	if (nfft > 2) {
		mp_mul_rcsqu_nt_in(nfft, &d1[1], ip[1], &w[ip[0]]);
	}
	d1[nfft + 1] *= d1[nfft + 1];
	cdft(nfft, -1, &d1[1], ip, w);
}


void mp_mul_d2i(fft_int n, fft_int radix, fft_int nfft, fft_float din[], fft_int out[])
{
	fft_int j, carry, carry1, carry2, shift, ndata;
	fft_float x, scale, d1_radix, d1_radix2, pow_radix, topdgt;

	scale = ((fft_float)2) / nfft;
	d1_radix = ((fft_float)1) / radix;
	d1_radix2 = d1_radix * d1_radix;
	topdgt = din[nfft + 1];
	x = topdgt < 0 ? -topdgt : topdgt;
	shift = x + FC_HALF >= radix ? 1 : 0;
	/* ---- correction of cyclic convolution of din[1] ---- */
	x *= nfft * FC_HALF;
	din[nfft + 1] = din[1] - x;
	din[1] = x;
	/* ---- output of digits ---- */
	ndata = n;
	if (n > nfft + 1 + shift) {
		ndata = nfft + 1 + shift;
		for (j = n + 1; j > ndata + 1; j--) {
			out[j] = 0;
		}
	}
	x = 0;
	pow_radix = 1;
	for (j = ndata + 1 - shift; j <= nfft + 1; j++) {
		x += pow_radix * din[j];
		pow_radix *= d1_radix;
		if (pow_radix < FFT_FLOAT_EPSILON) {
			break;
		}
	}
	x = d1_radix2 * (scale * x + FC_HALF);
	carry2 = ((fft_int)x) - 1;
	carry = (fft_int)(radix * (x - carry2) + FC_HALF);
	for (j = ndata; j > 1; j--) {
		x = d1_radix2 * (scale * din[j - shift] + carry + FC_HALF);
		carry = carry2;
		carry2 = ((fft_int)x) - 1;
		x = radix * (x - carry2);
		carry1 = (fft_int)x;
		out[j + 1] = (fft_int)(radix * (x - carry1));
		carry += carry1;
	}
	x = carry + ((fft_float)radix) * carry2 + FC_HALF;
	if (shift == 0) {
		x += scale * din[1];
	}
	carry = (fft_int)(d1_radix * x);
	out[2] = (fft_int)(x - ((fft_float)radix) * carry);
	if (carry > 0) {
		for (j = n + 1; j > 2; j--) {
			out[j] = out[j - 1];
		}
		out[2] = carry;
		shift++;
	}
	/* ---- output of exp, sgn ---- */
	x = din[0] + shift + FC_HALF;
	shift = ((fft_int)x) - 1;
	out[1] = shift + ((fft_int)(x - shift));
	out[0] = topdgt > FC_HALF ? 1 : -1;
	if (out[2] == 0) {
		out[0] = 0;
		out[1] = 0;
	}
}


fft_float mp_mul_d2i_test(fft_int radix, fft_int nfft, fft_float din[])
{
	fft_int j, carry, carry1, carry2;
	fft_float x, scale, d1_radix, d1_radix2, err;

	scale = ((fft_float)2) / nfft;
	d1_radix = ((fft_float)1) / radix;
	d1_radix2 = d1_radix * d1_radix;
	/* ---- correction of cyclic convolution of din[1] ---- */
	x = din[nfft + 1] * nfft * FC_HALF;
	if (x < 0) {
		x = -x;
	}
	din[nfft + 1] = din[1] - x;
	/* ---- check of digits ---- */
	err = 0;
	carry = 0;
	carry2 = 0;
	for (j = nfft + 1; j > 1; j--) {
		x = d1_radix2 * (scale * din[j] + carry + FC_HALF);
		carry = carry2;
		carry2 = ((fft_int)x) - 1;
		x = radix * (x - carry2);
		carry1 = (fft_int)x;
		x = radix * (x - carry1);
		carry += carry1;
		x = x - FC_HALF - ((fft_int)x);
		if (x > err) {
			err = x;
		}
		else if (-x > err) {
			err = -x;
		}
	}
	return err;
}


/* -------- mp_mul child^2 routines (mix RFFT routines) -------- */


void mp_mul_rcmul(fft_int n, fft_float* a, fft_float* b, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, ajr, aji, akr, aki, bjr, bji, bkr, bki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data a[] into RFFT data ---- */
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		ajr = a[j] - yr;
		aji = a[j + 1] - yi;
		akr = a[k] + yr;
		aki = a[k + 1] - yi;
		a[j] = ajr;
		a[j + 1] = aji;
		a[k] = akr;
		a[k + 1] = aki;
		/* ---- transform CFFT data b[] into RFFT data ---- */
		xr = b[j] - b[k];
		xi = b[j + 1] + b[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = b[j] - yr;
		xi = b[j + 1] - yi;
		yr = b[k] + yr;
		yi = b[k + 1] - yi;
		/* ---- cmul ---- */
		bjr = ajr * xr - aji * xi;
		bji = ajr * xi + aji * xr;
		bkr = akr * yr - aki * yi;
		bki = akr * yi + aki * yr;
		/* ---- transform RFFT data bxx into CFFT data ---- */
		xr = bjr - bkr;
		xi = bji + bki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		b[j] = bjr - yr;
		b[j + 1] = bji - yi;
		b[k] = bkr + yr;
		b[k + 1] = bki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	yr = b[m];
	yi = b[m + 1];
	b[m] = xr * yr - xi * yi;
	b[m + 1] = xr * yi + xi * yr;
}


void mp_mul_rcmul_nt_in1(fft_int n, fft_float* a, fft_float* b, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, bjr, bji, bkr, bki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data b[] into RFFT data ---- */
		xr = b[j] - b[k];
		xi = b[j + 1] + b[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = b[j] - yr;
		xi = b[j + 1] - yi;
		yr = b[k] + yr;
		yi = b[k + 1] - yi;
		/* ---- cmul ---- */
		bjr = a[j] * xr - a[j + 1] * xi;
		bji = a[j] * xi + a[j + 1] * xr;
		bkr = a[k] * yr - a[k + 1] * yi;
		bki = a[k] * yi + a[k + 1] * yr;
		/* ---- transform RFFT data bxx into CFFT data ---- */
		xr = bjr - bkr;
		xi = bji + bki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		b[j] = bjr - yr;
		b[j + 1] = bji - yi;
		b[k] = bkr + yr;
		b[k + 1] = bki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	yr = b[m];
	yi = b[m + 1];
	b[m] = xr * yr - xi * yi;
	b[m + 1] = xr * yi + xi * yr;
}


void mp_mul_rcmul_nt_in2(fft_int n, fft_float* a, fft_float* b, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, bjr, bji, bkr, bki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data a[] into RFFT data ---- */
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = a[j] - yr;
		xi = a[j + 1] - yi;
		yr = a[k] + yr;
		yi = a[k + 1] - yi;
		a[j] = xr;
		a[j + 1] = xi;
		a[k] = yr;
		a[k + 1] = yi;
		/* ---- cmul ---- */
		bjr = b[j] * xr - b[j + 1] * xi;
		bji = b[j] * xi + b[j + 1] * xr;
		bkr = b[k] * yr - b[k + 1] * yi;
		bki = b[k] * yi + b[k + 1] * yr;
		/* ---- transform RFFT data bxx into CFFT data ---- */
		xr = bjr - bkr;
		xi = bji + bki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		b[j] = bjr - yr;
		b[j + 1] = bji - yi;
		b[k] = bkr + yr;
		b[k + 1] = bki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	yr = b[m];
	yi = b[m + 1];
	b[m] = xr * yr - xi * yi;
	b[m + 1] = xr * yi + xi * yr;
}


void mp_mul_rcmul_nt_out(fft_int n, fft_float* a, fft_float* b, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, ajr, aji, akr, aki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data a[] into RFFT data ---- */
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		ajr = a[j] - yr;
		aji = a[j + 1] - yi;
		akr = a[k] + yr;
		aki = a[k + 1] - yi;
		a[j] = ajr;
		a[j + 1] = aji;
		a[k] = akr;
		a[k + 1] = aki;
		/* ---- transform CFFT data b[] into RFFT data ---- */
		xr = b[j] - b[k];
		xi = b[j + 1] + b[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = b[j] - yr;
		xi = b[j + 1] - yi;
		yr = b[k] + yr;
		yi = b[k + 1] - yi;
		/* ---- cmul ---- */
		b[j] = ajr * xr - aji * xi;
		b[j + 1] = ajr * xi + aji * xr;
		b[k] = akr * yr - aki * yi;
		b[k + 1] = akr * yi + aki * yr;
	}
	xr = a[m];
	xi = a[m + 1];
	yr = b[m];
	yi = b[m + 1];
	b[m] = xr * yr - xi * yi;
	b[m + 1] = xr * yi + xi * yr;
}


void mp_mul_rcmul_nt_in1_add(fft_int n, fft_float* a, fft_float* b, fft_float* badd,
	fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, bjr, bji, bkr, bki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data b[] into RFFT data ---- */
		xr = b[j] - b[k];
		xi = b[j + 1] + b[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = b[j] - yr;
		xi = b[j + 1] - yi;
		yr = b[k] + yr;
		yi = b[k + 1] - yi;
		/* ---- cmul + add ---- */
		bjr = badd[j] + (a[j] * xr - a[j + 1] * xi);
		bji = badd[j + 1] + (a[j] * xi + a[j + 1] * xr);
		bkr = badd[k] + (a[k] * yr - a[k + 1] * yi);
		bki = badd[k + 1] + (a[k] * yi + a[k + 1] * yr);
		/* ---- transform RFFT data bxx into CFFT data ---- */
		xr = bjr - bkr;
		xi = bji + bki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		badd[j] = bjr - yr;
		badd[j + 1] = bji - yi;
		badd[k] = bkr + yr;
		badd[k + 1] = bki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	yr = b[m];
	yi = b[m + 1];
	badd[m] += xr * yr - xi * yi;
	badd[m + 1] += xr * yi + xi * yr;
}


void mp_mul_rcsqu(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, ajr, aji, akr, aki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data a[] into RFFT data ---- */
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = a[j] - yr;
		xi = a[j + 1] - yi;
		yr = a[k] + yr;
		yi = a[k + 1] - yi;
		/* ---- csqu ---- */
		ajr = xr * xr - xi * xi;
		aji = 2 * xr * xi;
		akr = yr * yr - yi * yi;
		aki = 2 * yr * yi;
		/* ---- transform RFFT data axx into CFFT data ---- */
		xr = ajr - akr;
		xi = aji + aki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		a[j] = ajr - yr;
		a[j + 1] = aji - yi;
		a[k] = akr + yr;
		a[k + 1] = aki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	a[m] = xr * xr - xi * xi;
	a[m + 1] = 2 * xr * xi;
}


void mp_mul_rcsqu_save(fft_int n, fft_float* a, fft_float* b, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, ajr, aji, akr, aki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- transform CFFT data a[] into RFFT data ---- */
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		xr = a[j] - yr;
		xi = a[j + 1] - yi;
		yr = a[k] + yr;
		yi = a[k + 1] - yi;
		a[j] = xr;
		a[j + 1] = xi;
		a[k] = yr;
		a[k + 1] = yi;
		/* ---- csqu ---- */
		ajr = xr * xr - xi * xi;
		aji = 2 * xr * xi;
		akr = yr * yr - yi * yi;
		aki = 2 * yr * yi;
		/* ---- transform RFFT data axx into CFFT data ---- */
		xr = ajr - akr;
		xi = aji + aki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		b[j] = ajr - yr;
		b[j + 1] = aji - yi;
		b[k] = akr + yr;
		b[k + 1] = aki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	b[m] = xr * xr - xi * xi;
	b[m + 1] = 2 * xr * xi;
}


void mp_mul_rcsqu_nt_in(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi, ajr, aji, akr, aki;

	ks = (nc << 2) / n;
	kk = 0;
	m = n >> 1;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		/* ---- csqu ---- */
		xr = a[j];
		xi = a[j + 1];
		yr = a[k];
		yi = a[k + 1];
		ajr = xr * xr - xi * xi;
		aji = 2 * xr * xi;
		akr = yr * yr - yi * yi;
		aki = 2 * yr * yi;
		/* ---- transform RFFT data axx into CFFT data ---- */
		xr = ajr - akr;
		xi = aji + aki;
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		a[j] = ajr - yr;
		a[j + 1] = aji - yi;
		a[k] = akr + yr;
		a[k + 1] = aki - yi;
	}
	xr = a[m];
	xi = a[m + 1];
	a[m] = xr * xr - xi * xi;
	a[m + 1] = 2 * xr * xi;
}


/* -------- mp_inv routines -------- */


fft_int mp_inv(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[])
{
	fft_int mp_get_nfft_init(fft_int radix, fft_int nfft_max);
	void mp_inv_init(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	fft_int mp_inv_newton(fft_int n, fft_int radix, fft_int in[], fft_int inout[],
		fft_int tmp1[], fft_int tmp2[], fft_int nfft, fft_float tmp1fft[],
		fft_float tmp2fft[], fft_int ip[], fft_float w[]);
	fft_int n_nwt, nfft_nwt, thr, prc;

	if (in[0] == 0) {
		return -1;
	}
	nfft_nwt = mp_get_nfft_init(radix, nfft);
	n_nwt = nfft_nwt + 2;
	if (n_nwt > n) {
		n_nwt = n;
	}
	mp_inv_init(n_nwt, radix, in, out);
	thr = 8;
	do {
		n_nwt = nfft_nwt + 2;
		if (n_nwt > n) {
			n_nwt = n;
		}
		prc = mp_inv_newton(n_nwt, radix, in, out,
			tmp1, tmp2, nfft_nwt, tmp1fft, tmp2fft, ip, w);
		if (thr * nfft_nwt >= nfft) {
			thr = 0;
			if (2 * prc <= n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		else {
			if (3 * prc < n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		nfft_nwt <<= 1;
	} while (nfft_nwt <= nfft);
	return 0;
}


fft_int mp_sqrt(fft_int n, fft_int radix, fft_int in[], fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[])
{
	void mp_load_0(fft_int n, fft_int radix, fft_int out[]);
	fft_int mp_get_nfft_init(fft_int radix, fft_int nfft_max);
	void mp_sqrt_init(fft_int n, fft_int radix, fft_int in[], fft_int out[], fft_int out_rev[]);
	fft_int mp_sqrt_newton(fft_int n, fft_int radix, fft_int in[], fft_int inout[],
		fft_int inout_rev[], fft_int tmp[], fft_int nfft, fft_float tmp1fft[],
		fft_float tmp2fft[], fft_int ip[], fft_float w[], fft_int * n_tmp1fft);
	fft_int n_nwt, nfft_nwt, thr, prc, n_tmp1fft;

	if (in[0] < 0) {
		return -1;
	}
	else if (in[0] == 0) {
		mp_load_0(n, radix, out);
		return 0;
	}
	nfft_nwt = mp_get_nfft_init(radix, nfft);
	n_nwt = nfft_nwt + 2;
	if (n_nwt > n) {
		n_nwt = n;
	}
	mp_sqrt_init(n_nwt, radix, in, out, tmp1);
	n_tmp1fft = 0;
	thr = 8;
	do {
		n_nwt = nfft_nwt + 2;
		if (n_nwt > n) {
			n_nwt = n;
		}
		prc = mp_sqrt_newton(n_nwt, radix, in, out,
			tmp1, tmp2, nfft_nwt, tmp1fft, tmp2fft,
			ip, w, &n_tmp1fft);
		if (thr * nfft_nwt >= nfft) {
			thr = 0;
			if (2 * prc <= n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		else {
			if (3 * prc < n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		nfft_nwt <<= 1;
	} while (nfft_nwt <= nfft);
	return 0;
}


fft_int mp_invisqrt(fft_int n, fft_int radix, fft_int in, fft_int out[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft,
	fft_float tmp1fft[], fft_float tmp2fft[], fft_int ip[], fft_float w[])
{
	fft_int mp_get_nfft_init(fft_int radix, fft_int nfft_max);
	void mp_invisqrt_init(fft_int n, fft_int radix, fft_int in, fft_int out[]);
	fft_int mp_invisqrt_newton(fft_int n, fft_int radix, fft_int in, fft_int inout[],
		fft_int tmp1[], fft_int tmp2[], fft_int nfft, fft_float tmp1fft[],
		fft_float tmp2fft[], fft_int ip[], fft_float w[]);
	fft_int n_nwt, nfft_nwt, thr, prc;

	if (in <= 0) {
		return -1;
	}
	nfft_nwt = mp_get_nfft_init(radix, nfft);
	n_nwt = nfft_nwt + 2;
	if (n_nwt > n) {
		n_nwt = n;
	}
	mp_invisqrt_init(n_nwt, radix, in, out);
	thr = 8;
	do {
		n_nwt = nfft_nwt + 2;
		if (n_nwt > n) {
			n_nwt = n;
		}
		prc = mp_invisqrt_newton(n_nwt, radix, in, out,
			tmp1, tmp2, nfft_nwt, tmp1fft, tmp2fft, ip, w);
		if (thr * nfft_nwt >= nfft) {
			thr = 0;
			if (2 * prc <= n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		else {
			if (3 * prc < n_nwt - 2) {
				nfft_nwt >>= 1;
			}
		}
		nfft_nwt <<= 1;
	} while (nfft_nwt <= nfft);
	return 0;
}


/* -------- mp_inv child routines -------- */


fft_int mp_get_nfft_init(fft_int radix, fft_int nfft_max)
{
	fft_int nfft_init;
	fft_float r;

	r = radix;
	nfft_init = 1;
	do {
		r *= r;
		nfft_init <<= 1;
	} while (FFT_FLOAT_EPSILON * r < 1 && nfft_init < nfft_max);
	return nfft_init;
}


void mp_inv_init(fft_int n, fft_int radix, fft_int in[], fft_int out[])
{
	void mp_unexp_d2mp(fft_int n, fft_int radix, fft_float din, fft_int out[]);
	fft_float mp_unexp_mp2d(fft_int n, fft_int radix, fft_int in[]);
	fft_int outexp;
	fft_float din;

	out[0] = in[0];
	outexp = -in[1];
	din = ((fft_float)1) / mp_unexp_mp2d(n, radix, &in[2]);
	while (din < 1) {
		din *= radix;
		outexp--;
	}
	out[1] = outexp;
	mp_unexp_d2mp(n, radix, din, &out[2]);
}


void mp_sqrt_init(fft_int n, fft_int radix, fft_int in[], fft_int out[], fft_int out_rev[])
{
	void mp_unexp_d2mp(fft_int n, fft_int radix, fft_float din, fft_int out[]);
	fft_float mp_unexp_mp2d(fft_int n, fft_int radix, fft_int in[]);
	fft_int outexp;
	fft_float din;

	out[0] = 1;
	out_rev[0] = 1;
	outexp = in[1];
	din = mp_unexp_mp2d(n, radix, &in[2]);
	if (outexp % 2 != 0) {
		din *= radix;
		outexp--;
	}
	outexp /= 2;
	din = fft_sqrt(din);
	if (din < 1) {
		din *= radix;
		outexp--;
	}
	out[1] = outexp;
	mp_unexp_d2mp(n, radix, din, &out[2]);
	outexp = -outexp;
	din = ((fft_float)1) / din;
	while (din < 1) {
		din *= radix;
		outexp--;
	}
	out_rev[1] = outexp;
	mp_unexp_d2mp(n, radix, din, &out_rev[2]);
}


void mp_invisqrt_init(fft_int n, fft_int radix, fft_int in, fft_int out[])
{
	void mp_unexp_d2mp(fft_int n, fft_int radix, fft_float din, fft_int out[]);
	fft_int outexp;
	fft_float dout;

	out[0] = 1;
	outexp = 0;
	dout = fft_sqrt(((fft_float)1) / in);
	while (dout < 1) {
		dout *= radix;
		outexp--;
	}
	out[1] = outexp;
	mp_unexp_d2mp(n, radix, dout, &out[2]);
}


void mp_unexp_d2mp(fft_int n, fft_int radix, fft_float din, fft_int out[])
{
	fft_int j, x;

	for (j = 0; j < n; j++) {
		x = (fft_int)din;
		if (x >= radix) {
			x = radix - 1;
			din = radix;
		}
		din = radix * (din - x);
		out[j] = x;
	}
}


fft_float mp_unexp_mp2d(fft_int n, fft_int radix, fft_int in[])
{
	fft_int j;
	fft_float d1_radix, dout;

	d1_radix = ((fft_float)1) / radix;
	dout = 0;
	for (j = n - 1; j >= 0; j--) {
		dout = d1_radix * dout + in[j];
	}
	return dout;
}


fft_int mp_inv_newton(fft_int n, fft_int radix, fft_int in[], fft_int inout[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft, fft_float tmp1fft[],
	fft_float tmp2fft[], fft_int ip[], fft_float w[])
{
	void mp_load_1(fft_int n, fft_int radix, fft_int out[]);
	void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_mulh(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
		fft_int nfft, fft_float in1fft[], fft_float outfft[],
		fft_int ip[], fft_float w[]);
	void mp_mulh_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[],
		fft_int shift, fft_int in2[], fft_int out[], fft_int nfft, fft_float outfft[],
		fft_int ip[], fft_float w[]);
	fft_int n_h, shift, prc;

	shift = (nfft >> 1) + 1;
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp1 = inout * (upper) in (half to normal precision) ---- */
	mp_round(n, radix, shift, inout);
	mp_mulh(n, radix, inout, in, tmp1,
		nfft, tmp1fft, tmp2fft, ip, w);
	/* ---- tmp2 = 1 - tmp1 ---- */
	mp_load_1(n, radix, tmp2);
	mp_sub(n, radix, tmp2, tmp1, tmp2);
	/* ---- tmp2 -= inout * (lower) in (half precision) ---- */
	mp_mulh_use_in1fft(n, radix, tmp1fft, shift, in, tmp1,
		nfft, tmp2fft, ip, w);
	mp_sub(n_h, radix, tmp2, tmp1, tmp2);
	/* ---- get precision ---- */
	prc = -tmp2[1];
	if (tmp2[0] == 0) {
		prc = nfft + 1;
	}
	/* ---- tmp2 *= inout (half precision) ---- */
	mp_mulh_use_in1fft(n_h, radix, tmp1fft, 0, tmp2, tmp2,
		nfft, tmp2fft, ip, w);
	/* ---- inout += tmp2 ---- */
	if (tmp2[0] != 0) {
		mp_add(n, radix, inout, tmp2, inout);
	}
	return prc;
}


fft_int mp_sqrt_newton(fft_int n, fft_int radix, fft_int in[], fft_int inout[],
	fft_int inout_rev[], fft_int tmp[], fft_int nfft, fft_float tmp1fft[],
	fft_float tmp2fft[], fft_int ip[], fft_float w[], fft_int* n_tmp1fft)
{
	void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_idiv_2(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_mulh(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[],
		fft_int nfft, fft_float in1fft[], fft_float outfft[],
		fft_int ip[], fft_float w[]);
	void mp_squh(fft_int n, fft_int radix, fft_int in[], fft_int out[],
		fft_int nfft, fft_float outfft[], fft_int ip[], fft_float w[]);
	void mp_squh_use_in1fft(fft_int n, fft_int radix, fft_float inoutfft[], fft_int out[],
		fft_int nfft, fft_int ip[], fft_float w[]);
	fft_int n_h, nfft_h, shift, prc;

	nfft_h = nfft >> 1;
	shift = nfft_h + 1;
	if (nfft_h < 2) {
		nfft_h = 2;
	}
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp = inout_rev^2 (1/4 to half precision) ---- */
	mp_round(n_h, radix, (nfft_h >> 1) + 1, inout_rev);
	if (*n_tmp1fft != nfft_h) {
		mp_squh(n_h, radix, inout_rev, tmp,
			nfft_h, tmp1fft, ip, w);
	}
	else {
		mp_squh_use_in1fft(n_h, radix, tmp1fft, tmp,
			nfft_h, ip, w);
	}
	/* ---- tmp = inout_rev - inout * tmp (half precision) ---- */
	mp_round(n, radix, shift, inout);
	mp_mulh(n_h, radix, inout, tmp, tmp,
		nfft, tmp1fft, tmp2fft, ip, w);
	mp_sub(n_h, radix, inout_rev, tmp, tmp);
	/* ---- inout_rev += tmp ---- */
	mp_add(n_h, radix, inout_rev, tmp, inout_rev);
	/* ---- tmp = in - inout^2 (half to normal precision) ---- */
	mp_squh_use_in1fft(n, radix, tmp1fft, tmp,
		nfft, ip, w);
	mp_sub(n, radix, in, tmp, tmp);
	/* ---- get precision ---- */
	prc = in[1] - tmp[1];
	if (in[2] > tmp[2]) {
		prc++;
	}
	if (tmp[0] == 0) {
		prc = nfft + 1;
	}
	/* ---- tmp = tmp * inout_rev / 2 (half precision) ---- */
	mp_round(n_h, radix, shift, inout_rev);
	mp_mulh(n_h, radix, inout_rev, tmp, tmp,
		nfft, tmp1fft, tmp2fft, ip, w);
	*n_tmp1fft = nfft;
	mp_idiv_2(n_h, radix, tmp, tmp);
	/* ---- inout += tmp ---- */
	if (tmp[0] != 0) {
		mp_add(n, radix, inout, tmp, inout);
	}
	return prc;
}


fft_int mp_invisqrt_newton(fft_int n, fft_int radix, fft_int in, fft_int inout[],
	fft_int tmp1[], fft_int tmp2[], fft_int nfft, fft_float tmp1fft[],
	fft_float tmp2fft[], fft_int ip[], fft_float w[])
{
	void mp_load_1(fft_int n, fft_int radix, fft_int out[]);
	void mp_round(fft_int n, fft_int radix, fft_int m, fft_int inout[]);
	void mp_add(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_sub(fft_int n, fft_int radix, fft_int in1[], fft_int in2[], fft_int out[]);
	void mp_imul(fft_int n, fft_int radix, fft_int in1[], fft_int in2, fft_int out[]);
	void mp_idiv_2(fft_int n, fft_int radix, fft_int in[], fft_int out[]);
	void mp_squh_save_infft(fft_int n, fft_int radix, fft_int in[], fft_int out[],
		fft_int nfft, fft_float infft[], fft_float outfft[],
		fft_int ip[], fft_float w[]);
	void mp_mulh_use_in1fft(fft_int n, fft_int radix, fft_float in1fft[],
		fft_int shift, fft_int in2[], fft_int out[], fft_int nfft, fft_float outfft[],
		fft_int ip[], fft_float w[]);
	fft_int n_h, shift, prc;

	shift = (nfft >> 1) + 1;
	n_h = n / 2 + 1;
	if (n_h < n - shift) {
		n_h = n - shift;
	}
	/* ---- tmp1 = in * inout^2 (half to normal precision) ---- */
	mp_round(n, radix, shift, inout);
	mp_squh_save_infft(n, radix, inout, tmp1,
		nfft, tmp1fft, tmp2fft, ip, w);
	mp_imul(n, radix, tmp1, in, tmp1);
	/* ---- tmp2 = 1 - tmp1 ---- */
	mp_load_1(n, radix, tmp2);
	mp_sub(n, radix, tmp2, tmp1, tmp2);
	/* ---- get precision ---- */
	prc = -tmp2[1];
	if (tmp2[0] == 0) {
		prc = nfft + 1;
	}
	/* ---- tmp2 *= inout / 2 (half precision) ---- */
	mp_mulh_use_in1fft(n_h, radix, tmp1fft, 0, tmp2, tmp2,
		nfft, tmp2fft, ip, w);
	mp_idiv_2(n_h, radix, tmp2, tmp2);
	/* ---- inout += tmp2 ---- */
	if (tmp2[0] != 0) {
		mp_add(n, radix, inout, tmp2, inout);
	}
	return prc;
}


/* -------- mp_io routines -------- */


void mp_sprintf(fft_int n, fft_int log10_radix, fft_int in[], char out[])
{
	fft_int j, k, x, y, outexp, shift;

	if (in[0] < 0) {
		*out++ = '-';
	}
	x = in[2];
	shift = log10_radix;
	for (k = log10_radix; k > 0; k--) {
		y = x % 10;
		x /= 10;
		out[k] = '0' + y;
		if (y != 0) {
			shift = k;
		}
	}
	out[0] = out[shift];
	out[1] = '.';
	for (k = 1; k <= log10_radix - shift; k++) {
		out[k + 1] = out[k + shift];
	}
	outexp = log10_radix - shift;
	out += outexp + 2;
	for (j = 3; j <= n + 1; j++) {
		x = in[j];
		for (k = log10_radix - 1; k >= 0; k--) {
			y = x % 10;
			x /= 10;
			out[k] = '0' + y;
		}
		out += log10_radix;
	}
	*out++ = 'e';
	outexp += log10_radix * in[1];
	sprintf(out, "%.0f", (double)outexp);
}


void mp_sscanf(fft_int n, fft_int log10_radix, const char in[], fft_int out[])
{
	const char* s;
	fft_int j, x, outexp, outexp_mod;
	double t;

	while (*in == ' ') {
		in++;
	}
	out[0] = 1;
	if (*in == '-') {
		out[0] = -1;
		in++;
	}
	else if (*in == '+') {
		in++;
	}
	while (*in == ' ' || *in == '0') {
		in++;
	}
	outexp = 0;
	for (s = in; *s != '\0'; s++) {
		if (*s == 'e' || *s == 'E' || *s == 'd' || *s == 'D') {
			if (sscanf(++s, "%lg", &t) != 1) {
				outexp = 0;
			}
			outexp = (fft_int)t;
			break;
		}
	}
	if (*in == '.') {
		do {
			outexp--;
			while (*++in == ' ');
		} while (*in == '0' && *in != '\0');
	}
	else if (*in != '\0') {
		s = in;
		while (*++s == ' ');
		while (*s >= '0' && *s <= '9' && *s != '\0') {
			outexp++;
			while (*++s == ' ');
		}
	}
	x = outexp / log10_radix;
	outexp_mod = outexp - log10_radix * x;
	if (outexp_mod < 0) {
		x--;
		outexp_mod += log10_radix;
	}
	out[1] = x;
	x = 0;
	j = 2;
	for (s = in; *s != '\0'; s++) {
		if (*s == '.' || *s == ' ') {
			continue;
		}
		if (*s < '0' || *s > '9') {
			break;
		}
		x = 10 * x + (*s - '0');
		if (--outexp_mod < 0) {
			if (j > n + 1) {
				break;
			}
			out[j++] = x;
			x = 0;
			outexp_mod = log10_radix - 1;
		}
	}
	while (outexp_mod-- >= 0) {
		x *= 10;
	}
	while (j <= n + 1) {
		out[j++] = x;
		x = 0;
	}
	if (out[2] == 0) {
		out[0] = 0;
		out[1] = 0;
	}
}


void mp_fprintf(fft_int n, fft_int log10_radix, fft_int in[], FILE* fout)
{
	fft_int j, k, x, y, outexp, shift;
	char out[256];

	if (in[0] < 0) {
		putc('-', fout);
	}
	x = in[2];
	shift = log10_radix;
	for (k = log10_radix; k > 0; k--) {
		y = x % 10;
		x /= 10;
		out[k] = '0' + y;
		if (y != 0) {
			shift = k;
		}
	}
	putc(out[shift], fout);
	putc('.', fout);
	for (k = 1; k <= log10_radix - shift; k++) {
		putc(out[k + shift], fout);
	}
	outexp = log10_radix - shift;
	for (j = 3; j <= n + 1; j++) {
		x = in[j];
		for (k = log10_radix - 1; k >= 0; k--) {
			y = x % 10;
			x /= 10;
			out[k] = '0' + y;
		}
		for (k = 0; k < log10_radix; k++) {
			putc(out[k], fout);
		}
	}
	putc('e', fout);
	outexp += log10_radix * in[1];
	sprintf(out, "%.0f", (double)outexp);
	for (k = 0; out[k] != '\0'; k++) {
		putc(out[k], fout);
	}
}


fft_int mp_chksum(fft_int n, fft_int in[])
{
	fft_int j, sum;

	sum = 0;
	for (j = 0; j <= n + 1; j++) {
		sum ^= in[j];
	}
	return sum;
}


/* **************** MANUAL ****************
Fast Fourier/Cosine/Sine Transform
	dimension   :one
	data length :power of 2
	decimation  :frequency
	radix       :split-radix
	data        :inplace
	table       :use
functions
	cdft: Complex Discrete Fourier Transform
	rdft: Real Discrete Fourier Transform
	ddct: Discrete Cosine Transform
	ddst: Discrete Sine Transform
	dfct: Cosine Transform of RDFT (Real Symmetric DFT)
	dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
function prototypes
	void cdft(fft_int, fft_int, fft_float *, fft_int *, fft_float *);
	void rdft(fft_int, fft_int, fft_float *, fft_int *, fft_float *);
	void ddct(fft_int, fft_int, fft_float *, fft_int *, fft_float *);
	void ddst(fft_int, fft_int, fft_float *, fft_int *, fft_float *);
	void dfct(fft_int, fft_float *, fft_float *, fft_int *, fft_float *);
	void dfst(fft_int, fft_float *, fft_float *, fft_int *, fft_float *);
macro definitions
	USE_CDFT_PTHREADS : default=not defined
		CDFT_THREADS_BEGIN_N  : must be >= 512, default=8192
		CDFT_4THREADS_BEGIN_N : must be >= 512, default=65536
	USE_CDFT_WINTHREADS : default=not defined
		CDFT_THREADS_BEGIN_N  : must be >= 512, default=32768
		CDFT_4THREADS_BEGIN_N : must be >= 512, default=524288


-------- Complex DFT (Discrete Fourier Transform) --------
	[definition]
		<case1>
			X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
		<case2>
			X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
		(notes: sum_j=0^n-1 is a summation from j=0 to n-1)
	[usage]
		<case1>
			ip[0] = 0; // first time only
			cdft(2*n, 1, a, ip, w);
		<case2>
			ip[0] = 0; // first time only
			cdft(2*n, -1, a, ip, w);
	[parameters]
		2*n            :data length (fft_int)
						n >= 1, n = power of 2
		a[0...2*n-1]   :input/output data (fft_float *)
						input data
							a[2*j] = Re(x[j]),
							a[2*j+1] = Im(x[j]), 0<=j<n
						output data
							a[2*k] = Re(X[k]),
							a[2*k+1] = Im(X[k]), 0<=k<n
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n/2-1]   :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			cdft(2*n, -1, a, ip, w);
		is
			cdft(2*n, 1, a, ip, w);
			for (j = 0; j <= 2 * n - 1; j++) {
				a[j] *= 1.0 / n;
			}
		.


-------- Real DFT / Inverse of Real DFT --------
	[definition]
		<case1> RDFT
			R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
			I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
		<case2> IRDFT (excluding scale)
			a[k] = (R[0] + R[n/2]*cos(pi*k))/2 +
				   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) +
				   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
	[usage]
		<case1>
			ip[0] = 0; // first time only
			rdft(n, 1, a, ip, w);
		<case2>
			ip[0] = 0; // first time only
			rdft(n, -1, a, ip, w);
	[parameters]
		n              :data length (fft_int)
						n >= 2, n = power of 2
		a[0...n-1]     :input/output data (fft_float *)
						<case1>
							output data
								a[2*k] = R[k], 0<=k<n/2
								a[2*k+1] = I[k], 0<k<n/2
								a[1] = R[n/2]
						<case2>
							input data
								a[2*j] = R[j], 0<=j<n/2
								a[2*j+1] = I[j], 0<j<n/2
								a[1] = R[n/2]
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n/2)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n/2+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n/2-1]   :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			rdft(n, 1, a, ip, w);
		is
			rdft(n, -1, a, ip, w);
			for (j = 0; j <= n - 1; j++) {
				a[j] *= 2.0 / n;
			}
		.


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
	[definition]
		<case1> IDCT (excluding scale)
			C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
		<case2> DCT
			C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
	[usage]
		<case1>
			ip[0] = 0; // first time only
			ddct(n, 1, a, ip, w);
		<case2>
			ip[0] = 0; // first time only
			ddct(n, -1, a, ip, w);
	[parameters]
		n              :data length (fft_int)
						n >= 2, n = power of 2
		a[0...n-1]     :input/output data (fft_float *)
						output data
							a[k] = C[k], 0<=k<n
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n/2)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n/2+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n*5/4-1] :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			ddct(n, -1, a, ip, w);
		is
			a[0] *= 0.5;
			ddct(n, 1, a, ip, w);
			for (j = 0; j <= n - 1; j++) {
				a[j] *= 2.0 / n;
			}
		.


-------- DST (Discrete Sine Transform) / Inverse of DST --------
	[definition]
		<case1> IDST (excluding scale)
			S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
		<case2> DST
			S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
	[usage]
		<case1>
			ip[0] = 0; // first time only
			ddst(n, 1, a, ip, w);
		<case2>
			ip[0] = 0; // first time only
			ddst(n, -1, a, ip, w);
	[parameters]
		n              :data length (fft_int)
						n >= 2, n = power of 2
		a[0...n-1]     :input/output data (fft_float *)
						<case1>
							input data
								a[j] = A[j], 0<j<n
								a[0] = A[n]
							output data
								a[k] = S[k], 0<=k<n
						<case2>
							output data
								a[k] = S[k], 0<k<n
								a[0] = S[n]
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n/2)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n/2+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n*5/4-1] :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			ddst(n, -1, a, ip, w);
		is
			a[0] *= 0.5;
			ddst(n, 1, a, ip, w);
			for (j = 0; j <= n - 1; j++) {
				a[j] *= 2.0 / n;
			}
		.


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
	[definition]
		C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
	[usage]
		ip[0] = 0; // first time only
		dfct(n, a, t, ip, w);
	[parameters]
		n              :data length - 1 (fft_int)
						n >= 2, n = power of 2
		a[0...n]       :input/output data (fft_float *)
						output data
							a[k] = C[k], 0<=k<=n
		t[0...n/2]     :work area (fft_float *)
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n/4)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n/4+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n*5/8-1] :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			a[0] *= 0.5;
			a[n] *= 0.5;
			dfct(n, a, t, ip, w);
		is
			a[0] *= 0.5;
			a[n] *= 0.5;
			dfct(n, a, t, ip, w);
			for (j = 0; j <= n; j++) {
				a[j] *= 2.0 / n;
			}
		.


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
	[definition]
		S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
	[usage]
		ip[0] = 0; // first time only
		dfst(n, a, t, ip, w);
	[parameters]
		n              :data length + 1 (fft_int)
						n >= 2, n = power of 2
		a[0...n-1]     :input/output data (fft_float *)
						output data
							a[k] = S[k], 0<k<n
						(a[0] is used for work area)
		t[0...n/2-1]   :work area (fft_float *)
		ip[0...*]      :work area for bit reversal (fft_int *)
						length of ip >= 2+sqrt(n/4)
						strictly,
						length of ip >=
							2+(1<<(int)(log(n/4+0.5)/log(2))/2).
						ip[0],ip[1] are pointers of the cos/sin table.
		w[0...n*5/8-1] :cos/sin table (fft_float *)
						w[],ip[] are initialized if ip[0] == 0.
	[remark]
		Inverse of
			dfst(n, a, t, ip, w);
		is
			dfst(n, a, t, ip, w);
			for (j = 1; j <= n - 1; j++) {
				a[j] *= 2.0 / n;
			}
		.


Appendix :
	The cos/sin table is recalculated when the larger table required.
	w[] and ip[] are compatible with all routines.
	**************** END MANUAL ****************  */


void cdft(fft_int n, fft_int isgn, fft_float* a, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void cftbsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	fft_int nw;

	nw = ip[0];
	if (n > (nw << 2)) {
		nw = n >> 2;
		makewt(nw, ip, w);
	}
	if (isgn >= 0) {
		cftfsub(n, a, ip, nw, w);
	}
	else {
		cftbsub(n, a, ip, nw, w);
	}
}


void rdft(fft_int n, fft_int isgn, fft_float* a, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void makect(fft_int nc, fft_int * ip, fft_float * c);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void cftbsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void rftfsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void rftbsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_int nw, nc;
	fft_float xi;

	nw = ip[0];
	if (n > (nw << 2)) {
		nw = n >> 2;
		makewt(nw, ip, w);
	}
	nc = ip[1];
	if (n > (nc << 2)) {
		nc = n >> 2;
		makect(nc, ip, w + nw);
	}
	if (isgn >= 0) {
		if (n > 4) {
			cftfsub(n, a, ip, nw, w);
			rftfsub(n, a, nc, w + nw);
		}
		else if (n == 4) {
			cftfsub(n, a, ip, nw, w);
		}
		xi = a[0] - a[1];
		a[0] += a[1];
		a[1] = xi;
	}
	else {
		a[1] = FC_HALF * (a[0] - a[1]);
		a[0] -= a[1];
		if (n > 4) {
			rftbsub(n, a, nc, w + nw);
			cftbsub(n, a, ip, nw, w);
		}
		else if (n == 4) {
			cftbsub(n, a, ip, nw, w);
		}
	}
}


void ddct(fft_int n, fft_int isgn, fft_float* a, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void makect(fft_int nc, fft_int * ip, fft_float * c);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void cftbsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void rftfsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void rftbsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void dctsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_int j, nw, nc;
	fft_float xr;

	nw = ip[0];
	if (n > (nw << 2)) {
		nw = n >> 2;
		makewt(nw, ip, w);
	}
	nc = ip[1];
	if (n > nc) {
		nc = n;
		makect(nc, ip, w + nw);
	}
	if (isgn < 0) {
		xr = a[n - 1];
		for (j = n - 2; j >= 2; j -= 2) {
			a[j + 1] = a[j] - a[j - 1];
			a[j] += a[j - 1];
		}
		a[1] = a[0] - xr;
		a[0] += xr;
		if (n > 4) {
			rftbsub(n, a, nc, w + nw);
			cftbsub(n, a, ip, nw, w);
		}
		else if (n == 4) {
			cftbsub(n, a, ip, nw, w);
		}
	}
	dctsub(n, a, nc, w + nw);
	if (isgn >= 0) {
		if (n > 4) {
			cftfsub(n, a, ip, nw, w);
			rftfsub(n, a, nc, w + nw);
		}
		else if (n == 4) {
			cftfsub(n, a, ip, nw, w);
		}
		xr = a[0] - a[1];
		a[0] += a[1];
		for (j = 2; j < n; j += 2) {
			a[j - 1] = a[j] - a[j + 1];
			a[j] += a[j + 1];
		}
		a[n - 1] = xr;
	}
}


void ddst(fft_int n, fft_int isgn, fft_float* a, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void makect(fft_int nc, fft_int * ip, fft_float * c);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void cftbsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void rftfsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void rftbsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void dstsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_int j, nw, nc;
	fft_float xr;

	nw = ip[0];
	if (n > (nw << 2)) {
		nw = n >> 2;
		makewt(nw, ip, w);
	}
	nc = ip[1];
	if (n > nc) {
		nc = n;
		makect(nc, ip, w + nw);
	}
	if (isgn < 0) {
		xr = a[n - 1];
		for (j = n - 2; j >= 2; j -= 2) {
			a[j + 1] = -a[j] - a[j - 1];
			a[j] -= a[j - 1];
		}
		a[1] = a[0] + xr;
		a[0] -= xr;
		if (n > 4) {
			rftbsub(n, a, nc, w + nw);
			cftbsub(n, a, ip, nw, w);
		}
		else if (n == 4) {
			cftbsub(n, a, ip, nw, w);
		}
	}
	dstsub(n, a, nc, w + nw);
	if (isgn >= 0) {
		if (n > 4) {
			cftfsub(n, a, ip, nw, w);
			rftfsub(n, a, nc, w + nw);
		}
		else if (n == 4) {
			cftfsub(n, a, ip, nw, w);
		}
		xr = a[0] - a[1];
		a[0] += a[1];
		for (j = 2; j < n; j += 2) {
			a[j - 1] = -a[j] - a[j + 1];
			a[j] -= a[j + 1];
		}
		a[n - 1] = -xr;
	}
}


void dfct(fft_int n, fft_float* a, fft_float* t, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void makect(fft_int nc, fft_int * ip, fft_float * c);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void rftfsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void dctsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_int j, k, l, m, mh, nw, nc;
	fft_float xr, xi, yr, yi;

	nw = ip[0];
	if (n > (nw << 3)) {
		nw = n >> 3;
		makewt(nw, ip, w);
	}
	nc = ip[1];
	if (n > (nc << 1)) {
		nc = n >> 1;
		makect(nc, ip, w + nw);
	}
	m = n >> 1;
	yi = a[m];
	xi = a[0] + a[n];
	a[0] -= a[n];
	t[0] = xi - yi;
	t[m] = xi + yi;
	if (n > 2) {
		mh = m >> 1;
		for (j = 1; j < mh; j++) {
			k = m - j;
			xr = a[j] - a[n - j];
			xi = a[j] + a[n - j];
			yr = a[k] - a[n - k];
			yi = a[k] + a[n - k];
			a[j] = xr;
			a[k] = yr;
			t[j] = xi - yi;
			t[k] = xi + yi;
		}
		t[mh] = a[mh] + a[n - mh];
		a[mh] -= a[n - mh];
		dctsub(m, a, nc, w + nw);
		if (m > 4) {
			cftfsub(m, a, ip, nw, w);
			rftfsub(m, a, nc, w + nw);
		}
		else if (m == 4) {
			cftfsub(m, a, ip, nw, w);
		}
		a[n - 1] = a[0] - a[1];
		a[1] = a[0] + a[1];
		for (j = m - 2; j >= 2; j -= 2) {
			a[2 * j + 1] = a[j] + a[j + 1];
			a[2 * j - 1] = a[j] - a[j + 1];
		}
		l = 2;
		m = mh;
		while (m >= 2) {
			dctsub(m, t, nc, w + nw);
			if (m > 4) {
				cftfsub(m, t, ip, nw, w);
				rftfsub(m, t, nc, w + nw);
			}
			else if (m == 4) {
				cftfsub(m, t, ip, nw, w);
			}
			a[n - l] = t[0] - t[1];
			a[l] = t[0] + t[1];
			k = 0;
			for (j = 2; j < m; j += 2) {
				k += l << 2;
				a[k - l] = t[j] - t[j + 1];
				a[k + l] = t[j] + t[j + 1];
			}
			l <<= 1;
			mh = m >> 1;
			for (j = 0; j < mh; j++) {
				k = m - j;
				t[j] = t[m + k] - t[m + j];
				t[k] = t[m + k] + t[m + j];
			}
			t[mh] = t[m + mh];
			m = mh;
		}
		a[l] = t[0];
		a[n] = t[2] - t[1];
		a[0] = t[2] + t[1];
	}
	else {
		a[1] = a[0];
		a[2] = t[0];
		a[0] = t[1];
	}
}


void dfst(fft_int n, fft_float* a, fft_float* t, fft_int* ip, fft_float* w)
{
	void makewt(fft_int nw, fft_int * ip, fft_float * w);
	void makect(fft_int nc, fft_int * ip, fft_float * c);
	void cftfsub(fft_int n, fft_float * a, fft_int * ip, fft_int nw, fft_float * w);
	void rftfsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	void dstsub(fft_int n, fft_float * a, fft_int nc, fft_float * c);
	fft_int j, k, l, m, mh, nw, nc;
	fft_float xr, xi, yr, yi;

	nw = ip[0];
	if (n > (nw << 3)) {
		nw = n >> 3;
		makewt(nw, ip, w);
	}
	nc = ip[1];
	if (n > (nc << 1)) {
		nc = n >> 1;
		makect(nc, ip, w + nw);
	}
	if (n > 2) {
		m = n >> 1;
		mh = m >> 1;
		for (j = 1; j < mh; j++) {
			k = m - j;
			xr = a[j] + a[n - j];
			xi = a[j] - a[n - j];
			yr = a[k] + a[n - k];
			yi = a[k] - a[n - k];
			a[j] = xr;
			a[k] = yr;
			t[j] = xi + yi;
			t[k] = xi - yi;
		}
		t[0] = a[mh] - a[n - mh];
		a[mh] += a[n - mh];
		a[0] = a[m];
		dstsub(m, a, nc, w + nw);
		if (m > 4) {
			cftfsub(m, a, ip, nw, w);
			rftfsub(m, a, nc, w + nw);
		}
		else if (m == 4) {
			cftfsub(m, a, ip, nw, w);
		}
		a[n - 1] = a[1] - a[0];
		a[1] = a[0] + a[1];
		for (j = m - 2; j >= 2; j -= 2) {
			a[2 * j + 1] = a[j] - a[j + 1];
			a[2 * j - 1] = -a[j] - a[j + 1];
		}
		l = 2;
		m = mh;
		while (m >= 2) {
			dstsub(m, t, nc, w + nw);
			if (m > 4) {
				cftfsub(m, t, ip, nw, w);
				rftfsub(m, t, nc, w + nw);
			}
			else if (m == 4) {
				cftfsub(m, t, ip, nw, w);
			}
			a[n - l] = t[1] - t[0];
			a[l] = t[0] + t[1];
			k = 0;
			for (j = 2; j < m; j += 2) {
				k += l << 2;
				a[k - l] = -t[j] - t[j + 1];
				a[k + l] = t[j] - t[j + 1];
			}
			l <<= 1;
			mh = m >> 1;
			for (j = 1; j < mh; j++) {
				k = m - j;
				t[j] = t[m + k] + t[m + j];
				t[k] = t[m + k] - t[m + j];
			}
			t[0] = t[m + mh];
			m = mh;
		}
		a[l] = t[0];
	}
	a[0] = 0;
}


/* -------- initializing routines -------- */


#include <math.h>

void makewt(fft_int nw, fft_int* ip, fft_float* w)
{
	void makeipt(fft_int nw, fft_int * ip);
	fft_int j, nwh, nw0, nw1;
	fft_float delta, wn4r, wk1r, wk1i, wk3r, wk3i;

	ip[0] = nw;
	ip[1] = 1;
	if (nw > 2) {
		nwh = nw >> 1;
		delta = fft_atan(1.0) / nwh;
		wn4r = fft_cos(delta * nwh);
		w[0] = 1;
		w[1] = wn4r;
		if (nwh == 4) {
			w[2] = fft_cos(delta * 2);
			w[3] = fft_sin(delta * 2);
		}
		else if (nwh > 4) {
			makeipt(nw, ip);
			w[2] = FC_HALF / fft_cos(delta * 2);
			w[3] = FC_HALF / fft_cos(delta * 6);
			for (j = 4; j < nwh; j += 4) {
				w[j] = fft_cos(delta * j);
				w[j + 1] = fft_sin(delta * j);
				w[j + 2] = fft_cos(3 * delta * j);
				w[j + 3] = -fft_sin(3 * delta * j);
			}
		}
		nw0 = 0;
		while (nwh > 2) {
			nw1 = nw0 + nwh;
			nwh >>= 1;
			w[nw1] = 1;
			w[nw1 + 1] = wn4r;
			if (nwh == 4) {
				wk1r = w[nw0 + 4];
				wk1i = w[nw0 + 5];
				w[nw1 + 2] = wk1r;
				w[nw1 + 3] = wk1i;
			}
			else if (nwh > 4) {
				wk1r = w[nw0 + 4];
				wk3r = w[nw0 + 6];
				w[nw1 + 2] = FC_HALF / wk1r;
				w[nw1 + 3] = FC_HALF / wk3r;
				for (j = 4; j < nwh; j += 4) {
					wk1r = w[nw0 + 2 * j];
					wk1i = w[nw0 + 2 * j + 1];
					wk3r = w[nw0 + 2 * j + 2];
					wk3i = w[nw0 + 2 * j + 3];
					w[nw1 + j] = wk1r;
					w[nw1 + j + 1] = wk1i;
					w[nw1 + j + 2] = wk3r;
					w[nw1 + j + 3] = wk3i;
				}
			}
			nw0 = nw1;
		}
	}
}


void makeipt(fft_int nw, fft_int* ip)
{
	fft_int j, l, m, m2, p, q;

	ip[2] = 0;
	ip[3] = 16;
	m = 2;
	for (l = nw; l > 32; l >>= 2) {
		m2 = m << 1;
		q = m2 << 3;
		for (j = m; j < m2; j++) {
			p = ip[j] << 2;
			ip[m + j] = p;
			ip[m2 + j] = p + q;
		}
		m = m2;
	}
}


void makect(fft_int nc, fft_int* ip, fft_float* c)
{
	fft_int j, nch;
	fft_float delta;

	ip[1] = nc;
	if (nc > 1) {
		nch = nc >> 1;
		delta = fft_atan(1.0) / nch;
		c[0] = fft_cos(delta * nch);
		c[nch] = FC_HALF * c[0];
		for (j = 1; j < nch; j++) {
			c[j] = FC_HALF * fft_cos(delta * j);
			c[nc - j] = FC_HALF * fft_sin(delta * j);
		}
	}
}


/* -------- child routines -------- */


void cftfsub(fft_int n, fft_float* a, fft_int* ip, fft_int nw, fft_float* w)
{
	void bitrv2(fft_int n, fft_int * ip, fft_float * a);
	void bitrv216(fft_float * a);
	void bitrv208(fft_float * a);
	void cftf1st(fft_int n, fft_float * a, fft_float * w);
	void cftrec4(fft_int n, fft_float * a, fft_int nw, fft_float * w);
	void cftleaf(fft_int n, fft_int isplt, fft_float * a, fft_int nw, fft_float * w);
	void cftfx41(fft_int n, fft_float * a, fft_int nw, fft_float * w);
	void cftf161(fft_float * a, fft_float * w);
	void cftf081(fft_float * a, fft_float * w);
	void cftf040(fft_float * a);
	void cftx020(fft_float * a);
#ifdef USE_CDFT_THREADS
	void cftrec4_th(fft_int n, fft_float * a, fft_int nw, fft_float * w);
#endif /* USE_CDFT_THREADS */

	if (n > 8) {
		if (n > 32) {
			cftf1st(n, a, &w[nw - (n >> 2)]);
#ifdef USE_CDFT_THREADS
			if (n > CDFT_THREADS_BEGIN_N) {
				cftrec4_th(n, a, nw, w);
			}
			else
#endif /* USE_CDFT_THREADS */
				if (n > 512) {
					cftrec4(n, a, nw, w);
				}
				else if (n > 128) {
					cftleaf(n, 1, a, nw, w);
				}
				else {
					cftfx41(n, a, nw, w);
				}
			bitrv2(n, ip, a);
		}
		else if (n == 32) {
			cftf161(a, &w[nw - 8]);
			bitrv216(a);
		}
		else {
			cftf081(a, w);
			bitrv208(a);
		}
	}
	else if (n == 8) {
		cftf040(a);
	}
	else if (n == 4) {
		cftx020(a);
	}
}


void cftbsub(fft_int n, fft_float* a, fft_int* ip, fft_int nw, fft_float* w)
{
	void bitrv2conj(fft_int n, fft_int * ip, fft_float * a);
	void bitrv216neg(fft_float * a);
	void bitrv208neg(fft_float * a);
	void cftb1st(fft_int n, fft_float * a, fft_float * w);
	void cftrec4(fft_int n, fft_float * a, fft_int nw, fft_float * w);
	void cftleaf(fft_int n, fft_int isplt, fft_float * a, fft_int nw, fft_float * w);
	void cftfx41(fft_int n, fft_float * a, fft_int nw, fft_float * w);
	void cftf161(fft_float * a, fft_float * w);
	void cftf081(fft_float * a, fft_float * w);
	void cftb040(fft_float * a);
	void cftx020(fft_float * a);
#ifdef USE_CDFT_THREADS
	void cftrec4_th(fft_int n, fft_float * a, fft_int nw, fft_float * w);
#endif /* USE_CDFT_THREADS */

	if (n > 8) {
		if (n > 32) {
			cftb1st(n, a, &w[nw - (n >> 2)]);
#ifdef USE_CDFT_THREADS
			if (n > CDFT_THREADS_BEGIN_N) {
				cftrec4_th(n, a, nw, w);
			}
			else
#endif /* USE_CDFT_THREADS */
				if (n > 512) {
					cftrec4(n, a, nw, w);
				}
				else if (n > 128) {
					cftleaf(n, 1, a, nw, w);
				}
				else {
					cftfx41(n, a, nw, w);
				}
			bitrv2conj(n, ip, a);
		}
		else if (n == 32) {
			cftf161(a, &w[nw - 8]);
			bitrv216neg(a);
		}
		else {
			cftf081(a, w);
			bitrv208neg(a);
		}
	}
	else if (n == 8) {
		cftb040(a);
	}
	else if (n == 4) {
		cftx020(a);
	}
}


void bitrv2(fft_int n, fft_int* ip, fft_float* a)
{
	fft_int j, j1, k, k1, l, m, nh, nm;
	fft_float xr, xi, yr, yi;

	m = 1;
	for (l = n >> 2; l > 8; l >>= 2) {
		m <<= 1;
	}
	nh = n >> 1;
	nm = 4 * m;
	if (l == 8) {
		for (k = 0; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 4 * j + 2 * ip[m + k];
				k1 = 4 * k + 2 * ip[m + j];
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 -= nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nh;
				k1 += 2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 += nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += 2;
				k1 += nh;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 -= nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nh;
				k1 -= 2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 += nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
			k1 = 4 * k + 2 * ip[m + k];
			j1 = k1 + 2;
			k1 += nh;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nm;
			k1 += 2 * nm;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nm;
			k1 -= nm;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 -= 2;
			k1 -= nh;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nh + 2;
			k1 += nh + 2;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 -= nh - nm;
			k1 += 2 * nm - 2;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
		}
	}
	else {
		for (k = 0; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 4 * j + ip[m + k];
				k1 = 4 * k + ip[m + j];
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nh;
				k1 += 2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += 2;
				k1 += nh;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nh;
				k1 -= 2;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= nm;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
			k1 = 4 * k + ip[m + k];
			j1 = k1 + 2;
			k1 += nh;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nm;
			k1 += nm;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
		}
	}
}


void bitrv2conj(fft_int n, fft_int* ip, fft_float* a)
{
	fft_int j, j1, k, k1, l, m, nh, nm;
	fft_float xr, xi, yr, yi;

	m = 1;
	for (l = n >> 2; l > 8; l >>= 2) {
		m <<= 1;
	}
	nh = n >> 1;
	nm = 4 * m;
	if (l == 8) {
		for (k = 0; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 4 * j + 2 * ip[m + k];
				k1 = 4 * k + 2 * ip[m + j];
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 -= nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nh;
				k1 += 2;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 += nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += 2;
				k1 += nh;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 -= nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nh;
				k1 -= 2;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 += nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= 2 * nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
			k1 = 4 * k + 2 * ip[m + k];
			j1 = k1 + 2;
			k1 += nh;
			a[j1 - 1] = -a[j1 - 1];
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			a[k1 + 3] = -a[k1 + 3];
			j1 += nm;
			k1 += 2 * nm;
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nm;
			k1 -= nm;
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 -= 2;
			k1 -= nh;
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 += nh + 2;
			k1 += nh + 2;
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			j1 -= nh - nm;
			k1 += 2 * nm - 2;
			a[j1 - 1] = -a[j1 - 1];
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			a[k1 + 3] = -a[k1 + 3];
		}
	}
	else {
		for (k = 0; k < m; k++) {
			for (j = 0; j < k; j++) {
				j1 = 4 * j + ip[m + k];
				k1 = 4 * k + ip[m + j];
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nh;
				k1 += 2;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += 2;
				k1 += nh;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += nm;
				k1 += nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nh;
				k1 -= 2;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 -= nm;
				k1 -= nm;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
			}
			k1 = 4 * k + ip[m + k];
			j1 = k1 + 2;
			k1 += nh;
			a[j1 - 1] = -a[j1 - 1];
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			a[k1 + 3] = -a[k1 + 3];
			j1 += nm;
			k1 += nm;
			a[j1 - 1] = -a[j1 - 1];
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			a[k1 + 3] = -a[k1 + 3];
		}
	}
}


void bitrv216(fft_float* a)
{
	fft_float x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i,
		x5r, x5i, x7r, x7i, x8r, x8i, x10r, x10i,
		x11r, x11i, x12r, x12i, x13r, x13i, x14r, x14i;

	x1r = a[2];
	x1i = a[3];
	x2r = a[4];
	x2i = a[5];
	x3r = a[6];
	x3i = a[7];
	x4r = a[8];
	x4i = a[9];
	x5r = a[10];
	x5i = a[11];
	x7r = a[14];
	x7i = a[15];
	x8r = a[16];
	x8i = a[17];
	x10r = a[20];
	x10i = a[21];
	x11r = a[22];
	x11i = a[23];
	x12r = a[24];
	x12i = a[25];
	x13r = a[26];
	x13i = a[27];
	x14r = a[28];
	x14i = a[29];
	a[2] = x8r;
	a[3] = x8i;
	a[4] = x4r;
	a[5] = x4i;
	a[6] = x12r;
	a[7] = x12i;
	a[8] = x2r;
	a[9] = x2i;
	a[10] = x10r;
	a[11] = x10i;
	a[14] = x14r;
	a[15] = x14i;
	a[16] = x1r;
	a[17] = x1i;
	a[20] = x5r;
	a[21] = x5i;
	a[22] = x13r;
	a[23] = x13i;
	a[24] = x3r;
	a[25] = x3i;
	a[26] = x11r;
	a[27] = x11i;
	a[28] = x7r;
	a[29] = x7i;
}


void bitrv216neg(fft_float* a)
{
	fft_float x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i,
		x5r, x5i, x6r, x6i, x7r, x7i, x8r, x8i,
		x9r, x9i, x10r, x10i, x11r, x11i, x12r, x12i,
		x13r, x13i, x14r, x14i, x15r, x15i;

	x1r = a[2];
	x1i = a[3];
	x2r = a[4];
	x2i = a[5];
	x3r = a[6];
	x3i = a[7];
	x4r = a[8];
	x4i = a[9];
	x5r = a[10];
	x5i = a[11];
	x6r = a[12];
	x6i = a[13];
	x7r = a[14];
	x7i = a[15];
	x8r = a[16];
	x8i = a[17];
	x9r = a[18];
	x9i = a[19];
	x10r = a[20];
	x10i = a[21];
	x11r = a[22];
	x11i = a[23];
	x12r = a[24];
	x12i = a[25];
	x13r = a[26];
	x13i = a[27];
	x14r = a[28];
	x14i = a[29];
	x15r = a[30];
	x15i = a[31];
	a[2] = x15r;
	a[3] = x15i;
	a[4] = x7r;
	a[5] = x7i;
	a[6] = x11r;
	a[7] = x11i;
	a[8] = x3r;
	a[9] = x3i;
	a[10] = x13r;
	a[11] = x13i;
	a[12] = x5r;
	a[13] = x5i;
	a[14] = x9r;
	a[15] = x9i;
	a[16] = x1r;
	a[17] = x1i;
	a[18] = x14r;
	a[19] = x14i;
	a[20] = x6r;
	a[21] = x6i;
	a[22] = x10r;
	a[23] = x10i;
	a[24] = x2r;
	a[25] = x2i;
	a[26] = x12r;
	a[27] = x12i;
	a[28] = x4r;
	a[29] = x4i;
	a[30] = x8r;
	a[31] = x8i;
}


void bitrv208(fft_float* a)
{
	fft_float x1r, x1i, x3r, x3i, x4r, x4i, x6r, x6i;

	x1r = a[2];
	x1i = a[3];
	x3r = a[6];
	x3i = a[7];
	x4r = a[8];
	x4i = a[9];
	x6r = a[12];
	x6i = a[13];
	a[2] = x4r;
	a[3] = x4i;
	a[6] = x6r;
	a[7] = x6i;
	a[8] = x1r;
	a[9] = x1i;
	a[12] = x3r;
	a[13] = x3i;
}


void bitrv208neg(fft_float* a)
{
	fft_float x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i,
		x5r, x5i, x6r, x6i, x7r, x7i;

	x1r = a[2];
	x1i = a[3];
	x2r = a[4];
	x2i = a[5];
	x3r = a[6];
	x3i = a[7];
	x4r = a[8];
	x4i = a[9];
	x5r = a[10];
	x5i = a[11];
	x6r = a[12];
	x6i = a[13];
	x7r = a[14];
	x7i = a[15];
	a[2] = x7r;
	a[3] = x7i;
	a[4] = x3r;
	a[5] = x3i;
	a[6] = x5r;
	a[7] = x5i;
	a[8] = x1r;
	a[9] = x1i;
	a[10] = x6r;
	a[11] = x6i;
	a[12] = x2r;
	a[13] = x2i;
	a[14] = x4r;
	a[15] = x4i;
}


void cftf1st(fft_int n, fft_float* a, fft_float* w)
{
	fft_int j, j0, j1, j2, j3, k, m, mh;
	fft_float wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i,
		wd1r, wd1i, wd3r, wd3i;
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i;

	mh = n >> 3;
	m = 2 * mh;
	j1 = m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[0] + a[j2];
	x0i = a[1] + a[j2 + 1];
	x1r = a[0] - a[j2];
	x1i = a[1] - a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i - x2i;
	a[j2] = x1r - x3i;
	a[j2 + 1] = x1i + x3r;
	a[j3] = x1r + x3i;
	a[j3 + 1] = x1i - x3r;
	wn4r = w[1];
	csc1 = w[2];
	csc3 = w[3];
	wd1r = 1;
	wd1i = 0;
	wd3r = 1;
	wd3i = 0;
	k = 0;
	for (j = 2; j < mh - 2; j += 4) {
		k += 4;
		wk1r = csc1 * (wd1r + w[k]);
		wk1i = csc1 * (wd1i + w[k + 1]);
		wk3r = csc3 * (wd3r + w[k + 2]);
		wk3i = csc3 * (wd3i + w[k + 3]);
		wd1r = w[k];
		wd1i = w[k + 1];
		wd3r = w[k + 2];
		wd3i = w[k + 3];
		j1 = j + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j] + a[j2];
		x0i = a[j + 1] + a[j2 + 1];
		x1r = a[j] - a[j2];
		x1i = a[j + 1] - a[j2 + 1];
		y0r = a[j + 2] + a[j2 + 2];
		y0i = a[j + 3] + a[j2 + 3];
		y1r = a[j + 2] - a[j2 + 2];
		y1i = a[j + 3] - a[j2 + 3];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		y2r = a[j1 + 2] + a[j3 + 2];
		y2i = a[j1 + 3] + a[j3 + 3];
		y3r = a[j1 + 2] - a[j3 + 2];
		y3i = a[j1 + 3] - a[j3 + 3];
		a[j] = x0r + x2r;
		a[j + 1] = x0i + x2i;
		a[j + 2] = y0r + y2r;
		a[j + 3] = y0i + y2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i - x2i;
		a[j1 + 2] = y0r - y2r;
		a[j1 + 3] = y0i - y2i;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j2] = wk1r * x0r - wk1i * x0i;
		a[j2 + 1] = wk1r * x0i + wk1i * x0r;
		x0r = y1r - y3i;
		x0i = y1i + y3r;
		a[j2 + 2] = wd1r * x0r - wd1i * x0i;
		a[j2 + 3] = wd1r * x0i + wd1i * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j3] = wk3r * x0r + wk3i * x0i;
		a[j3 + 1] = wk3r * x0i - wk3i * x0r;
		x0r = y1r + y3i;
		x0i = y1i - y3r;
		a[j3 + 2] = wd3r * x0r + wd3i * x0i;
		a[j3 + 3] = wd3r * x0i - wd3i * x0r;
		j0 = m - j;
		j1 = j0 + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j0] + a[j2];
		x0i = a[j0 + 1] + a[j2 + 1];
		x1r = a[j0] - a[j2];
		x1i = a[j0 + 1] - a[j2 + 1];
		y0r = a[j0 - 2] + a[j2 - 2];
		y0i = a[j0 - 1] + a[j2 - 1];
		y1r = a[j0 - 2] - a[j2 - 2];
		y1i = a[j0 - 1] - a[j2 - 1];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		y2r = a[j1 - 2] + a[j3 - 2];
		y2i = a[j1 - 1] + a[j3 - 1];
		y3r = a[j1 - 2] - a[j3 - 2];
		y3i = a[j1 - 1] - a[j3 - 1];
		a[j0] = x0r + x2r;
		a[j0 + 1] = x0i + x2i;
		a[j0 - 2] = y0r + y2r;
		a[j0 - 1] = y0i + y2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i - x2i;
		a[j1 - 2] = y0r - y2r;
		a[j1 - 1] = y0i - y2i;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j2] = wk1i * x0r - wk1r * x0i;
		a[j2 + 1] = wk1i * x0i + wk1r * x0r;
		x0r = y1r - y3i;
		x0i = y1i + y3r;
		a[j2 - 2] = wd1i * x0r - wd1r * x0i;
		a[j2 - 1] = wd1i * x0i + wd1r * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j3] = wk3i * x0r + wk3r * x0i;
		a[j3 + 1] = wk3i * x0i - wk3r * x0r;
		x0r = y1r + y3i;
		x0i = y1i - y3r;
		a[j3 - 2] = wd3i * x0r + wd3r * x0i;
		a[j3 - 1] = wd3i * x0i - wd3r * x0r;
	}
	wk1r = csc1 * (wd1r + wn4r);
	wk1i = csc1 * (wd1i + wn4r);
	wk3r = csc3 * (wd3r - wn4r);
	wk3i = csc3 * (wd3i - wn4r);
	j0 = mh;
	j1 = j0 + m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[j0 - 2] + a[j2 - 2];
	x0i = a[j0 - 1] + a[j2 - 1];
	x1r = a[j0 - 2] - a[j2 - 2];
	x1i = a[j0 - 1] - a[j2 - 1];
	x2r = a[j1 - 2] + a[j3 - 2];
	x2i = a[j1 - 1] + a[j3 - 1];
	x3r = a[j1 - 2] - a[j3 - 2];
	x3i = a[j1 - 1] - a[j3 - 1];
	a[j0 - 2] = x0r + x2r;
	a[j0 - 1] = x0i + x2i;
	a[j1 - 2] = x0r - x2r;
	a[j1 - 1] = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[j2 - 2] = wk1r * x0r - wk1i * x0i;
	a[j2 - 1] = wk1r * x0i + wk1i * x0r;
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	a[j3 - 2] = wk3r * x0r + wk3i * x0i;
	a[j3 - 1] = wk3r * x0i - wk3i * x0r;
	x0r = a[j0] + a[j2];
	x0i = a[j0 + 1] + a[j2 + 1];
	x1r = a[j0] - a[j2];
	x1i = a[j0 + 1] - a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[j0] = x0r + x2r;
	a[j0 + 1] = x0i + x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[j2] = wn4r * (x0r - x0i);
	a[j2 + 1] = wn4r * (x0i + x0r);
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	a[j3] = -wn4r * (x0r + x0i);
	a[j3 + 1] = -wn4r * (x0i - x0r);
	x0r = a[j0 + 2] + a[j2 + 2];
	x0i = a[j0 + 3] + a[j2 + 3];
	x1r = a[j0 + 2] - a[j2 + 2];
	x1i = a[j0 + 3] - a[j2 + 3];
	x2r = a[j1 + 2] + a[j3 + 2];
	x2i = a[j1 + 3] + a[j3 + 3];
	x3r = a[j1 + 2] - a[j3 + 2];
	x3i = a[j1 + 3] - a[j3 + 3];
	a[j0 + 2] = x0r + x2r;
	a[j0 + 3] = x0i + x2i;
	a[j1 + 2] = x0r - x2r;
	a[j1 + 3] = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[j2 + 2] = wk1i * x0r - wk1r * x0i;
	a[j2 + 3] = wk1i * x0i + wk1r * x0r;
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	a[j3 + 2] = wk3i * x0r + wk3r * x0i;
	a[j3 + 3] = wk3i * x0i - wk3r * x0r;
}


void cftb1st(fft_int n, fft_float* a, fft_float* w)
{
	fft_int j, j0, j1, j2, j3, k, m, mh;
	fft_float wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i,
		wd1r, wd1i, wd3r, wd3i;
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i;

	mh = n >> 3;
	m = 2 * mh;
	j1 = m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[0] + a[j2];
	x0i = -a[1] - a[j2 + 1];
	x1r = a[0] - a[j2];
	x1i = -a[1] + a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[0] = x0r + x2r;
	a[1] = x0i - x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i + x2i;
	a[j2] = x1r + x3i;
	a[j2 + 1] = x1i + x3r;
	a[j3] = x1r - x3i;
	a[j3 + 1] = x1i - x3r;
	wn4r = w[1];
	csc1 = w[2];
	csc3 = w[3];
	wd1r = 1;
	wd1i = 0;
	wd3r = 1;
	wd3i = 0;
	k = 0;
	for (j = 2; j < mh - 2; j += 4) {
		k += 4;
		wk1r = csc1 * (wd1r + w[k]);
		wk1i = csc1 * (wd1i + w[k + 1]);
		wk3r = csc3 * (wd3r + w[k + 2]);
		wk3i = csc3 * (wd3i + w[k + 3]);
		wd1r = w[k];
		wd1i = w[k + 1];
		wd3r = w[k + 2];
		wd3i = w[k + 3];
		j1 = j + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j] + a[j2];
		x0i = -a[j + 1] - a[j2 + 1];
		x1r = a[j] - a[j2];
		x1i = -a[j + 1] + a[j2 + 1];
		y0r = a[j + 2] + a[j2 + 2];
		y0i = -a[j + 3] - a[j2 + 3];
		y1r = a[j + 2] - a[j2 + 2];
		y1i = -a[j + 3] + a[j2 + 3];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		y2r = a[j1 + 2] + a[j3 + 2];
		y2i = a[j1 + 3] + a[j3 + 3];
		y3r = a[j1 + 2] - a[j3 + 2];
		y3i = a[j1 + 3] - a[j3 + 3];
		a[j] = x0r + x2r;
		a[j + 1] = x0i - x2i;
		a[j + 2] = y0r + y2r;
		a[j + 3] = y0i - y2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i + x2i;
		a[j1 + 2] = y0r - y2r;
		a[j1 + 3] = y0i + y2i;
		x0r = x1r + x3i;
		x0i = x1i + x3r;
		a[j2] = wk1r * x0r - wk1i * x0i;
		a[j2 + 1] = wk1r * x0i + wk1i * x0r;
		x0r = y1r + y3i;
		x0i = y1i + y3r;
		a[j2 + 2] = wd1r * x0r - wd1i * x0i;
		a[j2 + 3] = wd1r * x0i + wd1i * x0r;
		x0r = x1r - x3i;
		x0i = x1i - x3r;
		a[j3] = wk3r * x0r + wk3i * x0i;
		a[j3 + 1] = wk3r * x0i - wk3i * x0r;
		x0r = y1r - y3i;
		x0i = y1i - y3r;
		a[j3 + 2] = wd3r * x0r + wd3i * x0i;
		a[j3 + 3] = wd3r * x0i - wd3i * x0r;
		j0 = m - j;
		j1 = j0 + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j0] + a[j2];
		x0i = -a[j0 + 1] - a[j2 + 1];
		x1r = a[j0] - a[j2];
		x1i = -a[j0 + 1] + a[j2 + 1];
		y0r = a[j0 - 2] + a[j2 - 2];
		y0i = -a[j0 - 1] - a[j2 - 1];
		y1r = a[j0 - 2] - a[j2 - 2];
		y1i = -a[j0 - 1] + a[j2 - 1];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		y2r = a[j1 - 2] + a[j3 - 2];
		y2i = a[j1 - 1] + a[j3 - 1];
		y3r = a[j1 - 2] - a[j3 - 2];
		y3i = a[j1 - 1] - a[j3 - 1];
		a[j0] = x0r + x2r;
		a[j0 + 1] = x0i - x2i;
		a[j0 - 2] = y0r + y2r;
		a[j0 - 1] = y0i - y2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i + x2i;
		a[j1 - 2] = y0r - y2r;
		a[j1 - 1] = y0i + y2i;
		x0r = x1r + x3i;
		x0i = x1i + x3r;
		a[j2] = wk1i * x0r - wk1r * x0i;
		a[j2 + 1] = wk1i * x0i + wk1r * x0r;
		x0r = y1r + y3i;
		x0i = y1i + y3r;
		a[j2 - 2] = wd1i * x0r - wd1r * x0i;
		a[j2 - 1] = wd1i * x0i + wd1r * x0r;
		x0r = x1r - x3i;
		x0i = x1i - x3r;
		a[j3] = wk3i * x0r + wk3r * x0i;
		a[j3 + 1] = wk3i * x0i - wk3r * x0r;
		x0r = y1r - y3i;
		x0i = y1i - y3r;
		a[j3 - 2] = wd3i * x0r + wd3r * x0i;
		a[j3 - 1] = wd3i * x0i - wd3r * x0r;
	}
	wk1r = csc1 * (wd1r + wn4r);
	wk1i = csc1 * (wd1i + wn4r);
	wk3r = csc3 * (wd3r - wn4r);
	wk3i = csc3 * (wd3i - wn4r);
	j0 = mh;
	j1 = j0 + m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[j0 - 2] + a[j2 - 2];
	x0i = -a[j0 - 1] - a[j2 - 1];
	x1r = a[j0 - 2] - a[j2 - 2];
	x1i = -a[j0 - 1] + a[j2 - 1];
	x2r = a[j1 - 2] + a[j3 - 2];
	x2i = a[j1 - 1] + a[j3 - 1];
	x3r = a[j1 - 2] - a[j3 - 2];
	x3i = a[j1 - 1] - a[j3 - 1];
	a[j0 - 2] = x0r + x2r;
	a[j0 - 1] = x0i - x2i;
	a[j1 - 2] = x0r - x2r;
	a[j1 - 1] = x0i + x2i;
	x0r = x1r + x3i;
	x0i = x1i + x3r;
	a[j2 - 2] = wk1r * x0r - wk1i * x0i;
	a[j2 - 1] = wk1r * x0i + wk1i * x0r;
	x0r = x1r - x3i;
	x0i = x1i - x3r;
	a[j3 - 2] = wk3r * x0r + wk3i * x0i;
	a[j3 - 1] = wk3r * x0i - wk3i * x0r;
	x0r = a[j0] + a[j2];
	x0i = -a[j0 + 1] - a[j2 + 1];
	x1r = a[j0] - a[j2];
	x1i = -a[j0 + 1] + a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[j0] = x0r + x2r;
	a[j0 + 1] = x0i - x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i + x2i;
	x0r = x1r + x3i;
	x0i = x1i + x3r;
	a[j2] = wn4r * (x0r - x0i);
	a[j2 + 1] = wn4r * (x0i + x0r);
	x0r = x1r - x3i;
	x0i = x1i - x3r;
	a[j3] = -wn4r * (x0r + x0i);
	a[j3 + 1] = -wn4r * (x0i - x0r);
	x0r = a[j0 + 2] + a[j2 + 2];
	x0i = -a[j0 + 3] - a[j2 + 3];
	x1r = a[j0 + 2] - a[j2 + 2];
	x1i = -a[j0 + 3] + a[j2 + 3];
	x2r = a[j1 + 2] + a[j3 + 2];
	x2i = a[j1 + 3] + a[j3 + 3];
	x3r = a[j1 + 2] - a[j3 + 2];
	x3i = a[j1 + 3] - a[j3 + 3];
	a[j0 + 2] = x0r + x2r;
	a[j0 + 3] = x0i - x2i;
	a[j1 + 2] = x0r - x2r;
	a[j1 + 3] = x0i + x2i;
	x0r = x1r + x3i;
	x0i = x1i + x3r;
	a[j2 + 2] = wk1i * x0r - wk1r * x0i;
	a[j2 + 3] = wk1i * x0i + wk1r * x0r;
	x0r = x1r - x3i;
	x0i = x1i - x3r;
	a[j3 + 2] = wk3i * x0r + wk3r * x0i;
	a[j3 + 3] = wk3i * x0i - wk3r * x0r;
}


#ifdef USE_CDFT_THREADS
struct cdft_arg_st {
	fft_int n0;
	fft_int n;
	fft_float* a;
	fft_int nw;
	fft_float* w;
};
typedef struct cdft_arg_st cdft_arg_t;


void cftrec4_th(fft_int n, fft_float* a, fft_int nw, fft_float* w)
{
	void* cftrec1_th(void* p);
	void* cftrec2_th(void* p);
	fft_int i, idiv4, m, nthread;
	cdft_thread_t th[4];
	cdft_arg_t ag[4];

	nthread = 2;
	idiv4 = 0;
	m = n >> 1;
	if (n > CDFT_4THREADS_BEGIN_N) {
		nthread = 4;
		idiv4 = 1;
		m >>= 1;
	}
	for (i = 0; i < nthread; i++) {
		ag[i].n0 = n;
		ag[i].n = m;
		ag[i].a = &a[i * m];
		ag[i].nw = nw;
		ag[i].w = w;
		if (i != idiv4) {
			cdft_thread_create(&th[i], cftrec1_th, &ag[i]);
		}
		else {
			cdft_thread_create(&th[i], cftrec2_th, &ag[i]);
		}
	}
	for (i = 0; i < nthread; i++) {
		cdft_thread_wait(th[i]);
	}
}


void* cftrec1_th(void* p)
{
	fft_int cfttree(fft_int n, fft_int j, fft_int k, fft_float * a, fft_int nw, fft_float * w);
	void cftleaf(fft_int n, fft_int isplt, fft_float * a, fft_int nw, fft_float * w);
	void cftmdl1(fft_int n, fft_float * a, fft_float * w);
	fft_int isplt, j, k, m, n, n0, nw;
	fft_float* a, * w;

	n0 = ((cdft_arg_t*)p)->n0;
	n = ((cdft_arg_t*)p)->n;
	a = ((cdft_arg_t*)p)->a;
	nw = ((cdft_arg_t*)p)->nw;
	w = ((cdft_arg_t*)p)->w;
	m = n0;
	while (m > 512) {
		m >>= 2;
		cftmdl1(m, &a[n - m], &w[nw - (m >> 1)]);
	}
	cftleaf(m, 1, &a[n - m], nw, w);
	k = 0;
	for (j = n - m; j > 0; j -= m) {
		k++;
		isplt = cfttree(m, j, k, a, nw, w);
		cftleaf(m, isplt, &a[j - m], nw, w);
	}
	return (void*)0;
}


void* cftrec2_th(void* p)
{
	fft_int cfttree(fft_int n, fft_int j, fft_int k, fft_float * a, fft_int nw, fft_float * w);
	void cftleaf(fft_int n, fft_int isplt, fft_float * a, fft_int nw, fft_float * w);
	void cftmdl2(fft_int n, fft_float * a, fft_float * w);
	fft_int isplt, j, k, m, n, n0, nw;
	fft_float* a, * w;

	n0 = ((cdft_arg_t*)p)->n0;
	n = ((cdft_arg_t*)p)->n;
	a = ((cdft_arg_t*)p)->a;
	nw = ((cdft_arg_t*)p)->nw;
	w = ((cdft_arg_t*)p)->w;
	k = 1;
	m = n0;
	while (m > 512) {
		m >>= 2;
		k <<= 2;
		cftmdl2(m, &a[n - m], &w[nw - m]);
	}
	cftleaf(m, 0, &a[n - m], nw, w);
	k >>= 1;
	for (j = n - m; j > 0; j -= m) {
		k++;
		isplt = cfttree(m, j, k, a, nw, w);
		cftleaf(m, isplt, &a[j - m], nw, w);
	}
	return (void*)0;
}
#endif /* USE_CDFT_THREADS */


void cftrec4(fft_int n, fft_float* a, fft_int nw, fft_float* w)
{
	fft_int cfttree(fft_int n, fft_int j, fft_int k, fft_float * a, fft_int nw, fft_float * w);
	void cftleaf(fft_int n, fft_int isplt, fft_float * a, fft_int nw, fft_float * w);
	void cftmdl1(fft_int n, fft_float * a, fft_float * w);
	fft_int isplt, j, k, m;

	m = n;
	while (m > 512) {
		m >>= 2;
		cftmdl1(m, &a[n - m], &w[nw - (m >> 1)]);
	}
	cftleaf(m, 1, &a[n - m], nw, w);
	k = 0;
	for (j = n - m; j > 0; j -= m) {
		k++;
		isplt = cfttree(m, j, k, a, nw, w);
		cftleaf(m, isplt, &a[j - m], nw, w);
	}
}


fft_int cfttree(fft_int n, fft_int j, fft_int k, fft_float* a, fft_int nw, fft_float* w)
{
	void cftmdl1(fft_int n, fft_float * a, fft_float * w);
	void cftmdl2(fft_int n, fft_float * a, fft_float * w);
	fft_int i, isplt, m;

	if ((k & 3) != 0) {
		isplt = k & 1;
		if (isplt != 0) {
			cftmdl1(n, &a[j - n], &w[nw - (n >> 1)]);
		}
		else {
			cftmdl2(n, &a[j - n], &w[nw - n]);
		}
	}
	else {
		m = n;
		for (i = k; (i & 3) == 0; i >>= 2) {
			m <<= 2;
		}
		isplt = i & 1;
		if (isplt != 0) {
			while (m > 128) {
				cftmdl1(m, &a[j - m], &w[nw - (m >> 1)]);
				m >>= 2;
			}
		}
		else {
			while (m > 128) {
				cftmdl2(m, &a[j - m], &w[nw - m]);
				m >>= 2;
			}
		}
	}
	return isplt;
}


void cftleaf(fft_int n, fft_int isplt, fft_float* a, fft_int nw, fft_float* w)
{
	void cftmdl1(fft_int n, fft_float * a, fft_float * w);
	void cftmdl2(fft_int n, fft_float * a, fft_float * w);
	void cftf161(fft_float * a, fft_float * w);
	void cftf162(fft_float * a, fft_float * w);
	void cftf081(fft_float * a, fft_float * w);
	void cftf082(fft_float * a, fft_float * w);

	if (n == 512) {
		cftmdl1(128, a, &w[nw - 64]);
		cftf161(a, &w[nw - 8]);
		cftf162(&a[32], &w[nw - 32]);
		cftf161(&a[64], &w[nw - 8]);
		cftf161(&a[96], &w[nw - 8]);
		cftmdl2(128, &a[128], &w[nw - 128]);
		cftf161(&a[128], &w[nw - 8]);
		cftf162(&a[160], &w[nw - 32]);
		cftf161(&a[192], &w[nw - 8]);
		cftf162(&a[224], &w[nw - 32]);
		cftmdl1(128, &a[256], &w[nw - 64]);
		cftf161(&a[256], &w[nw - 8]);
		cftf162(&a[288], &w[nw - 32]);
		cftf161(&a[320], &w[nw - 8]);
		cftf161(&a[352], &w[nw - 8]);
		if (isplt != 0) {
			cftmdl1(128, &a[384], &w[nw - 64]);
			cftf161(&a[480], &w[nw - 8]);
		}
		else {
			cftmdl2(128, &a[384], &w[nw - 128]);
			cftf162(&a[480], &w[nw - 32]);
		}
		cftf161(&a[384], &w[nw - 8]);
		cftf162(&a[416], &w[nw - 32]);
		cftf161(&a[448], &w[nw - 8]);
	}
	else {
		cftmdl1(64, a, &w[nw - 32]);
		cftf081(a, &w[nw - 8]);
		cftf082(&a[16], &w[nw - 8]);
		cftf081(&a[32], &w[nw - 8]);
		cftf081(&a[48], &w[nw - 8]);
		cftmdl2(64, &a[64], &w[nw - 64]);
		cftf081(&a[64], &w[nw - 8]);
		cftf082(&a[80], &w[nw - 8]);
		cftf081(&a[96], &w[nw - 8]);
		cftf082(&a[112], &w[nw - 8]);
		cftmdl1(64, &a[128], &w[nw - 32]);
		cftf081(&a[128], &w[nw - 8]);
		cftf082(&a[144], &w[nw - 8]);
		cftf081(&a[160], &w[nw - 8]);
		cftf081(&a[176], &w[nw - 8]);
		if (isplt != 0) {
			cftmdl1(64, &a[192], &w[nw - 32]);
			cftf081(&a[240], &w[nw - 8]);
		}
		else {
			cftmdl2(64, &a[192], &w[nw - 64]);
			cftf082(&a[240], &w[nw - 8]);
		}
		cftf081(&a[192], &w[nw - 8]);
		cftf082(&a[208], &w[nw - 8]);
		cftf081(&a[224], &w[nw - 8]);
	}
}


void cftmdl1(fft_int n, fft_float* a, fft_float* w)
{
	fft_int j, j0, j1, j2, j3, k, m, mh;
	fft_float wn4r, wk1r, wk1i, wk3r, wk3i;
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	mh = n >> 3;
	m = 2 * mh;
	j1 = m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[0] + a[j2];
	x0i = a[1] + a[j2 + 1];
	x1r = a[0] - a[j2];
	x1i = a[1] - a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i - x2i;
	a[j2] = x1r - x3i;
	a[j2 + 1] = x1i + x3r;
	a[j3] = x1r + x3i;
	a[j3 + 1] = x1i - x3r;
	wn4r = w[1];
	k = 0;
	for (j = 2; j < mh; j += 2) {
		k += 4;
		wk1r = w[k];
		wk1i = w[k + 1];
		wk3r = w[k + 2];
		wk3i = w[k + 3];
		j1 = j + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j] + a[j2];
		x0i = a[j + 1] + a[j2 + 1];
		x1r = a[j] - a[j2];
		x1i = a[j + 1] - a[j2 + 1];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		a[j] = x0r + x2r;
		a[j + 1] = x0i + x2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i - x2i;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j2] = wk1r * x0r - wk1i * x0i;
		a[j2 + 1] = wk1r * x0i + wk1i * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j3] = wk3r * x0r + wk3i * x0i;
		a[j3 + 1] = wk3r * x0i - wk3i * x0r;
		j0 = m - j;
		j1 = j0 + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j0] + a[j2];
		x0i = a[j0 + 1] + a[j2 + 1];
		x1r = a[j0] - a[j2];
		x1i = a[j0 + 1] - a[j2 + 1];
		x2r = a[j1] + a[j3];
		x2i = a[j1 + 1] + a[j3 + 1];
		x3r = a[j1] - a[j3];
		x3i = a[j1 + 1] - a[j3 + 1];
		a[j0] = x0r + x2r;
		a[j0 + 1] = x0i + x2i;
		a[j1] = x0r - x2r;
		a[j1 + 1] = x0i - x2i;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j2] = wk1i * x0r - wk1r * x0i;
		a[j2 + 1] = wk1i * x0i + wk1r * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j3] = wk3i * x0r + wk3r * x0i;
		a[j3 + 1] = wk3i * x0i - wk3r * x0r;
	}
	j0 = mh;
	j1 = j0 + m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[j0] + a[j2];
	x0i = a[j0 + 1] + a[j2 + 1];
	x1r = a[j0] - a[j2];
	x1i = a[j0 + 1] - a[j2 + 1];
	x2r = a[j1] + a[j3];
	x2i = a[j1 + 1] + a[j3 + 1];
	x3r = a[j1] - a[j3];
	x3i = a[j1 + 1] - a[j3 + 1];
	a[j0] = x0r + x2r;
	a[j0 + 1] = x0i + x2i;
	a[j1] = x0r - x2r;
	a[j1 + 1] = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[j2] = wn4r * (x0r - x0i);
	a[j2 + 1] = wn4r * (x0i + x0r);
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	a[j3] = -wn4r * (x0r + x0i);
	a[j3 + 1] = -wn4r * (x0i - x0r);
}


void cftmdl2(fft_int n, fft_float* a, fft_float* w)
{
	fft_int j, j0, j1, j2, j3, k, kr, m, mh;
	fft_float wn4r, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i;
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y2r, y2i;

	mh = n >> 3;
	m = 2 * mh;
	wn4r = w[1];
	j1 = m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[0] - a[j2 + 1];
	x0i = a[1] + a[j2];
	x1r = a[0] + a[j2 + 1];
	x1i = a[1] - a[j2];
	x2r = a[j1] - a[j3 + 1];
	x2i = a[j1 + 1] + a[j3];
	x3r = a[j1] + a[j3 + 1];
	x3i = a[j1 + 1] - a[j3];
	y0r = wn4r * (x2r - x2i);
	y0i = wn4r * (x2i + x2r);
	a[0] = x0r + y0r;
	a[1] = x0i + y0i;
	a[j1] = x0r - y0r;
	a[j1 + 1] = x0i - y0i;
	y0r = wn4r * (x3r - x3i);
	y0i = wn4r * (x3i + x3r);
	a[j2] = x1r - y0i;
	a[j2 + 1] = x1i + y0r;
	a[j3] = x1r + y0i;
	a[j3 + 1] = x1i - y0r;
	k = 0;
	kr = 2 * m;
	for (j = 2; j < mh; j += 2) {
		k += 4;
		wk1r = w[k];
		wk1i = w[k + 1];
		wk3r = w[k + 2];
		wk3i = w[k + 3];
		kr -= 4;
		wd1i = w[kr];
		wd1r = w[kr + 1];
		wd3i = w[kr + 2];
		wd3r = w[kr + 3];
		j1 = j + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j] - a[j2 + 1];
		x0i = a[j + 1] + a[j2];
		x1r = a[j] + a[j2 + 1];
		x1i = a[j + 1] - a[j2];
		x2r = a[j1] - a[j3 + 1];
		x2i = a[j1 + 1] + a[j3];
		x3r = a[j1] + a[j3 + 1];
		x3i = a[j1 + 1] - a[j3];
		y0r = wk1r * x0r - wk1i * x0i;
		y0i = wk1r * x0i + wk1i * x0r;
		y2r = wd1r * x2r - wd1i * x2i;
		y2i = wd1r * x2i + wd1i * x2r;
		a[j] = y0r + y2r;
		a[j + 1] = y0i + y2i;
		a[j1] = y0r - y2r;
		a[j1 + 1] = y0i - y2i;
		y0r = wk3r * x1r + wk3i * x1i;
		y0i = wk3r * x1i - wk3i * x1r;
		y2r = wd3r * x3r + wd3i * x3i;
		y2i = wd3r * x3i - wd3i * x3r;
		a[j2] = y0r + y2r;
		a[j2 + 1] = y0i + y2i;
		a[j3] = y0r - y2r;
		a[j3 + 1] = y0i - y2i;
		j0 = m - j;
		j1 = j0 + m;
		j2 = j1 + m;
		j3 = j2 + m;
		x0r = a[j0] - a[j2 + 1];
		x0i = a[j0 + 1] + a[j2];
		x1r = a[j0] + a[j2 + 1];
		x1i = a[j0 + 1] - a[j2];
		x2r = a[j1] - a[j3 + 1];
		x2i = a[j1 + 1] + a[j3];
		x3r = a[j1] + a[j3 + 1];
		x3i = a[j1 + 1] - a[j3];
		y0r = wd1i * x0r - wd1r * x0i;
		y0i = wd1i * x0i + wd1r * x0r;
		y2r = wk1i * x2r - wk1r * x2i;
		y2i = wk1i * x2i + wk1r * x2r;
		a[j0] = y0r + y2r;
		a[j0 + 1] = y0i + y2i;
		a[j1] = y0r - y2r;
		a[j1 + 1] = y0i - y2i;
		y0r = wd3i * x1r + wd3r * x1i;
		y0i = wd3i * x1i - wd3r * x1r;
		y2r = wk3i * x3r + wk3r * x3i;
		y2i = wk3i * x3i - wk3r * x3r;
		a[j2] = y0r + y2r;
		a[j2 + 1] = y0i + y2i;
		a[j3] = y0r - y2r;
		a[j3 + 1] = y0i - y2i;
	}
	wk1r = w[m];
	wk1i = w[m + 1];
	j0 = mh;
	j1 = j0 + m;
	j2 = j1 + m;
	j3 = j2 + m;
	x0r = a[j0] - a[j2 + 1];
	x0i = a[j0 + 1] + a[j2];
	x1r = a[j0] + a[j2 + 1];
	x1i = a[j0 + 1] - a[j2];
	x2r = a[j1] - a[j3 + 1];
	x2i = a[j1 + 1] + a[j3];
	x3r = a[j1] + a[j3 + 1];
	x3i = a[j1 + 1] - a[j3];
	y0r = wk1r * x0r - wk1i * x0i;
	y0i = wk1r * x0i + wk1i * x0r;
	y2r = wk1i * x2r - wk1r * x2i;
	y2i = wk1i * x2i + wk1r * x2r;
	a[j0] = y0r + y2r;
	a[j0 + 1] = y0i + y2i;
	a[j1] = y0r - y2r;
	a[j1 + 1] = y0i - y2i;
	y0r = wk1i * x1r - wk1r * x1i;
	y0i = wk1i * x1i + wk1r * x1r;
	y2r = wk1r * x3r - wk1i * x3i;
	y2i = wk1r * x3i + wk1i * x3r;
	a[j2] = y0r - y2r;
	a[j2 + 1] = y0i - y2i;
	a[j3] = y0r + y2r;
	a[j3 + 1] = y0i + y2i;
}


void cftfx41(fft_int n, fft_float* a, fft_int nw, fft_float* w)
{
	void cftf161(fft_float * a, fft_float * w);
	void cftf162(fft_float * a, fft_float * w);
	void cftf081(fft_float * a, fft_float * w);
	void cftf082(fft_float * a, fft_float * w);

	if (n == 128) {
		cftf161(a, &w[nw - 8]);
		cftf162(&a[32], &w[nw - 32]);
		cftf161(&a[64], &w[nw - 8]);
		cftf161(&a[96], &w[nw - 8]);
	}
	else {
		cftf081(a, &w[nw - 8]);
		cftf082(&a[16], &w[nw - 8]);
		cftf081(&a[32], &w[nw - 8]);
		cftf081(&a[48], &w[nw - 8]);
	}
}


void cftf161(fft_float* a, fft_float* w)
{
	fft_float wn4r, wk1r, wk1i,
		x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
		y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i,
		y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i,
		y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i;

	wn4r = w[1];
	wk1r = w[2];
	wk1i = w[3];
	x0r = a[0] + a[16];
	x0i = a[1] + a[17];
	x1r = a[0] - a[16];
	x1i = a[1] - a[17];
	x2r = a[8] + a[24];
	x2i = a[9] + a[25];
	x3r = a[8] - a[24];
	x3i = a[9] - a[25];
	y0r = x0r + x2r;
	y0i = x0i + x2i;
	y4r = x0r - x2r;
	y4i = x0i - x2i;
	y8r = x1r - x3i;
	y8i = x1i + x3r;
	y12r = x1r + x3i;
	y12i = x1i - x3r;
	x0r = a[2] + a[18];
	x0i = a[3] + a[19];
	x1r = a[2] - a[18];
	x1i = a[3] - a[19];
	x2r = a[10] + a[26];
	x2i = a[11] + a[27];
	x3r = a[10] - a[26];
	x3i = a[11] - a[27];
	y1r = x0r + x2r;
	y1i = x0i + x2i;
	y5r = x0r - x2r;
	y5i = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	y9r = wk1r * x0r - wk1i * x0i;
	y9i = wk1r * x0i + wk1i * x0r;
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	y13r = wk1i * x0r - wk1r * x0i;
	y13i = wk1i * x0i + wk1r * x0r;
	x0r = a[4] + a[20];
	x0i = a[5] + a[21];
	x1r = a[4] - a[20];
	x1i = a[5] - a[21];
	x2r = a[12] + a[28];
	x2i = a[13] + a[29];
	x3r = a[12] - a[28];
	x3i = a[13] - a[29];
	y2r = x0r + x2r;
	y2i = x0i + x2i;
	y6r = x0r - x2r;
	y6i = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	y10r = wn4r * (x0r - x0i);
	y10i = wn4r * (x0i + x0r);
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	y14r = wn4r * (x0r + x0i);
	y14i = wn4r * (x0i - x0r);
	x0r = a[6] + a[22];
	x0i = a[7] + a[23];
	x1r = a[6] - a[22];
	x1i = a[7] - a[23];
	x2r = a[14] + a[30];
	x2i = a[15] + a[31];
	x3r = a[14] - a[30];
	x3i = a[15] - a[31];
	y3r = x0r + x2r;
	y3i = x0i + x2i;
	y7r = x0r - x2r;
	y7i = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	y11r = wk1i * x0r - wk1r * x0i;
	y11i = wk1i * x0i + wk1r * x0r;
	x0r = x1r + x3i;
	x0i = x1i - x3r;
	y15r = wk1r * x0r - wk1i * x0i;
	y15i = wk1r * x0i + wk1i * x0r;
	x0r = y12r - y14r;
	x0i = y12i - y14i;
	x1r = y12r + y14r;
	x1i = y12i + y14i;
	x2r = y13r - y15r;
	x2i = y13i - y15i;
	x3r = y13r + y15r;
	x3i = y13i + y15i;
	a[24] = x0r + x2r;
	a[25] = x0i + x2i;
	a[26] = x0r - x2r;
	a[27] = x0i - x2i;
	a[28] = x1r - x3i;
	a[29] = x1i + x3r;
	a[30] = x1r + x3i;
	a[31] = x1i - x3r;
	x0r = y8r + y10r;
	x0i = y8i + y10i;
	x1r = y8r - y10r;
	x1i = y8i - y10i;
	x2r = y9r + y11r;
	x2i = y9i + y11i;
	x3r = y9r - y11r;
	x3i = y9i - y11i;
	a[16] = x0r + x2r;
	a[17] = x0i + x2i;
	a[18] = x0r - x2r;
	a[19] = x0i - x2i;
	a[20] = x1r - x3i;
	a[21] = x1i + x3r;
	a[22] = x1r + x3i;
	a[23] = x1i - x3r;
	x0r = y5r - y7i;
	x0i = y5i + y7r;
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	x0r = y5r + y7i;
	x0i = y5i - y7r;
	x3r = wn4r * (x0r - x0i);
	x3i = wn4r * (x0i + x0r);
	x0r = y4r - y6i;
	x0i = y4i + y6r;
	x1r = y4r + y6i;
	x1i = y4i - y6r;
	a[8] = x0r + x2r;
	a[9] = x0i + x2i;
	a[10] = x0r - x2r;
	a[11] = x0i - x2i;
	a[12] = x1r - x3i;
	a[13] = x1i + x3r;
	a[14] = x1r + x3i;
	a[15] = x1i - x3r;
	x0r = y0r + y2r;
	x0i = y0i + y2i;
	x1r = y0r - y2r;
	x1i = y0i - y2i;
	x2r = y1r + y3r;
	x2i = y1i + y3i;
	x3r = y1r - y3r;
	x3i = y1i - y3i;
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[2] = x0r - x2r;
	a[3] = x0i - x2i;
	a[4] = x1r - x3i;
	a[5] = x1i + x3r;
	a[6] = x1r + x3i;
	a[7] = x1i - x3r;
}


void cftf162(fft_float* a, fft_float* w)
{
	fft_float wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i,
		x0r, x0i, x1r, x1i, x2r, x2i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
		y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i,
		y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i,
		y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i;

	wn4r = w[1];
	wk1r = w[4];
	wk1i = w[5];
	wk3r = w[6];
	wk3i = -w[7];
	wk2r = w[8];
	wk2i = w[9];
	x1r = a[0] - a[17];
	x1i = a[1] + a[16];
	x0r = a[8] - a[25];
	x0i = a[9] + a[24];
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	y0r = x1r + x2r;
	y0i = x1i + x2i;
	y4r = x1r - x2r;
	y4i = x1i - x2i;
	x1r = a[0] + a[17];
	x1i = a[1] - a[16];
	x0r = a[8] + a[25];
	x0i = a[9] - a[24];
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	y8r = x1r - x2i;
	y8i = x1i + x2r;
	y12r = x1r + x2i;
	y12i = x1i - x2r;
	x0r = a[2] - a[19];
	x0i = a[3] + a[18];
	x1r = wk1r * x0r - wk1i * x0i;
	x1i = wk1r * x0i + wk1i * x0r;
	x0r = a[10] - a[27];
	x0i = a[11] + a[26];
	x2r = wk3i * x0r - wk3r * x0i;
	x2i = wk3i * x0i + wk3r * x0r;
	y1r = x1r + x2r;
	y1i = x1i + x2i;
	y5r = x1r - x2r;
	y5i = x1i - x2i;
	x0r = a[2] + a[19];
	x0i = a[3] - a[18];
	x1r = wk3r * x0r - wk3i * x0i;
	x1i = wk3r * x0i + wk3i * x0r;
	x0r = a[10] + a[27];
	x0i = a[11] - a[26];
	x2r = wk1r * x0r + wk1i * x0i;
	x2i = wk1r * x0i - wk1i * x0r;
	y9r = x1r - x2r;
	y9i = x1i - x2i;
	y13r = x1r + x2r;
	y13i = x1i + x2i;
	x0r = a[4] - a[21];
	x0i = a[5] + a[20];
	x1r = wk2r * x0r - wk2i * x0i;
	x1i = wk2r * x0i + wk2i * x0r;
	x0r = a[12] - a[29];
	x0i = a[13] + a[28];
	x2r = wk2i * x0r - wk2r * x0i;
	x2i = wk2i * x0i + wk2r * x0r;
	y2r = x1r + x2r;
	y2i = x1i + x2i;
	y6r = x1r - x2r;
	y6i = x1i - x2i;
	x0r = a[4] + a[21];
	x0i = a[5] - a[20];
	x1r = wk2i * x0r - wk2r * x0i;
	x1i = wk2i * x0i + wk2r * x0r;
	x0r = a[12] + a[29];
	x0i = a[13] - a[28];
	x2r = wk2r * x0r - wk2i * x0i;
	x2i = wk2r * x0i + wk2i * x0r;
	y10r = x1r - x2r;
	y10i = x1i - x2i;
	y14r = x1r + x2r;
	y14i = x1i + x2i;
	x0r = a[6] - a[23];
	x0i = a[7] + a[22];
	x1r = wk3r * x0r - wk3i * x0i;
	x1i = wk3r * x0i + wk3i * x0r;
	x0r = a[14] - a[31];
	x0i = a[15] + a[30];
	x2r = wk1i * x0r - wk1r * x0i;
	x2i = wk1i * x0i + wk1r * x0r;
	y3r = x1r + x2r;
	y3i = x1i + x2i;
	y7r = x1r - x2r;
	y7i = x1i - x2i;
	x0r = a[6] + a[23];
	x0i = a[7] - a[22];
	x1r = wk1i * x0r + wk1r * x0i;
	x1i = wk1i * x0i - wk1r * x0r;
	x0r = a[14] + a[31];
	x0i = a[15] - a[30];
	x2r = wk3i * x0r - wk3r * x0i;
	x2i = wk3i * x0i + wk3r * x0r;
	y11r = x1r + x2r;
	y11i = x1i + x2i;
	y15r = x1r - x2r;
	y15i = x1i - x2i;
	x1r = y0r + y2r;
	x1i = y0i + y2i;
	x2r = y1r + y3r;
	x2i = y1i + y3i;
	a[0] = x1r + x2r;
	a[1] = x1i + x2i;
	a[2] = x1r - x2r;
	a[3] = x1i - x2i;
	x1r = y0r - y2r;
	x1i = y0i - y2i;
	x2r = y1r - y3r;
	x2i = y1i - y3i;
	a[4] = x1r - x2i;
	a[5] = x1i + x2r;
	a[6] = x1r + x2i;
	a[7] = x1i - x2r;
	x1r = y4r - y6i;
	x1i = y4i + y6r;
	x0r = y5r - y7i;
	x0i = y5i + y7r;
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	a[8] = x1r + x2r;
	a[9] = x1i + x2i;
	a[10] = x1r - x2r;
	a[11] = x1i - x2i;
	x1r = y4r + y6i;
	x1i = y4i - y6r;
	x0r = y5r + y7i;
	x0i = y5i - y7r;
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	a[12] = x1r - x2i;
	a[13] = x1i + x2r;
	a[14] = x1r + x2i;
	a[15] = x1i - x2r;
	x1r = y8r + y10r;
	x1i = y8i + y10i;
	x2r = y9r - y11r;
	x2i = y9i - y11i;
	a[16] = x1r + x2r;
	a[17] = x1i + x2i;
	a[18] = x1r - x2r;
	a[19] = x1i - x2i;
	x1r = y8r - y10r;
	x1i = y8i - y10i;
	x2r = y9r + y11r;
	x2i = y9i + y11i;
	a[20] = x1r - x2i;
	a[21] = x1i + x2r;
	a[22] = x1r + x2i;
	a[23] = x1i - x2r;
	x1r = y12r - y14i;
	x1i = y12i + y14r;
	x0r = y13r + y15i;
	x0i = y13i - y15r;
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	a[24] = x1r + x2r;
	a[25] = x1i + x2i;
	a[26] = x1r - x2r;
	a[27] = x1i - x2i;
	x1r = y12r + y14i;
	x1i = y12i - y14r;
	x0r = y13r - y15i;
	x0i = y13i + y15r;
	x2r = wn4r * (x0r - x0i);
	x2i = wn4r * (x0i + x0r);
	a[28] = x1r - x2i;
	a[29] = x1i + x2r;
	a[30] = x1r + x2i;
	a[31] = x1i - x2r;
}


void cftf081(fft_float* a, fft_float* w)
{
	fft_float wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
		y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i;

	wn4r = w[1];
	x0r = a[0] + a[8];
	x0i = a[1] + a[9];
	x1r = a[0] - a[8];
	x1i = a[1] - a[9];
	x2r = a[4] + a[12];
	x2i = a[5] + a[13];
	x3r = a[4] - a[12];
	x3i = a[5] - a[13];
	y0r = x0r + x2r;
	y0i = x0i + x2i;
	y2r = x0r - x2r;
	y2i = x0i - x2i;
	y1r = x1r - x3i;
	y1i = x1i + x3r;
	y3r = x1r + x3i;
	y3i = x1i - x3r;
	x0r = a[2] + a[10];
	x0i = a[3] + a[11];
	x1r = a[2] - a[10];
	x1i = a[3] - a[11];
	x2r = a[6] + a[14];
	x2i = a[7] + a[15];
	x3r = a[6] - a[14];
	x3i = a[7] - a[15];
	y4r = x0r + x2r;
	y4i = x0i + x2i;
	y6r = x0r - x2r;
	y6i = x0i - x2i;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	x2r = x1r + x3i;
	x2i = x1i - x3r;
	y5r = wn4r * (x0r - x0i);
	y5i = wn4r * (x0r + x0i);
	y7r = wn4r * (x2r - x2i);
	y7i = wn4r * (x2r + x2i);
	a[8] = y1r + y5r;
	a[9] = y1i + y5i;
	a[10] = y1r - y5r;
	a[11] = y1i - y5i;
	a[12] = y3r - y7i;
	a[13] = y3i + y7r;
	a[14] = y3r + y7i;
	a[15] = y3i - y7r;
	a[0] = y0r + y4r;
	a[1] = y0i + y4i;
	a[2] = y0r - y4r;
	a[3] = y0i - y4i;
	a[4] = y2r - y6i;
	a[5] = y2i + y6r;
	a[6] = y2r + y6i;
	a[7] = y2i - y6r;
}


void cftf082(fft_float* a, fft_float* w)
{
	fft_float wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i,
		y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
		y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i;

	wn4r = w[1];
	wk1r = w[2];
	wk1i = w[3];
	y0r = a[0] - a[9];
	y0i = a[1] + a[8];
	y1r = a[0] + a[9];
	y1i = a[1] - a[8];
	x0r = a[4] - a[13];
	x0i = a[5] + a[12];
	y2r = wn4r * (x0r - x0i);
	y2i = wn4r * (x0i + x0r);
	x0r = a[4] + a[13];
	x0i = a[5] - a[12];
	y3r = wn4r * (x0r - x0i);
	y3i = wn4r * (x0i + x0r);
	x0r = a[2] - a[11];
	x0i = a[3] + a[10];
	y4r = wk1r * x0r - wk1i * x0i;
	y4i = wk1r * x0i + wk1i * x0r;
	x0r = a[2] + a[11];
	x0i = a[3] - a[10];
	y5r = wk1i * x0r - wk1r * x0i;
	y5i = wk1i * x0i + wk1r * x0r;
	x0r = a[6] - a[15];
	x0i = a[7] + a[14];
	y6r = wk1i * x0r - wk1r * x0i;
	y6i = wk1i * x0i + wk1r * x0r;
	x0r = a[6] + a[15];
	x0i = a[7] - a[14];
	y7r = wk1r * x0r - wk1i * x0i;
	y7i = wk1r * x0i + wk1i * x0r;
	x0r = y0r + y2r;
	x0i = y0i + y2i;
	x1r = y4r + y6r;
	x1i = y4i + y6i;
	a[0] = x0r + x1r;
	a[1] = x0i + x1i;
	a[2] = x0r - x1r;
	a[3] = x0i - x1i;
	x0r = y0r - y2r;
	x0i = y0i - y2i;
	x1r = y4r - y6r;
	x1i = y4i - y6i;
	a[4] = x0r - x1i;
	a[5] = x0i + x1r;
	a[6] = x0r + x1i;
	a[7] = x0i - x1r;
	x0r = y1r - y3i;
	x0i = y1i + y3r;
	x1r = y5r - y7r;
	x1i = y5i - y7i;
	a[8] = x0r + x1r;
	a[9] = x0i + x1i;
	a[10] = x0r - x1r;
	a[11] = x0i - x1i;
	x0r = y1r + y3i;
	x0i = y1i - y3r;
	x1r = y5r + y7r;
	x1i = y5i + y7i;
	a[12] = x0r - x1i;
	a[13] = x0i + x1r;
	a[14] = x0r + x1i;
	a[15] = x0i - x1r;
}


void cftf040(fft_float* a)
{
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	x0r = a[0] + a[4];
	x0i = a[1] + a[5];
	x1r = a[0] - a[4];
	x1i = a[1] - a[5];
	x2r = a[2] + a[6];
	x2i = a[3] + a[7];
	x3r = a[2] - a[6];
	x3i = a[3] - a[7];
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[2] = x1r - x3i;
	a[3] = x1i + x3r;
	a[4] = x0r - x2r;
	a[5] = x0i - x2i;
	a[6] = x1r + x3i;
	a[7] = x1i - x3r;
}


void cftb040(fft_float* a)
{
	fft_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	x0r = a[0] + a[4];
	x0i = a[1] + a[5];
	x1r = a[0] - a[4];
	x1i = a[1] - a[5];
	x2r = a[2] + a[6];
	x2i = a[3] + a[7];
	x3r = a[2] - a[6];
	x3i = a[3] - a[7];
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[2] = x1r + x3i;
	a[3] = x1i - x3r;
	a[4] = x0r - x2r;
	a[5] = x0i - x2i;
	a[6] = x1r - x3i;
	a[7] = x1i + x3r;
}


void cftx020(fft_float* a)
{
	fft_float x0r, x0i;

	x0r = a[0] - a[2];
	x0i = a[1] - a[3];
	a[0] += a[2];
	a[1] += a[3];
	a[2] = x0r;
	a[3] = x0i;
}


void rftfsub(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi;

	m = n >> 1;
	ks = 2 * nc / m;
	kk = 0;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		a[j] -= yr;
		a[j + 1] -= yi;
		a[k] += yr;
		a[k + 1] -= yi;
	}
}


void rftbsub(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr, xi, yr, yi;

	m = n >> 1;
	ks = 2 * nc / m;
	kk = 0;
	for (j = 2; j < m; j += 2) {
		k = n - j;
		kk += ks;
		wkr = FC_HALF - c[nc - kk];
		wki = c[kk];
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		a[j] -= yr;
		a[j + 1] -= yi;
		a[k] += yr;
		a[k + 1] -= yi;
	}
}


void dctsub(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr;

	m = n >> 1;
	ks = nc / n;
	kk = 0;
	for (j = 1; j < m; j++) {
		k = n - j;
		kk += ks;
		wkr = c[kk] - c[nc - kk];
		wki = c[kk] + c[nc - kk];
		xr = wki * a[j] - wkr * a[k];
		a[j] = wkr * a[j] + wki * a[k];
		a[k] = xr;
	}
	a[m] *= c[0];
}


void dstsub(fft_int n, fft_float* a, fft_int nc, fft_float* c)
{
	fft_int j, k, kk, ks, m;
	fft_float wkr, wki, xr;

	m = n >> 1;
	ks = nc / n;
	kk = 0;
	for (j = 1; j < m; j++) {
		k = n - j;
		kk += ks;
		wkr = c[kk] - c[nc - kk];
		wki = c[kk] + c[nc - kk];
		xr = wki * a[k] - wkr * a[j];
		a[k] = wkr * a[k] + wki * a[j];
		a[j] = xr;
	}
	a[m] *= c[0];
}

