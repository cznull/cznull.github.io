
#include <iostream>
#include <opencv2/opencv.hpp>
#include <Windows.h>
#include <math.h>
#include <string>
#include <io.h>
#include <stdio.h>
#include <stdlib.h>


char cmap[20] = "+-sincox.0123456789";
int pos[19] = { 14,5,25,39,40,23,24,27,7,6,9,8,22,32,20,47,12,11,18 };

int distance(cv::Mat& a, cv::Mat& b) {
	int x = 0;
	int m;
	for (int i = 0; i < a.cols * a.rows; i++) {
		m = (int)(a.data[i * 3]) - (int)(b.data[i * 3]);
		x += m * m;
	}
	return x;
}

int readpara(const char* s, double* x) {
	int cur = 0;
	sscanf_s(s, "%lf", x);
	for (int i = 0; i < 10; i++) {
		cur++;
		for (; cur < 166; cur++) {
			if (!(s[cur] == '.' || ('0' <= s[cur] && s[cur] <= '9'))) {//|| s[cur] == '+' || s[cur] == '-' 
				break;
			}
		}
		if (i) {
			cur += 5;
		}
		sscanf_s(s + cur, "%lf", x + 2 + i);
	}
	return 0;
}

double gety(double* para, int n, double x) {
	double y = 0;
	for (int i = 0; i < n; i++) {
		y += para[i * 2] * cos(x * i);
		y += para[i * 2 + 1] * sin(x * i);
	}
	return y;
}

int drawimg(cv::Mat& img, double* para, int n) {
	memset(img.data, 255, 512 * 512 * 3);
	int y;
	for (int x = 0; x < 512; x++) {
		y = 256 - 51.2 * gety(para, n, (x + 0.5) / 512 * 6.2832);
		if (y < 0) {
			y = 0;
		}
		if (y > 511) {
			y = 511;
		}

		for (int w = -3; w < 4; w++) {
			for (int h = -3; h < 4; h++) {
				if (0 <= x + w && x + w <= 511 && 0 <= y + h && y + h <= 511) {
					if (w * w + h * h < 6) {
						img.data[((y + h) * 512 + x + w) * 3] = 0;
						img.data[((y + h) * 512 + x + w) * 3 + 1] = 0;
						img.data[((y + h) * 512 + x + w) * 3 + 2] = 0;

					}
				}
			}
		}
	}
	return 0;
}

int main()
{
	const int w = 1920, h = 1080;
	cv::VideoCapture vc("badfx.mp4");
	cv::VideoWriter writer;
	writer.open("demo23.mp4", cv::VideoWriter::fourcc('a', 'v', 'c', '1'), 30.0, cv::Size(512, 512));
	cv::Mat frame(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::Mat c(23, 11, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::Mat f(512, 512, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::Mat white, black;
	cv::Mat cimg[19];
	std::vector<char> str(167, 0);
	double para[12] = { 0.0 };

	white = cv::imread("white.png");
	black = cv::imread("black.png");
	int i;
	i = 0;
	for (i = 0; i < 481; i++) {
		vc >> frame;
	}
	for (; i < 7055; i++) {
		vc >> frame;
		if (i == 481) {
			for (int j = 0; j < 19; j++) {
				cimg[j] = cv::Mat(23, 11, CV_8UC3, cv::Scalar(0, 0, 0));
				frame(cv::Rect(pos[j] * 11.52, 40, 11, 23)).copyTo(cimg[j]);
				char s[128];
				sprintf_s(s, "i_%d.png", j);
				cv::imwrite(s, cimg[j]);
			}
		}
		for (int i = 0; i < w * h * 3; i++) {
			float a = (frame.data[i] - black.data[i]) * 255.0 / (white.data[i] - black.data[i]);
			if (a < 0) {
				a = 0;
			}
			if (a > 255) {
				a = 255;
			}
			frame.data[i] = a;
		}
		for (int y = 0; y < 1; y++) {
			for (int j = 0; j < 167; j++) {
				if (j == 166) {
					frame(cv::Rect(j * 11.52, 40 + y * 22.76, 8, 23)).copyTo(c);
				}
				else {
					frame(cv::Rect(j * 11.52, 40 + y * 22.76, 11, 23)).copyTo(c);
				}
				int min = 255 * 255 * 11 * 13, minindex = 0;
				int d;
				for (int k = 0; k < 19; k++) {
					d = distance(c, cimg[k]);
					if (d < min) {
						min = d;
						minindex = k;
					}
				}
				if (j > 4) {
					str[j - 5] = cmap[minindex];
				}
			}
		}
		readpara(str.data(), para);

		for (int i = 0; i < 12; i++) {
			printf("%f,", para[i]);
		}
		printf("%d\n", i);
		drawimg(f, para, 6);
		writer << f;
	}
	return 0;
}