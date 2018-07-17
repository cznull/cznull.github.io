#include"stdafx.h"
#include "vec.h"

template <int n>
double dot(doublen<n> x1, doublen<n> x2) {
	double sum = 0;
	int i;
	for (i = 0; i < n; i++) {
		sum += x1.x[i] * x2.x[i];
	}
	return sum;
}

template <int n>
double norm(doublen<n> x1) {
	double sum = 0;
	int i;
	for (i = 0; i < n; i++) {
		sum += x1.x[i] * x1.x[i];
	}
	return sum;
}

template <int n>
doublen<n> operator-(doublen<n> x1, doublen<n> x2) {
	doublen<n> x;
	int i;
	for (i = 0; i < n; i++) {
		x.x[i] = x1.x[i] - x1.x[i];
	}
	return x;
}

template <int n>
doublen<n> operator+(doublen<n> x1, doublen<n> x2) {
	doublen<n> x;
	int i;
	for (i = 0; i < n; i++) {
		x.x[i] = x1.x[i] + x1.x[i];
	}
	return x;
}

double dot(double2 x1, double2 x2) {
	return x1.x*x2.x + x1.y*x2.y;
}

double norm(double2 x1) {
	return x1.x*x1.x + x1.y*x1.y;
}

double vol(double2 x1, double2 x2) {
	return x1.x*x2.y - x1.y*x2.x;
}

double2 operator-(double2 x1, double2 x2) {
	return { x1.x - x2.x,x1.y - x2.y };
}

double2 operator+(double2 x1, double2 x2) {
	return { x1.x + x2.x,x1.y + x2.y };
}

double2 operator*(double s, double2 x) {
	return { s*x.x,s*x.y };
}

double2 operator*(double2 x, double s) {
	return { s*x.x,s*x.y };
}

double dot(double3 x1, double3 x2) {
	return x1.x*x2.x + x1.y*x2.y + x1.z*x2.z;
}

double norm(double3 x1) {
	return x1.x*x1.x + x1.y*x1.y + x1.z*x1.z;
}

double vol(double3 x1, double3 x2, double3 x3) {
	return x1.x*(x2.y*x3.z - x2.z*x3.y) + x1.y*(x2.z*x3.x - x2.x*x3.z) + x1.z*(x2.x*x3.y - x2.y*x3.x);
}

double3 cross(double3 x1, double3 x2) {
	return { x1.y*x2.z - x1.z*x2.y ,x1.z*x2.x - x1.x*x2.z ,x1.x*x2.y - x1.y*x2.x };
}

double3 operator-(double3 x1, double3 x2) {
	return { x1.x - x2.x,x1.y - x2.y ,x1.z - x2.z };
}

double3 operator+(double3 x1, double3 x2) {
	return { x1.x + x2.x,x1.y + x2.y ,x1.z + x2.z };
}

double3 operator*(double3 x, double s) {
	return { x.x*s,x.y*s,x.z*s };
}

float dot(float3 x1, float3 x2) {
	return x1.x*x2.x + x1.y*x2.y + x1.z*x2.z;
}

float norm(float3 x1) {
	return x1.x*x1.x + x1.y*x1.y + x1.z*x1.z;
}

float vol(float3 x1, float3 x2, float3 x3) {
	return x1.x*(x2.y*x3.z - x2.z*x3.y) + x1.y*(x2.z*x3.x - x2.x*x3.z) + x1.z*(x2.x*x3.y - x2.y*x3.x);
}

float3 cross(float3 x1, float3 x2) {
	return { x1.y*x2.z - x1.z*x2.y ,x1.z*x2.x - x1.x*x2.z ,x1.x*x2.y - x1.y*x2.x };
}

float3 operator-(float3 x1, float3 x2) {
	return { x1.x - x2.x,x1.y - x2.y ,x1.z - x2.z };
}

float3 operator+(float3 x1, float3 x2) {
	return { x1.x + x2.x,x1.y + x2.y ,x1.z + x2.z };
}