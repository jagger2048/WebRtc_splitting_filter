#pragma once
#include"three_band_filter_bank.h"
#include<iostream>
using namespace std;

/* Test case */
int three_band_filter_cpp(const float* in, size_t length, float* const* out) {
	ThreeBandFilterBank T1(480);

	T1.Analysis(in, length, out);


}