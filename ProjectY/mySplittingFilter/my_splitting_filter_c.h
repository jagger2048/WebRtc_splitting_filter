#pragma once
#ifndef MY_SPLITTING_FILTER_H
#define MY_SPLITTING_FILTER_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>	// memset

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include "two_band_QMF_filter.h"		// two-band_qmf
#include "my_three_band_filter.h"

	//#define kTwoBandStateSize 6
	//typedef struct mSplittingFilter{
	typedef struct {
		size_t sample_rate_;							// 32000,48000
														//static const int kTwoBandStateSize = 6;
														//int32_t two_band_analysis_state[2][kTwoBandStateSize];
														//int32_t two_band_synthesis_state[2][kTwoBandStateSize];
		int32_t two_band_analysis_state[2][6];
		int32_t two_band_synthesis_state[2][6];

		ThreeBandFilter* three_band_filter_48k;
		float data_f32[480];
		int16_t data_s16[480];
		float *three_band_f32[3];
		int16_t *three_band_s16[3];
	}mSplittingFilter;

	mSplittingFilter* SplittingFilter_Create(size_t sample_rate);
	int SplittingFilter_Analysis_s16(mSplittingFilter *handles, const int16_t *data, int16_t* const* bands);
	int SplittingFilter_Synthesis_s16(mSplittingFilter *handles, const int16_t* const* bands, int16_t *data);
	int SplittingFilter_Destory(mSplittingFilter *handles);
	void f32_to_s16(const float* pIn, size_t sampleCount, int16_t* pOut);
	void s16_to_f32(const int16_t* pIn, size_t sampleCount, float* pOut);





#ifdef __cplusplus
}
#endif // __cpluseplus

#endif // MY_SPLITTING_FILTER_H




//#pragma once
//#ifdef MY_SPLITTING_FILTER_H
//#ifdef __cplusplus
//extern "c" {
//#endif // __cplusplus
//
//
//
//
//
//
//
//#ifdef __cpluseplus
//}
//#endif // __cpluseplus
//
//#endif // MY_SPLITTING_FILTER_H