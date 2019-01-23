#pragma once
#ifndef SPLITTING_FILTER_H
#define SPLITTING_FILTER_H


#include <stdlib.h>
#include <stdint.h>
#include <string.h>	// memset
#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.141592654
#endif // !M_PI
//const size_t kNumBands = 3;
//const size_t kSparsity = 4;
//const size_t kNumCoeffs = 4;
#define kNumBands 3
#define kSparsity 4
#define kNumCoeffs 4

#ifndef WEBRTC_SPL_SCALEDIFF32
// C + the 32 most significant bits of A * B
#define WEBRTC_SPL_SCALEDIFF32(A, B, C) \
  (C + (B >> 16) * A + (((uint32_t)(B & 0x0000FFFF) * A) >> 16))
#endif // !WEBRTC_SPL_SCALEDIFF32

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

	//=========Three-band filter=========//

	typedef struct mSparseFIRFilter
	{
		size_t sparsity_;
		size_t offset_;

		float * nonzero_coeffs_;
		size_t nozero_coeffs_len_;
		float * state_;
		size_t state_len_;
	}mSparseFIRFilter;
	int  SparseFIRFilter_Init(
		mSparseFIRFilter * handles,
		const float* nonzero_coeffs,
		size_t num_nonzero_coeffs,
		size_t sparsity,
		size_t offset);
	void SparseFIRFilter_Filter(mSparseFIRFilter *handles, const float* in, size_t length, float* out);

	int SparseFIRFilter_Destory(mSparseFIRFilter *handles);



	// Factors to take into account when choosing |kNumCoeffs|:
	//   1. Higher |kNumCoeffs|, means faster transition, which ensures less
	//      aliasing. This is especially important when there is non-linear
	//      processing between the splitting and merging.
	//   2. The delay that this filter bank introduces is
	//      |kNumBands| * |kSparsity| * |kNumCoeffs| / 2, so it increases linearly
	//      with |kNumCoeffs|.
	//   3. The computation complexity also increases linearly with |kNumCoeffs|.


	// The Matlab code to generate these |kLowpassCoeffs| is:
	//
	// N = kNumBands * kSparsity * kNumCoeffs - 1;
	// h = fir1(N, 1 / (2 * kNumBands), kaiser(N + 1, 3.5));
	// reshape(h, kNumBands * kSparsity, kNumCoeffs);
	//
	// Because the total bandwidth of the lower and higher band is double the middle
	// one (because of the spectrum parity), the low-pass prototype is half the
	// bandwidth of 1 / (2 * |kNumBands|) and is then shifted with cosine modulation
	// to the right places.
	// A Kaiser window is used because of its flexibility and the alpha is set to
	// 3.5, since that sets a stop band attenuation of 40dB ensuring a fast
	// transition.

	extern float kLowpassCoeffs[kNumBands * kSparsity][kNumCoeffs];
	//= {
	//	{ -0.00047749f, -0.00496888f, +0.16547118f, +0.00425496f },
	//{ -0.00173287f, -0.01585778f, +0.14989004f, +0.00994113f },
	//{ -0.00304815f, -0.02536082f, +0.12154542f, +0.01157993f },
	//{ -0.00383509f, -0.02982767f, +0.08543175f, +0.00983212f },
	//{ -0.00346946f, -0.02587886f, +0.04760441f, +0.00607594f },
	//{ -0.00154717f, -0.01136076f, +0.01387458f, +0.00186353f },
	//{ +0.00186353f, +0.01387458f, -0.01136076f, -0.00154717f },
	//{ +0.00607594f, +0.04760441f, -0.02587886f, -0.00346946f },
	//{ +0.00983212f, +0.08543175f, -0.02982767f, -0.00383509f },
	//{ +0.01157993f, +0.12154542f, -0.02536082f, -0.00304815f },
	//{ +0.00994113f, +0.14989004f, -0.01585778f, -0.00173287f },
	//{ +0.00425496f, +0.16547118f, -0.00496888f, -0.00047749f } };

	typedef struct
	{
		float* in_buffer_;
		size_t buffer_len_;
		float* out_buffer_;
		float dct_modulation_[kSparsity*kNumBands][kNumBands];	// [kNumBands]
		mSparseFIRFilter *pAnalysis[kNumBands * kSparsity];
		mSparseFIRFilter *pSynthesis[kNumBands * kSparsity];
	}ThreeBandFilter;

	int ThreeBandFilter_Init(ThreeBandFilter *handles, size_t length);
	ThreeBandFilter* ThreeBandFilter_Create(size_t length);
	// Downsamples |in| into |out|, taking one every |kNumbands| starting from
	// |offset|. |split_length| is the |out| length. |in| has to be at least
	// |kNumBands| * |split_length| long.
	void mDownsample(const float* in,
		size_t split_length,
		size_t offset,
		float* out);

	// Upsamples |in| into |out|, scaling by |kNumBands| and accumulating it every
	// |kNumBands| starting from |offset|. |split_length| is the |in| length. |out|
	// has to be at least |kNumBands| * |split_length| long.
	void mUpsample(const float* in, size_t split_length, size_t offset, float* out);
	// Modulates |in| by |dct_modulation_| and accumulates it in each of the
	// |kNumBands| bands of |out|. |offset| is the index in the period of the
	// cosines used for modulation. |split_length| is the length of |in| and each
	// band of |out|.

	void ThreeBandFilter_DownModulate(
		ThreeBandFilter *handles,
		const float* in,
		size_t split_length,
		size_t offset,
		float* const* out);

	// Modulates each of the |kNumBands| bands of |in| by |dct_modulation_| and
	// accumulates them in |out|. |out| is cleared before starting to accumulate.
	// |offset| is the index in the period of the cosines used for modulation.
	// |split_length| is the length of each band of |in| and |out|.
	void ThreeBandFilter_UpModulate(
		ThreeBandFilter *handles,
		const float* const* in,
		size_t split_length,
		size_t offset,
		float* out);


	void ThreeBandFilter_Analysis(ThreeBandFilter *handles, const float* in,
		size_t length,
		float* const* out);

	void ThreeBandFilter_Synthesis(ThreeBandFilter *handles, const float* const* in,
		size_t split_length,
		float* out);

	int ThreeBandFilter_Destory(ThreeBandFilter *handles);

//=========Two-band QMF filter=========//

	void WebRtcSpl_AllPassQMF(int32_t* in_data, size_t data_length,
		int32_t* out_data, const uint16_t* filter_coefficients,
		int32_t* filter_state);

	void WebRtcSpl_AnalysisQMF(const int16_t* in_data,
		size_t in_data_length,
		int16_t* low_band,
		int16_t* high_band,
		int32_t* filter_state1,
		int32_t* filter_state2);

	void WebRtcSpl_SynthesisQMF(const int16_t* low_band,
		const int16_t* high_band,
		size_t band_length,
		int16_t* out_data,
		int32_t* filter_state1,
		int32_t* filter_state2);

	//=========Splitting filter=========//
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
	

	//================SPL library=============//
	extern const int8_t kWebRtcSpl_CountLeadingZeros32_Table[64];

	// Don't call this directly except in tests!
	static __inline int WebRtcSpl_CountLeadingZeros32_NotBuiltin(uint32_t n) {
		// Normalize n by rounding up to the nearest number that is a sequence of 0
		// bits followed by a sequence of 1 bits. This number has the same number of
		// leading zeros as the original n. There are exactly 33 such values.
		n |= n >> 1;
		n |= n >> 2;
		n |= n >> 4;
		n |= n >> 8;
		n |= n >> 16;

		// Multiply the modified n with a constant selected (by exhaustive search)
		// such that each of the 33 possible values of n give a product whose 6 most
		// significant bits are unique. Then look up the answer in the table.
		return kWebRtcSpl_CountLeadingZeros32_Table[(n * 0x8c0b2891) >> 26];
}

	// Don't call this directly except in tests!
	static __inline int WebRtcSpl_CountLeadingZeros64_NotBuiltin(uint64_t n) {
		const int leading_zeros = n >> 32 == 0 ? 32 : 0;
		return leading_zeros + WebRtcSpl_CountLeadingZeros32_NotBuiltin(
			(uint32_t)(n >> (32 - leading_zeros)));
	}

	// Returns the number of leading zero bits in the argument.
	static __inline int WebRtcSpl_CountLeadingZeros32(uint32_t n) {
#ifdef __GNUC__
		RTC_COMPILE_ASSERT(sizeof(unsigned int) == sizeof(uint32_t));
		return n == 0 ? 32 : __builtin_clz(n);
#else
		return WebRtcSpl_CountLeadingZeros32_NotBuiltin(n);
#endif
	}

	// Returns the number of leading zero bits in the argument.
	static __inline int WebRtcSpl_CountLeadingZeros64(uint64_t n) {
#ifdef __GNUC__
		RTC_COMPILE_ASSERT(sizeof(unsigned long long) == sizeof(uint64_t));  // NOLINT
		return n == 0 ? 64 : __builtin_clzll(n);
#else
		return WebRtcSpl_CountLeadingZeros64_NotBuiltin(n);
#endif
	}

#ifdef WEBRTC_ARCH_ARM_V7
#include "common_audio/signal_processing/include/spl_inl_armv7.h"
#else

#if !defined(MIPS_DSP_R1_LE)
	static __inline int16_t WebRtcSpl_SatW32ToW16(int32_t value32) {
		int16_t out16 = (int16_t)value32;

		if (value32 > 32767)
			out16 = 32767;
		else if (value32 < -32768)
			out16 = -32768;

		return out16;
	}

	static __inline int32_t WebRtcSpl_AddSatW32(int32_t a, int32_t b) {
		// Do the addition in unsigned numbers, since signed overflow is undefined
		// behavior.
		const int32_t sum = (int32_t)((uint32_t)a + (uint32_t)b);

		// a + b can't overflow if a and b have different signs. If they have the
		// same sign, a + b also has the same sign iff it didn't overflow.
		if ((a < 0) == (b < 0) && (a < 0) != (sum < 0)) {
			// The direction of the overflow is obvious from the sign of a + b.
			return sum < 0 ? INT32_MAX : INT32_MIN;
		}
		return sum;
	}

	static __inline int32_t WebRtcSpl_SubSatW32(int32_t a, int32_t b) {
		// Do the subtraction in unsigned numbers, since signed overflow is undefined
		// behavior.
		const int32_t diff = (int32_t)((uint32_t)a - (uint32_t)b);

		// a - b can't overflow if a and b have the same sign. If they have different
		// signs, a - b has the same sign as a iff it didn't overflow.
		if ((a < 0) != (b < 0) && (a < 0) != (diff < 0)) {
			// The direction of the overflow is obvious from the sign of a - b.
			return diff < 0 ? INT32_MAX : INT32_MIN;
		}
		return diff;
	}

	static __inline int16_t WebRtcSpl_AddSatW16(int16_t a, int16_t b) {
		return WebRtcSpl_SatW32ToW16((int32_t)a + (int32_t)b);
	}

	static __inline int16_t WebRtcSpl_SubSatW16(int16_t var1, int16_t var2) {
		return WebRtcSpl_SatW32ToW16((int32_t)var1 - (int32_t)var2);
	}
#endif  // #if !defined(MIPS_DSP_R1_LE)

#if !defined(MIPS32_LE)
	static __inline int16_t WebRtcSpl_GetSizeInBits(uint32_t n) {
		return 32 - WebRtcSpl_CountLeadingZeros32(n);
	}

	// Return the number of steps a can be left-shifted without overflow,
	// or 0 if a == 0.
	static __inline int16_t WebRtcSpl_NormW32(int32_t a) {
		return a == 0 ? 0 : WebRtcSpl_CountLeadingZeros32(a < 0 ? ~a : a) - 1;
	}

	// Return the number of steps a can be left-shifted without overflow,
	// or 0 if a == 0.
	static __inline int16_t WebRtcSpl_NormU32(uint32_t a) {
		return a == 0 ? 0 : WebRtcSpl_CountLeadingZeros32(a);
	}

	// Return the number of steps a can be left-shifted without overflow,
	// or 0 if a == 0.
	static __inline int16_t WebRtcSpl_NormW16(int16_t a) {
		const int32_t a32 = a;
		return a == 0 ? 0 : WebRtcSpl_CountLeadingZeros32(a < 0 ? ~a32 : a32) - 17;
	}

	static __inline int32_t WebRtc_MulAccumW16(int16_t a, int16_t b, int32_t c) {
		return (a * b + c);
	}
#endif  // #if !defined(MIPS32_LE)

#endif  // WEBRTC_ARCH_ARM_V7

#ifdef __cplusplus
}
#endif // __cpluseplus

#endif // !SPLITTING_FILTER_H
