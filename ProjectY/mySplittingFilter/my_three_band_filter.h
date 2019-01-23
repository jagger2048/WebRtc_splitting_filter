#pragma once
#ifndef MY_THREE_BAND_FILTER_H
#define MY_THREE_BAND_FILTER_H
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
#ifndef M_PI
#define M_PI 3.141592654
#endif // !M_PI
	//const size_t kNumBands = 3;
	//const size_t kSparsity = 4;
	//const size_t kNumCoeffs = 4;
#define kNumBands 3
#define kSparsity 4
#define kNumCoeffs 4

// ===
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

	// ==
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


#ifdef __cplusplus
}
#endif // __cplusplus

#endif // !MY_THREE_BAND_FILTER_H

