#include"splitting_filter.h"
// ====== splitting filter ======//
mSplittingFilter * SplittingFilter_Create(size_t sample_rate) {

	mSplittingFilter*handles = (mSplittingFilter*)malloc(sizeof(mSplittingFilter));
	if (handles != NULL)
	{
		handles->sample_rate_ = sample_rate;
		if (handles->sample_rate_ == 32000)
		{
			// Two band mode
			memset(handles->two_band_analysis_state, 0, sizeof(handles->two_band_analysis_state));
			memset(handles->two_band_synthesis_state, 0, sizeof(handles->two_band_synthesis_state));
		}
		else if (handles->sample_rate_ == 48000)
		{
			// Three band mode.
			handles->three_band_filter_48k = ThreeBandFilter_Create(480);
			if (handles->three_band_filter_48k == NULL)
			{
				printf("Three-band splitting filter initialize fail \n");
				free(handles);
				return NULL;
			}
			memset(handles->data_f32, 0, sizeof(handles->data_f32));
			memset(handles->data_s16, 0, sizeof(handles->data_s16));
			for (size_t kBands = 0; kBands < 3; kBands++)
			{
				handles->three_band_f32[kBands] = (float*)malloc(160 * sizeof(float));
				handles->three_band_s16[kBands] = (int16_t*)malloc(160 * sizeof(int16_t));
			}
		}
		else
		{
			printf("Only support 32khz or 48 khz\n");
			free(handles);
			return NULL;
		}
	}

	return handles;
}

int SplittingFilter_Analysis_s16(mSplittingFilter * handles, const int16_t * data, int16_t * const * bands) {
	if (handles->sample_rate_ == 32000)
	{
		WebRtcSpl_AnalysisQMF(
			data,
			320,
			bands[0], bands[1],
			handles->two_band_analysis_state[0], handles->two_band_analysis_state[1]);
	}
	else
	{

		s16_to_f32(data, 480, handles->data_f32);
		ThreeBandFilter_Analysis(handles->three_band_filter_48k, handles->data_f32, 480, handles->three_band_f32);
		for (size_t kBands = 0; kBands < 3; kBands++)
		{
			f32_to_s16(handles->three_band_f32[kBands], 160, bands[kBands]);
		}
	}
	return 0;
}

int SplittingFilter_Synthesis_s16(mSplittingFilter * handles, const int16_t * const * bands, int16_t * data) {
	if (handles->sample_rate_ == 32000)
	{
		WebRtcSpl_SynthesisQMF(
			bands[0], bands[1],
			160,
			data,
			handles->two_band_synthesis_state[0], handles->two_band_synthesis_state[1]);
	}
	else {
		for (size_t kBands = 0; kBands < 3; kBands++)
		{
			s16_to_f32(bands[kBands], 160, handles->three_band_f32[kBands]);
		}
		ThreeBandFilter_Synthesis(handles->three_band_filter_48k, handles->three_band_f32, 160, handles->data_f32);
		f32_to_s16(handles->data_f32, 480, data);
	}
	return 0;
}

int SplittingFilter_Destory(mSplittingFilter * handles) {
	if (handles != NULL)
	{
		if (handles->sample_rate_ == 48000)
		{
			if (handles->three_band_filter_48k != NULL)
			{
				free(handles->three_band_filter_48k);
				for (size_t kBands = 0; kBands < 3; kBands++)
				{
					free(handles->three_band_f32[kBands]);
					free(handles->three_band_s16[kBands]);
				}
				return 0;
			}

		}
	}
	return -1;
}
void f32_to_s16(const float * pIn, size_t sampleCount, int16_t * pOut)
{
	int r;
	for (size_t i = 0; i < sampleCount; ++i) {
		float x = pIn[i];
		float c;
		c = ((x < -1) ? -1 : ((x > 1) ? 1 : x));
		c = c + 1;
		r = (int)(c * 32767.5f);
		r = r - 32768;
		pOut[i] = (short)r;
	}
}

void s16_to_f32(const int16_t * pIn, size_t sampleCount, float * pOut)
{
	if (pOut == NULL || pIn == NULL) {
		return;
	}

	for (size_t i = 0; i < sampleCount; ++i) {
		*pOut++ = pIn[i] / 32768.0f;
	}
}


//====== Three-band filter======//
float kLowpassCoeffs[kNumBands * kSparsity][kNumCoeffs] = {
	{ -0.00047749f, -0.00496888f, +0.16547118f, +0.00425496f },
{ -0.00173287f, -0.01585778f, +0.14989004f, +0.00994113f },
{ -0.00304815f, -0.02536082f, +0.12154542f, +0.01157993f },
{ -0.00383509f, -0.02982767f, +0.08543175f, +0.00983212f },
{ -0.00346946f, -0.02587886f, +0.04760441f, +0.00607594f },
{ -0.00154717f, -0.01136076f, +0.01387458f, +0.00186353f },
{ +0.00186353f, +0.01387458f, -0.01136076f, -0.00154717f },
{ +0.00607594f, +0.04760441f, -0.02587886f, -0.00346946f },
{ +0.00983212f, +0.08543175f, -0.02982767f, -0.00383509f },
{ +0.01157993f, +0.12154542f, -0.02536082f, -0.00304815f },
{ +0.00994113f, +0.14989004f, -0.01585778f, -0.00173287f },
{ +0.00425496f, +0.16547118f, -0.00496888f, -0.00047749f } };

int ThreeBandFilter_Init(ThreeBandFilter * handles, size_t length) {

	handles->buffer_len_ = length / kNumBands;										// 480 / 3 = 160
	handles->in_buffer_ = (float *)malloc(sizeof(float) * handles->buffer_len_);
	handles->out_buffer_ = (float*)malloc(sizeof(float)*handles->buffer_len_);

	// sparsity FIR LPF init
	for (size_t i = 0; i < kSparsity; ++i) {
		for (size_t j = 0; j < kNumBands; ++j) {
			SparseFIRFilter_Init(handles->pAnalysis[i * kNumBands + j], kLowpassCoeffs[i * kNumBands + j], kNumCoeffs, kSparsity, i);
			SparseFIRFilter_Init(handles->pSynthesis[i * kNumBands + j], kLowpassCoeffs[i * kNumBands + j], kNumCoeffs, kSparsity, i);
		}
	}

	// 12 * 3
	for (size_t i = 0; i < kSparsity*kNumBands; ++i) {
		for (size_t j = 0; j < kNumBands; ++j) {
			handles->dct_modulation_[i][j] =
				2.f * cos(2.f * M_PI * i * (2.f * j + 1.f) / (float)(kSparsity*kNumBands));
		}
	}
	return 0;
}

ThreeBandFilter * ThreeBandFilter_Create(size_t length) {

	// create and malloc memory
	ThreeBandFilter *handles = (ThreeBandFilter*)malloc(sizeof(ThreeBandFilter));
	if (handles != NULL)
	{
		for (size_t i = 0; i < kSparsity*kNumBands; i++)
		{
			handles->pAnalysis[i] = (mSparseFIRFilter *)malloc(sizeof(mSparseFIRFilter));
			handles->pSynthesis[i] = (mSparseFIRFilter *)malloc(sizeof(mSparseFIRFilter));
		}
		ThreeBandFilter_Init(handles, length);
	}

	return handles;
}

// Downsamples |in| into |out|, taking one every |kNumbands| starting from
// |offset|. |split_length| is the |out| length. |in| has to be at least
// |kNumBands| * |split_length| long.
void mDownsample(const float * in, size_t split_length, size_t offset, float * out) {
	for (size_t i = 0; i < split_length; ++i) {
		out[i] = in[kNumBands * i + offset];
	}
}

// Upsamples |in| into |out|, scaling by |kNumBands| and accumulating it every
// |kNumBands| starting from |offset|. |split_length| is the |in| length. |out|
// has to be at least |kNumBands| * |split_length| long.
void mUpsample(const float * in, size_t split_length, size_t offset, float * out) {
	for (size_t i = 0; i < split_length; ++i) {
		out[kNumBands * i + offset] += kNumBands * in[i];
	}
}

void ThreeBandFilter_DownModulate(ThreeBandFilter * handles, const float * in, size_t split_length, size_t offset, float * const * out) {
	for (size_t i = 0; i < kNumBands; ++i) {
		for (size_t j = 0; j < split_length; ++j) {
			out[i][j] += handles->dct_modulation_[offset][i] * in[j];
		}
	}
}

// Modulates each of the |kNumBands| bands of |in| by |dct_modulation_| and
// accumulates them in |out|. |out| is cleared before starting to accumulate.
// |offset| is the index in the period of the cosines used for modulation.
// |split_length| is the length of each band of |in| and |out|.
void ThreeBandFilter_UpModulate(ThreeBandFilter * handles, const float * const * in, size_t split_length, size_t offset, float * out) {
	memset(out, 0, split_length * sizeof(*out));
	for (size_t i = 0; i < kNumBands; ++i) {
		for (size_t j = 0; j < split_length; ++j) {
			out[j] += handles->dct_modulation_[offset][i] * in[i][j];
		}
	}
}

void ThreeBandFilter_Analysis(ThreeBandFilter * handles, const float * in, size_t length, float * const * out) {
	//RTC_CHECK_EQ(in_buffer_.size(), rtc::CheckedDivExact(length, kNumBands));
	for (size_t i = 0; i < kNumBands; ++i) {
		memset(out[i], 0, handles->buffer_len_ * sizeof(*out[i]));
	}
	for (size_t i = 0; i < kNumBands; ++i) {
		mDownsample(in, handles->buffer_len_, kNumBands - i - 1, handles->in_buffer_);
		for (size_t j = 0; j < kSparsity; ++j) {
			const size_t offset = i + j * kNumBands;
			SparseFIRFilter_Filter(
				handles->pAnalysis[offset],
				handles->in_buffer_,
				handles->buffer_len_,
				handles->out_buffer_);
			ThreeBandFilter_DownModulate(handles, handles->out_buffer_, handles->buffer_len_, offset, out);
		}
	}
}

void ThreeBandFilter_Synthesis(ThreeBandFilter * handles, const float * const * in, size_t split_length, float * out) {
	//RTC_CHECK_EQ(in_buffer_.size(), split_length);
	memset(out, 0, kNumBands * handles->buffer_len_ * sizeof(*out));
	for (size_t i = 0; i < kNumBands; ++i) {
		for (size_t j = 0; j < kSparsity; ++j) {
			const size_t offset = i + j * kNumBands;
			ThreeBandFilter_UpModulate(handles, in, handles->buffer_len_, offset, handles->in_buffer_);
			SparseFIRFilter_Filter(
				handles->pSynthesis[offset],
				handles->in_buffer_,
				handles->buffer_len_,
				handles->out_buffer_);
			mUpsample(handles->out_buffer_, handles->buffer_len_, i, out);
		}
	}
}

int ThreeBandFilter_Destory(ThreeBandFilter * handles) {

	if (handles != NULL)
	{
		free(handles->in_buffer_);
		free(handles->out_buffer_);
		for (size_t i = 0; i < kNumBands * kSparsity; i++)
		{
			SparseFIRFilter_Destory(handles->pAnalysis[i]);
			SparseFIRFilter_Destory(handles->pSynthesis[i]);
		}
		return 0;
	}
	return -1;
}

//==
int SparseFIRFilter_Init(mSparseFIRFilter * handles, const float * nonzero_coeffs, size_t num_nonzero_coeffs, size_t sparsity, size_t offset) {

	if (num_nonzero_coeffs<1 || sparsity < 1)
	{
		return -1;
	}
	handles->sparsity_ = sparsity;
	handles->offset_ = offset;

	handles->nozero_coeffs_len_ = num_nonzero_coeffs;
	handles->nonzero_coeffs_ = (float*)malloc(sizeof(float)*handles->nozero_coeffs_len_);
	memmove(handles->nonzero_coeffs_, nonzero_coeffs, handles->nozero_coeffs_len_ * sizeof(*nonzero_coeffs));

	handles->state_len_ = handles->sparsity_ * (num_nonzero_coeffs - 1) + handles->offset_;
	handles->state_ = (float *)malloc(sizeof(float)*(handles->state_len_));
	memset(handles->state_, 0, sizeof(float)*(handles->state_len_));


	return 0;
}

void SparseFIRFilter_Filter(mSparseFIRFilter * handles, const float * in, size_t length, float * out) {

	for (size_t i = 0; i < length; ++i) {
		out[i] = 0.f;
		size_t j;
		for (j = 0; i >= j * handles->sparsity_ + handles->offset_ && j < handles->nozero_coeffs_len_;
			++j) {
			out[i] += in[i - j * handles->sparsity_ - handles->offset_] * handles->nonzero_coeffs_[j];
		}
		for (; j < handles->nozero_coeffs_len_; ++j) {
			out[i] += handles->state_[i + (handles->nozero_coeffs_len_ - j - 1) * handles->sparsity_] *
				handles->nonzero_coeffs_[j];
		}
	}

	// Update current state.
	if (handles->state_len_ > 0u) {
		if (length >= handles->state_len_) {
			memcpy(handles->state_, &in[length - handles->state_len_],
				handles->state_len_ * sizeof(*in));
		}
		else {
			memmove(handles->state_, handles->state_ + length,
				(handles->state_len_ - length) * sizeof(handles->state_[0]));
			memcpy(handles->state_ + handles->state_len_ - length, in, length * sizeof(*in));
		}
	}
}

int SparseFIRFilter_Destory(mSparseFIRFilter * handles) {

	if (handles != NULL)
	{
		if (handles->nonzero_coeffs_ != NULL)
		{
			free(handles->nonzero_coeffs_);
		}
		else
		{
			return -1;
		}
		if (handles->state_ != NULL)
		{
			free(handles->state_);
		}
		else
		{
			return -1;
		}
		return 0;
	}
	return -1;
}


//====== Two-band QMF filter ======//

// Maximum number of samples in a low/high-band frame.
enum
{
	kMaxBandFrameLength = 320  // 10 ms at 64 kHz.
};

// QMF filter coefficients in Q16.
static const uint16_t WebRtcSpl_kAllPassFilter1[3] = { 6418, 36982, 57261 };
static const uint16_t WebRtcSpl_kAllPassFilter2[3] = { 21333, 49062, 63010 };

///////////////////////////////////////////////////////////////////////////////////////////////
// WebRtcSpl_AllPassQMF(...)
//
// Allpass filter used by the analysis and synthesis parts of the QMF filter.
//
// Input:
//    - in_data             : Input data sequence (Q10)
//    - data_length         : Length of data sequence (>2)
//    - filter_coefficients : Filter coefficients (length 3, Q16)
//
// Input & Output:
//    - filter_state        : Filter state (length 6, Q10).
//
// Output:
//    - out_data            : Output data sequence (Q10), length equal to
//                            |data_length|
//

void WebRtcSpl_AllPassQMF(int32_t* in_data, size_t data_length,
	int32_t* out_data, const uint16_t* filter_coefficients,
	int32_t* filter_state)
{
	// The procedure is to filter the input with three first order all pass filters
	// (cascade operations).
	//
	//         a_3 + q^-1    a_2 + q^-1    a_1 + q^-1
	// y[n] =  -----------   -----------   -----------   x[n]
	//         1 + a_3q^-1   1 + a_2q^-1   1 + a_1q^-1
	//
	// The input vector |filter_coefficients| includes these three filter coefficients.
	// The filter state contains the in_data state, in_data[-1], followed by
	// the out_data state, out_data[-1]. This is repeated for each cascade.
	// The first cascade filter will filter the |in_data| and store the output in
	// |out_data|. The second will the take the |out_data| as input and make an
	// intermediate storage in |in_data|, to save memory. The third, and final, cascade
	// filter operation takes the |in_data| (which is the output from the previous cascade
	// filter) and store the output in |out_data|.
	// Note that the input vector values are changed during the process.
	size_t k;
	int32_t diff;
	// First all-pass cascade; filter from in_data to out_data.

	// Let y_i[n] indicate the output of cascade filter i (with filter coefficient a_i) at
	// vector position n. Then the final output will be y[n] = y_3[n]

	// First loop, use the states stored in memory.
	// "diff" should be safe from wrap around since max values are 2^25
	// diff = (x[0] - y_1[-1])
	diff = WebRtcSpl_SubSatW32(in_data[0], filter_state[1]);
	// y_1[0] =  x[-1] + a_1 * (x[0] - y_1[-1])
	out_data[0] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[0], diff, filter_state[0]);

	// For the remaining loops, use previous values.
	for (k = 1; k < data_length; k++)
	{
		// diff = (x[n] - y_1[n-1])
		diff = WebRtcSpl_SubSatW32(in_data[k], out_data[k - 1]);
		// y_1[n] =  x[n-1] + a_1 * (x[n] - y_1[n-1])
		out_data[k] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[0], diff, in_data[k - 1]);
	}

	// Update states.
	filter_state[0] = in_data[data_length - 1]; // x[N-1], becomes x[-1] next time
	filter_state[1] = out_data[data_length - 1]; // y_1[N-1], becomes y_1[-1] next time

												 // Second all-pass cascade; filter from out_data to in_data.
												 // diff = (y_1[0] - y_2[-1])
	diff = WebRtcSpl_SubSatW32(out_data[0], filter_state[3]);
	// y_2[0] =  y_1[-1] + a_2 * (y_1[0] - y_2[-1])
	in_data[0] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[1], diff, filter_state[2]);
	for (k = 1; k < data_length; k++)
	{
		// diff = (y_1[n] - y_2[n-1])
		diff = WebRtcSpl_SubSatW32(out_data[k], in_data[k - 1]);
		// y_2[0] =  y_1[-1] + a_2 * (y_1[0] - y_2[-1])
		in_data[k] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[1], diff, out_data[k - 1]);
	}

	filter_state[2] = out_data[data_length - 1]; // y_1[N-1], becomes y_1[-1] next time
	filter_state[3] = in_data[data_length - 1]; // y_2[N-1], becomes y_2[-1] next time

												// Third all-pass cascade; filter from in_data to out_data.
												// diff = (y_2[0] - y[-1])
	diff = WebRtcSpl_SubSatW32(in_data[0], filter_state[5]);
	// y[0] =  y_2[-1] + a_3 * (y_2[0] - y[-1])
	out_data[0] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[2], diff, filter_state[4]);
	for (k = 1; k < data_length; k++)
	{
		// diff = (y_2[n] - y[n-1])
		diff = WebRtcSpl_SubSatW32(in_data[k], out_data[k - 1]);
		// y[n] =  y_2[n-1] + a_3 * (y_2[n] - y[n-1])
		out_data[k] = WEBRTC_SPL_SCALEDIFF32(filter_coefficients[2], diff, in_data[k - 1]);
	}
	filter_state[4] = in_data[data_length - 1]; // y_2[N-1], becomes y_2[-1] next time
	filter_state[5] = out_data[data_length - 1]; // y[N-1], becomes y[-1] next time
}

void WebRtcSpl_AnalysisQMF(const int16_t* in_data, size_t in_data_length,
	int16_t* low_band, int16_t* high_band,
	int32_t* filter_state1, int32_t* filter_state2)
{
	size_t i;
	int16_t k;
	int32_t tmp;
	int32_t half_in1[kMaxBandFrameLength];
	int32_t half_in2[kMaxBandFrameLength];
	int32_t filter1[kMaxBandFrameLength];
	int32_t filter2[kMaxBandFrameLength];
	const size_t band_length = in_data_length / 2;
	//RTC_DCHECK_EQ(0, in_data_length % 2);
	//RTC_DCHECK_LE(band_length, kMaxBandFrameLength);

	// Split even and odd samples. Also shift them to Q10.
	for (i = 0, k = 0; i < band_length; i++, k += 2)
	{
		half_in2[i] = ((int32_t)in_data[k]) * (1 << 10);
		half_in1[i] = ((int32_t)in_data[k + 1]) * (1 << 10);
	}

	// All pass filter even and odd samples, independently.
	WebRtcSpl_AllPassQMF(half_in1, band_length, filter1,
		WebRtcSpl_kAllPassFilter1, filter_state1);
	WebRtcSpl_AllPassQMF(half_in2, band_length, filter2,
		WebRtcSpl_kAllPassFilter2, filter_state2);

	// Take the sum and difference of filtered version of odd and even
	// branches to get upper & lower band.
	for (i = 0; i < band_length; i++)
	{
		tmp = (filter1[i] + filter2[i] + 1024) >> 11;
		low_band[i] = WebRtcSpl_SatW32ToW16(tmp);

		tmp = (filter1[i] - filter2[i] + 1024) >> 11;
		high_band[i] = WebRtcSpl_SatW32ToW16(tmp);
	}
}

void WebRtcSpl_SynthesisQMF(const int16_t* low_band, const int16_t* high_band,
	size_t band_length, int16_t* out_data,
	int32_t* filter_state1, int32_t* filter_state2)
{
	int32_t tmp;
	int32_t half_in1[kMaxBandFrameLength];
	int32_t half_in2[kMaxBandFrameLength];
	int32_t filter1[kMaxBandFrameLength];
	int32_t filter2[kMaxBandFrameLength];
	size_t i;
	int16_t k;
	//RTC_DCHECK_LE(band_length, kMaxBandFrameLength);

	// Obtain the sum and difference channels out of upper and lower-band channels.
	// Also shift to Q10 domain.
	for (i = 0; i < band_length; i++)
	{
		tmp = (int32_t)low_band[i] + (int32_t)high_band[i];
		half_in1[i] = tmp * (1 << 10);
		tmp = (int32_t)low_band[i] - (int32_t)high_band[i];
		half_in2[i] = tmp * (1 << 10);
	}

	// all-pass filter the sum and difference channels
	WebRtcSpl_AllPassQMF(half_in1, band_length, filter1,
		WebRtcSpl_kAllPassFilter2, filter_state1);
	WebRtcSpl_AllPassQMF(half_in2, band_length, filter2,
		WebRtcSpl_kAllPassFilter1, filter_state2);

	// The filtered signals are even and odd samples of the output. Combine
	// them. The signals are Q10 should shift them back to Q0 and take care of
	// saturation.
	for (i = 0, k = 0; i < band_length; i++)
	{
		tmp = (filter2[i] + 512) >> 10;
		out_data[k++] = WebRtcSpl_SatW32ToW16(tmp);

		tmp = (filter1[i] + 512) >> 10;
		out_data[k++] = WebRtcSpl_SatW32ToW16(tmp);
	}

}

//====== spl_inl ======//
// Table used by WebRtcSpl_CountLeadingZeros32_NotBuiltin. For each uint32_t n
// that's a sequence of 0 bits followed by a sequence of 1 bits, the entry at
// index (n * 0x8c0b2891) >> 26 in this table gives the number of zero bits in
// n.
const int8_t kWebRtcSpl_CountLeadingZeros32_Table[64] = {
	32, 8,  17, -1, -1, 14, -1, -1, -1, 20, -1, -1, -1, 28, -1, 18,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,  26, 25, 24,
	4,  11, 23, 31, 3,  7,  10, 16, 22, 30, -1, -1, 2,  6,  13, 9,
	-1, 15, -1, 21, -1, 29, 19, -1, -1, -1, -1, -1, 1,  27, 5,  12,
};
