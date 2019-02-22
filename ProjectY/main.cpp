//
#define DR_WAV_IMPLEMENTATION	// for the dr_wav.h
#include "wavfile.h"
#include "splitting_filter.h"	
#include <stdio.h>
#include <stdlib.h>
int main()
{
	// Read the wav file.
	WAV *wavfile = wavfile_read("sweep-32k.wav");
	size_t sample_rate = wavfile->sampleRate;
	size_t len = wavfile->totalPCMFrameCount;
	mSplittingFilter* sf = (mSplittingFilter*)SplittingFilter_Create(sample_rate);
	//
	int16_t bands_data[3][160]{};
	int16_t *pBands[3] = { bands_data[0],bands_data[1],bands_data[2] };
	int16_t *pOut = (int16_t*)malloc(len*sizeof(int16_t));

	unsigned int nFrames = sample_rate / 100;

	int16_t *pCurIn = wavfile->pDataS16[0];
	int16_t *pCurOut = pOut;

	for (size_t n = 0; n < len / nFrames; n++)
	{
		SplittingFilter_Analysis_s16(sf, pCurIn, pBands);
		// your band processing here,e.g. AGC, noise suppression.

		SplittingFilter_Synthesis_s16(sf, pBands, pCurOut);

		pCurIn += nFrames;
		pCurOut += nFrames;
	}

	wavfile_write_s16("output 333.wav", &pOut, len, 1, 32000);

	SplittingFilter_Destory(sf);
	wavfile_destory(wavfile);
	free(pOut);
	return 0;
}