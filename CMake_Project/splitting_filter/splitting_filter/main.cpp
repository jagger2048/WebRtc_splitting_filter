//
#define DR_WAV_IMPLEMENTATION	// for the dr_wav.h
#include "wavfile.h"
//#include "splitting_filter.h"	
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char * argv[])
{
	// Read the wav file. 
	// test the wavfile library in c mode.
	wav wavfile;
	wavread("C:\\workspace\\Repository\\ProjectY\\CMake_Project\\splitting_filter\\DEBUG\\sweep-32k.wav", &wavfile);
	//FILE *fp;
	//fopen_s(&fp, "C:\\workspace\\Repository\\ProjectY\\CMake_Project\\splitting_filter\\DEBUG\\sweep-32k.wav", "rb");


	wav_destory(&wavfile);
	printf("destoried/n");
	return 0;
}