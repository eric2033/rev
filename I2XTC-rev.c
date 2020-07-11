//Eric Zhang - Interactive Individualized Crosstalk Cancellation
//Portaudio code by M. Farbood
//FFTW base implementation by O. Nieto

#include <stdlib.h>
#include <stdio.h>
#include <portaudio.h>
#include <sndfile.h>
#include <string.h>
#include <ncurses.h>
#include "inputlib.h"
#include <math.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define FRAMES_PER_BUFFER 1024
#define MONO 1
#define STEREO 2
#define DELAY_INCREMENT 1
#define AMPLITUDE_INCREMENT 0.1

//data struct
typedef struct
{
	float sampleRate;
	int delay;
	float amplitude_scale;
	SNDFILE *infile;
	SNDFILE *H30L;
	SNDFILE *H30R;
	SNDFILE *H330L;
	SNDFILE *H330R;
	SF_INFO infileinfo;
	SF_INFO H30Linfo;
	SF_INFO H30Rinfo;
	SF_INFO H330Linfo;
	SF_INFO H330Rinfo;
} paData;

//complex data struct
typedef struct
{
	double a;
	double b;
} Complex;

//global variables
unsigned flags = FFTW_MEASURE;

//wisdom file
const char *filename = "I2XTC.wis";

//global readcount
int readcount = 0;

//check for the start of the signal to read HRTFs and zero OLA buffer
bool start = TRUE;

//global pointers for prev arrays for OLA and HRTFs
float *leftOLA;
float *rightOLA;
float *H30LHRTF;
float *H30RHRTF;
float *H330LHRTF;
float *H330RHRTF;

//global fftw_complex vars
fftw_complex *fftw_H30L_in, *fftw_H30R_in, *fftw_H330L_in, *fftw_H330R_in, *fftw_H30L_out, *fftw_H30R_out, *fftw_H330L_out, *fftw_H330R_out;

//callback function declaration
static int paCallback(const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData);

//complex arithmetic declaration
Complex divide(Complex x, Complex y);
Complex multiply(Complex x, Complex y);
Complex reciprocal(Complex x);

int main(int argc, char **argv)
{
	paData data;
	PaStream *stream;
	PaError err;
	PaStreamParameters outputParams;

	//check input arguments
	if (argc != 6)
	{
		printf("\nPlease input five audio files as input arguments: \n");
		printf("	Usage: I2XTC <sigfile> <H30L> <H30R> <H330L> <H330R>\n");
		exit(1);
	}

	//open files
	if ((data.infile = sf_open(argv[1], SFM_READ, &data.infileinfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[1]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	if ((data.H30L = sf_open(argv[2], SFM_READ, &data.H30Linfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[2]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	if ((data.H30R = sf_open(argv[3], SFM_READ, &data.H30Rinfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[3]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	if ((data.H330L = sf_open(argv[4], SFM_READ, &data.H330Linfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[4]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	if ((data.H330R = sf_open(argv[5], SFM_READ, &data.H330Rinfo)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[5]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	printf("Signal audio file: Frames: %d Channels: %d Sample Rate: %d\n",
		(int)data.infileinfo.frames, data.infileinfo.channels, data.infileinfo.samplerate);

	//init struct parameters
	data.sampleRate = data.infileinfo.samplerate;
	data.delay = 0;
	data.amplitude_scale = 1;

	//initialize and start PortAudio
	err = Pa_Initialize();
	if (err != paNoError)
	{
		printf("PortAudio error: %s\n", Pa_GetErrorText(err));
		printf("\nExiting.\n");
		exit(1);
	}

	//output stream parameters
	outputParams.device = Pa_GetDefaultOutputDevice();
	outputParams.channelCount = STEREO;
	outputParams.sampleFormat = paFloat32;
	outputParams.suggestedLatency =
		Pa_GetDeviceInfo(outputParams.device)->defaultLowOutputLatency;
	outputParams.hostApiSpecificStreamInfo = NULL;

	err = Pa_OpenStream(&stream,
		NULL,
		&outputParams,
		data.sampleRate,
		FRAMES_PER_BUFFER/2,
		paNoFlag,
		paCallback,
		&data);

	if (err != paNoError)
	{
		printf("PortAudio error: open stream: %s\n", Pa_GetErrorText(err));
		exit(2);
	}

	//start audio stream
	err = Pa_StartStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: start stream: %s\n", Pa_GetErrorText(err));
		exit(3);
	}

	//interactive character input
	initscr();
	cbreak();
	noecho();

	char ch;
	ch = '\0';

	mvprintw(0, 0, "Delay: %d fftbins Amplitude: %.2f x \n \n[a/z] increases/decreases effect\n"\
		"[k/m] increases/decreases amplitude scaling \n"\
		"[q] to quit \n", data.delay, data.amplitude_scale);

	while (ch != 'q')
	{
		ch = getch();
		switch(ch){
			case 'k':
			{
				data.amplitude_scale += AMPLITUDE_INCREMENT;
				break;
			}
			case 'm':
			{
				data.amplitude_scale -= AMPLITUDE_INCREMENT;
				if (data.amplitude_scale < 0){
					data.amplitude_scale = 0;
				}
				break;
			}
			case 'a':
			{
				data.delay += DELAY_INCREMENT;
				break;
			}
			case 'z':
			{
				data.delay -= DELAY_INCREMENT;
				break;
			}
		}
		mvprintw(0, 0, "Delay: %d samples  Amplitude: %.2fx \n ", data.delay, data.amplitude_scale);
	}

	//end curses
	endwin();

	//close files
	sf_close(data.infile);
	sf_close(data.H30L);
	sf_close(data.H30R);
	sf_close(data.H330L);
	sf_close(data.H330R);

	//free global pointers
	free(leftOLA);
	free(rightOLA);
	free(H30LHRTF);
	free(H30RHRTF);
	free(H330LHRTF);
	free(H330RHRTF);
	fftw_free(fftw_H30L_in);
	fftw_free(fftw_H30R_in);
	fftw_free(fftw_H330L_in);
	fftw_free(fftw_H330R_in);
	fftw_free(fftw_H30L_out);
	fftw_free(fftw_H30R_out);
	fftw_free(fftw_H330L_out);
	fftw_free(fftw_H330R_out);

	//stop stream
	err = Pa_StopStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: stop stream: %s\n", Pa_GetErrorText(err));
	}

	//close stream
	err = Pa_CloseStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: close stream: %s\n", Pa_GetErrorText(err));
	}

	//terminate PortAudio
	err = Pa_Terminate();
	if (err != paNoError)
	{
		printf("PortAudio error: terminate: %s\n", Pa_GetErrorText(err));
	}

	return 0;
}

//callback function
static int paCallback( const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData)
{
	//Cast input params
	paData *data = (paData*)userData;
	float *outbuff = (float *)outputBuffer;
	int readcount = 0;
	float normL = 0;
	float normR = 0;

	//Complex numbers
	Complex x, y, w, z, o, p, q, r, s, t, u ,v, m, n;
	Complex one, negone;
	one.a = 1.f;
	one.b = 0.f;
	negone.a = -1;
	negone.b = 0;

	//allocate buffers and memset
	float *buffer = (float *)fftw_malloc(sizeof(float)*framesPerBuffer*STEREO);
	float *H30Lshift = (float *)fftw_malloc(sizeof(float)*data->H30Linfo.frames*MONO);
	float *H30Rshift = (float *)fftw_malloc(sizeof(float)*data->H30Rinfo.frames*MONO);
	float *H330Lshift = (float *)fftw_malloc(sizeof(float)*data->H330Linfo.frames*MONO);
	float *H330Rshift = (float *)fftw_malloc(sizeof(float)*data->H330Rinfo.frames*MONO);
	memset(H30Lshift, 0.f, sizeof(float)*data->H30Linfo.frames*MONO);
	memset(H30Rshift, 0.f, sizeof(float)*data->H30Rinfo.frames*MONO);
	memset(H330Lshift, 0.f, sizeof(float)*data->H330Linfo.frames*MONO);
	memset(H330Rshift, 0.f, sizeof(float)*data->H330Rinfo.frames*MONO);

	//read HRTFs and set starting OLA to zero
	if (start == TRUE)
	{
		H30LHRTF = (float *)fftw_malloc(sizeof(float)*data->H30Linfo.frames*MONO);
		H30RHRTF = (float *)fftw_malloc(sizeof(float)*data->H30Rinfo.frames*MONO);
		H330LHRTF = (float *)fftw_malloc(sizeof(float)*data->H330Linfo.frames*MONO);
		H330RHRTF = (float *)fftw_malloc(sizeof(float)*data->H330Rinfo.frames*MONO);
		sf_readf_float(data->H30L, H30LHRTF, data->H30Linfo.frames);
		sf_readf_float(data->H30R, H30RHRTF, data->H30Rinfo.frames);
		sf_readf_float(data->H330L, H330LHRTF, data->H330Linfo.frames);
		sf_readf_float(data->H330R, H330RHRTF, data->H330Rinfo.frames);
		leftOLA = (float *)fftw_malloc(sizeof(float)*framesPerBuffer*MONO);
		rightOLA = (float *)fftw_malloc(sizeof(float)*framesPerBuffer*MONO);
		memset(leftOLA, 0.f, sizeof(float)*framesPerBuffer);
		memset(rightOLA, 0.f, sizeof(float)*framesPerBuffer);
	}

	//1s substituted for HRTFs
	//for(int i = 0; i < data->H30Linfo.frames; i++)
	//{
	//	H30LHRTF[i] = 1.f;
	//	H30RHRTF[i] = 1.f;
	//	H330LHRTF[i] = 1.f;
	//	H330RHRTF[i] = 1.f;
	//}

	//TODO: add delay and scale amplitude of HRTFs
	if (data->delay < 0)
	{
		for(int i = 0; i < data->H30Linfo.frames+data->delay; i++)
		{
			H30Lshift[i] = H30LHRTF[i-data->delay] * (float)data->amplitude_scale;
			H30Rshift[i] = H30RHRTF[i-data->delay] * (float)data->amplitude_scale;
			H330Lshift[i] = H330LHRTF[i-data->delay] * (float)data->amplitude_scale;
			H330Rshift[i] = H330RHRTF[i-data->delay] * (float)data->amplitude_scale;
		}
		for(int i = data->H30Linfo.frames+data->delay; i < data->H30Linfo.frames; i++)
		{
			H30Lshift[i] = 0.f;
			H30Rshift[i] = 0.f;
			H330Lshift[i] = 0.f;
			H330Rshift[i] = 0.f;
		}
	}
	else if (data->delay > 0)
	{
		for(int i = 0; i < data->H30Linfo.frames-data->delay; i++)
		{
			H30Lshift[i+data->delay] = H30LHRTF[i] * (float)data->amplitude_scale;
			H30Rshift[i+data->delay] = H30RHRTF[i] * (float)data->amplitude_scale;
			H330Lshift[i+data->delay] = H330LHRTF[i] * (float)data->amplitude_scale;
			H330Rshift[i+data->delay] = H330RHRTF[i] * (float)data->amplitude_scale;
		}
		for(int i = 0; i < data->delay; i++)
		{
			H30Lshift[i] = 0.f;
			H30Rshift[i] = 0.f;
			H330Lshift[i] = 0.f;
			H330Rshift[i] = 0.f;
		}
	}
	else if (data->delay == 0)
	{
		for(int i = 0; i < data->H30Linfo.frames; i++)
		{
			H30Lshift[i] = H30LHRTF[i] * (float)data->amplitude_scale;
			H30Rshift[i] = H30RHRTF[i] * (float)data->amplitude_scale;
			H330Lshift[i] = H330LHRTF[i] * (float)data->amplitude_scale;
			H330Rshift[i] = H330RHRTF[i] * (float)data->amplitude_scale;
		}
	}

	//Read input file to buffer (interleaved) and rewind if necessary
	readcount = sf_readf_float(data->infile, buffer,framesPerBuffer);
	if (readcount < framesPerBuffer)
	{
		sf_seek(data->infile, 0, SEEK_SET);
		sf_readf_float(data->infile, buffer+readcount, framesPerBuffer-readcount);
	}

	//fftw variables
	fftw_complex *mult_fft_l, *mult_fft_r, *in_fft_l, *out_fft_l, *in_fft_r, *out_fft_r, *crosstalk_IR_L, *crosstalk_IR_R, *left_ch, *right_ch, *proc_IR_L, *proc_IR_R, *eq_IR_L, *eq_IR_R, *cross_cancel_L, *cross_cancel_R;
	fftw_plan plan_l, plan_r, plan_H30L, plan_H30R, plan_H330L, plan_H330R;

	//allocate arrays and memset
	in_fft_l = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	in_fft_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	out_fft_l = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	out_fft_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	//mult_fft_l = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*((framesPerBuffer+data->H30Linfo.frames-1)/2));
	//mult_fft_r = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*((framesPerBuffer+data->H30Linfo.frames-1)/2));
	memset(in_fft_l, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(in_fft_r, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(out_fft_l, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(out_fft_r, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	//memset(mult_fft_l, 0.f, sizeof(fftw_complex)*((framesPerBuffer+data->H30Linfo.frames-1)/2));
	//memset(mult_fft_r, 0.f, sizeof(fftw_complex)*((framesPerBuffer+data->H30Linfo.frames-1)/2));
	crosstalk_IR_L = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	crosstalk_IR_R = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	cross_cancel_L = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	cross_cancel_R = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	left_ch = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	right_ch = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	proc_IR_L = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	proc_IR_R = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	eq_IR_L = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	eq_IR_R = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(crosstalk_IR_L, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(crosstalk_IR_R, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(left_ch, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(right_ch, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(proc_IR_L, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(proc_IR_R, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(eq_IR_L, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(eq_IR_R, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	// memset(leftOLA, 0.f, sizeof(float)*framesPerBuffer);
	// memset(rightOLA, 0.f, sizeof(float)*framesPerBuffer);

	//allocate HRTF arrays
	if (start == TRUE)
	{
	fftw_H30L_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H30R_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H330L_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H330R_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H30L_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H30R_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H330L_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	fftw_H330R_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	start = FALSE;
	}

	//memset each HRTF
	memset(fftw_H30L_in, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(fftw_H30R_in, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(fftw_H330L_in, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));
	memset(fftw_H330R_in, 0.f, sizeof(fftw_complex)*(framesPerBuffer+data->H30Linfo.frames-1));

	//deinterleave input stereo file
	for (int i = 0; i < framesPerBuffer; i++)
	{
		in_fft_l[i][0] = buffer[2*i];
		in_fft_l[i][1] = 0.f;
		in_fft_r[i][0] = buffer[2*i+1];
		in_fft_r[i][1] = 0.f;
		out_fft_l[i][0] = 0.f;
		out_fft_l[i][1] = 0.f;
		out_fft_l[i][0] = 0.f;
		out_fft_l[i][1] = 0.f;
	}

	//memset in_fft_r and in_fft_l manually
	for (int i = framesPerBuffer; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	{
		in_fft_l[i][0] = 0.f;
		in_fft_l[i][1] = 0.f;
		in_fft_r[i][0] = 0.f;
		in_fft_r[i][1] = 0.f;
		out_fft_l[i][0] = 0.f;
		out_fft_l[i][1] = 0.f;
		out_fft_l[i][0] = 0.f;
		out_fft_l[i][1] = 0.f;
	}

	//fill HRTF complex buffers
	for (int i = 0; i < data->H30Linfo.frames; i++)
	{
		fftw_H30L_in[i][0] = H30Lshift[i];
		fftw_H30L_in[i][1] = 0.f;
		fftw_H30R_in[i][0] = H30Rshift[i];
		fftw_H30R_in[i][1] = 0.f;
		fftw_H330L_in[i][0] = H330Lshift[i];
		fftw_H330L_in[i][1] = 0.f;
		fftw_H330R_in[i][0] = H330Rshift[i];
		fftw_H330R_in[i][1] = 0.f;
	}

	//memset HRTF complex buffers manually
	for (int i = data->H30Linfo.frames; i < framesPerBuffer+data->H30Linfo.frames-1 ; i++)
	{
		fftw_H30L_in[i][0] = 0.f;
		fftw_H30L_in[i][1] = 0.f;
		fftw_H30R_in[i][0] = 0.f;
		fftw_H30R_in[i][1] = 0.f;
		fftw_H330L_in[i][0] = 0.f;
		fftw_H330L_in[i][1] = 0.f;
		fftw_H330R_in[i][0] = 0.f;
		fftw_H330R_in[i][1] = 0.f;
	}

	//import wisdom file
	if (fftw_import_wisdom_from_filename(filename) == 0)
	{
		printf("Error importing wisdom file.");
		exit(1);
	}

	//take FFTs of buffered audio
	plan_l = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, in_fft_l, out_fft_l, FFTW_FORWARD, flags);
	plan_r = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, in_fft_r, out_fft_r, FFTW_FORWARD, flags);
	fftw_execute(plan_l);
	fftw_execute(plan_r);

	//take FFTs of HRTFS
	plan_H30L = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, fftw_H30L_in, fftw_H30L_out, FFTW_FORWARD, flags);
	plan_H30R = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, fftw_H30R_in, fftw_H30R_out, FFTW_FORWARD, flags);
	plan_H330L = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, fftw_H330L_in, fftw_H330L_out, FFTW_FORWARD, flags);
	plan_H330R = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, fftw_H330R_in, fftw_H330R_out, FFTW_FORWARD, flags);
	fftw_execute(plan_H30L);
	fftw_execute(plan_H30R);
	fftw_execute(plan_H330L);
	fftw_execute(plan_H330R);

	//export wisdom file
	//if (fftw_export_wisdom_to_filename(filename) == 0)
	//{
	//	printf("Error exporting wisdom file.");
	//	exit(1);
	//}

	//destroy plans
	fftw_destroy_plan(plan_l);
	fftw_destroy_plan(plan_r);
	fftw_destroy_plan(plan_H30L);
	fftw_destroy_plan(plan_H30R);
	fftw_destroy_plan(plan_H330L);
	fftw_destroy_plan(plan_H330R);

	//test
	for (int i = 0; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	{
		w.a = out_fft_r[i][0];
		w.b = out_fft_r[i][1];
		z.a = out_fft_l[i][0];
		z.b = out_fft_l[i][1];
		x.a = fftw_H330R_out[i][0];
		x.b = fftw_H330R_out[i][1];
		y.a = fftw_H30L_out[i][0];
		y.b = fftw_H30L_out[i][1];
		o = multiply(w,x);
		p = multiply(z,y);
		out_fft_l[i][0] = o.a;
		out_fft_l[i][1] = o.b;
		out_fft_r[i][0] = p.a;
		out_fft_r[i][1] = p.b;
	}

	//TODO process HRTFs, convolve HRTFS with audio, OLA and output, save OLA segment
	// for (int i = 0; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	// {
	// 	//crosstalk_IR
	// 	x.a = fftw_H330R_out[i][0];
	// 	x.b = fftw_H330R_out[i][1];
	// 	y.a = fftw_H330L_out[i][0];
	// 	y.b = fftw_H330L_out[i][1];
	// 	w.a = fftw_H30L_out[i][0];
	// 	w.b = fftw_H30L_out[i][1];
	// 	z.a = fftw_H30R_out[i][0];
	// 	z.b = fftw_H30R_out[i][1];
	// 	o = divide(x,y);
	// 	p = divide(w,z);
	// 	crosstalk_IR_L[i][0] = -o.a;
	// 	crosstalk_IR_L[i][1] = -o.b;
	// 	crosstalk_IR_R[i][0] = -p.a;
	// 	crosstalk_IR_R[i][1] = -p.b;

	// 	//proc_IR
	// 	q = multiply(o,o);
	// 	r = multiply(p,p);
	// 	s.a = one.a - q.a;
	// 	s.b = one.b - q.b;
	// 	t.a = one.a - r.a;
	// 	t.b = one.b - r.b;
	// 	u = reciprocal(s);
	// 	v = reciprocal(t);

	// 	//cross_cancel
	// 	x.a = crosstalk_IR_R[i][0];
	// 	x.b = crosstalk_IR_R[i][1];
	// 	y.a = crosstalk_IR_L[i][0];
	// 	y.b = crosstalk_IR_L[i][1];
	// 	w.a = out_fft_r[i][0];
	// 	w.b = out_fft_r[i][1];
	// 	z.a = out_fft_l[i][0];
	// 	z.b = out_fft_l[i][1];
	// 	o = multiply(w,x);
	// 	p = multiply(z,y);

	// 	//add cross_cancel
	// 	x.a = out_fft_l[i][0]+o.a;
	// 	x.b = out_fft_l[i][1]+o.b;
	// 	y.a = out_fft_r[i][0]+p.a;
	// 	y.b = out_fft_r[i][1]+p.b;

	// 	//convolve proc_IR
	// 	z = multiply(x,u);
	// 	w = multiply(y,v);

	// 	//equalize and convolve final step
	// 	m.a = fftw_H330L_out[i][0];
	// 	m.b = fftw_H330L_out[i][1];
	// 	n.a = fftw_H30R_out[i][0];
	// 	n.b = fftw_H30R_out[i][1];
	// 	o = reciprocal(m);
	// 	p = reciprocal(n);
	// 	q = multiply(z,o);
	// 	r = multiply(w,p);

	// 	//to ifft
	// 	out_fft_l[i][0] = q.a;
	// 	out_fft_l[i][1] = q.b;
	// 	out_fft_r[i][0] = r.a;
	// 	out_fft_r[i][1] = r.b;
	// }

	// //memset out_fft
	// for (int i = (framesPerBuffer+data->H30Linfo.frames-1)/2+1; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	// {
	// 	out_fft_l[i][0] = 0.f;
	// 	out_fft_l[i][1] = 0.f;
	// 	out_fft_r[i][0] = 0.f;
	// 	out_fft_r[i][1] = 0.f;
	// }

	//ifft plans
	plan_l = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, out_fft_l, in_fft_l, FFTW_BACKWARD, flags);
	plan_r = fftw_plan_dft_1d(framesPerBuffer+data->H30Linfo.frames-1, out_fft_r, in_fft_r, FFTW_BACKWARD, flags);
	fftw_execute(plan_l);
	fftw_execute(plan_r);

	//find normalization factors
	// for (int i = 0; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	// {
	// 	if (fabs(in_fft_l[i][0]) > normL)
	// 	{
	// 		normL = fabs(in_fft_l[i][0]);
	// 	}
	// 	if (fabs(in_fft_r[i][0]) > normR)
	// 	{
	// 		normR = fabs(in_fft_r[i][0]);
	// 	}
	// }

	//normalize and copy to output buffer
	for (int i = 0; i < framesPerBuffer; i++)
	{
		outbuff[2*i] = ((float)in_fft_l[i][0]/(float)(framesPerBuffer+data->H30Linfo.frames-1) + leftOLA[i])/2.0f;
		outbuff[2*i+1] = ((float)in_fft_r[i][0]/(float)(framesPerBuffer+data->H30Linfo.frames-1) + rightOLA[i])/2.0f;
	}

	//save to OLA buffer
	for (int i = framesPerBuffer; i < framesPerBuffer+data->H30Linfo.frames-1; i++)
	{
		leftOLA[i - framesPerBuffer] = (float)in_fft_l[i][0]/(float)(framesPerBuffer+data->H30Linfo.frames-1);
		rightOLA[i - framesPerBuffer] = (float)in_fft_r[i][0]/(float)(framesPerBuffer+data->H30Linfo.frames-1);
	}

	//free vars
	fftw_destroy_plan(plan_l);
	fftw_destroy_plan(plan_r);
	fftw_free(buffer);
	fftw_free(H30Rshift);
	fftw_free(H30Lshift);
	fftw_free(H330Rshift);
	fftw_free(H330Lshift);
	fftw_free(in_fft_r);
	fftw_free(in_fft_l);
	fftw_free(out_fft_r);
	fftw_free(out_fft_l);
	fftw_free(right_ch);
	fftw_free(left_ch);
	fftw_free(crosstalk_IR_R);
	fftw_free(crosstalk_IR_L);
	fftw_free(cross_cancel_R);
	fftw_free(cross_cancel_L);
	fftw_free(proc_IR_L);
	fftw_free(proc_IR_R);
	fftw_free(eq_IR_R);
	fftw_free(eq_IR_L);
	//fftw_free(mult_fft_l);
	//fftw_free(mult_fft_r);

	return paContinue;
}

Complex divide(Complex x, Complex y)
{
	Complex z;
	z.a = (x.a*y.a + x.b*y.b)/(y.a*y.a + y.b*y.b);
	z.b = (x.b*y.a - x.a*y.b)/(y.a*y.a + y.b*y.b);
	return z;
}

Complex multiply(Complex x, Complex y)
{
	Complex z;
	z.a = x.a*y.a - x.b*y.b;
	z.b = x.a*y.b + x.b*y.a;
	return z;
}

Complex reciprocal(Complex x)
{
	Complex z;
	z.a = x.a/(x.a*x.a + x.b*x.b);
	z.b = -x.b/(x.a*x.a + x.b*x.b);
	return z;
}
