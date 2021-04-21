#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <complex.h>
#include <string.h>

#define Form ("WAVE")
#define Chunk_ID ("RIFF")
#define SubC_1_ID ("fmt ")
#define SubC_2_ID ("data")
#define Sub_C_1_size (16)
#define Audio_format (1)
#define Num_cha (1)
#define Bits_Per_Sample (16)


// header file'S format
typedef struct header{
	char ChunkID[4];
	int ChunkSize;
	char Format[4];
	char Subchunk1ID[4];
	int Subchunk1Size;
	short AudioFormat;
	short NumChannels;
	int SampleRate;
	int ByteRate;
	short BlockAlign;
	short BitsPerSample;
	char SubChunk2ID[4];
	int SubChunk2Size;


}wav_header;


//this struct is used to do complex calculation
typedef struct {
	double real;
	double imag;
}comp;


// set up the header file of wave file
void header_filling(int Sample_Rate, int Total_Sample, FILE* wav_file){
	wav_header cos_header;

	for (int i = 0; i <4 ; i++)
		cos_header.ChunkID[i] = Chunk_ID[i];
	for(int i = 0; i <4; i++)
		cos_header.Format[i] = Form[i];
	for(int i = 0; i <4; i++)
		cos_header.Subchunk1ID[i] = SubC_1_ID[i];
	for(int i=0;i<4;i++)
		cos_header.SubChunk2ID[i] = SubC_2_ID[i];

	cos_header.Subchunk1Size = Sub_C_1_size;
	cos_header.AudioFormat = Audio_format;
	cos_header.NumChannels = Num_cha;
	cos_header.SampleRate = Sample_Rate;
	cos_header.ByteRate = (cos_header.SampleRate * cos_header.NumChannels) *2;
	cos_header.BlockAlign = cos_header.NumChannels * Bits_Per_Sample/8;
	cos_header.BitsPerSample = Bits_Per_Sample;
	cos_header.SubChunk2Size = Total_Sample * Num_cha * 2;
	cos_header.ChunkSize = 20 + cos_header.SubChunk2Size + Sub_C_1_size;
	fwrite(&cos_header, sizeof(wav_header), 1, wav_file);

}


// takes in freq, total sample to make cosine waves with different frequency and sampling rate
void Cos_wave_generating(double Total_Sample, short *temp_storing, short Freq_Hz)
{
	double Sample =0;
	for (int Sample_time = 0; Sample_time < Total_Sample; Sample_time++){
		Sample = Sample_time / Total_Sample ;
		temp_storing[Sample_time] = 127 + 10000 * cos(2 * acos(-1) * Freq_Hz * Sample);
	}
}

//dynamically allocation
short * allocating_memory(unsigned int Total_Sample){
	return (short*)malloc(sizeof(short) * Total_Sample);
}


//dynamically allocation of float numbers
float *allocating_memory_float(unsigned int Total_Sample){
	return (float*)malloc(sizeof(float)* Total_Sample);
}


//DFT function with rectangular function window
void calculating_DFT_rec(int M ,int K,float *windowing_result, unsigned int Total_Sample, FILE* cos_file, float win){
	int N = K;
	int last_n = 0;
	int last_boundary = K;
	comp sum;
	comp result[M][K];
	comp buffer[N];
	int jump = Total_Sample / M; //every other value of M will be a new range of n through adding this jump value

	int time = 0;
	int total_win = win * Total_Sample;   // how many points there are for windowing function
	float ** amplitude;  // used to store final results
	amplitude = (float**)malloc(sizeof(float *) * M);
	for (int n = 0 ; n < M; n++)
		amplitude[n] = (float *)malloc(sizeof(float)*K);

	sum.real=0;
	sum.imag=0;  // initialize

	for ( int m = 0; m < M; m++){
		for(int k = 0 ; k < K; k++){
			for(int n = last_n; n < last_boundary; n++, time++){
				if (time < total_win){
					buffer[time].real =cos((2 * acos(-1) * k * time) / N) * windowing_result[n];
					buffer[time].imag = 0-sin((2* acos(-1)* k * time)/ N) * windowing_result[n];

					sum.real += buffer[time].real; //store the sum of calculation
					sum.imag += buffer[time].imag;
				}
				else{
					sum.real +=0;
					sum.imag +=0;
				}
		   	}
		result[m][k].real = sum.real;
		result[m][k].imag = sum.imag;
		sum.real = 0;
		sum.imag = 0;
		time = 0;
		}

	last_n += jump;
	last_boundary += jump;

	if(last_boundary >= Total_Sample ){
		last_boundary = Total_Sample;    //if n is over the total sample point make it back to total_sample // ex: 8024- 8000
		}
	}

	for(int i = 0 ; i < M; i++){
		for(int j = 0 ; j< K ;j++){
			if(result[i][j].real==0 && result[i][j].imag==0)
				amplitude[i][j]=0;
			else
				amplitude[i][j] = 20 * log10(sqrt(result[i][j].real * result[i][j].real + result[i][j].imag * result[i][j].imag));   //doing the log calculation
		}
	}

	for(int i = 0; i < M ; i++){
		for(int j = 0; j < K ; j++){
		fprintf(cos_file, "%f ", amplitude[i][j]);
		}
		fprintf(cos_file, "\n");  // store in all the data
	}

}

// except for the windowing part it's all the same
void calculating_DFT_ham(int M ,int K,float *windowing_result, unsigned int Total_Sample, FILE* cos_file, float win){
	printf("M= %d K=%d ", M, K);
	int N = K;
	int last_n = 0;
	int last_boundary = K;
	comp sum;
	comp result[M][K];
	comp buffer[N];
	int jump = Total_Sample / M;
	int time = 0;
	int total_win = win * Total_Sample;
	float ** amplitude;
	amplitude = (float**)malloc(sizeof(float *) * M);
	for (int n = 0 ; n < M; n++)
		amplitude[n] = (float *)malloc(sizeof(float)*K);

	sum.real=0;
	sum.imag=0;

	for ( int m = 0; m < M; m++){
		for(int k = 0 ; k < K; k++){
			for(int n = last_n; n < last_boundary; n++, time++){
				if (time < total_win){
					//this is the hamming windowing function , if it is over p , just make the value to be 0
					buffer[time].real =cos((2 * acos(-1) * k * time) / N) * windowing_result[n] * (0.54 - 0.46 *cos(2 * acos(-1)*time)/(total_win-1));
					buffer[time].imag = 0-sin((2* acos(-1)* k * time)/ N) * windowing_result[n] * (0.54 - 0.46 *cos(2 * acos(-1)*time)/(total_win-1));

					sum.real += buffer[time].real;
					sum.imag += buffer[time].imag;
				}
				else{
					sum.real +=0;
					sum.imag +=0;
				}
		   	}
		result[m][k].real = sum.real;
		result[m][k].imag = sum.imag;
		sum.real = 0;
		sum.imag = 0;
		time = 0;
		}

	last_n += jump;
	last_boundary += jump;

	if(last_boundary >= Total_Sample ){
		last_boundary = Total_Sample;
		}
	}

	for(int i = 0 ; i < M; i++){
		for(int j = 0 ; j< K ;j++){
			if(result[i][j].real==0 && result[i][j].imag==0)
				amplitude[i][j]=0;
			else
				amplitude[i][j] = 20 * log10(sqrt(result[i][j].real * result[i][j].real + result[i][j].imag * result[i][j].imag));
		}
	}

	for(int i = 0; i < M ; i++){
		for(int j = 0; j < K ; j++){
		fprintf(cos_file, "%f ", amplitude[i][j]);
		}
		fprintf(cos_file, "\n");
	}

}


// this was used to do vowel calculation since i coundn't make it successfully so i left it here for your reference
// it takes in one more integer (sample rate ) to do more calculation for longer period of wave file
void calculating_DFT_rec_vowel(int M ,int K,float *windowing_result, unsigned int Total_Sample, FILE* cos_file, float win, int Sample_Rate){
	int N = K;
	int last_n = 0;
	int last_boundary = K;
	comp sum;
	comp result[M][K];
	comp buffer[N];
	int jump = Total_Sample / M;

	int time = 0;
	int total_win = win * Sample_Rate;
	float ** amplitude;
	amplitude = (float**)malloc(sizeof(float *) * M);
	for (int n = 0 ; n < M; n++)
		amplitude[n] = (float *)malloc(sizeof(float)*K);

	sum.real=0;
	sum.imag=0;

	for ( int m = 0; m < M; m++){
		for(int k = 0 ; k < K; k++){
			for(int n = last_n; n < last_boundary; n++, time++){
				if (time < total_win){
					buffer[time].real =cos((2 * acos(-1) * k * time) / N) * windowing_result[n] ;
					buffer[time].imag = 0-sin((2* acos(-1)* k * time)/ N) * windowing_result[n];

					sum.real += buffer[time].real;
					sum.imag += buffer[time].imag;
				}
				else{
					sum.real +=0;
					sum.imag +=0;
				}
		   	}
		result[m][k].real = sum.real;
		result[m][k].imag = sum.imag;
		sum.real = 0;
		sum.imag = 0;
		time = 0;
		}

	last_n += jump;
	last_boundary += jump;

	if(last_boundary >= Total_Sample ){
		last_boundary = Total_Sample;
		}
	}

	for(int i = 0 ; i < M; i++){
		for(int j = 0 ; j< K ;j++){
			if(result[i][j].real==0 && result[i][j].imag==0)
				amplitude[i][j]=0;
			else
				amplitude[i][j] = 20 * log10(sqrt(result[i][j].real * result[i][j].real + result[i][j].imag * result[i][j].imag));
		}
	}

	for(int i = 0; i < M ; i++){
		for(int j = 0; j < K ; j++){
		fprintf(cos_file, "%f ", amplitude[i][j]);
		}
		fprintf(cos_file, "\n");
	}

}


// this is the same except for the windowing part
void calculating_DFT_ham_vowel(int M ,int K,float *windowing_result, unsigned int Total_Sample, FILE* cos_file, float win, int Sample_Rate){
	printf("M= %d K=%d ", M, K);
	int N = K;
	int last_n = 0;
	int last_boundary = K;
	comp sum;
	comp result[M][K];
	comp buffer[N];
	int jump = Total_Sample / M;
	int time = 0;
	int total_win = win * Sample_Rate;
	float ** amplitude;
	amplitude = (float**)malloc(sizeof(float *) * M);
	for (int n = 0 ; n < M; n++)
		amplitude[n] = (float *)malloc(sizeof(float)*K);

	sum.real=0;
	sum.imag=0;

	for ( int m = 0; m < M; m++){
		for(int k = 0 ; k < K; k++){
			for(int n = last_n; n < last_boundary; n++, time++){
				if (time < total_win){
					buffer[time].real =cos((2 * acos(-1) * k * time) / N) * windowing_result[n] * (0.54 - 0.46 *cos(2 * acos(-1)*time)/(total_win-1));
					buffer[time].imag = 0-sin((2* acos(-1)* k * time)/ N) * windowing_result[n] * (0.54 - 0.46 *cos(2 * acos(-1)*time)/(total_win-1));

					sum.real += buffer[time].real;
					sum.imag += buffer[time].imag;
				}
				else{
					sum.real +=0;
					sum.imag +=0;
				}
		   	}
		result[m][k].real = sum.real;
		result[m][k].imag = sum.imag;
		sum.real = 0;
		sum.imag = 0;
		time = 0;
		}

	last_n += jump;
	last_boundary += jump;

	if(last_boundary >= Total_Sample ){
		last_boundary = Total_Sample;
		}
	}

	for(int i = 0 ; i < M; i++){
		for(int j = 0 ; j< K ;j++){
			if(result[i][j].real==0 && result[i][j].imag==0)
				amplitude[i][j]=0;
			else
				amplitude[i][j] = 20 * log10(sqrt(result[i][j].real * result[i][j].real + result[i][j].imag * result[i][j].imag));
		}
	}

	for(int i = 0; i < M ; i++){
		for(int j = 0; j < K ; j++){
		fprintf(cos_file, "%f ", amplitude[i][j]);
		}
		fprintf(cos_file, "\n");
	}

}

	int main(){
// some basic values for all the rest of the calculation
		short T = 1;
		short Freq_Hz_0 = 50;
		short Freq_Hz_1 = 55;
		short Freq_Hz_2 = 200;
		short Freq_Hz_3 = 220;
		float DFT_window_size_0 = 0.008;
		float DFT_window_size_1 = 0.032;
		float frame_interval_0 = 0.005;
		float frame_interval_1 = 0.01;
		int Sample_Rate_0 = 8000;
		int Sample_Rate_1 = 16000;

		int setting_0_M = 1/frame_interval_0 ;
		int setting_0_K_8k = DFT_window_size_0 * Sample_Rate_0;
		int setting_0_K_16k = DFT_window_size_0 * Sample_Rate_1;

		int setting_1_M = 1/frame_interval_1;
		int setting_1_K_8k = DFT_window_size_1 * Sample_Rate_0;
		int setting_1_K_16k = DFT_window_size_1 * Sample_Rate_1;

		int32_t Total_Sample_0 = T * Sample_Rate_0;
		int32_t Total_Sample_1 = T * Sample_Rate_1;


// fopening all the txt file in here
		FILE * cos_8k_50_set1;
		cos_8k_50_set1 = fopen("cos_050Hz-8k.{set1}.txt", "w+");
		FILE * cos_8k_50_set2;
		cos_8k_50_set2 = fopen("cos_050Hz-8k.{set2}.txt", "w+");
		FILE * cos_8k_50_set3;
		cos_8k_50_set3 = fopen("cos_050Hz-8k.{set3}.txt", "w+");
		FILE * cos_8k_50_set4;
		cos_8k_50_set4 = fopen("cos_050Hz-8k.{set4}.txt", "w+");



		FILE * cos_8k_55_set1;
		cos_8k_55_set1 = fopen("cos_055Hz-8k.{set1}.txt", "w+");
		FILE * cos_8k_55_set2;
		cos_8k_55_set2 = fopen("cos_055Hz-8k.{set2}.txt", "w+");
		FILE * cos_8k_55_set3;
		cos_8k_55_set3 = fopen("cos_055Hz-8k.{set3}.txt", "w+");
		FILE * cos_8k_55_set4;
		cos_8k_55_set4 = fopen("cos_055Hz-8k.{set4}.txt", "w+");

		FILE * cos_16k_50_set1;
		cos_16k_50_set1 = fopen("cos_050Hz-16k.{set1}.txt", "w+");
		FILE * cos_16k_50_set2;
		cos_16k_50_set2 = fopen("cos_050Hz-16k.{set2}.txt", "w+");
		FILE * cos_16k_50_set3;
		cos_16k_50_set3 = fopen("cos_050Hz-16k.{set3}.txt", "w+");
		FILE * cos_16k_50_set4;
		cos_16k_50_set4 = fopen("cos_050Hz-16k.{set4}.txt", "w+");

		FILE * cos_16k_55_set1;
		cos_16k_55_set1 = fopen("cos_055Hz-16k.{set1}.txt", "w+");
		FILE * cos_16k_55_set2;
		cos_16k_55_set2 = fopen("cos_055Hz-16k.{set2}.txt", "w+");
		FILE * cos_16k_55_set3;
		cos_16k_55_set3 = fopen("cos_055Hz-16k.{set3}.txt", "w+");
		FILE * cos_16k_55_set4;
		cos_16k_55_set4 = fopen("cos_055Hz-16k.{set4}.txt", "w+");

		FILE * cos_8k_200_set1;
		cos_8k_200_set1 = fopen("cos_200Hz-8k.{set1}.txt", "w+");
		FILE * cos_8k_200_set2;
		cos_8k_200_set2 = fopen("cos_200Hz-8k.{set2}.txt", "w+");
		FILE * cos_8k_200_set3;
		cos_8k_200_set3 = fopen("cos_200Hz-8k.{set3}.txt", "w+");
		FILE * cos_8k_200_set4;
		cos_8k_200_set4 = fopen("cos_200Hz-8k.{set4}.txt", "w+");

		FILE * cos_8k_220_set1;
		cos_8k_220_set1 = fopen("cos_220Hz-8k.{set1}.txt", "w+");
		FILE * cos_8k_220_set2;
		cos_8k_220_set2 = fopen("cos_220Hz-8k.{set2}.txt", "w+");
		FILE * cos_8k_220_set3;
		cos_8k_220_set3 = fopen("cos_220Hz-8k.{set3}.txt", "w+");
		FILE * cos_8k_220_set4;
		cos_8k_220_set4 = fopen("cos_220Hz-8k.{set4}.txt", "w+");

		FILE * cos_16k_200_set1;
		cos_16k_200_set1 = fopen("cos_200Hz-16k.{set1}.txt", "w+");
		FILE * cos_16k_200_set2;
		cos_16k_200_set2 = fopen("cos_200Hz-16k.{set2}.txt", "w+");
		FILE * cos_16k_200_set3;
		cos_16k_200_set3 = fopen("cos_200Hz-16k.{set3}.txt", "w+");
		FILE * cos_16k_200_set4;
		cos_16k_200_set4 = fopen("cos_200Hz-16k.{set4}.txt", "w+");

		FILE * cos_16k_220_set1;
		cos_16k_220_set1 = fopen("cos_220Hz-16k.{set1}.txt", "w+");
		FILE * cos_16k_220_set2;
		cos_16k_220_set2 = fopen("cos_220Hz-16k.{set2}.txt", "w+");
		FILE * cos_16k_220_set3;
		cos_16k_220_set3 = fopen("cos_220Hz-16k.{set3}.txt", "w+");
		FILE * cos_16k_220_set4;
		cos_16k_220_set4 = fopen("cos_220Hz-16k.{set4}.txt", "w+");


// some buffers used to temporarily store value from fread function
		float *windowing_result_8k = NULL;
		windowing_result_8k = allocating_memory_float(Total_Sample_0);
		float *windowing_result_16k_50Hz = NULL;
		windowing_result_16k_50Hz = allocating_memory_float(Total_Sample_1);

		short *buffer_8k = NULL;
		buffer_8k = allocating_memory(Total_Sample_0);
		short *buffer_16k = NULL;
		buffer_16k = allocating_memory(Total_Sample_1);
		short *temp_storing = NULL;
		temp_storing = allocating_memory(Total_Sample_0);


// this is where the major calculation begins
		FILE* cos050Hz_8k;
		cos050Hz_8k = fopen("cos_050Hz-8k.wav","wb+");
		Cos_wave_generating(Total_Sample_0, temp_storing, Freq_Hz_0); //generating cos wave with the corresponding Hertz and sampling rate
		header_filling(Sample_Rate_0, Total_Sample_0, cos050Hz_8k); // write the header data with it's corresponding variable
		fwrite(temp_storing, sizeof(short), Total_Sample_0, cos050Hz_8k); //write the data into the wave file with the size of short
		fseek(cos050Hz_8k , 44, SEEK_SET); //seek to the part where real data is stored
		fread(buffer_8k, sizeof(short), Total_Sample_0, cos050Hz_8k); // store in value to the buffer
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_50_set1,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_ham(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_50_set2,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_50_set3,0.02);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_50_set4,0.02);
		fclose(cos050Hz_8k);

// the comments for all the rest will be the same
		FILE* cos055Hz_8k;
		cos055Hz_8k = fopen("cos_055Hz-8k.wav","wb+");
		Cos_wave_generating(Total_Sample_0, temp_storing, Freq_Hz_1);
		header_filling(Sample_Rate_0, Total_Sample_0, cos055Hz_8k);
		fwrite(temp_storing, sizeof(short), Total_Sample_0, cos055Hz_8k);
		fseek(cos055Hz_8k , 44, SEEK_SET);
		fread(buffer_8k, sizeof(short), Total_Sample_0, cos055Hz_8k);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_55_set1,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_55_set2,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_55_set3,0.02);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_55_set4,0.02);
		fclose(cos055Hz_8k);

		FILE* cos200Hz_8k;
		cos200Hz_8k = fopen("cos_200Hz-8k.wav","wb+");
		Cos_wave_generating(Total_Sample_0, temp_storing, Freq_Hz_2);
		header_filling(Sample_Rate_0, Total_Sample_0, cos200Hz_8k);
		fwrite(temp_storing, sizeof(short), Total_Sample_0, cos200Hz_8k);
		fseek(cos200Hz_8k , 44, SEEK_SET);
		fread(buffer_8k, sizeof(short), Total_Sample_0, cos200Hz_8k);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i];
		calculating_DFT_rec(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_200_set1,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i];
		calculating_DFT_ham(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_200_set2,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_200_set3,0.02);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_200_set4,0.02);
		fclose(cos200Hz_8k);


		FILE* cos220Hz_8k;
		cos220Hz_8k = fopen("cos_220Hz-8k.wav","wb+");
		Cos_wave_generating(Total_Sample_0, temp_storing, Freq_Hz_3);
		header_filling(Sample_Rate_0, Total_Sample_0, cos220Hz_8k);
		fwrite(temp_storing, sizeof(short), Total_Sample_0, cos220Hz_8k);
		fseek(cos220Hz_8k , 44, SEEK_SET);
		fread(buffer_8k, sizeof(short), Total_Sample_0, cos220Hz_8k);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_220_set1,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i];
		calculating_DFT_ham(setting_0_M, setting_0_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_220_set2,0.005);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i];
		calculating_DFT_rec(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_220_set3,0.02);
		for (int i = 0 ; i <Total_Sample_0; i++)
			windowing_result_8k[i] = buffer_8k[i];
		calculating_DFT_ham(setting_1_M, setting_1_K_8k, windowing_result_8k, Total_Sample_0, cos_8k_220_set4,0.02);
		fclose(cos220Hz_8k);

		short *temp_storing_1 = NULL;
		temp_storing_1 = allocating_memory(Total_Sample_1);


		FILE* cos050Hz_16k;
		cos050Hz_16k = fopen("cos_050Hz-16k.wav","wb+");
		Cos_wave_generating(Total_Sample_1, temp_storing_1, Freq_Hz_0);
		header_filling(Sample_Rate_1, Total_Sample_1, cos050Hz_16k);
		fwrite(temp_storing_1, sizeof(short), Total_Sample_1, cos050Hz_16k);
		fseek(cos050Hz_16k , 44, SEEK_SET);
		fread(buffer_16k, sizeof(short), Total_Sample_1, cos050Hz_16k);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_50_set1,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_50_set2,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_50_set3,0.02);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_50_set4,0.02);
		fclose(cos050Hz_16k);


		FILE* cos055Hz_16k;
		cos055Hz_16k = fopen("cos_055Hz-16k.wav", "wb+");
		Cos_wave_generating(Total_Sample_1, temp_storing_1, Freq_Hz_1);
		header_filling(Sample_Rate_1, Total_Sample_1, cos055Hz_16k);
		fwrite(temp_storing_1, sizeof(short), Total_Sample_1, cos055Hz_16k);
		fseek(cos055Hz_16k , 44, SEEK_SET);
		fread(buffer_16k, sizeof(short), Total_Sample_1, cos055Hz_16k);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_55_set1,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_55_set2,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_55_set3,0.02);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_55_set4,0.02);
		fclose(cos055Hz_16k);

		FILE* cos200Hz_16k;
		cos200Hz_16k = fopen("cos_200Hz-16k.wav", "wb+");
		Cos_wave_generating(Total_Sample_1, temp_storing_1, Freq_Hz_2);
		header_filling(Sample_Rate_1, Total_Sample_1, cos200Hz_16k);
		fwrite(temp_storing_1, sizeof(short), Total_Sample_1, cos200Hz_16k);
		fseek(cos200Hz_16k , 44, SEEK_SET);
		fread(buffer_16k, sizeof(short), Total_Sample_1, cos200Hz_16k);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_200_set1,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_200_set2,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_rec(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_200_set3,0.02);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_200_set4,0.02);
		fclose(cos200Hz_16k);

		FILE* cos220Hz_16k;
		cos220Hz_16k = fopen("cos_220Hz-16k.wav", "wb+");
		Cos_wave_generating(Total_Sample_1, temp_storing_1, Freq_Hz_3);
		header_filling(Sample_Rate_1, Total_Sample_1, cos220Hz_16k);
		fwrite(temp_storing_1, sizeof(short), Total_Sample_1, cos220Hz_16k);
		fseek(cos220Hz_16k , 44, SEEK_SET);
		fread(buffer_16k, sizeof(short), Total_Sample_1, cos220Hz_16k);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i];
		calculating_DFT_rec(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_220_set1,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i] ;
		calculating_DFT_ham(setting_0_M, setting_0_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_220_set2,0.005);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i];
		calculating_DFT_rec(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_220_set3,0.02);
		for (int i = 0 ; i < Total_Sample_1; i++)
			windowing_result_16k_50Hz[i] = buffer_16k[i];
		calculating_DFT_ham(setting_1_M, setting_1_K_16k, windowing_result_16k_50Hz, Total_Sample_1, cos_16k_220_set4,0.02);
		fclose(cos220Hz_16k);



// this files were used to do matlab plotting but then i realized i didn't really need this to run matlab
// and i decided just to leave it there

		FILE * x;  //200
		x = fopen("x{format1}.txt", "w+");
		for (int i = 0 ; i < setting_0_M; i ++){
			fprintf(x,"%d", i+1);
			fprintf(x,"\n");
		}
		fclose(x);

		FILE * y; // 64
		y = fopen("y{format1}.txt", "w+");
		for (int i = 0 ; i < setting_0_K_8k; i ++){
			fprintf(y,"%d", i+1);
			fprintf(y,"\n");
		}
		fclose(y);

		FILE * x_1; //100
		x_1 = fopen("x{format2}.txt", "w+");
		for (int i = 0 ; i < setting_1_M; i ++){
			fprintf(x_1,"%d", i+1);
			fprintf(x_1,"\n");
		}
		fclose(x_1);

		FILE * y_1; //256
		y_1 = fopen("y{format2}.txt", "w+");
		for (int i = 0 ; i < setting_1_K_8k; i ++){
			fprintf(y_1,"%d", i+1);
			fprintf(y_1,"\n");
		}
		fclose(y_1);


		FILE * y_2;
		y_2 = fopen("y{format3}.txt", "w+");
		for (int i = 0 ; i < setting_0_K_16k; i ++){
			fprintf(y_2,"%d", i+1);
			fprintf(y_2,"\n");
		}
		fclose(y_2);


		FILE * y_3;
		y_3 = fopen("y{format4}.txt", "w+");
		for (int i = 0 ; i < setting_1_K_16k; i ++){
			fprintf(y_3,"%d", i+1);
			fprintf(y_3,"\n");
		}
		fclose(y_3);
		printf("program ends");
	}


