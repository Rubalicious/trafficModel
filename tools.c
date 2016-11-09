#include <math.h>
#include <stdlib.h>
#include "tools.h"

void get_positions(int num, Car *cars, double pos[]){
	Car *ptr;
	int i;
	for (ptr = cars, i = 0; (ptr-cars) < num; ptr++, i++){
		pos[i] = ptr->x;
	}
}

double gauss_rand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

double uniform_rand() { return ((double)rand()/(double)RAND_MAX); }

double min(int num, double arr[]){
	int i;
	double m = arr[0];
	for(i = 0; i < num; i++){
		if (m <= arr[i])
		{
			m = arr[i];
		}
	}		
	return m;
}

double max(int num, double arr[]){
	int i;
	double m = arr[0];
	for(i = 0; i < num; i++){
		if (m >= arr[i])
		{
			m = arr[i];
		}
	}		
	return m;
}

int is_monotonic(int num, double arr[]){
	int i;
	int is_inc = 1, is_dec = 1; // assume are increasing
	for (i = 0; i < num-1; ++i)
	{
		if( arr[i] < arr[i+1] ){	is_dec = 0;		}
		if( arr[i] > arr[i+1] ){	is_inc = 0;		}
	}
	if( is_inc || is_dec ){		return 1;	}
	else{	return 0;	}
}

int index_of(double val, int num, double arr[]){
	int i;
	for (i = 0; i < num; ++i)
	{
		if (arr[i] == val)
		{
			return i;	// return the index
		}
	}
	return -1;	// not in the array
}