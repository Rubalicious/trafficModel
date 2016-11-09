// tools.h 
#ifndef TOOLS_H_   /* Include guard */
#define TOOLS_H_

typedef struct Car
{	
	// DRIVER PARAMETERS
	double car_size, tao, alpha, a_max, a_min, v_steep;
	int v_opt;
	double x, v; // position, velocity
	double x1, v1; // next position, next velocity
	struct Car *next;
} Car;

// #####################################################
// #########               TOOLS              ##########
// #####################################################

// copies to pos all current positions of num cars
void get_positions(int num, struct Car *cars, double pos[]);

// Returns the smallest value for a given array
double min(int num, double arr[]);

// Returns the smallest value for a given array
double max(int num, double arr[]);

// will check if a given array of values is monotonically inc/dec
int is_monotonic(int num, double arr[]);

// returns index of a given value in an array
int index_of(double val, int num, double arr[]);

// gaussian random variable generator with mu = 0, sigma = 1
double gauss_rand();

// uniform random variable generator with range (0,1]
double uniform_rand();

#endif // TOOLS_H_