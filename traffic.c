/*
* File: traffic.c
* Author: Ruby Abrams
* Purpose: This will run traffic simulations.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "tools.h"
#include <gsl/gsl_types.h>
#include "gsl/rng/gsl_rng.h"
#include "gsl/randist/gsl_randist.h"
#include "gsl/specfunc/gsl_sf_trig.h"
#include "gsl/complex/gsl_complex.h"
#include "gsl/complex/gsl_complex_math.h"

// #################################################################	
// 							GLOBAL VARIABLES
// #################################################################
double Alpha, A_max, A_min, Car_size, Dt, Tao, V_opt, V_steep, v_std; 	
int distance;

// initializing random number generator
const gsl_rng_type * T;
gsl_rng * r;

// #####################################################
// #########        Simulation Methods        ##########
// #####################################################

// Gives all cars their unique properties
void initialize_cars(int n, Car *cars, gsl_rng * r);

// Gives all cars their unique properties
void initialize_cars_bug1(Car *cars, gsl_rng * r);

void initialize_cars_two_car_experiment(Car *cars, gsl_rng *r);

// This method establishes random initial conditions for num cars
void setup_Initial_Conditions(int num, Car *cars);

// Initial conditions to fix a bug in code for 2 cars
void setup_Initial_Conditions_bug1(Car *cars);

// two car experiment
void setup_Initial_Conditions_two_car_experiment(Car *cars);

// This will implement wrap-around boundary conditions
void incorporate_boundary_conditions(int num, Car *cars);

// collision time of car1 from car2
double collision_time(Car *car1, Car *car2);

// update position of a given car
// double update_position(double x, double v);
double update_position(Car *car);

// velocity function
// double update_velocity(double v, double a);
double update_velocity(Car *car, double acc);

// acceleration function
double update_acceleration(double t_collision, Car *car);

// Collision detection: returns 1 if any pairs of cars have collided and 0 otherwise
int collision_detection(int num, double pos_1[], double time_);

int mod (int a, int b)
{
   if(b < 0) //you can check for b == 0 separately and do what you want
     return mod(a, -b);   
   int ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}


int main()
{	

	// #################################################################	
	// 							GLOBAL VARIABLES
	// #################################################################
	Dt = .001; 				// time step ~ approximately 1/1000 of a second
	v_std = 0.1;			// standard deviation of velocities
	distance = 10000; 		// length of road in [m]

	// DRIVER PARAMETERS
	Car_size = 4;			// determines how close cars get to each other [m]
	Tao = 3000*Dt; 			// breaking reaction time ~ approximately 3 seconds
	Alpha = 10;				// driver's breaking coefficient [m/s^2]
	A_max = 2.5;			// max acceleration of a vehicle [m/s^2]
	A_min = -A_max;			// min acceleration of a vehicle [m/s^2]
	V_opt = 27; 			// speed limit [m/s]
	V_steep = Alpha/10;		// steepness of region where acceleration changes from max to min
	
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T); // instantiating random number generator

	// #####################################################
	// #########            Simulation            ##########
	// #####################################################

	int num = 2;	// number of cars simulated

	Car *cars = NULL;
	if(!cars) cars = (Car *)malloc(num*sizeof(Car));
	initialize_cars(num, cars, r);	// gives each car, its properties;
	// initialize_cars_bug1(cars, r);
	initialize_cars_two_car_experiment(cars, r);

	// Spreads cars normally in equally sized intervals 
	// for each car along 1000 m road.
	// setup_Initial_Conditions(num, cars);
	setup_Initial_Conditions_two_car_experiment(cars);
	// setup_Initial_Conditions_bug1(cars);

	double time_ = 0; // time variable to track duration of simulation in seconds
	// double avg_vel = 0; // this will track the average velocity of all cars
	int count = 0, increment = 10; // keeps track of the iteration of the while loop
	double t_collision, acc; // collision time	
	// double pos[num]; // records positions

	// Car *robot = malloc(sizeof(Car));
	// robot->x = 0;
	// robot->v = 0;
	while(time_ < 10*60){ // 3600 sec = 60 min = 1 hr

		if(mod(count, increment) == 0){
			printf("%.3f", time_);	
		}
		
		Car *car1, *car2;
		// // for each car, move the car
		for (car1 = cars, car2 = cars+1; (car1-cars) < num-1 ; car1++, car2++){ // 
			t_collision = collision_time(car1, car2);
			acc = update_acceleration(t_collision, car1);
			car1->v1 = update_velocity(car1, acc);
			car1->x1 = update_position(car1);

			// displaying the position and velocity of the ith car
			if(mod(count, increment) == 0){
				printf(" %.3f %.3f", car1->x, car1->v);	
			}
			
		}

		// side cases
		// car1 now points to first Car
		t_collision = collision_time(car1, cars); // compare last car to first car
		acc = update_acceleration(t_collision, car1);
		car1->v1 = update_velocity(car1, acc);
		car1->x1 = update_position(car1);
		if(mod(count, increment) == 0){
			printf(" %.3f %.3f\n", car1->x, car1->v);	
		}

		// incorporate boundary conditions: making a loop

		incorporate_boundary_conditions(num, cars);

		// COLLISION DETECTION
		// if(cars[0].x >= cars[1].x - .5*cars[1].car_size){
		// 	break;
		// }
		
		// get_positions(num, cars, pos);

		// // Collision detection
		// int collision_occurred = collision_detection(num, pos, time_);
		// if (collision_occurred)
		// {
		// 	break;
		// }

		// experiment: after traveling two hours,
		// place a red traffic light in the middle of the road.

		// // find the index of the car with smallest position
		// double least_pos = min(num, pos);
		// // find the index of the car with largest position
		// double max_pos = max(num, pos);
		// Car *minptr, *maxptr;
		// for(minptr = cars; minptr->x != least_pos; minptr++);
		// for(maxptr = cars; maxptr->x != max_pos; maxptr++);
		// if(time_ == 4){	
		// 	// adding a robot
		// 	maxptr->next = robot;
		// 	robot->next = minptr;

		// }

		time_ += Dt;	// incriment time step
		count++;		// incriment iteration of loop

		Car *ptr;
		for (ptr = cars; (ptr-cars) < num; ptr++){
			ptr->x = ptr->x1;		// updating current position array
			ptr->v = ptr->v1;		// updating current velocity array
		}
	}
	fprintf(stderr, "        Summary      \n");
	fprintf(stderr, "=====================\n");
	fprintf(stderr, "time taken:      %.3f\n", time_);
	fprintf(stderr, "iteration count:   %d\n", count);
	// printf("Total time: %.2f\n", time_ );
	gsl_rng_free(r);
	free(cars);
	return 0;
}

void initialize_cars(int n, Car *cars, gsl_rng * r){

	time_t mytime = time(NULL); /* GET TIME FOR RANDOM SEED */
 	gsl_rng_set(r,(unsigned long int)mytime);

	Car *car, *nextc;
	for (car = cars, nextc = car+1; (car-cars) < n; car++){
		car->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
		car->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
		car->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
		car->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
		car->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
		car->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
		car->a_min = -car->a_max;
		car->next = nextc; // cars point to the next car
	}
	// edge case: last car points to first car
	car->next = cars; // this is a potential memory leak
}


void setup_Initial_Conditions(int num, Car *cars){
	printf("time");
	// size of interval is distance/(num cars)
	// mean at each i is size of interval over 2 + i*size of interval 
	double interval_size = (double)distance/(double)num;

	Car *ptr;
	int i; // indexing the cars
	for (ptr = cars, i = 0; (ptr-cars) < num; ptr++, i++){
		printf(" x%d v%d", i+1, i+1);
		ptr->x = .06*(interval_size/2 +i*interval_size); // (interval_size/6)*gauss_rand() +
		ptr->v = 2*i;
		ptr->x1 = ptr->v1 = 0;
	}
	printf("\n");
}

void initialize_cars_bug1(Car *cars, gsl_rng * r){
	
	time_t mytime = time(NULL); /* GET TIME FOR RANDOM SEED */
 	gsl_rng_set(r,(unsigned long int)mytime);

	// setting up two cars alone
	Car *car1, *car2, *car3;
	car1 = cars;
	car2 = car1+1;
	car3 = car1+2;

	// first car
	car1->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
	car1->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
	car1->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
	car1->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
	car1->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
	car1->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
	car1->a_min = -car1->a_max;
	car1->next = car2;

	// second car
	car2->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
	car2->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
	car2->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
	car2->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
	car2->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
	car2->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
	car2->a_min = -car2->a_max;
	car2->next = car3;

	// third car
	car3->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
	car3->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
	car3->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
	car3->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
	car3->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
	car3->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
	car3->a_min = -car3->a_max;
	car3->next = car1;
}

void initialize_cars_two_car_experiment(Car *cars, gsl_rng * r){
	
	time_t mytime = time(NULL); /* GET TIME FOR RANDOM SEED */
 	gsl_rng_set(r,(unsigned long int)mytime);

	// setting up two cars alone
	Car *car1, *car2;
	car1 = cars;
	car2 = car1+1;

	// first car
	car1->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
	car1->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
	car1->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
	car1->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
	car1->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
	car1->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
	car1->a_min = -car1->a_max;
	car1->next = car2;

	// second car
	car2->car_size = (int)3*gsl_rng_uniform(r) + Car_size;
	car2->tao = gsl_ran_gaussian(r, .05*Tao) + Tao;
	car2->v_opt = gsl_ran_gaussian(r, .1*V_opt) + V_opt;
	car2->v_steep = gsl_ran_gaussian(r, .05*V_steep) + V_steep;
	car2->alpha = gsl_ran_gaussian(r, .05*Alpha) + Alpha;
	car2->a_max = gsl_ran_gaussian(r, .1*A_max) + A_max;
	car2->a_min = -car2->a_max;
	car2->next = car1;
}

// getting 2 cars set up correctly
void setup_Initial_Conditions_bug1(Car *cars){
	printf("time x1 v1 x2 v2 x3 v3\n");
	// set up two cars
	Car *ptr = cars;
	ptr->x = 0;
	ptr->v = 20;

	// next car
	ptr++;
	ptr->x = 50;
	ptr->v = 10;

	// next car
	ptr++;
	ptr->x = 100;
	ptr->v = 5;
}

// getting 2 cars set up correctly
void setup_Initial_Conditions_two_car_experiment(Car *cars){
	printf("time x1 v1 x2 v2\n");
	// set up two cars
	Car *ptr = cars;
	ptr->x = 0;
	ptr->v = 20;

	// next car
	ptr++;
	ptr->x = 50;
	ptr->v = 10;
}

void incorporate_boundary_conditions(int num, Car *cars){
	Car *ptr;
	for (ptr = cars; (ptr-cars) < num; ptr++){
		if(ptr->x1 >= distance){	ptr->x1 -= distance;	 }
	}
}

double collision_time(Car *ptr, Car *ptr1){
	// if (ptr->x == ptr1->x && ptr->v == ptr1->v) // single car case
	// {
	// 	return 0;
	// }
	return ((ptr1->x)-(ptr->x)-(ptr->car_size))/((ptr->v)-(ptr1->v));
}

double update_acceleration(double t, Car *car){
	if (0 < t && t < car->tao){
		return (car->alpha)*(1-(car->tao)/t);
	}
	else{
		double shift = ((car->a_max)+(car->a_min))/2;
		double scale = 0.5*((car->a_max)-(car->a_min));
		double theta = ((car->v_opt)-(car->v))/(car->v_steep);
		
		gsl_complex z = gsl_complex_rect(theta, 0);
		double acc = shift + scale*GSL_REAL(gsl_complex_tanh(z));
		return acc;
	}
}

double update_velocity(Car *car, double a){
	return car->v + Dt*a;
}

double update_position(Car *car){
	return car->x + Dt*(car->v);
}

// NEEDS TO BE REDESIGNED
// Has been redesigned for a circular track
int collision_detection(int num, double pos_1[], double time_){
	double m = min(num, pos_1);
	int min_index = index_of(m, num, pos_1);
	double buf[num];
	memcpy(buf, pos_1, num*sizeof(double));
	// copying everything after min_index to front of buf array
	int i;
	for (i = min_index; i < num; ++i)
	{
		buf[i-min_index] = pos_1[i];
	}
	// copying everything from start to min_index to back of array
	for(i = 0; i < min_index; ++i){
		buf[(num - min_index) +i] = pos_1[i];
	} 
	if (is_monotonic(num, buf)){	return 1;	}
	else{	return 0;	}
}

// My code
int gcd(int a, int b){
	int r = a%b;
	while(a%b != 0){
		r = a%b;
		//int p = (a-r)/b;
		a = b;
		b = r;
	}
	return r;
}

/*
*
* SOURCE CODE FOUND ON
* http://www.geeksforgeeks.org/array-rotation/ 
*
**/

/*Function to left rotate arr[] of siz n by d*/
void leftRotate(int arr[], int d, int n)
{
  int i, j, k, temp;
  for (i = 0; i < gcd(d, n); i++)
  {
    /* move i-th values of blocks */
    temp = arr[i];
    j = i;
    while(1)
    {
      k = j + d;
      if (k >= n)
        k = k - n;
      if (k == i)
        break;
      arr[j] = arr[k];
      j = k;
    }
    arr[j] = temp;
  }
}

/* 
* 
* SOURCE CODE FOUND ON
* http://www.geeksforgeeks.org/a-program-to-check-if-strings-are-rotations-of-each-other-or-not/
*
*/

/* Function checks if passed strings (str1 and str2)
   are rotations of each other */
int areRotations(char *str1, char *str2)
{
  int size1   = strlen(str1);
  int size2   = strlen(str2);
  char *temp;
  void *ptr;
 
  /* Check if sizes of two strings are same */
  if (size1 != size2)
     return 0;
 
  /* Create a temp string with value str1.str1 */
  temp   = (char *)malloc(sizeof(char)*(size1*2 + 1));
  temp[0] = '\0';
  strcat(temp, str1);
  strcat(temp, str1);
 
  /* Now check if str2 is a substring of temp */
  ptr = strstr(temp, str2);
 
  free(temp); // Free dynamically allocated memory
 
  /* strstr returns NULL if the second string is NOT a
    substring of first string */
  if (ptr != NULL)
    return 1;
  else
    return 0;
}

/*
*   OLD CODE
*
**/
// , *car3;
// 		car1 = cars;
// 		car2 = car1->next;
// 		car3 = car2->next;

// 		// update car 1
// 		t_collision = collision_time(car1, car2);
// 		acc = update_acceleration(t_collision, car1);
// 		car1->v1 = update_velocity(car1, acc);
// 		car1->x1 = update_position(car1);
// 		// fprintf(stderr, "Car1\n");
// 		// fprintf(stderr, "collision time: %.3f\n", t_collision );
// 		// fprintf(stderr, "acc: %.3f\n", acc );

// 		// update car 2
// 		t_collision = collision_time(car2, car3);
// 		acc = update_acceleration(t_collision, car2);
// 		car2->v1 = update_velocity(car2, acc);
// 		car2->x1 = update_position(car2);
// 		// fprintf(stderr, "Car2\n");
// 		// fprintf(stderr, "collision time: %.3f\n", t_collision );
// 		// fprintf(stderr, "acc: %.3f\n", acc );

// 		// update car 3
// 		t_collision = collision_time(car3, car1);
// 		acc = update_acceleration(t_collision, car3);
// 		car3->v1 = update_velocity(car3, acc);
// 		car3->x1 = update_position(car3);
// 		// fprintf(stderr, "Car3\n");
// 		// fprintf(stderr, "collision time: %.3f\n", t_collision );
// 		// fprintf(stderr, "acc: %.3f\n", acc );

// 		// printf(" %.2f %.2f %.2f %.2f %.2f %.2f\n", car1->x, car1->v, car2->x, car2->v, car3->x, car3->v);

