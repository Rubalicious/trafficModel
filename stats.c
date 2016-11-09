/*
* File: stats.c
* Author: Ruby Abrams
* Purpose: This will read the output of the traffic
* model and calculate some statistics.
*/

#include <stdio.h>
#include <stdlib.h>

int main()
{	

	// double averages[]; 
	size_t sz = 1024;
	char buf[sz];
	int count = 0;
	
	// read first line to find number of cars
	int bytes_read = getline(&buf, &sz, stdin);
	printf("%d\n", bytes_read); 

	while(fgets(buf, sz, stdin)){
		count++;
	}
	printf("%d\n", count);
	return 0;
}
