all: traffic model

traffic: traffic.c tools.c
	gcc -g -Wall -I/usr/local/include traffic.c tools.c -o traffic -lgsl -lgslcblas -lm

model: model.py
	./traffic > results; python model.py