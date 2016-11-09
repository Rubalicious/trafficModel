import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
from turtle import *
import turtle
import pandas as pd 
from time import ctime

# # to initialize random velocities for cars
# def random_IC(sample_num, mew, sigma):
# 	return np.random.normal(mew, sigma, sample_num)

# use turtle graphics for visual representation
def animate(t, x):
	# first we draw the road
	wn = turtle.Screen()
	bob = turtle.Turtle()
	bob.goto(300, 0)
	bob.goto(-300, 0)
	bob.penup()

	bob.shape("circle")
	# now we draw our car
	for dt, dx in zip(t, x):
		bob.goto(600*dt-300, 10)
		bob.stamp()


# read from a file and plot results
def main():

	file_name = "results"

	df = pd.read_csv(file_name, sep=' ') # for the time being

	times = df['time'].values
	pos_cols = np.array([col_name for col_name in df.columns if 'x' in col_name])
	vel_cols = np.array([col_name for col_name in df.columns if 'v' in col_name])

	for col in vel_cols:
		plt.plot(times, df[col].values)

	speedlimit = np.array([40 for Dt in times])
	# plt.plot(times, speedlimit, 'r--')
		
	# animate(t, x)

	print "Initial Conditions:\n%s"%df.loc[0]

	plt.title("Velocity")
	plt.xlabel("time in seconds")
	plt.ylabel("velocity in [m/s]")
	# plt.legend({'car 1', 'car 2', 'car 3'}, loc=4)
	# plt.show()
	plt.figure()
	# plt.savefig('pictures/velocities_for_%d_cars_%s.jpg'%(len(pos_cols), ctime().replace(' ', '_').replace(':','-')))

	for col in pos_cols:
		plt.plot(times, df[col].values)

	plt.title("Position")
	plt.xlabel("time in seconds")
	plt.ylabel("position in [m]")
	# plt.legend({'car 1', 'car 2', 'car 3'}, loc=4)
	# plt.figure()
	# plt.savefig('pictures/positions_for_%d_cars_%s.jpg'%(len(pos_cols), ctime().replace(' ', '_').replace(':', '-')))

	# Plot of Initial Conditions
	# x0 = df.loc[0][1::2]
	# v0 = df.loc[0][2::2]
	# plt.plot(range(len(x0)), x0)
	# plt.show()
	# plt.plot(range(len(v0)), v0)
	plt.show()


if __name__ == '__main__':
	main()