#!/usr/bin/env python

""" fuzzy K-means clustering script for BE 562. Inputs: a text file with clusters
and the number of clusters desired. Output: cluster coordinates """

import sys
import csv
import random
import numpy
import matplotlib.pyplot as plt

#get and format the data points for processing
def get_points(fname):
	with open(fname) as f:
		reader = csv.reader(f, delimiter="\t")
		raw = list(reader) # Create list of lists, each one being a row in the original file.

	#format all values
	points = []
	for i in range(len(raw)):
		points.append([])
		for j in range(len(raw[i])-1): #remove the last attribute
			points[i].append(float(raw[i][j]))

	#set the points as actual coordinates
	newpoints = []
	for point in points:
		newpoints.append(numpy.array(point))
	return newpoints

#get the dimensions of the data
def get_dimensions(points):
	return len(points[1])

#randomly initialize the centers
def initialize_clusters(points, k):
	return random.sample(points,k)

#assign the points to clusters
def assign_clusters_fuzzy(points,mu):
	clusters = []
	all_distances = [] #saves distances per point between point and centroid of each cluster
	pointcount = 0
	for point in points:
		all_distances.append([])
		distances = []
		for center in mu:
			distance = numpy.linalg.norm(point-center)
			distances.append(distance) #calculate all the distances from each point to each cluster
			all_distances[pointcount].append(distance)
		clusters.append(distances.index(min(distances)))
		pointcount += 1

	#get the probabilities of each point per center
	probabilities = []
	min_distance = min(min(all_distances))
	max_distance = max(max(all_distances))
	for point in all_distances:
		probabilities.append([])
		for distance in point:
			prob = (distance-min_distance)/(max_distance-min_distance)
			probabilities[all_distances.index(point)].append(prob)

	#scale all the points by their probabilities
	scaledpoints = []
	pointcount = 0
	for point in points:
		scaledpoints.append([])
		centercount = 0

		for center in mu:
			scaledpoints[pointcount].append((point[centercount]*probabilities[pointcount][centercount])+0.001)
			centercount += 1
		pointcount += 1

	#format
	newscaledpoints = []
	for scaledpoint in scaledpoints:
		newscaledpoints.append(numpy.array(scaledpoint))

	d = newscaledpoints #update d to the new points
	return clusters

#sort the points into their assigned clusters
def sort_by_clusters(points,clusters):
	sorted_by_cluster = []
	for i in range(K): #number of clusters
		sorted_by_cluster.append([])
		for j in range(len(clusters)):
			if clusters[j] == i:
				sorted_by_cluster[i].append(points[j])
	return sorted_by_cluster

#calculate the new centroids
def redefine_centers(points,clusters):
	sorted_by_cluster = sort_by_clusters(points,clusters)
	new_mu = []
	for cluster in sorted_by_cluster:
		clusternum = sorted_by_cluster.index(cluster)
		new_mu.append([])
		for m in range(M): #for each dimension
			m_values = [] ##x values, y values,...,n-dim values
			for point in cluster:
				m_values.append(point[m]) #get the m-dim value
			m_average = numpy.mean(m_values) #and take the average
			new_mu[clusternum].append(m_average)
	
	#format and return the new centroids
	format_new_mu = []
	for mu in new_mu:
		format_new_mu.append(numpy.array(mu))
	return format_new_mu

#test if convergence has been reached 
def converged(oldmu,currentmu,current_iter,max_iter):
	if current_iter > 1: #allow for at least 1 full iteration
		if current_iter < max_iter: #if the maximum number of iterations has not been reached
		#compare the difference of one value in the old and new centroid coordinates
			if abs(currentmu[0][0] - oldmu[0][0]) < 0.000000005:
				return True
			else:
				return False
		else:
			return True

#format the output for console
def format_output(rawmu,rawpoints):
	for k in range(K):
		print 'Cluster',str(k),'mean vector:',rawmu[k]
	print
	print 'Data points and cluster assignments:'
	for i in range(len(rawpoints)):
		print d[i],'--> Cluster',rawpoints[i]

#plot if input data is 2D
def plot(f_mu,f_clusters):
	if M == 2 and K <= 7:
		print 'plotting...'
		colors = ['red','blue','green','magenta','cyan','yellow','black']
		sorted_by_cluster = sort_by_clusters(d,f_clusters)

		for cluster in sorted_by_cluster:
			xvals = []
			yvals = []

			for i in range(len(cluster)):
	 			xvals.append(cluster[i][0])
	 			yvals.append(cluster[i][1])

	 		current_cluster_index = sorted_by_cluster.index(cluster)
			plt.scatter(xvals,yvals,color=colors[current_cluster_index])
			plt.scatter(f_mu[current_cluster_index][0],f_mu[current_cluster_index][1],color = colors[current_cluster_index],marker = '*',s=200)
		plt.show()
	elif M == 2 and K > 7:
		print 'plot only available for 7 clusters or fewer.'
	else:
		print 'plot only available for 2 dimensions.'

#run the kmeans algorithm
def run_Kmeans(d,k):
	current_iteration = 0
	max_iterations = 10

	initial_mu = initialize_clusters(d,K)
	initial_clusters = assign_clusters_fuzzy(d,initial_mu)

	current_mu = redefine_centers(d,initial_clusters)
	old_mu = current_mu
	current_clusters = initial_clusters

	while not converged(old_mu,current_mu,current_iteration,max_iterations):
		old_mu = current_mu
		current_mu = redefine_centers(d,current_clusters)
		current_clusters = assign_clusters_fuzzy(d,current_mu)
		current_iteration += 1
	format_output(current_mu,current_clusters)
	return current_mu, current_clusters

#Execute
file_name = sys.argv[1]
K = int(sys.argv[2])
d = get_points(file_name)
M = get_dimensions(d)






