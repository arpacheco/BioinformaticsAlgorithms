#!/usr/bin/env python

""" K-means clustering script for BE 562. Inputs: a text file with clusters
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
def assign_clusters(points,mu):
	clusters = []
	for point in points:
		distances = []
		for center in mu:
			distances.append(numpy.linalg.norm(point-center)) #calculate all the distances from each point to each cluster
		clusters.append(distances.index(min(distances)))
	return clusters #in order of the points provided originally

#assign the points to clusters using weights
def assign_clusters_fuzzy(points,mu):
	clusters = []
	all_distances = [] #saves distances per point between point and centroid of each cluster
	pointcount = 0
	for point in points:
		distances = []
		for center in mu:
			distance = numpy.linalg.norm(point-center)
			distances.append(distance) #calculate all the distances from each point to each cluster
		clusters.append(distances.index(min(distances)))
		all_distances.append(min(distances))
		pointcount += 1

	#get the probabilities of each point per center, using a scaled (0-1) distance metric to generate a probability of a point being near its most likely cluster
	probabilities = []
	min_distance = min(all_distances)
	max_distance = max(all_distances)
	for distance in all_distances:
		probabilities.append(1-(distance-min_distance)/(max_distance-min_distance))

	#scale all the points by their probabilities
	scaledpoints = []
	pointcount = 0
	for point in points:
		scaledpoints.append([])

		for dimension in range(M):
			scaledpoints[pointcount].append((point[dimension]*probabilities[pointcount]))
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
	for i in range(len(sorted_by_cluster)):
		new_mu.append([])
		for m in range(M): #for each dimension
			m_values = [] ##x values, y values,...,n-dim values
			for point in sorted_by_cluster[i]:
				m_values.append(point[m]) #get the m-dim value
			m_average = numpy.mean(m_values) #and take the average
			new_mu[i].append(m_average)
	
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
		print 'Cluster',str(k),'mean vector:','\t',rawmu[k]
	print
	#print
	#print 'Data points and cluster assignments:'
	#for i in range(len(rawpoints)):
	#	print d[i],'--> Cluster',rawpoints[i]

#plot if input data is 2D
def plot(f_mu,f_clusters):
	if M == 2 and K <= 7:
		print 'plotting...'
		colors = ['red','blue','green','yellow','cyan','magenta','black']
		sorted_by_cluster = sort_by_clusters(d,f_clusters)


		#plot "correct" clusters from original distribution
		startpoint = 0
		for k in range(3): #number of clusters originally prescribed
			color = colors[k]
			startpoint = k*50
			
			xvals = []
			yvals = []
			for i in range(startpoint,startpoint+50):
				xvals.append(d[i][0])
				yvals.append(d[i][1])
			plt.scatter(xvals,yvals,color=color,marker = 'o',s=50)

		#plot predicted clusters
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

#make histogram for large values of K
def histogram(f_clusters):
	print 'plotting histogram...'
	clusternums = []
	clustercounts = []
	for i in range(K):
		clustercount = 0
		for j in range(len(f_clusters)):
			if f_clusters[j] == i:
				clustercount += 1
		clustercounts.append(clustercount)
		clusternums.append(i)

	plt.bar(clusternums,clustercounts)
	plt.xlabel('Cluster Number')
	plt.ylabel('Number of points in cluster')
	plt.show()

#run the kmeans algorithm
def run_Kmeans(d,k):
	current_iteration = 0
	max_iterations = 10

	initial_mu = initialize_clusters(d,K)
	initial_clusters = assign_clusters(d,initial_mu)

	current_mu = redefine_centers(d,initial_clusters)
	old_mu = current_mu
	current_clusters = initial_clusters

	while not converged(old_mu,current_mu,current_iteration,max_iterations):
		old_mu = current_mu
		current_mu = redefine_centers(d,current_clusters)
		current_clusters = assign_clusters(d,current_mu)
		current_iteration += 1
	format_output(current_mu,current_clusters)
	return current_mu, current_clusters

def run_fuzzy_Kmeans(d,k):
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

#For part e, report on the cluster of the most affected gene
def find_max(d,clusters):
	pointmaxes = []
	for point in d:
		pointmaxes.append(max(point))
	clusternum = clusters[pointmaxes.index(max(pointmaxes))]

	others = []
	for i in range(len(clusters)):
		if clusters[i] == clusternum:
			others.append(i)
	print 'Maximum:', max(pointmaxes), 'Index:',pointmaxes.index(max(pointmaxes)), 'Cluster:',clusternum
	print 'Other genes in cluster '+str(clusternum)+': '+str(others)


#Execute
file_name = sys.argv[1]
K = int(sys.argv[2])
d = get_points(file_name)
M = get_dimensions(d)

final_mu, final_clusters = run_Kmeans(d,K)

