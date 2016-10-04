# exe 2, assignment 2

from pylab import *
from random import *

#The gaussian PDF
def PDF(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) * exp(-0.5*(x-mean)**2/sigma**2)

N = 10#1e5 #nr of random walks

#--------------------------------------
#random walks
epsilon = 0.5 # changes S_c and will decide how fast the random walk converges
realization_time = 0
maximum_number_of_iterations = 500

delta_density = zeros(maximum_number_of_iterations) #array to store the density values
S_c = zeros(maximum_number_of_iterations)
final_delta_density = zeros(N) #storing the density values after realization to make histogram

for k in range(N-1):
	#first values:
	mean = 0 #defined to zero
	sigma_old = 0.5*1e-4 #< given in the exe.
	S_c[0] = sqrt(sqrt(pi/sigma_old))
	delta_density[0] = gauss(mean, sigma_old) #drawing a random nr from Gaussian distribution w. a mean, and sigma is the standard deviation.

	for i in range(maximum_number_of_iterations):
		S_c[i+1] = S_c[i] - epsilon
		sigma_new = pi/S_c[i+1]**4 
		print S_c[i]
		beta = gauss(mean, sqrt(sigma_new**2 - sigma_old**2))
		delta_density[i+1] = delta_density[i] + beta

		sigma_old = sigma_new #updating sigma
		print realization_time
		realization_time += 1

		if S_c[i] < 1:
			print "Realization happened! realization_time= ", realization_time
			print "S_c= ", S_c[i]
			final_delta_density[k] = delta_density[i+1]
			break
		if realization_time > maximum_number_of_iterations:
			print "Maximum nr of iterations reached!!! (", maximum_number_of_iterations, ")"
			final_delta_density[k] = delta_density[i+1]
			break

#	figure(1)
	plot(S_c[0:realization_time], delta_density[0:realization_time])
	plot(S_c[0], delta_density[0], 'bo-') #plotting first point blue

plot(S_c[0:realization_time+1], 1 + zeros(realization_time+1), 'g')
xlabel('S_c')
ylabel('delta_density')
#axis([0.5,S_c[0],-5,5])
title('Random walk of the density versus S_c')
show()


#--------------------------------------
#plotting the PDF:

density = linspace(-0.0002,0.0002,1000) # picking small density pertubations since the sigma is so small
#figure(2)
#plot(density, PDF(density, mean, sigma_old))
#show()



