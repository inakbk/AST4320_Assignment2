# exe 2, assignment 2

from pylab import *
from random import *
import sys

#The gaussian PDF
def PDF(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) * exp(-0.5*(x-mean)**2/sigma**2)

delta_crit = 1
def PDF_nc(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) *( exp(-0.5*(x-mean)**2/sigma**2) - exp(-0.5*(2*delta_crit - x)**2/sigma**2) ) 

#--------------------------------------
#random walks

N = 1000 #1e5 nr of random walks! :)

epsilon = 0.99 # changes S_c and will decide how fast the random walk converges (small=too fast convergence)
maximum_number_of_iterations = 500
abort_program = 0

final_delta_density = zeros(N) #storing the density values after realization to make histogram
threshold_delta_density = []
mean = 0 #defined to zero
sigma_start = 0.5*1e-4 #< given in the exe.
S_cstart = sqrt(sqrt(pi/sigma_start))

figure(1)
for k in range(N):
	#values reset:
	delta_density = zeros(maximum_number_of_iterations) #array to store the density values to plot random walk
	S_c = zeros(maximum_number_of_iterations)
	realization_time = 0

	sigma_old = sigma_start
	S_c[0] = S_cstart

	#should draw new for each walk or same starting for each walk?
	delta_density[0] = gauss(mean, sigma_old) #drawing a random nr from Gaussian distribution w. a mean, and sigma is the standard deviation.

	#one random walk / realization:
	for i in range(maximum_number_of_iterations):
		S_c[i+1] = S_c[i] * epsilon
		sigma_new = pi/S_c[i+1]**4 
		beta = gauss(mean, sqrt(sigma_new**2 - sigma_old**2))
		delta_density[i+1] = delta_density[i] + beta

		sigma_old = sigma_new #resetting sigma
		realization_time += 1

		if S_c[i] < 1:
			#print "Realization happened! realization_time= ", realization_time
			#print "S_c= ", S_c[i]
			final_delta_density[k] = delta_density[i+1]
			if final_delta_density[k] < 1:
				threshold_delta_density.append(final_delta_density[k])
			#print "density= ", final_delta_density[k]
			#print "-------"
			break
		if realization_time +1 >= maximum_number_of_iterations:
			print "Maximum nr of iterations reached!!! (", maximum_number_of_iterations, ")"
			final_delta_density[k] = delta_density[i+1]
			print "density= ", final_delta_density[k]
			print "-------"
			abort_program += 1
			break

	plot(S_c[0:realization_time+1], delta_density[0:realization_time+1])
	if abort_program >= 10:
		print "Program aborted! Maximum number of iterations exceeded in the loop 10 times!"
		sys.exit()

print "last realization_time: ", realization_time
print "S_c at last realization: ", S_c[i]

loglog(S_c[0:realization_time], delta_crit + zeros(realization_time), 'ro') #plotting delta_crit
xlabel('S_c')
ylabel('delta_density')
axis([0.5,S_c[0],-5,5])
title('Random walk of the density versus S_c')
#show()

#--------------------------------------
#plotting the normalized histogram and PDF for all densities:
figure(2)
hist(final_delta_density, normed=1, bins=100)

density_PDF = linspace(-15, 15, 1000) 
sigma_PDF = sqrt(pi) # since S_c=1 at the end of the chain?
plot(density_PDF, PDF(density_PDF, mean, sigma_PDF), 'r')

#--------------------------------------
#plotting the normalized histogram and PDF for all densities under delta_crit:
figure(3)
threshold_delta_density = array(threshold_delta_density)
hist(threshold_delta_density, normed=1, bins=800)

#other PDF:
plot(density_PDF, PDF_nc(density_PDF, mean, sigma_PDF), 'r')

show()



