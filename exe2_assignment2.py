# exe 2, assignment 2

from pylab import *
import math
import random
import sys

#The gaussian PDF
def PDF(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) * exp(-0.5*(x-mean)**2/sigma**2)

delta_crit = 1
def PDF_nc(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) *( exp(-0.5*(x-mean)**2/sigma**2) - exp(-0.5*(2*delta_crit - x)**2/sigma**2) ) 

#--------------------------------------
#random walks

N = 10000 #1e5 nr of random walks! :)

epsilon = 0.99 # changes S_c and will decide how fast the random walk converges (small=too fast convergence)
maximum_number_of_iterations = 500
abort_program = 0

final_delta_density = zeros(N) #storing the density values after realization to make histogram
threshold_delta_density = []
mean = 0 #defined to zero
sigma_start = sqrt(0.5*1e-4) #< given in the exe.
S_cstart = (pi/sigma_start**2)**(1./4.)

figure(1)
for k in range(N):
	#values reset:
	delta_density = zeros(maximum_number_of_iterations) #array to store the density values to plot random walk
	S_c = zeros(maximum_number_of_iterations)
	realization_time = 0

	sigma_old = sigma_start
	S_c[0] = S_cstart

	#should draw new for each walk or same starting for each walk?
	delta_density[0] = random.gauss(mean, sigma_old) #drawing a random nr from Gaussian distribution w. a mean, and sigma is the standard deviation.

	#one random walk / realization:
	for i in range(maximum_number_of_iterations):
		S_c[i+1] = S_c[i] * epsilon
		sigma_new = sqrt(pi/S_c[i+1]**4)
		beta = random.gauss(mean, sqrt(sigma_new**2 - sigma_old**2))
		delta_density[i+1] = delta_density[i] + beta
		#print sigma_old
		#print sigma_new
		#exit()

		sigma_old = sigma_new #resetting sigma
		realization_time += 1

		if S_c[i+1] < 1:
			#print "Realization happened! realization_time= ", realization_time
			#print "S_c= ", S_c[i+1]
			#print "sigma at realization minus sqrt(pi): ", sqrt(pi/S_c[i+1]**4)-sqrt(pi)
			final_delta_density[k] = delta_density[i+1]
			if max(delta_density) < 1:
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

	loglog(S_c[0:realization_time+1], delta_density[0:realization_time+1])
	loglog(S_c[0], delta_density[0], 'ro') #plotting starting point
	#plot(S_c[realization_time], delta_density[realization_time], 'ko') #plotting ending point
	
	if abort_program >= 10:
		print "Program aborted! Maximum number of iterations exceeded in the loop 10 times!"
		sys.exit()

print "last realization_time: ", realization_time
print "S_c at last realization: ", S_c[i+1]

loglog(S_c[0:realization_time], delta_crit + zeros(realization_time), 'ro') #plotting delta_crit
xlabel('S_c')
ylabel('delta_density')
axis([0.5,S_c[0],-5,5])
title('Random walk of the density versus S_c')

#--------------------------------------
#plotting the normalized histogram and PDF for all walkers at realization:
figure(2)
hist(final_delta_density, normed=1, bins=42, color='b')

density_PDF = linspace(-9, 9, 1000) 
sigma_PDF = sqrt(pi) # since S_c=1 at the end of the chain
plot(density_PDF, PDF(density_PDF, mean, sigma_PDF), 'r', linewidth=4)
xlabel('delta_density')
ylabel('# counts / probability')
legend(['Numerical', "Analytical"])
title('Histogram of the density for $S_C = 1$ together with the analytical \n solution given.')

#--------------------------------------
#plotting the normalized histogram and PDF for all densities under delta_crit:
figure(3)
threshold_delta_density = array(threshold_delta_density)
hist(threshold_delta_density, normed=True, bins=28, color='g')

#the other PDF:
raw = PDF_nc(density_PDF, mean, sigma_PDF)
#setting values over delta_crit to zero
for i in range(len(density_PDF)):
	if density_PDF[i] >= delta_crit:
		raw[i] = 0
normalization_constant = math.erf(1./(sqrt(2)*sigma_PDF)) #calculated analytically

plot(density_PDF, raw/normalization_constant, 'r', linewidth=4)
xlabel('delta_density < delta_crit')
ylabel('# counts / probability')
legend(['Numerical', "Analytical"])
title('Histogram of the density pertubations that never crossed $\delta_{crit}$ \n for $S_C = 1$ together with the analytical solution given.')


show()
