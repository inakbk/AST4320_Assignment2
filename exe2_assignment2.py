# exe 2, assignment 2

from pylab import *
from random import *

#The gaussian PDF
def PDF(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) * exp(-0.5*(x-mean)**2/sigma**2)

N = 1e5 #nr of random walks
delta_density = zeros(N) #array to store the density values
S_c = zeros(N)

#first values:
mean = 0 #defined to zero
sigma_old = 0.5*1e-4 #< given in the exe.
S_c[0] = sqrt(sqrt(pi/sigma_old))
delta_density[0] = gauss(mean, sigma_old) #drawing a random nr from Gaussian distribution w. a mean, and sigma is the standard deviation.
print delta_density

#--------------------------------------
#first plotting the PDF:

density = linspace(-0.0002,0.0002,1000) # picking small density pertubations since the sigma is so small
figure(1)
plot(density, PDF(density, mean, sigma_old))
#show()

#--------------------------------------
#random walks
epsilon = 1e-1 # changes S_c and will decide how fast the random walk converges

S_c[1] = S_c[0] - epsilon
sigma_new = pi/S_c[1]**4 #sigma_2(S_2)
beta = gauss(mean, sqrt(sigma_new**2 - sigma_old**2))
delta_density[1] = delta_density[0] + beta

sigma_old = sigma_new #updating sigma



figure(2)
plot(delta_density[0:2], S_c[0:2], 'ro')
show()

















