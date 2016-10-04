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
sigma = 0.5*1e-4 #< given in the exe.
S_c[0] = sqrt(sqrt(pi/sigma))
delta_density[0] = gauss(mean, sigma) #drawing a random nr from Gaussian distribution w. a mean, and sigma is the standard deviation.
#print delta



#--------------------------------------
#first plotting the PDF:

density = linspace(-0.0002,0.0002,1000) # picking small density pertubations since the sigma is so small
plot(density, PDF(density, mean, sigma))
#show()

#--------------------------------------
#random walks





















