# exe 2, assignment 2

from pylab import *
from random import *

#The gaussian PDF
def PDF(x, mean, sigma): 
	return (1./(sigma*sqrt(2*pi))) * exp(-0.5*(x-mean)**2/sigma**2)

#first values:
mean = 0 #defined to zero
sigma = 0.5*1e-4 #< given in the exe.
S_c = sqrt(sqrt(pi/sigma))

#--------------------------------------
#first plotting the PDF:

density = linspace(-0.0002,0.0002,1000) # picking small density pertubations since the sigma is so small
plot(density, PDF(density, mean, sigma))
#show()

#--------------------------------------
#random walks


N = 1e5 #nr of random walks



delta = gauss(mean, sigma) #gauss(mu, sigma), Gaussian distribution. mu is the mean, and sigma is the standard deviation.
#print delta















