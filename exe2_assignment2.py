# exe 2, assignment 2

from pylab import *
from random import *

N = 1e5 #nr of random walks

#first values
mean = 0 #defined to zero
sigma = 0.5*1e-4 #< given in the exe.
S_c = sqrt(sqrt(pi/sigma))

delta = gauss(mean, sigma) #gauss(mu, sigma), Gaussian distribution. mu is the mean, and sigma is the standard deviation.
print delta

