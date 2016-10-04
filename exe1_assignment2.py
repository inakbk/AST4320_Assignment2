# exe 1, assignment 2

from pylab import *
from scipy.interpolate import UnivariateSpline

R = 1. #something
N = 1000

def W(k):
	return (1./(pi*k))*sin(2*pi*R*k)

k = linspace(-10,10,N)

W_max = max(W(k))
W_min = min(W(k))

#finding the FWHM
k2 = linspace(-1,1,N)
r = UnivariateSpline(k2, W(k2))
r1, r2 = r.roots()/2.

plot(k, W(k))
plot(r1, W(r1), 'ro')
plot(r2, W(r2), 'ro')
xlabel('$k$', fontsize=22)
ylabel('$W^*(k)$', fontsize=22)
title('The fourier transform of the smoothing function W(k) in 1D')
legend(['$W^*(k)$', 'FWHM'], fontsize=22)
show()


