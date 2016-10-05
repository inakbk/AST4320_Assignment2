# exe 1, assignment 2

from pylab import *
from scipy.interpolate import UnivariateSpline

R = 1. #something
N = 1000

def W(k):
	return (1./(pi*k))*sin(2*pi*R*k)

k = linspace(-10,10,N)

#finding the FWHM
k2 = linspace(-1,1,N)
r = UnivariateSpline(k2, W(k2))
print r.roots()



r1, r2 = r.roots()/2. 

W_max = max(W(k2))
HM = W_max/2 #half maximum

for i in range(N):
	if W(k[i]) > HM:
		print i
		r1, r2 = k[i]
		break

plot(k, W(k))
plot(r1, HM, 'ro')
plot(r2, HM, 'ro')
xlabel('$k$', fontsize=22)
ylabel('$W^*(k)$', fontsize=22)
title('The fourier transform of the smoothing function W(k) in 1D')
legend(['$W^*(k)$', 'FWHM'], fontsize=22)
show()


