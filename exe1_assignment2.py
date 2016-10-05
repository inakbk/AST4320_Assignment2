# exe 1, assignment 2

from pylab import *
from scipy.interpolate import UnivariateSpline

R = 1. #something
N = 10000

def W(k):
	return (1./(pi*k))*sin(2*pi*R*k)

k = linspace(-10,10,N)

#finding the FWHM
k2 = linspace(-1,1,N)
W_max = max(W(k2))
HM = W_max/2 #half maximum

#finding the widt very approximately
for i in range(N):
	if W(k[i]) > HM:
		print i
		r1 = k[i-1]
		break
for j in range(i,N):
	if W(k[j]) < HM:
		print j
		r2 = k[j]
		break

plot(k, W(k))
plot(r1, HM, 'ro')
plot(r2, HM, 'ro')
k3 = linspace(r1,r2,N)
plot(k3, W(r1 + zeros(N)))
xlabel('$k$', fontsize=22)
ylabel('$W^*(k)$', fontsize=22)
title('The fourier transform of the smoothing function W(k) in 1D', fontsize=18)
legend(['$W^*(k)$', 'FWHM'], fontsize=22)
show()


