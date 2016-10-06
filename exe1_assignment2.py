# exe 1, assignment 2

from pylab import *

R = 1. #something
N = 10000

def W(k):
	return (1./(pi*k))*sin(2*pi*R*k)

k = linspace(-2,2,N)

#finding the FWHM
k2 = linspace(-1,1,N)
W_max = max(W(k2))
HM = W_max/2 #half maximum

#finding the widt (very approximately because it is more fun to do it yourself than googling a fancy function)
for i in range(N):
	if W(k[i]) > HM:
		r1 = k[i-1]
		break
for j in range(i,N):
	if W(k[j]) < HM:
		r2 = k[j]
		break

FWHM = r2 - r1
print 'FWHM= ', FWHM

plot(k, W(k))
#plot(r1, HM, 'ro')
#plot(r2, HM, 'ro')
k3 = linspace(r1,r2,N) #plotting the FWHM
plot(k3, W(r1 + zeros(N)), 'r', linewidth=4)
xlabel('$k$', fontsize=22)
ylabel('$W^*(k)$', fontsize=22)
title('The fourier transform of the smoothing function W(k) in 1D', fontsize=18)
legend(['$W^*(k)$', 'FWHM'], fontsize=22)
show()

"""
Output:

"""
