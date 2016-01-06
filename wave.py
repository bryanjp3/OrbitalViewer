
import math 
import numpy as np
from scipy.special import sph_harm
from scipy.special import factorial

import matplotlib.pyplot as plt

def fac(n):
	return factorial(n)
"""
computes associated Laguerre polynomial coefficients for l and n
and returns the polynomial function of input x
http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html
"""
def associatedLaguerre(l, n, x):
	coeffs = []
	for k1 in range(n-l-1+1):
		c = (-1)**k1*fac(n+l)**2/fac(n-l-k1-1)/fac(2*l+1+k1)/fac(k1)
		coeffs.append(c)
		
	out = [0.0 for z in range(len(x))] # output
	for m in range(len(coeffs)):
		#print m, " ", coeffs[m], "jgfk"
		out = out + coeffs[m]*(x**m)
	return out

# define constants
eps0 = 8.85e-14
pi = 3.141592658979323
hbar = 6.626e-34
a = .529e-10

# https://en.wikipedia.org/wiki/Quantum_number
# principal quantum number: n is a postive number
#   angular quantum number: l = 0,1,2,...,n-1
#  magnetic quantum number: m = -l, -l+1,...,l-1,l

l = 3
m = 1
n = 6

# angle dependence of psi
theta = np.linspace(0.0, 2*np.pi, 100)
phi = theta
s = sph_harm(m, l, theta*0.0, phi)


# radial dependence of psi
R1 = np.sqrt((2.0/n/a)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3 ))

r = np.linspace(0.0, 4*a, 10)
R2 = np.exp(-r/n/a)*(2*r/n/a)**l*associatedLaguerre(l, n, 2*r/n/a)

#test = 1/np.sqrt(2)*a**(-1.5)*(1-.5*r/a)*np.exp(-r/2/a)
#test = 8.0/27.0/np.sqrt(6)*a**(-1.5)*(1-r/a/6)*(r/a)*np.exp(-r/3/a)

def psi(r1, theta1, phi1):
	
	R2 = np.exp(-r1/n/a)*(2*r1/n/a)**l*associatedLaguerre(l, n, 2*r1/n/a)
	s = sph_harm(m, l, theta1, phi1)
	return R1*R2*s

# Screen coordinates: x0, y0, x1, y1.
S = (-25.0*a, -25.0*a, 25.0*a, 25.0*a)

# Loop through all pixels.
w = 500
h = 500

r = []
theta = []
coords = []
pixels = []
for i, x in enumerate(np.linspace(S[0], S[2], w)):
    for j, y in enumerate(np.linspace(S[1], S[3], h)):
		#print (x,y)
		pixels.append((i,j))
		rad = np.sqrt(x**2 + y**2)
		angle = math.atan2(x,y)
		#print (rad, angle*180/pi)
		coords.append((rad, angle))
		r.append(rad)
		theta.append(angle)

r = np.array(r)
theta = np.array(theta)
phi = np.array(theta)

wave = psi(r, theta*0.0, phi)
wave = np.abs(wave)*np.abs(wave)
wave = wave/np.sum(wave)

image = np.zeros((w,h,3))

for c in range(len(pixels)):
	pixel = pixels[c]
	x = pixel[0]
	y = pixel[1]
	image[x,y,:] = np.array([1.0, 1.0, 1.0])*wave[c]/np.amax(wave)
	

plt.imsave('fig2.png', image)
plt.show()








