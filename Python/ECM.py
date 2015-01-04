"""
My implementation of the 'elementary' ECM factorization method.
"""

import random
from fractions import gcd

def dec2bin(n):
	"""
	Converts a decimal number to a binary array.
	
	eg: dec2bin(9) = [1,0,0,1]
	"""
	
	bin = []
	
	if (n == 0):
		return [0]
		
	while (n > 0):
		bin.insert(0, n % 2)
		n = n >> 1
	return bin


def ext_euc(a, n):
	
    """
    Does the extended Euclidean algorithm to find a triple: (g, x, y), 
    such that ax + ny = g = gcd(a,n). 
    """
    if a == 0:
        return (n, 0, 1)
    else:
        g, y, x = ext_euc(n % a, a)
        return (g, x - (n // a) * y, y)
        
        
def check_mod_inverse(a,n):
	"""
	Attempts to compute the modular inverse of a mod n. Since n is not
	necessarily prime, this process will hopefully uncover a non-trivial
	factor of n.
	
	If a factor, p, is found, it will return the pair [p, True]
	otherwise it will return the inverse [a^-1, False]
	"""
	
	g,x,y = ext_euc(a,n)
	
	# Non-trivial factor has been found
	if g != 1:
		return [g, True]
	
	return [x % n, False]
	
	
	
def add(E,P1,P2,n):
	"""
	Adds the affine coordinates P1 and P2 on the affine curve E_a,b(Z_n)
	"""
	
	m = 0
	a,b = E
	X1, Y1, Z1 = P1
	X2, Y2, Z2 = P2
	
	#If either point is the identity
	if (Z1 == 0):
		return P2
	
	if (Z2 == 0):
		return P1
	
	#if the x coordinates are the same
	if (X1 == X2):
		
		#If they are inverses
		if ((Y1 + Y2) % n == 0):
			return (0,1,0)
		
		#Here is where the check happens
		find_fac = (check_mod_inverse(2 * Y1,n))
		
		if(find_fac[1]):
			return 'Factor Found! : ' + str(find_fac[0])
		
		m = (3 * (X1 ** 2) + a) * find_fac[0]
		
	else:
		#Here is where the check happens
		find_fac = (check_mod_inverse(2 * Y1,n))
		
		if(find_fac[1]):
			return 'Factor Found! : ' + str(find_fac[0])
			
		m = (Y2 - Y1) * find_fac[0]
		
	X3 = m ** 2 - X1 - X2
	
	return (X3 % n, (m * (X1 - X3) - Y1) % n,1)
	
	
def double(E,P,n):
	"""
	Doubles the affine coordinate P = (x,y,1) on E defined over Z_n
	"""
	
	return add(E,P,P,n)
	

def neg(E,P,n):
	"""
	Returns the inverse of the affine point P
	"""
	x,y,z = P
	return (x,(-1)*y % n,z)
	
	
def subtract(E,P1,P2,n):
	"""
	Returns the difference of the affine points P1 - P2
	"""
	return add(E,P1,neg(E,P2,n),n)
	

def mult_triv(E,P,k,n):
	"""
	Does the multiplication on E trivially by repeatedly adding P to 
	itself
	"""
	if (k == 0):
		return (0,1,0)
		
	Q = P
	for x in xrange(k):
		Q = add(E,Q,P,n)
		if (type(Q) == type('str')):
			return Q
			break
		
		
	return Q


def multiply(E,P,k,n):
	"""
	Multiplies the point P by k \in Z on E defined over n using a
	multiplication ladder while simultaneously checking to see of the 
	process of checking for modular inverses has broken down or not.
	
	If it has, it will return the a factor of n
	"""
	if (k == 0):
		return (0,1,0)
	Q = P
	
	# Binary representation of M = 3n and N, each B-bits, N padded on the
	# left
	M = dec2bin(3 * n)
	N = dec2bin(n)
	B = len(M)
	while B != len(N):
		N.insert(0,0)
		
	# Comparing bits and doing the addition/subtraction, if at any point 
	# a factor is found, we return the string representation of the factor
	for j in range(B-2,0,-1):
		Q = double(E,Q,n)
		if (type(Q) == type('str')):
			return Q
		
		if ((M[j],N[j]) == (1,0)):
			Q = add(E,Q,P,n)
			if (type(Q) == type('str')):
				return Q
		
		if ((M[j],N[j]) == (0,1)):
			Q = subtract(E,Q,P,n)
			if (type(Q) == type('str')):
				return Q
	
	return Q
				
		
	
	
def find_curve(n):
	
	""" 
	Finds an elliptic curve and a point with the appropriate 
	properties to factor the integer n.
	
	Returns: (A,B) that describe a curve E: y^2 = x^3 + A x + B and 
	the point (x,y) on E in the form of a list [(A,B),(x,y)]
	"""
	
	# Choose random x,y,a in  [0,n-1]
	x = random.randint(0,n)
	y = random.randint(0,n)
	a = random.randint(0,n)
	
	# Define b, g
	b = ((y ** 2) % n - (x ** 3) % n  - (a*x) % n)
	
	# The Discriminant
	disc = (4 * (a ** 3)) + (27 * (b **2))
	
	g = gcd(disc,n)
	
	# g cannot be n
	if (g == n):
		return find_curve(n)
	
	# if g is not one, we have found a non-trivial factor
	if (g > 1):
		return 'Factor Found! : ' + str(g)
		
	return [(a,b),(x,y,1)]
	
	
def factor(n,B):
	"""
	Factors n using the basic ECM-algorithm. B is the stage 1 limit
	"""
	
	# Ensure gcd(n,6) = 1
	if (n % 2 == 0):
		return "Factor found! : " + str(2)
		
	if (n % 3 == 0):
		return "Factor found! : " + str(3)
		
	# Find a suitable Elliptic Curve
	
	temp = find_curve(n)
	if (type(temp) == type('str')):
		return temp
		
	E,P = temp
	
	# A list of the first 100000 prime numbers
	f = open('primes.txt')
	primes = []	
	while len(primes) < B:
		line = f.readline()
		q = line.split()
		for x in q:
			primes.append(int(x))
			
	f.close()		
	for p_i in primes:
		a = 1
		while pow(p_i,a) <= B:
			a += 1
			
		
		for j in xrange(a):
			
			P = mult_triv(E,P,p_i,n)
			
			if (type(P) == type('str')):
				return P
				break
	
	return "Factorization Failed... Increase B?"
	

if __name__ == '__main__':
	print('Enter n: ')
	n = input()
	print('Enter Threshold: ')
	B = input()
	print factor(n,B)
	

	
	
