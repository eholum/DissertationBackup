from fractions import gcd
from fractions import Fraction
from math import sqrt
from math import floor
import sys
import random


#Trivial Primality Test
def is_prime(n):
	
	for x in range(2, int(floor(sqrt(n))) + 1):
		
		if n % x == 0:
			return False
			break
	else:
		return True
		

#Converts an integer to its binary representation
def toBinary(n):
	bin = []
	
	if (n == 0):
		return [0]
		
	while (n > 0):
		bin.insert(0, n % 2)
		n = n >> 1
	return bin


#Helper method for Miller rabin
def mr_test(a, n):
	
  	b = toBinary(n - 1)
  	d = 1
  	for i in xrange(len(b) - 1, -1, -1):
  		x = d
  		d = (d * d) % n
  		if d == 1 and x != 1 and x != n - 1:
  			return True # Composite
  		if b[i] == 1:
  			d = (d * a) % n
  		if d != 1:
  			return True # Composite
  	
	return False # Prime


#Miller Rabin Primaility Test s is the number of tests to make. Probability of
#Error is exactly 2^(-s)
def mr_prime(n, s=75):

  for j in xrange(1, s + 1):
    a = random.randint(1, n - 1)
    if (mr_test(a, n)):
      return False # n is complex
  return True # n is prime
  
  
  
#Returns a list of the prime factors of an integer n
def factor(n):
	if mr_prime(n):
		return [(n,1)]
	else:
		factors = []
		for x in range(2,n/2 + 1):
			power = 0
			while(n%x == 0):
				power += 1
				n = n/x
				if(n == 1):
					factors.append((x,power))
					return factors		
			
			if not(power == 0):
				factors.append((x,power))
		return factors	

	
	
#Returns the modular inverse of the integer n mod p. If there is no such 
#integer, returns -1
def ext_euc(a, n):
	
    """
    Returns a triple (g, x, y), such that ax + by = g = gcd(a,b).
    Assumes a, b >= 0, and that at least one of them is > 0.
    Bounds on output values: |x|, |y| <= max(a, b).
    """
    if a == 0:
        return (n, 0, 1)
    else:
        g, y, x = ext_euc(n % a, a)
        return (g, x - (n // a) * y, y)
        

#Returns the inverse of a mod p:
def mod_inverse(a,p):
	g,x,y = ext_euc(a,p)
	if g != 1:
		return None
	
	return x % p

    
	
#Returns the Legendre symbol of an integer x mod p
def legendre(a, p):
        if a%p == 0: return 0
        a = a%p
        x, y, L = a, p, 1
        while 1:
            if x > (y >> 1):
                x = y - x
                if y & 3 == 3: L = -L
            while x & 3 == 0:
                x = x >> 2
            if x & 1 == 0:
                x = x >> 1
                if y & 7 == 3 or y & 7 == 5: L = -L
            if x == 1: return L
            if x & 3 == 3 and y & 3 == 3: L = -L
            x, y = y % x, x


#Returns all of the quadratic residues mod p
def quad_res(p):
	l = []	
	for x in range(p):
		if(legendre(x,p) != -1):
			l.append(x)		
	return l


#Returns all of the quadratic residues and their roots mod p
def quad_roots(p):
	l = []	
	for x in range(p):
		if(legendre(x,p) != -1):
			r = modular_sqrt(x,p)
			l.append((x,r))		
	return l
	
           
#Finds a solution to the equation x^2 = a mod p
#Assuming the riemann hypothesis this will run in polynomial time
def modular_sqrt(a, p):
	if (is_prime(p) == False):
		raise NameError('p not prime!!!')
	
	if (p == 2):
		return a % 2
		
	if legendre(a,p) == -1:
		return -1
		
	a = a % p

	#Check Simple Cases
	check = p%8
	if (check == 3 or check ==7):
		return int(pow(a,(p+1)/4,p))
	
	elif (check == 5):
		x = pow(a,((p+3)/8),p)
		c = x**2 % p
		if (c != a % p):
			temp = pow(2,(p-1)/4,p)
			x = x*temp % p
		return int(x)
		
	# check = 1 mod p
	else:
		list = range(2,p)
		
		d = random.randint(2,p)
		while (legendre(d,p) != -1):
			d = random.randint(2,p)

		s = 0				
		while (p-1)%(pow(2,s)) == 0:
			s += 1		
		s = s - 1
				
		t = (p-1)/(pow(2,s))
		
		A = pow(a,t,p) 
		D = pow(d,t,p) 
		m = 0
		
		#Find the number m such that AD^m = 1 mod p
		for i in range(s):
			temp = (A)*(pow(D,m))
			if (pow(temp,pow(2,(s-1-i)))) % p == p-1:
				m = m + pow(2,i)
			
		root = ((pow(a,(t+1)/2,p))*(pow(D,(m/2),p))) % p
		return int(root)		
	
	







