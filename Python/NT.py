import math
import cmath

"""
Used solely for the examination. Unhelpful otherwise.
"""
def croot(n):
	"""
	returns f-p approx for cube root of n
	"""
	return math.pow(n,1.0/3.0)
	
	
def trips(M):
	"""
	returns a list of all triples (a,b,c) with M-1 <= a,b,c <= M-1
	"""
	l = []
	
	for a in range(-M + 1,M):
		for b in range(-M+1,M):
			for c in range(-M+1,M):
				l.append((a,b,c))
				
	l.remove((0,0,0))
				
	return l


def moduli(P,d):
	"""
	For B = P[0] + P[1]t + P[2]t, computes and returns:
	[ B , B' , B'' ]
	"""
	
	w = (1.0 + cmath.sqrt(-3))/2.0
	t = croot(d)
	
	B = 1.0 * P[0] + t * P[1] + P[2] * pow(t,2)
	B1 = P[0] + (P[1] * w * t) + (P[2] * pow(w,2) * pow(t,2))
	B2 = B1.conjugate()
	
	return [B,B1,B2]


def find_all(M,d):
	"""
	Finds all triples B = (x,y,z) such that |B| <= \mu * M^-2
	"""
	
	L = trips(M)
	sat = []
	t = math.pow(d,1.0/3.0)
	check = (1 + t + math.pow(t,2.0)) * math.pow(M,-2)
	for P in L:
		B = P[0] + (P[1] * t) + (P[2] * math.pow(t,2.0))
		
		if B.__abs__() <= check:
			sat.append(P)
						
	return sat


def find_mod(M,d):
	"""
	For each B in L checks the if |B'| < mu M
	"""
	
	L = trips(M)
	sat = []
	t = math.pow(d,1.0/3.0)
	check = check = (1 + t + math.pow(t,2.0)) * M
	
	for P in L:
		B,B1, B2 = moduli(P,d)
		
		if (B1.__abs__()) <= check:
			sat.append(P)
			
	return sat
	
	
def find_pairs(M,d):
	"""
	For each satisfying point, find a pair of elements (a,b,c),(x,y,z) in
	L_M such that B = a_1 - a_2
	"""
	
	L = find_mod(M,d)
	points = []
	
	for B in L:
		a1 = []
		a2 = []
		for P in B:
			x1, x2 = mod_sum(M,P)
			a1.append(x1)
			a2.append(x2)
			
		points.append([text(M,a1),text(M,a2)])
		
	return points
	

def intersection(L1,L2):
	"""
	returns the intersection of 2 lists
	"""
	L = []
	for P in L1:
		if P in L2:
			L.append(P)
			
	return L
	
	
def mod_sum(M,x):
	"""
	finds 2 integers a,b \in [1,M] such that
	a - b = x for x \in [1 - M, M - 1]
	"""
	if (x == 0):
		return [M,M]
	elif (x < 0):
		return [x % M, M]
	else:
		return [M, M - x]
		
		
def text(M,P):
	"""
	writes P = (a,b,c) as:
	P = ('M' - (M - a), 'M' - (M - b)...)
	"""
	if (P[0] == M):
		a = 'M - ' + str(M-P[0])
	else:
		a = 0
	if (P[1] == M):
		b = 'M - ' + str(M-P[0])
	else:
		b = 0
	if (P[2] == M):
		c = 'M - ' + str(M-P[0])
	else:
		c = 0
		
	return (a,b,c)