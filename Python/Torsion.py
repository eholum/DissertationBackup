import ModEC
from math import floor, sqrt, cos
from fractions import Fraction

"""
A module used to calculate the torsion subgroup of an elliptic curve of the 
form:
E: Y^2 = X^3 + aX + b
"""
def proportion(M):
	"""
	returns the number of trivial torsion curves and the number of
	total curves for (M,b),(a,M) ,|a|,|b| < M
	"""
	count = 0
	total = 0
	for a in range(M+1):
		total += 4
		if tors(a,M):
			count += 1
		if tors(-1*a,M):
			count += 1
		if tors(a,-1*M):
			count += 1
		if tors(-1*a,-1*M):
			count += 1
	for b in range(M+1):		
		total += 4
		if tors(M,b):
			count += 1
		if tors(M,-1*b):
			count += 1
		if tors(-1*M,b):
			count += 1
		if tors(-1*M,-1*b):
			count += 1
			
	return (total-count,total)
	
	
	
def find_tors(a,b):
	"""
	Finds the elements of torsion subgroup of the curve E_a,b, along
	with their orders
	"""
	ps = points(a,b)
	
	point = ['o']
	
	for P in ps:
		r = order(a,b,(P[0],P[1],1))
		
		if (r > 0):
			print P
			Q = P
			for x in xrange(r):
				print add(a,b,Q+(1,),P + (1,))
			point.append([P,r])
			print
	
	return point
	
def find_triv_tors(M):
	"""
	finds all (a,b) such that |a|,|b| <= M and E_ab has trivial torsion
	"""
	l = []
	
	for a in range(M):
		for b in range(M):
			if tors(a,b):
				l.append(ModEC.EC(a,b))
			if tors(-1*a,b):
				l.append(ModEC.EC(-1*a,b))
			if tors(a,-1*b):
				l.append(ModEC.EC(a,-1*b))
			if tors(-1*a,-1*b):
				l.append(ModEC.EC(-1*a,-1*b))
	
	return l

def tors(a,b):
	"""
	Returns True of the curve Y^2 = X^3 + aX + b has trivial torsion,
	False otherwise.
	"""
	
	if len(cubic(0,a,b)) > 0:
		return False
		
	ps = points(a,b)
	
	E = ModEC.EC(a,b)
	
	for P in ps:
		if order(a,b,(P[0],P[1],1)) > 0:
			return False
			
	return True

def order(a,b,P):
	"""
	Determines the order of P in E(Q). Returns -1 if the order
	is infinte (ie > 12).
	"""
	if P[2] == 0:
		return 1
	if P[1] == 0:
		return 2

	count = 1
	Q = P
	while (Q != (0,1,0)):	
		Q = add(a,b,P,Q)
	
		for x in Q:
			if x.denominator > 1:
				return -1
		count += 1
		if count >= 13:
			return -1
	
	return count




def add(A,B,P,Q):
	"""
	Adds the points P and Q on E over Q.
	"""
	x1,y1,z1 = (Fraction(P[0]),Fraction(P[1]),Fraction(P[2]))
	x2,y2,z2 = (Fraction(Q[0]),Fraction(Q[1]),Fraction(Q[2]))

	if (y1 == -1*y2):
		return (0,1,0)
		
	if not(x1 == x2):
		l = (y1 - y2)/(x1 - x2)
		m = (x1 * y2 - x2 * y1)/(x1 - x2)
		x3 = ((x1*(x2**2)) + ((x1**2)*x2) + (A*(x1 + x2)) + (2*B) - (2*y1*y2))
		x3 = x3/((x1 - x2)**2)
		y3 = (-1*l*x3) - m
		return (x3,y3,1)
		
	else:
		if (y1*-1 == y2):
			return (0,1,0)
		
		
		l = ((3*(x1**2)) + A)/(2*y1)
		m =((-1*(x1**3)) + A*(x1) + 2*B)/(2*y1)
		x3 = ((x1**4) - (2*A*(x1**2)) - (8*B*x1) + A**2)
		x3 = x3/(4*(y1**2))
		y3 = (-1*l*x3) - m
		return (x3,y3,1)
	

def on_e(a,b,x,y):
	"""
	returns True if the point x,y is on the curve, False otherwise
	"""
	lhs = y ** 2
	rhs = (x**3) + (a*x) + (b)
	
	return lhs == rhs
	
def points(a,b):
	"""
	Returns a list of all the points (x,y) on the curve E. Note this 
	is exactly determining the solutions to the equation:
	0 = X^3 + aX + (b - y^2)
	These points will be exactly the possible torsion points on the curve E
	"""
	ys = [0]+ posy(a,b)
	l = []
	
	for y in ys:
		temp = cubic(0,a,b-y ** 2)
		for x in temp:
			if not((x,y) in l):
				l.append((x,y))
				
	return l



def posy(a,b):
	"""
	returns all possible values for y as dictated by the Nagell Luts
	Theorem, excluding 0
	"""
	return square_factor((4 * pow(a,3)) + (27 * pow(b,2)))


def square_factor(n):
	n = abs(n)
	"""
	Returns a list of the square prime factors of n
	"""
	if is_prime(n):
		return [1]
	else:
		factors = [1]
		for x in range(2,n/2 + 1):
			if (n % (x ** 2)) == 0:
				factors.append(x)
				while(n % (x) == 0):
					n = n/x
					if(n == 1):
						return factors					
		return factors	


def is_prime(n):
	"""
	Trivially decides if an integer n is prime.
	"""
	for x in range(2, int(floor(sqrt(n))) + 1):
		
		if n % x == 0:
			return False
			break
	else:
		return True
		

def check(P):
	"""
	Returns the real integer elements of a tuple.
	"""
	l = []
	for x in P:
		if not(type(x) == type(0j)):
			if math.floor(x) == x:
				l.append(int(x))
				
	return l
		
"""
x^3 + ax^2 + bx + c = 0  (or ax^3 + bx^2 + cx + d = 0)
With substitution x = y-t and t = a/3, the cubic equation reduces to    
    y^3 + py + q = 0,
where p = b-3t^2 and q = c-bt+2t^3.  Then, one real root y1 = u+v can
be determined by solving 
    w^2 + qw - (p/3)^3 = 0
where w = u^3, v^3.  From Vieta's theorem,
    y1 + y2 + y3 = 0
    y1 y2 + y1 y3 + y2 y3 = p
    y1 y2 y3 = -q,
the other two (real or complex) roots can be obtained by solving
    y^2 + (y1)y + (p+y1^2) = 0
"""
def cubic(a, b, c, d=None):
    
    if d:			# (ax^3 + bx^2 + cx + d = 0)
	a, b, c = b / float(a), c / float(a), d / float(a)
    t = a / 3.0
    p, q = b - 3 * t**2, c - b * t + 2 * t**3
    u, v = quadratic(q, -(p/3.0)**3)
    if type(u) == type(0j):	# complex cubic root
	r, w = polar(u.real, u.imag)
	y1 = 2 * cbrt(r) * cos(w / 3.0)
    else:			# real root
        y1 = cbrt(u) + cbrt(v)
    y2, y3 = quadratic(y1, p + y1**2)
    return check((y1 - t, y2 - t, y3 - t))
    
    
def quadratic(a, b, c=None):
    import math, cmath
    if c:		# (ax^2 + bx + c = 0)
	a, b = b / float(a), c / float(a)
    t = a / 2.0
    r = t**2 - b
    if r >= 0:		# real roots
	y1 = math.sqrt(r)
    else:		# complex roots
	y1 = cmath.sqrt(r)
    y2 = -y1
    return y1 - t, y2 - t
    
def cbrt(x):
    from math import pow
    if x >= 0: 
	return pow(x, 1.0/3.0)
    else:
	return -pow(abs(x), 1.0/3.0)
	
"""
Convert from rectangular (x,y) to polar (r,w)
    r = sqrt(x^2 + y^2)
    w = arctan(y/x) = [-\pi,\pi] = [-180,180]
"""
def polar(x, y, deg=0):		# radian if deg=0; degree if deg=1
    from math import hypot, atan2, pi
    if deg:
	return hypot(x, y), 180.0 * atan2(y, x) / pi
    else:
	return hypot(x, y), atan2(y, x)

