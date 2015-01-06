from fractions import gcd
import ModEC
import Numbers

def get(M, P):
	"""
	Find all curves E: Y^2 = X^3 + aX + b such that a,b < M+1 and that
	gcd(Np, Nq) = 1 over all primes p,q in P, with p,q not dividing 2\Delta
	"""
	l = []
	
	for a in range(M+1):
		for b in range(M+1):
			disc = (4*(a**3)) + (2*(b**2))
			badred = Numbers.factor(disc*2)

			E = ModEC.EC(a,b)
			
			sizes = 0
			for p in P:
				if not(p in badred):
					Np = E.find_all_points(p)
					sizes = gcd(sizes,len(Np))
			
			if sizes == 1:
				l.append((a,b))
	return l


def doubs(M):
	l = []
	for a in range(M+1):
		for b in range(M+1):
			l.append((a,b))
			
	return l