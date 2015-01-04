"""
Codes an algorithm for constructing a rational point on an elliptic
curve defined over the rationals.
	
Author: Erik Holum
E-mail: EHolum@gmail.com
"""
import Numbers
from sage.schemes.elliptic_curves.heegner import *


def point(E,D,prec=50):
	"""
	Returns a point of infinite order on the curve E(Q) by using the
	heegner point construction for modular elliptic curves
	
	Currently it produces an estimate very close to the actual point,
	but i'm hoping to implement a reverse continued fractions method
	to get an estimate on the rational value of this point
	"""
	points = get_heegner_points(E,D,prec)
	
	P = 0
	
	for x in points:
		P = P + x[1]
		
	return P


def get_heegner_points(E,D,prec=50):
	"""
	Returns pairs of heegner points of discriminant D on X_0(N) as well as
	their mappings to heegner points on the curve E(C) under the canonical
	modular mapping 
	
	\phi: X_0(N) \to E(C)
	
	INPUT:
	E - An Elliptic Curve defined over the rationals
	D - The fundamental discriminant of an imaginary quadratic field
	prec - bit precision, am hoping to implement a continued fractions
	
	OUTPUT:
	a list of pairs (P,X) P \in X_0(N) X \on E(C), the list length is 
	dependent on the class number of the field Q(\sqrt(D))
	"""
	if not satisfies_heegner_hypothesis(E,D):
		raise ValueError('E and D must satisfy the heegner hypothesis')
	
	rho = E.weierstrass_p(prec)
	rho1 = rho.derivative()
	
	H = Heegner(E,prec)
	P = H.heegner_points(D)
	
	L = []
	for x in P:
		temp = x.map_to_curve(E).map_to_complex_numbers(prec)
		L.append((x,(rho(temp),rho1(temp))))
		
	return L
	
	
		
def next_disc(E,bound=-3):
	""" 
	Returns a suitable negative discriminant that is a quadratic
	residue mod 4N and coprime to 2N. This discriminant represents
	the fundamental discriminant of an imaginary quadratic field
	K = Q(\sqrt{d})
	
	INPUT:
	bound - returns first satisfying D > bound
	
	RETURN:
	D < bound that is a discriminant for a potential heegner point
	"""		
	bound -= 1
	while not satisfies_heegner_hypothesis(E,bound):
		bound -= 1
	return bound



class HeegnerPoint:
	
	def __init__(self, E,prec=50):
		""" 
		Constructs an object based on the Elliptic Curve E
		
		INPUT:
		E - An Elliptic curve
		prec - Precision for the weierstrass laurent series and complex
			maps
		"""
		self.E = E
		self.N = E.conductor()
		
		
		
	def __repr__(self):
		return 'Tools for constructing rational points on %s using heegner points'% self.E
		
		
	def curve(self):
		""" Returns the curve associated to this object.
		"""
		return self.E
		
		
	def heegner_points(self,D):
		"""
		Returns a list of all heegner points on X_0(N) of the chosen
		discriminant using the process described in Zagiers paper [Z]
		We consider points tau \in H of the form:
		
		..::Math::..
		tau = b + i\sqrt{d} / 2a
		
		with
		
		a,b in ZZ, a > 0, N|a, and b^2 \cong -d (mod 4a)
		
		If legendre(-d,N) == -1, then there is no such point, so we
		return [0]. Otherwise there are h such points, where h is 
		equal to the class number of K
		
		INPUT:
		N - Level of X_0(N)
		
		OUTPUT:
		List of Heegner points together with their maps onto the curve
		E
		""" 
		
		#If legendre(-d,N) == -1, there are no such points
		if (Numbers.legendre(D,self.N) == -1):
			return [0]
		else:
			if (D is None):
				raise(ValueError("Need to define discriminant first"))
			temp = heegner_points(self.N,D,1)
			H = []
			for p in temp:
				H.append(p)
			return H
		
		
		

		
	
