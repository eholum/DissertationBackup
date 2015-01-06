from sage.schemes.elliptic_curves.heegner import *
import sage.rings.number_field.number_field as number_field
from sage.schemes.elliptic_curves.constructor import *
from sage.symbolic.function import *
import sage.rings.all as rings
from sage.rings.all import (ZZ, GF, QQ, CDF, CC,
                            Integers, RealField, ComplexField, QuadraticField,
                            gcd, lcm, is_fundamental_discriminant)

"""
Module for representing and providing basic functions on heegner points. Including their
construction from modular curves.
"""
import pickle
from math import sqrt,log

# For saving
#	E = EllipticCurve([1,-1,0,-751055859,-7922219731979])
#	E = EllipticCurve([0,1,1,-4912150272,-132513750628709])
#	E = EllipticCurve([0,0,1,-5115523309,-140826120488927])
	
#E = EllipticCurve([0,0,1,-1,0])
#E = EllipticCurve([1, 0, 1, -362, 2615])
#E = EllipitcCurve([1,-1,0,-36502495762,-2684284892271276])
#E = EllipticCurve([0,1,1,-693591631266,-222333206671717949])

def find_large_height_SteinWatkins(n,c=0):
	"""
	Returns an elliptic curve with a heegner point height > n by iterating
	over the Large Stein-Watkins database starting with curves of conductor
	10^5(c) to 10^5(c+1)
	"""
	from sage.databases.stein_watkins import SteinWatkinsAllData
	D = SteinWatkinsAllData(c)
	C = D.next()
	height = 0
	while height < n:
		C = D.next()
		for curve in C.curves:
			E = EllipticCurve(curve[0])
			an_rank = E.analytic_rank()
			if (an_rank == 1):
				print curve[0]
				HP = HeegnerPointsOnCurve(E)
				height = HP.heegnerPointHeight()
				del HP,E
				print height
	return (E,height)
			
			
def find_large_height_Cremona(n,c=11):
	"""
	Attempts to find a curve with generator height > n in the 
	cremona database, starting with curves of conductor = 11
	"""
	from sage.databases.cremona import LargeCremonaDatabase
	C = LargeCremonaDatabase()
	E = None
	height = 0
	while height < n:
		for conductor in range(c,180000):
			curves = C.allcurves(conductor)
			for curve in curves:
				E = EllipticCurve(curves.get(curve)[0])
				an_rank = E.analytic_rank()
				if (an_rank == 1):
					print E
					HP = HeegnerPointsOnCurve(E)
					height = HP.heegnerPointHeight()
					print height
					del HP,E
	return (E,height)
		
	
		
	 
def sat_heegner_hypothesis(N,D,c=1):
	"""
	Checks to ensure that N and D (and c) satisfy the Heegner Hypothesis
	"""
	if D>=0:
		return False
	if D.gcd(c) != 1:
		return False	
	
	L = N.factor()
	for x in L:
		p = x[0]
		l = x[1]
		if (D%p == 0):
			if l > 1:
				return False
		elif D.kronecker(p) != 1:
			return False
			
	if not number_field.is_fundamental_discriminant(D):
		return False
	
	return True
	
	
def tuple_sat_heegner_hypothesis(f,N,D):
	"""
	Checks if a tuple f satisfies the heegner hypothesis with N
	"""
	A, B, C = f
        
	if B**2 - 4*A*C != D:
        	return False
        if A%N != 0:
        	return False	
        if gcd(A/N,gcd(B,C*N)) != 1:
        	return False
        	
	return True


class HeegnerPointsOnCurve:
	"""
	A class to store the set of all heegner points on X_0(N) with respect
	to an elliptic curve E.
	"""
	def forSage(self,D,l=None,prec=100):
		print str(self.discriminants())
		print "$ \\beta = " + str(self.betas(D))
		forms = self.forms(D)
		z = 0
		for f in forms:
			bqf = BinaryQF(f).reduced_form()
			p = HeegnerPoint(self.N,D,1,f)
			c = p.mapToComplexPlane(self.E,prec)
			z += c
			print "$" + str(f) + "$ & $" + str(p) + "$ & $" + str(c) +"$\\" + "\\"
		print "z = " + str(z)
		C = self.E.period_lattice()
		P = C.elliptic_exponential(z)
		print P
		if l is None:
			l = sqrt(self.computeGeneratorHeight(D)[1])
		
		omega_re, omega_im = C.basis()
		Etors = self.E.torsion_subgroup()
		m = l.gcd(Etors.exponent())
		
		if self.E.discriminant() > 0:
			barz = []
			for u in range(1,l*m):
				temp = (m*z.real() + u*omega_re) / (m*l)
				barz.append((u,l,temp))
				print temp
				temp += (omega_im)/2
				barz.append((u,l,temp))
				print temp
			return (z,P,m,l,barz)
		else:
			barz = []
			o = z.imag()/omega_re
			for u in range(1,l*m):
				temp = ((m*z.real() + u*omega_re) / (m*l)) + ((o*omega_re)/2)
				barz.append((u,l,temp))
				print temp
			return (z,P,m,l,barz)
	
	
	def range_over_u(self,z,l=None):
		"""
		range over u
		"""
		omega_re, omega_im = C.basis()
		Etors = self.E.torsion_subgroup()
		m = l.gcd(Etors.exponent())
		
		if self.E.discriminant() > 0:
			barz = []
			for u in range(1,l*m):
				temp = (m*z.real() + u*omega_re) / (m*l)
				barz.append((u,l,temp))
				print temp
				temp += (omega_im)/2
				barz.append((u,l,temp))
				print temp
			return barz
		else:
			barz = []
			o = z.imag()/omega_re
			for u in range(1,l*m):
				temp = ((m*z.real() + u*omega_re) / (m*l)) + ((o*omega_re)/2)
				barz.append((u,l,temp))
				print temp
			return barz
	

	def __init__(self,E):
		self.N = ZZ(E.conductor())
		self.E = E
	
	
	def __repr__(self):
		return "The set of all Heegner Points in X_0(%s) covering %s"%(self.N,self.E)
		
		
	def betas(self,D):
		"""
		Returns the allowable set of betas that are square roots of d mod 
		4n mod 2n.
		"""
		N = self.N
		R = Integers(4*N)
		m = 2*N
		return tuple(sorted( set([a%m for a in R(D).sqrt(all=True)]) ))
	
	
	def discriminants(self,n=10):
		"""
		Returns the first n integers <= -4 that satisfy the heegner
		hypothesis
		"""
		l = []
		test = ZZ(-4)
		while len(l) < n:
			if sat_heegner_hypothesis(self.N,test):
				l.append(test)
			test -= 1
		
		return l
	
	
	def getPoints(self,D,beta=None,c=1):
		"""
		Returns a list of heegner points such that (A,B,C) has 
		B congruent to Beta. If beta is None, we pick it
		"""
		
		beta = self.betas(D)[0]
		#If Beta is not a root of 2N mod 4N, raise error
		if (beta not in self.betas(D)):
			raise ValueError("%s^2 must be congruent to %s mod %s"%(Beta,2*self.N,4*self.N))
		
		forms = self.forms(D,beta,c)
		points = []
		
		for bqf in forms:
			points.append(HeegnerPoint(self.N,D,c,bqf))
			
		return sorted(points)
		
		
	def forms(self,D,beta=None,c=1):
		"""
		Returns the list of BQFs (not reduced) of binary quadratic forms
		associated to the set of Heegner points of discrimant D and 
		beta
		"""
		N = ZZ(self.N)
		D = ZZ(D)
		c = ZZ(c)
		
		if not sat_heegner_hypothesis(N,D):
			raise ValueError("N and D must satisfy Heegner Hypothesis")
		if c.gcd(N) != 1:
			raise ValueError("Conductor and N must be coprime")
			
		if beta is None:
			beta = self.betas(D)[0]
		
		# size of set of binary quadratic forms equivalence class representatives
		class_number = QuadraticField(D,'x').class_number()
		disc = ZZ(D*(c**2))
		U = []
		R = []
		b = ZZ(beta%(2*N))
		a = 1
		while len(U) < class_number:
			while c.gcd(a) != 1:
				a += 1
			
			C = ZZ((b**2 - disc)/(4*N)) 
			for x in Integers(a):
				if N*x**2 + b*x + C == 0:
					x = ZZ(x)
					Ap =a*N
					Bp = b+2*N*x
					Cp = ZZ(((b + 2*N*x)**2 - disc)/(4*a*N))
					f = (Ap,Bp,Cp)
					bqf = BinaryQF(f).reduced_form()
					if bqf not in U:
						U.append(bqf)
						R.append(f)
						if len(U) == class_number: 
							break
			a += 1
		return sorted(R)
		
		
	def heegnerPointHeight(self,D=None,prec=2):
		"""
		Computes the height of the Heegner point to a certain level 
		of precision using the Gross-Zagier Theorem.
		
		we have
		H(P) = L'(E,1) #E_tors^2 / 2 omega T 
		where T is the Tamagawa number of E
		
		We will compute the height to an accuracy of 10^-prec
		"""
		if D is None:
			D = self.discriminants(1)[0]
		RR = rings.RealField(53)
		IR = rings.RealIntervalField(53)
		E = self.E
		D = ZZ(D)

		if (E.lseries().L1_vanishes() and E.root_number() == 1):
			return IR(0,0)
		
		
		if not sat_heegner_hypothesis(self.N,D):
			raise ValueError("Discriminant %s and %s must satisfy the Heegner hypothesis"%(D,self.N))
		
		ED = E.quadratic_twist(D)
		#Area of the period lattice
		omega = RR(E.period_lattice().complex_area())
		#Number of units in Q(\sqrt(D)) and Q(sqrt(gcd(D,N)))
		WD = QuadraticField(D,'x').unit_group().order()
		WND = len((D.gcd(self.N)).factor())		
		
		#WHY THESE VALUES?
		precE = prec*sqrt(self.N) + 20
		precED = prec*sqrt(ED.conductor()) + 20
		
		MAX_ERROR = RR(0.00000001)
		temp = E.lseries().deriv_at1(precE)
		if (temp == 0):
			return 0
		LE,LE_Error = temp
		temp = ED.lseries().at1(precED)
		if (temp == 0):
			return 0
		LED,LED_Error = temp
		
		M_E = max(MAX_ERROR,LE_Error,LED_Error)
		
		height = RR((sqrt(abs(D))/(4*omega))*(LE*LED)*((WD/2)**2)*(pow(2,WND)))
		
		return IR(height-M_E,height+M_E)
		
		
	def computeGeneratorHeight(self,D=None,prec=2):
		""" 
		Computes the height of a generator using a combination of the
		Gross-Zagier theorem and the BSD conjecture
		"""
		if D is None:
			D = self.discriminants(1)[0]
		RR = rings.RealField(53)
		IR = rings.RealIntervalField(53)
		E = self.E
		D = ZZ(D)
		if not sat_heegner_hypothesis(self.N,D):
			raise ValueError("Discriminant %s and %s must satisfy the Heegner hypothesis"%(D,self.N))
		if (E.lseries().L1_vanishes() and E.root_number() == 1):
			return IR(0,0)
		
			
		MAX_ERROR = RR(0.00000001)
		ED = E.quadratic_twist(D)
		
		#Area of the period lattice
		period_lattice = E.period_lattice()
		omega = RR(period_lattice.complex_area())
		omega_re = period_lattice.basis()[0]
		
		#Number of units in Q(\sqrt(D)) and len(factor(sqrt(gcd(D,N))))
		WD = QuadraticField(D,'x').unit_group().order()
		WND = len((D.gcd(self.N)).factor())		
		
		#WHY THESE VALUES?
		precED = prec*sqrt(ED.conductor()) + 20
		precE = precE = prec*sqrt(self.N) + 20
		temp = E.lseries().deriv_at1(precE)
		if (temp == 0):
			return 0
		LE,LE_Error = temp
		temp = ED.lseries().at1(precED)
		if (temp == 0):
			return 0
		LED,LED_Error = temp
		
		#Tamagawa Product
		TP = 1
		for (p,l) in self.N.factor():
			TP *= E.tamagawa_number(p)
		TP *= E.sha().an_numerical()
		
		p1 = (omega_re/(4*omega))
		p2 = TP
		p3 = sqrt(abs(D))/((E.torsion_subgroup().order())**2)
		p4 =((WD/2)**2)*(pow(2,WND))
		
		l2 = RR(p1*p2*LED*p3*p4)
		height = RR((sqrt(abs(D))/(4*omega))*(LE*LED)*((WD/2)**2)*(pow(2,WND)))
		M_E = max(MAX_ERROR,LE_Error,LED_Error)
		
		return [IR((height/l2)-M_E,(height/l2)+M_E), sqrt(l2)]
		
	
	def getHeegnerPoint(self,D,beta=None,c=1,prec=10,verbose=False):
		"""
		Attempts to construct a rational point on an elliptic curve
		using the Heegner point construction. Returns the Heegner
		point as an element of C/Lambda
		"""
		if beta is None:
			beta = self.betas(D)[0]
		
		if verbose:
			print "Getting list of Heegner Points..."
		points = sorted(self.getPoints(D,beta,c))
		if verbose:
			print "... done (%s points)"%len(points)
		
		if verbose:
			print "Computing the height of the Heegner point..."
		height = self.heegnerPointHeight(D)
		
		precision = 2*height.upper()
		if verbose:
			print "...Heegner point has height %s+"%height.lower()
		
		M = ((log(10)**(-1*precision)))
		
		
		if verbose:
			print "Computing L-series to sufficient precision (%s terms)..."%nterms
		enumerated_an = list(enumerate(self.E.anlist(M)))[1:]
		if verbose:
			print "...done"

		if verbose:
			print "Mapping Heegner points to complex plane..."
		heegnerPointOnLattice = 0
		
		for p in points:	
			if verbose:
				print "Mapping %s to C/L..."%p
			z = CC(p.tau)
			q = (2*CC.gen()*CC.pi()*z).exp()
			
			complexpoint = 0
			for n,an in reversed(enumerated_an):
				complexpoint += (an/n)
				complexpoint *= q
			
			heegnerPointOnLattice += complexpoint
		
		if verbose:
			print "...done"
		
		#heegnerpoint = self.E.elliptic_exponential(heegnerPointOnLattice,False)
		
		return (heegnerPointOnLattice)#,heegnerpoint)
		

class HeegnerPoint:
	
	def __init__(self,N,D,c=1,f=None):
		""" 	
		Constructs a heegner point on H\Gamma_0(N) with discriminant D.
		and conductor c (default 1)
		"""
		self._N = ZZ(N) 
		self._D = ZZ(D)
		self._C = ZZ(c)
		
		if self._C.gcd(self._N) != 1:
			raise ValueError("Conductor and N must be coprime")
	
		if not sat_heegner_hypothesis(N,D):
			raise ValueError("N and D must satisfy Heegner Hypothesis")
		
		if f is None:
			A = N
			B = ZZ(Integers(4*N)(D*c**2).sqrt(extend=False) % (2*N))
			C = ZZ((B*B - D*c**2)/(4*A))
			f = (A,B,C)
			
		else:
			if not len(f) == 3:
				raise IOError("f must be a quadratic form")
			f = tuple(ZZ(x) for x in f)
			A,B,C = f[0],f[1],f[2]
			disc = (B**2) - (4*A*C)
			if disc != D*c**2:
				raise ValueError("f must have disc %s"%d)

		X = "sqrt_minus_%s"%-D
		self._QF = number_field.QuadraticField(D,X)
		self._f = f        
		d = self._QF.gen()*c
		self.tau = (-B + d)/(2*A)
		

	def __repr__(self):
		"""
		Return a String representation of this heegner point
		"""
		num = self.tau*self._f[0]*2
		den = 2*self._f[0]
		#Simplified representation
		string = (str(num/den)).replace("sqrt_minus_%s"%-self._D,"sqrt(%s)"%self._D)
		#Unreduced
		string = ("\\frac{" + str(num) + "}{" + str(den) + "}").replace("sqrt_minus_%s"%-self._D,"\\sqrt{%s}"%self._D)
		return str(string)
		
	def tau(self):
		"""
		Return the value of this element as viewed in the upper half
		plane modulo the action of \Gamma_0(N)
		"""
		return self.tau
		
	def bqf(self):
		"""
		returns the binary quadratic form attached to this curve
		"""
		return self._f
		
	def mapToComplexPlane(self,E,prec=100,anseq=None):
		"""
		maps this point to the elliptic curve E as viewed as a lattice
		C/L using the standard modular parametrisation phi: X_0(N) to
		E
		"""
		if (E.conductor() != self._N):
			raise ValueError("E must have conductor %s"%N)
		
		#M = (log(10,2)**(prec))/(2*3.15159265358979323846264*Function_imag_part(self.tau))
		CC = ComplexField(prec)	
		z = CC(self.tau)
		q = (2*CC.gen()*CC.pi()*z).exp()
		nterms = (-(prec+10)/q.abs().log2()).ceil()
		enumerated_an = list(enumerate(E.anlist(nterms)))[1:]
		point = 0
		for n, an in reversed(enumerated_an):
			point += an/n
			point *= q
		return point
		
		
			
	
	def disc(self):
		"""
		Return the discriminant of this heegner point
		"""
		return self._D
		
		
	def level(self):
		"""
		Return the level of this heegner point
		"""
		return self._N
		
		

		
		
