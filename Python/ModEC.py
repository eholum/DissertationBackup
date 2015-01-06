from random import randint
from random import choice
from Numbers import mod_inverse
from Numbers import legendre
from Numbers import modular_sqrt
from Numbers import factor
from fractions import gcd

"""
Module for operations on modular elliptic curves.
"""
class EC:
	
	
	def np(self,P):
		"""
		Determines the size of E(F_p) for each p in P.
		"""
		l = []
		for p in P:
			N = len(self.find_all_points(p))
			l.append(N)
			
		return l
	
	#Creates a curve with coefficients a,b,c where the curve is:
	#y^2 = x^3 + a x + b defined over n where n is in Z (usually prime)
	def __init__(self, a,b):
		self.C = [a,b]
		
		
	#Returns string representative of elliptic Curve
	def tostring(self):
		string = 'Y^2 = X^3' 
		
		if(self.C[0] < 0):
			string = string + ' - ' + str(-1*self.C[0]) + 'X'
		elif(self.C[0] > 0):
			string = string + ' + ' + str(self.C[0]) + 'X'
			
		if(self.C[1] < 0):
			string = string + ' - ' + str(-1*self.C[1])
		elif(self.C[1] > 0):
			string = string + ' + ' + str(self.C[1])
				
		return string
		
	
	#Returns true if the point P is on the curve mod n, false otherwise
	def on_curve(self,P,n):
		t = self.evaluate(P[0],n)
		lhs = modular_sqrt(t,n)
		
		if (lhs == P[1]):
			return True
		elif (mod_inverse(lhs,n) == P[1]):
			return True
		
		return False
		
	#Evaluates f(X) mod n
	def evaluate(self,X,n):
		a,b = self.C
		return (pow(X,3,n) + a * X + b) % n
	
	
	
	#Finds a point on the curve over F_n by choosing randomly
	def find_point(self,n):
		
		l = list(range(0,n))		
		x = choice(l)
		l.remove(x)
		
		t = self.evaluate(x,n)
		
		#while t is not a quadratic residue mod n, keep trying to find one
		while(legendre(t,n) == -1):
			x = choice(l)
			l.remove(x)
			t = self.evaluate(x,n)
		
		return (x,modular_sqrt(t,n))
		
		
	#Returns all integer points on the elliptic curve mod n
	def find_all_points(self,n):
		
		points = [(0,1,0)]
		
		for x in range(0,n):
			t = self.evaluate(x,n)
			if (legendre(t,n) != -1):
				points.append((x,modular_sqrt(t,n),1))
				
		for P in points:
			Q = P
			while Q != (0,1,0):
				Q = self.add(Q,P,n)
				if not(Q in points):
					
					points.append(Q)
				
		
		return points
		
	
	def add(self,P1,P2,n):
		"""
		Adds the affine coordinates P1 and P2 on the affine curve E_a,b(Z_n)
		"""
	
		m = 0
		a,b = self.C
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
			
			inv = mod_inverse(2 * Y1,n)
			m = (3 * (X1 ** 2) + a) * inv
		
		else:
			inv = mod_inverse((X2 - X1) % n,n)
			m = (Y2 - Y1) * inv
		
		X3 = m ** 2 - X1 - X2
	
		return (X3 % n, (m * (X1 - X3) - Y1) % n,1)
		
	
		#Returns the negative of a point
		def neg(self,P,n):
			return (P[0],-1*P[1],P[2])
		
	
	#Doubles the point
	def double(self,P,n):
		return self.add(P,P,n)
		
		
	#Subtracts two points
	def subtract(self,P1,P2,n):
		return self.add(P1,neg(P2,n),n)
		
	
	def order(self,P,n):
		"""
		Determines the order of P in E(F_n)
		"""
		
		count = 1
		Q = P
		while (Q != (0,1,0)):
			Q = self.add(Q,P,n)
			count += 1
			
		return count
		
	
	def group_table(self,n):
		"""
		prints group table of the = curve: E(F_n) suitable for 
		copying and pasting into latex as a tabular environment
		"""
		
		l = self.find_all_points(n)
		table = []

		for P in l:
			temp = []
			
			for Y in l:
				print Y
				temp.append(self.add(P,Y,n))
			table.append(temp)
		
		
		string = ''
		
		temp = l
		
		for P in temp:
			if P[2] == 0:
				string = string + "& $\id$"
			else:
				string = string + '& $(' + str(P[0]) + ', ' + str(P[1]) + ')$ '
		print string
		
		print "\hline"
		
		j = 0
		for line in table:
			if j == 0:
				string = "& $\id$"
			else:
				string = string + '& $(' + str(l[j][0]) + ', ' + str(l[j][1]) + ')$'
			temp = line
			for P in temp:
				if P[2] == 0:
					string = string + '&$\id$'
				else:
					string = string + '&$(' + str(P[0]) + ', ' + str(P[1]) + ')$'
			print string + "\\"
			string = ""
			print string
			j += 1


		
		
		
	def modular_sqrt(a, p):
		"""
		Finds the squareroot of a mod p. eg:
		Finds a solution to the equation x^2 \cong a mod p
		Assuming the riemann hypothesis this will run in polynomial time
		"""
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
			
		
	def legendre(a, p):
		"""
		Computes the Legendre Symbol of a mod p
		"""
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
			


		
