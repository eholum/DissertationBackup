"""
Module for operations on projected elliptic curves.
"""

def add_proj(E,P1,P2,n):
	"""
	adds the projective point P1 = (X,Y,Z) to P2 = (X,Y,Z)
	on the curve E: [A,B] -- E: y^2 = x^3 + A x + B in Z_n
	"""
	
	# To make the code slightly more legible
	X1, Y1, Z1 = P1
	X2, Y2, Z2 = P2
	a = E[0]
	b = E[1]
	
	# If either point is the identity
	if (Z1 == 0):
		return P2
		
	if (Z2 == 0):
		return P1
		
	U1 = (X2 * (Z1 ** 2)) % n
	U2 = (X1 * (Z2 ** 2)) % n
	
	S1 = (Y2 * (Z1 ** 3)) % n
	S2 = (Y1 * (Z2 ** 3)) % n
	
	W = (U1 - U2) % n
	R = (S1 - S2) % n
	
	# If the X coordinates match
	if (W == 0):
		if (R == 0):
			return double(E,P1,n)
		return (0,1,0)
	
	T = (U1 + U2) % n
	M = S1 + S2
	
	X3 = ((R ** 2) - (T * (W ** 2))) % n
	Y3 = (((1/2) * ((T * (W ** 2)) - (2 * X3)) * R) - (M * (W ** 3))) % n
	Z3 = (Z1 * Z2 * W) % n
	
	return (X3,Y3,Z3)
	
	
def double_proj(E,P1,n):
	"""
	Adds the projective point P1 to itself on E defined over Z_n
	"""
	X = P1[0]
	Y = P1[1]
	Z = P1[2]
	a = E[0]
	b = E[1]
	
	if (Y == 0 or Z == 0):
		return (0,1,0)
	
	M = ((3 * (X ** 2)) + (a * (Z ** 4))) % n
	S = (4 * X * (Y ** 2)) % n
	L = (M ** 2 - 2 ** S) % n 
	Q = ((M * (S - X ** 2)) - (8 * (Y ** 4))) % n
	R = 2 * Y * Z
	
	return (L,Q,R)
	
	
def subtract_proj(E,P1,P2,n):
	"""
	Subtracts the projective P2 from P1 on E defined over Z_n
	"""
	return add(E,P1,neg(P2),n)
	
	
def neg_proj(E,P,n):
	""" 
	Inverts the projective point P on E over Z_n, ie, returns P^-1
	"""
	return (P[0]%n,(-1*P[1])%n,P[2]%n)
