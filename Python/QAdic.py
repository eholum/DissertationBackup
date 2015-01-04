import fractions
import Numbers
import math

#Given rational number n and prime p, computes the highest power v such that
#p^n m = n such that p does not divide m
def vp(n,p):
	a = Numbers.factor(n.numerator)
	b = Numbers.factor(n.denominator)
	v = -1*(a.count(p) - b.count(p))
	return v

#Computes the lowest common multiple of 2 integers
def lcm(num1, num2):
    result = num1*num2/fractions.gcd(num1,num2)
    return result
    
  
#Determines if a rational number is a square in Q_p, note that this function 
#only works for numbers a such that |a|_p = 1. Returns 1 if yes, -1 if no,
#and 0 if it can't decide
def isSquare(n,p):
	v = vp(n,p)
	l = math.pow(p,v)
	if(l != 1):
		return 0
	else:
		if(p==2):
			if (n % 8 == 1):
				return 1
			else:
				return -1
		else:
			return Numbers.legendre(n,p)
			

