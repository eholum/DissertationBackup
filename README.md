DissertationBackup
==================
"Heegner Point Constructions for Modular Elliptic Curves"

Backup of subset of my masters thesis.
I seem to have lost some of it over the years, which is kind of a drag.

Note to self: Hard drives are physical objects that can be lost in the hands of careless
folks such as myself.

Missing some Magma as well as some clarity in the python. There's some fun stuff in here though,
like an implementation of factoring integers using elliptic curves, as well as a lot of random
operations on elliptic curves over various rings.

Abstract:

Let $E/\Q$ be an elliptic curve in minimal Weierstrass form
$$E: y^2 + a_1xy + a_3y = x^3 + a_2x^2 + a_4x + a_6.$$
Finding rational points on $E(\Q)$ has been a long standing problem, and currently 
there are several methods. The brute force method of looping over the value $a/d^2$ 
for $d = 1,2,3,...$ and $a = \pm 1,\pm 2, \pm 3,...$ is disappointingly slow, and 
other algorithms have been developed to speed up the process. In this paper, we examine 
one such algorithm which uses the modular parameterization of $E/\Q$ and so called 
Heegner points to construct a rational point on a curve of analytic (and hence algebraic) 
rank 1. 

The aim of this paper is precisely to describe each step in the process as well as 
enough of the under-lying mathematics for someone with limited experience working with 
modular curves to be able to understand the majority of the computation. We first 
present background on elliptic functions, complex lattices, and isomorphism classes of 
elliptic curves, then continue with a brief overview of modular forms and the modularity 
theorem for elliptic curves, and then define Heegner points in both an abstract and concrete 
manner. We then discuss the process of constructing a rational point on $E/\Q$, as well as 
provide several examples ranging in height complexity from less than 1 to over 1000. 