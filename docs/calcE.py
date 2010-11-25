from math import *
def calcE(beta,U,E0):
	E1=0
	E2=-U/2.
	E3=-U/2.
	E4=0
	r1=exp(-beta*(E1-E0))
	r2=exp(-beta*(E2-E0))
	r3=exp(-beta*(E3-E0))
	r4=exp(-beta*(E4-E0))
	Z=r1+r2+r3+r4
	avE_=(E1-E0)*r1+(E2-E0)*r2+(E3-E0)*r3+(E4-E0)*r4
	avE_/=Z
	return avE_+E0


def easyE(beta,U):
	return -U/2.*exp(beta*U/2.)/(1+exp(beta*U/2.))
	
