#==========================================================#
# Process: e+e- -> Z/gamma -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.16
# Version: 1.0
#==========================================================#

#import guojin's module
import param_card, run_card

import numpy as np

class Calc_Dsigma(object):

	def dsigma(self, costh):
		# CL and CR
		CV = -0.5 + 2 * param_card.sw2
		CA = -0.5
		# constants and functions that appear in the differential cross section
		kappa = np.sqrt(2) * param_card.G_Fermi * param_card.MZ**2 / (4 * np.pi * param_card.alpha)
		chi1 = kappa * run_card.hats * ( run_card.hats - param_card.MZ**2 ) / ((run_card.hats-param_card.MZ**2)**2 + param_card.GAMMAZ**2*param_card.MZ**2)
		chi2 = kappa**2 * run_card.hats**2 / ((run_card.hats-param_card.MZ**2)**2 + param_card.GAMMAZ**2*param_card.MZ**2)
		A0 = 1 + 2 * CV**2 * chi1 + (CA**2 + CV**2)**2 * chi2
		A1 = 4 * CA**2 * chi1 + 8 * CA**2 * CV**2 * chi2

		PREFAC = (2 * np.pi) * param_card.alpha**2 / (4 * run_card.hats) # 2 * pi comes from d-phi integral
		return  PREFAC * ( A0 * ( 1 + costh[0]**2 ) + A1 * costh[0] )

		# return (4*param_card.alpha**2*np.pi**2*(4*param_card.MZ**4*(-1 + param_card.sw2)**2*param_card.sw2**2 + 2*run_card.ECM**2*param_card.MZ**2*param_card.sw2*(-1 - 7*param_card.sw2 +
		# 	8*param_card.sw2**2) + run_card.ECM**4*(1 + 24*param_card.sw2**2) - 2*run_card.ECM**2*(2*param_card.MZ**2*(-1 + param_card.sw2)*param_card.sw2 +
		# 	run_card.ECM**2*(1 + 8*param_card.sw2**2))*costh + (4*param_card.MZ**4*(-1 + param_card.sw2)**2*param_card.sw2**2 + 2*run_card.ECM**2*param_card.MZ**2*param_card.sw2*(-1 - 7*param_card.sw2 + 8*param_card.sw2**2) + run_card.ECM**4*(1 + 24*param_card.sw2**2))*costh**2))/((-4*run_card.ECM**2 + param_card.MZ**2)**2*(-1 + param_card.sw2)**2*param_card.sw2**2)

		# return costh+1
