#==========================================================#
# Process: e+e- -> Z/gamma -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.16
# Version: 1.0
#==========================================================#

#import guojin's module
import calc_dsigma, param_card, run_card

#import some useful packages
import math, random

class Calc_Xsection(object):

	def __init__(self):
		self.dsig = calc_dsigma.Calc_Dsigma()

	def xsec(self):
		print "The e+e- center-of-momentum frame energy: ", run_card.ECM, "GeV"

		print '\n'

		random.seed(run_card.seed)

		#=====================Integrate the differential cross section, and find the weight_max in the meantime========================#
		sum_weight = 0
		sum_weight_sq = 0
		weight_max = 0
		costh_max = 0
		delta = 2

		for i in range(0, run_card.N):
			costh_i = -1 + random.random() * delta
			weight_i = self.dsig.dsigma(costh_i) * delta
			sum_weight = sum_weight + weight_i
			sum_weight_sq = sum_weight_sq + weight_i**2
			if weight_i > weight_max:
				weight_max = weight_i
				costh_max = costh_i

		sigma = sum_weight / run_card.N

		variance = sum_weight_sq / run_card.N - (sum_weight / run_card.N)**2
		error = math.sqrt(variance / run_card.N)
		#=====================Integrate the differential cross section, and find the weight_max in the meantime========================#

		# Show results and compare MC value with analytical value
		print 'Find the maximum value of dsigma = ', weight_max, 'at costh = ', costh_max

		print '\n'

		print 'Total Xection(MC value) =', sigma * param_card.pb_convert, '+-', error * param_card.pb_convert, 'pb'

		print 'Total Xection(Analytical value) =', self.dsig.dsigma(0) * param_card.pb_convert * 8/3., 'pb'

		print '\n'

		return weight_max
