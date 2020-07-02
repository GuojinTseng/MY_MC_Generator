# -*- coding: utf-8 -*-
#==========================================================#
# Process: e+e- -> Z/gamma -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.16
# Version: 1.0
#==========================================================#

#import guojin's module
import calc_dsigma, param_card, run_card

#import some useful packages
import math
import random
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import numpy as np

class Calc_Xsection(object):

	def __init__(self):
		self.dsig = calc_dsigma.Calc_Dsigma()

	def generate_histo(array, name):
		plt.hist(array, bins=1000)
		plt.title('Sampling Importance Resampling')
		plt.show()

	def ZGJ(self, regn, init, itmx):
		print "The e+e- center-of-momentum frame energy: ", run_card.ECM, "GeV"

		print '\n'
		iteration_interval = []
		total_interval = regn
		iteration_interval.append(total_interval)
		abswgt = []
		area = 0

		for it in range(itmx):
			xo = regn[0][0] # xo = 0.0
			xn = regn[0][1] # xn = 2.0
			tl_interval = []
			i = 0

			pp = [xo] # 用来记录调整后的区间
			subarea = 0
			num_points = 1000
			pp_diff = []

			for interval in total_interval: # 对每一个interval撒点，排序，记录权重，
				num_bin = 3
				i += 1

				num_points = (int(1000*(interval[1] - interval[0])/(regn[0][1] - regn[0][0]))+1)
				points = interval[0] + (interval[1] - interval[0]) * np.random.rand(num_points) #随机点的个数一定要足够大（大于num_bin），否则可能导致点不够分，最后使得区间数目达不到要求
				points.sort(axis = 0)

				wgts = abs(self.dsig.dsigma(points))
				sum_wgts = np.sum(wgts)
				cut = sum_wgts/num_bin
				# print '当前区间的cut：', cut, i

				k = 0
				sum_wgt = 0
				xo = interval[0]

				for point in points:
					k += 1
					sum_wgt += abs(self.dsig.dsigma(point)) # 先考虑fx大于零的情况
					if sum_wgt > cut:
						pp.append(point)
						pp_diff.append(pp[-1]-pp[-2])
						tl_interval.append([xo, point])
						xo = point
						sum_wgt = 0
						k = 0

				pp.append(interval[1])
				pp_diff.append(pp[-1]-pp[-2])
				tl_interval.append([xo, interval[1]])

			total_interval = tl_interval
			iteration_interval.append(total_interval)

			# print 'iteration: ', it + 1
			# generate_histo(pp, 'test') # 用来查看区间划分情况

		fxpp, sigma = step_function(pp, tl_interval, self.dsig.dsigma)
		tl_sigma = np.sqrt(np.sum(np.array(sigma)**2))
		expectation = fxpp * np.array(pp_diff)

		print 'Total Xection(MC value) =', np.sum(expectation) * param_card.pb_convert, '+-', tl_sigma * param_card.pb_convert, 'pb'
		print 'Total Xection(analytic value) =', 1060.93990308, 'pb'

		return iteration_interval[-1], np.max(fxpp)


def step_function(points, inter, f):
	pp = []
	sigma = []
	for point in points:
		for x in inter:
			if (x[0] <= point < x[1]):
				xx = x[0] + (x[1] - x[0]) * np.random.rand(1000)
				xx= np.hstack((xx, x[1]))
				xx= np.hstack((xx, x[0]))
				pp.append(np.sum(f(xx))/1002)
				sigma.append(np.sqrt(abs(np.sum(f(xx)*(x[1] - x[0]) * f(xx)*(x[1] - x[0]))/1002 - (np.sum(f(xx)*(x[1] - x[0]))/1002)**2)/1002))
				break
	return np.array(pp), sigma
