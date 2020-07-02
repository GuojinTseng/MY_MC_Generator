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
		ndim = len(regn)

		total_interval = [0 for i in range(ndim)]

		for di in range(ndim):
			pp = [[] for i in range(ndim)]
			pp_diff = [[] for i in range(ndim)]


			total_interval[di] = regn[di]
			area = 0

			for it in range(itmx):
				xo = regn[di][0][0] # xo = 0.0
				xn = regn[di][0][1] # xn = 2.0
				tl_interval = []
				i = 0

				subarea = 0
				num_points = 1000

				for interval in total_interval[di]: # 对每一个interval撒点，排序，记录权重，
					num_bin = 3
					i += 1

					num_points = (int(1000*(interval[1] - interval[0])/(regn[di][0][1] - regn[di][0][0]))+1)
					points = interval[0] + (interval[1] - interval[0]) * np.random.rand(num_points) #随机点的个数一定要足够大（大于num_bin），否则可能导致点不够分，最后使得区间数目达不到要求
					points.sort(axis = 0)

					yy = [[] for i in range(ndim)]

					for ddi in range(ndim): # 不同变量的区间不一样
						yy[ddi] = regn[ddi][0][0] + (regn[ddi][0][1] - regn[ddi][0][0]) * np.random.rand(10)

					del(yy[di])

					yyy = []


					for n in range(10): # 在剩下ndim-1个维度上生成随机数
						yyy.append([yy[ddi][n] for ddi in range(ndim-1)])

					# print yyy

					xxx = []
					wgts = [0 for ii in range(len(points))]
					for ii in range(len(points)):
						for n in range(10):
							# print 'yyy' , yyy[n], di
							yyy[n].insert(di,points[ii])
							# print yyy[n]
							xxx.append(yyy[n])
							# print points[ii],'ys'
							
							wgts[ii] += abs(self.dsig.dsigma(xxx[0]))

					sum_wgts = np.sum(wgts)
					cut = sum_wgts/num_bin
					# print '当前区间的cut：', cut, i

					k = 0
					sum_wgt = 0
					xo = interval[0]

					for point in points:
						sum_wgt += wgts[k] # 先考虑fx大于零的情况
						k += 1
						if sum_wgt > cut:
							tl_interval.append([xo, point])
							xo = point
							sum_wgt = 0
							k = 0

					tl_interval.append([xo, interval[1]])

				total_interval[di] = tl_interval


		fxpp, subgrid_area = step_function(total_interval, self.dsig.dsigma, ndim)
		expectation = np.sum(np.array(fxpp) * np.array(subgrid_area))

		print 'Total Xection(MC value) =', np.sum(expectation) * param_card.pb_convert, '+-', 'pb'
		print 'Total Xection(analytic value) =', 1060.93990308, 'pb'

		return expectation

def step_function(inter, fx, ndim):

	do = inter[0]
	for di in range(ndim - 1):
		do = itertools.product(do, inter[di + 1]) # 获得区间的笛卡尔积

	grids = do

	xx = [0 for i in range(ndim)]

	pp = []
	subgrid_area = []

	for grid in grids:

		if ndim > 1:
			inverse_grid = []
			for i in range(ndim-1):
				inverse_grid.append(grid[-1])
				grid = grid[0]
				if i == ndim-2:
					inverse_grid.append(grid)

			inverse_grid.reverse()
			grid  = inverse_grid
		
		else:
			gridgrid = [0]
			gridgrid[0] = grid
			grid = gridgrid

		area_grid = 1
		for di in range(ndim):
			xx[di] = grid[di][0] + (grid[di][1] - grid[di][0]) * np.random.rand(15)
			area_grid *= (grid[di][1] - grid[di][0])

		subgrid_area.append(area_grid)

		xxx = [0 for n in range(15)]

		for n in range(15):
			xxx[n] = [xx[di][n] for di in range(ndim)]

		pp.append(np.mean([fx(xxx[n]) for n in range(15)])) # 每一格点里的函数平均值已经得到，同时被收集到pp里面

	return pp, subgrid_area
