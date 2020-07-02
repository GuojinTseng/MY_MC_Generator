# -*- coding: utf-8 -*-

import math
import random
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import numpy as np
import itertools

def generate_histo(array, name):
	plt.hist(array, bins=1000)
	plt.title('Sampling Importance Resampling')
	plt.show()

def ZGJ(regn, fx, init, itmx, step_function):
	
	ndim = len(regn[0])/2

	total_interval = [0 for i in range(ndim)]

	for di in range(ndim):
		pp = [[] for i in range(ndim)]
		pp_diff = [[] for i in range(ndim)]
		

		total_interval[di] = [[regn[0][di], regn[0][ndim + di]]]
		area = 0

		for it in range(itmx):
			xo = regn[0][di] # xo = 0.0
			xn = regn[0][ndim + di] # xn = 2.0
			tl_interval = []
			i = 0

			subarea = 0
			num_points = 100

			for interval in total_interval[di]: # 对每一个interval撒点，排序，记录权重，
				num_bin = 3
				i += 1

				num_points = (int(100*(interval[1] - interval[0])/(regn[0][ndim + di] - regn[0][di]))+1)
				points = interval[0] + (interval[1] - interval[0]) * np.random.rand(num_points) #随机点的个数一定要足够大（大于num_bin），否则可能导致点不够分，最后使得区间数目达不到要求
				points.sort(axis = 0)

				wgts = abs(fx(points))
				sum_wgts = np.sum(wgts)
				cut = sum_wgts/num_bin
				# print '当前区间的cut：', cut, i

				k = 0
				sum_wgt = 0
				xo = interval[0]

				for point in points:
					k += 1
					sum_wgt += abs(fx([point])) # 先考虑fx大于零的情况
					if sum_wgt > cut:
						tl_interval.append([xo, point])
						xo = point
						sum_wgt = 0
						k = 0

				tl_interval.append([xo, interval[1]])

			total_interval[di] = tl_interval

			print 'iteration: ', it + 1, len(total_interval)

	fxpp, subgrid_area = step_function(total_interval, fx, ndim)
	expectation = np.sum(np.array(fxpp) * np.array(subgrid_area))

	return expectation

def fx(x):
	# return x**2+2
	# return np.exp((x-0.5)**2)
	return np.cos(x[0])+1
	# return 2


def step_function(inter, fx, ndim):
	do = inter[0]
	for di in range(ndim - 1):
		do = itertools.product(do, inter[di + 1])

	grids = do
	xx = [0 for i in range(ndim)]

	pp = []
	subgrid_area = []

	for grid in grids:
		area_grid = 1
		for di in range(ndim):
			xx[di] = grid[di][0] + (grid[di][1] - grid[di][0]) * np.random.rand(10)
			area_grid *= (grid[di][1] - grid[di][0])

		subgrid_area.append(area_grid)

		xxx = [0 for n in range(10)]

		for n in range(10):
			xxx[n] = [xx[di][n] for di in range(ndim)]

		pp.append(np.mean([fx(xxx[n]) for n in range(10)])) # 每一格点里的函数平均值已经得到，同时被收集到pp里面

	return pp, subgrid_area


################### To Begin ###################


if __name__ == "__main__":
	init = 1
	itmx = 2
	EE = ZGJ([[0.0, 0.0, 2.0, 1.0]], fx, init, itmx, step_function) # regen should be float
	print 'iteration: ', itmx, 'times'
	print 'integral: ', EE



################### To-Do List ###################
# 1. damp function
# 2. 处理一下函数值可能为负数的情况
# 3. the more iteration times, the lower speed
# 4. sigma estimate
# 5. alpha running
# 6. multi-dimension