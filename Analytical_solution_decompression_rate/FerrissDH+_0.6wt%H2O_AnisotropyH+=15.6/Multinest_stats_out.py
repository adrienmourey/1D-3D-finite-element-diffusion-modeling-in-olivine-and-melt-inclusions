#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:45:51 2017

@author: ejfm2
"""

# Script for extracting pymultinest data into a panda dataframe

from __future__ import absolute_import, unicode_literals, print_function

import numpy
from numpy import exp, log
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest
import pandas as pd


# Use Profile name as second argument then make extra column with profile name 
#- then append this file to a master file which can then be used for manipulation of whole dataset

os.environ['D'] = '2'

if len(sys.argv) < 2:
	sys.stderr.write("""SYNOPSIS: %s <output-root> 

	output-root: 	Where the output of a MultiNest run has been written to. 
	            	Example: chains/1-
%s""" % (sys.argv[0], __doc__))
	sys.exit(1)

prefix = sys.argv[1]
profile = sys.argv[2]
print('model "%s"' % prefix)
parameters = json.load(open(prefix + 'params.json'))
n_params = len(parameters)

#==============================================================================
# if len(sys.argv) == 3:
#         Nlim = int(sys.argv[2])
# else:
#==============================================================================
Nlim = n_params

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s = a.get_stats()

json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)

values = a.get_equal_weighted_posterior()

m = s['marginals'][0]

center = m['median']
low1, high1 = m['1sigma']
low2, high2 = m['2sigma']

xerr=numpy.transpose([[center - low1, high1 - center]])

xerrl = center - low1
xerrh = high1 - center

#'col1': [1, 2], 'col2': [3, 4]

#print(center)
#print(low1)
#print(high1)
#print(high2)
#print(low2)
#print(xerrl)
#print(xerrh)

df = pd.DataFrame(index=[0], data= {'Median': center, '1 sigma low': low1, '1 sigma high': high1, '2 sigma low': low2, '2 sigma high': high2, '1sigma error low': xerrl, '1sigma error high': xerrh})

df['Profile'] = profile

df.to_csv('{0}_pyMN_stats.csv'.format(profile), sep=',')

df2 = pd.DataFrame(index=[0])

# Create empty pandas dataframe here

for i in range(0,Nlim):
#==============================================================================
#     plt.subplot(Nlim, Nlim, i + 1)
#     plt.xlabel(parameters[i])
#     m = s['marginals'][i]
#     plt.xlim(m['5sigma'])
#==============================================================================
    m = s['marginals'][i]
    df2[parameters[i]] = m['median']
	
 
     # Append to empty parameters dataframe
     # Use parameters[i] as name and values[:,i] as array

# Create profile 
df2['Profile'] = profile

df2.to_csv('{0}_pyMN_medians.csv'.format(profile), sep=',')


# Plot up cumulative frequency curves for t
#t1 = df['t'].values
#T1 = df['T'].values
#P1 = df['P'].values
#Fe1 = df['fe_3'].values

#t2 = numpy.sort(t1)
#T2 = numpy.sort(T1)
#P2 = numpy.sort(P1)
#Fe2 = numpy.sort(Fe1)

#n = len(t1)

#f2 = numpy.array(range(n))/float(n)

#fig, axes = plt.subplots(nrows=2,ncols=2)
#axes[0,0].plot(t2, f2)
#axes[0,1].plot(T2, f2)
#axes[1,0].plot(P2, f2)
#axes[1,1].plot(Fe2, f2)
#axes[0,0].set_ylabel('Probability')
#axes[0,0].set_xlabel('t') 
#axes[0,1].set_ylabel('Probability')
#axes[0,1].set_xlabel('T') 
#axes[1,0].set_ylabel('Probability')
#axes[1,0].set_xlabel('P') 
#axes[1,1].set_ylabel('Probability')
#axes[1,1].set_xlabel('fe_3') 

#plt.savefig('{0}_CFD.pdf'.format(profile))


 
#==============================================================================
# 	oldax = plt.gca()
# 	x,w,patches = oldax.hist(values[:,i], bins=nbins, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
# 	oldax.set_ylim(0, x.max())
# 	
# 	newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
# 	p.plot_marginal(i, ls='-', color='blue', linewidth=3)
# 	newax.set_ylim(0, 1)
# 	
# 	ylim = newax.get_ylim()
# 	y = ylim[0] + 0.05*(ylim[1] - ylim[0])
# 	center = m['median']
# 	low1, high1 = m['1sigma']
# 	#print(center, low1, high1)
# 	newax.errorbar(x=center, y=y,
# 		xerr=numpy.transpose([[center - low1, high1 - center]]), 
# 		color='blue', linewidth=2, marker='s')
# 	oldax.set_yticks([])
# 	#newax.set_yticks([])
# 	newax.set_ylabel("Probability")
# 	ylim = oldax.get_ylim()
# 	newax.set_xlim(m['5sigma'])
# 	oldax.set_xlim(m['5sigma'])
#==============================================================================
	#plt.close()
 
 
# Save pandas dataframe here - differntiate between constant and Al conditions etc? 

# Could also try plotting up cumaltive frequency 
# Order values in array and then normalise to 1? Something like that look it up. 

# append in shell script
