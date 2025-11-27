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

df = pd.DataFrame(index=range(len(values[:,0])))

# Create empty pandas dataframe here

for i in range(0,Nlim):
#==============================================================================
#     plt.subplot(Nlim, Nlim, i + 1)
#     plt.xlabel(parameters[i])
#     m = s['marginals'][i]
#     plt.xlim(m['5sigma'])
#==============================================================================
    
    df[parameters[i]] = values[:, i]
	
 
     # Append to empty parameters dataframe
     # Use parameters[i] as name and values[:,i] as array

# Create profile 
df['Profile'] = profile

df.to_csv('{0}_pyMN_dataframe.csv'.format(profile), sep=',')


# Plot up cumulative frequency curves for t
t1 = df['dpdt'].values
T1 = df['T'].values
P1 = df['Kd'].values
Fe1 = df['MI_radius'].values


t2 = numpy.sort(t1)
T2 = numpy.sort(T1)
P2 = numpy.sort(P1)
Fe2 = numpy.sort(Fe1)

n = len(t1)

f2 = numpy.array(range(n))/float(n)

fig, axes = plt.subplots(nrows=2,ncols=2)
axes[0,0].plot(t2, f2)
axes[0,1].plot(T2, f2)
axes[1,0].plot(P2, f2)
axes[1,1].plot(Fe2, f2)
axes[0,0].set_ylabel('Probability')
axes[0,0].set_xlabel('dpdt (MPa/s)') 
axes[0,1].set_ylabel('Probability')
axes[0,1].set_xlabel('T (C)') 
axes[1,0].set_ylabel('Probability')
axes[1,0].set_xlabel('Kd') 
axes[1,1].set_ylabel('Probability')
axes[1,1].set_xlabel('MI_radius (um)') 

plt.savefig('{0}_CFD.pdf'.format(profile))


 
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
