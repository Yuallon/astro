#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 12:48:52 2019

@author: huangyulong
"""

#实现画Jarrett选择条件图

import matplotlib.pyplot as plt


#plt.vlines(2.2, 0.6, 1.7, colors = "b", linestyles = "dashed")
#plt.vlines(4.2, 0.6, 1.7, colors = "b", linestyles = "dashed")
#plt.hlines(1.7, 2.2, 4.2, colors = "b", linestyles = "dashed")

def frange(start, stop=None, step=None):

    if stop == None:
        stop = start + 0.0
        start = 0.0

    if step == None:
        step = 1.0

    while True:
        if step > 0 and start >= stop:
            break
        elif step < 0 and start <= stop:
            break
        yield ("%g" % start) # return float number
        start = start + step

print ("Printing float range")
floatList = frange(2.2, 4.21, 0.01)
xf = []
for num in floatList:
    xf.append(float(num))
#    print (num)
    
yf= [i*0.1 + 0.38 for i in xf]
plt.vlines(2.2, 0.6, 1.7, colors = "b", linestyles = "dashed")
plt.vlines(4.2, 0.8, 1.7, colors = "b", linestyles = "dashed")
plt.hlines(1.7, 2.2, 4.2, colors = "b", linestyles = "dashed")
plt.plot(xf,yf,'b--')

