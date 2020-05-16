#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 19:47:39 2019

@author: huangyulong
"""

#这个代码是用来作maximum likelihood estimation
#星等值和观测误差值
#mag=[19.832449 ,20.406380 ,20.455864, 20.286505, 20.504390 ,20.072280, 21.025521,20.698955]
#magerr=[0.13522495, 0.19616871 ,0.31400712, 0.19169995, 0.30148667 ,0.50314934 ,0.40691672,0.30421210]


import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import ascii
#from scipy import stats
#from scipy.stats import chisquare
import math

####Read WISE datas

#data = ascii.read('/Users/huangyulong/Desktop/Summary/Liupick_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/Sartori2015_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/Sartori_pick1_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/Sartori_pick2_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/Galaxypick_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/red_agn_neo.tbl.txt')
#data = ascii.read('/Users/huangyulong/Desktop/Summary/Reines2013_neo.tbl.txt')
data = ascii.read('/Users/huangyulong/Desktop/1111.tbl.txt')


#Good datas
data1 = data[data['qi_fact'] >= 1.0]
data2 = data1[data1['saa_sep'] > 5.0]
data3 = data2[data2['moon_masked'] == '00']
n1 = len(data3)

# Extract index, w1，w2, mjd, error
x1 = data3['cntr_01']
ra = data3['ra']
dec = data3['dec']
w_1 = data3['w1mpro']
w_2 = data3['w2mpro']
mjd = data3['mjd']
erro_1 = data3['w1sigmpro']
erro_2 = data3['w2sigmpro']

#Curve the masked into zero value 
w1 = w_1.filled(fill_value=0)
w2 = w_2.filled(fill_value=0)
e1 = erro_1.filled(fill_value=0)
e2 = erro_2.filled(fill_value=0)

df = pd.DataFrame([x1,w1,w2,mjd,e1,e2,ra,dec])
data4 = df.T

#排除掉0值
i = 0
data5 = []
while i < len(data4):
   if min(data4.iloc[i][:6])>0:
       data5.append(list(data4.iloc[i]))
   i += 1

data6 = pd.DataFrame(data5)

s = len(set(data6[0]))

n2 = 0

#给出每个源的坐标，中值，中值误差，MJD1
RA = []
DEC = []
Z = []
#magchange = []

#中值，中值误差
MAG1 = []
ERR1 = []
MAG2 = []
ERR2 = []
MJD1 = []

w1_list = []
w2_list = []
mjd_list = []

for i in set(data6[0]):
    data7 = data6[data6[0] == i]
    data7.columns = ['index','w1','w2','mjd','e1','e2','ra','dec']
    
    data8 = data7.sort_values('mjd')
    
    
    
    index = list(data7['index'])
    ra = list(data7['ra'])
    dec = list(data7['dec'])
    RA.append(ra[0])
    DEC.append(dec[0])
    
    w_11 = list(data8['w1'])
    w_21 = list(data8['w2'])
    mjd1 = list(data8['mjd'])
    e_11 = list(data8['e1']) 
    e_21 = list(data8['e2'])
    

    MJD = []
    w11 = []
    w21 = []
    e11 = []
    e21 = []
        
    i = 0
    j = 0
    while i <= len(mjd1) - 1: 
    
        MJD.append([])
        MJD[j].append(mjd1[i])
        w11.append([])
        w11[j].append(w_11[i])
        w21.append([])
        w21[j].append(w_21[i])
        e11.append([])
        e11[j].append(e_11[i])
        e21.append([])
        e21[j].append(e_21[i])
    
        if i == len(mjd1) - 1:break #此处是为了保证最后一个数据也被成功归为一组
    
    #例如：a = [1.0,1.5,2.0,2.1,4.3]，运行后保证4.3作为单独一个分组：[[1.0, 1.5, 2.0, 2.1], [4.3]]
    
        while mjd1[i+1] - mjd1[i] < 1:
            MJD[j].append(mjd1[i+1])
            w11[j].append(w_11[i+1])
            w21[j].append(w_21[i+1])
            e11[j].append(e_11[i+1])
            e21[j].append(e_21[i+1])
            i += 1 #此处i值是为了将满足条件的数据点放到一个列表里面
            if i > len(mjd1) - 2 :break #保证i+1不会超过列表a的长度
        j += 1
        i += 1 #此处i值是为了实现第一个while循环，且实现数据分段

#     
    # median and standard deviation    
    def pp(x):
        median = []
        std = []
#        for i in x:
        for i in range(len(x)):   
            if len(x[i]) > 1: #如果一个周期观测数据点为1，舍弃不用
                median.append(np.median(x[i]))
                std.append(np.nanstd(x[i]))
        return(median,std)
        
#    #对于同一时间段内的多次观测，要求每次观测数据点大于1
#    MJD1 = [i for i in MJD if len(i) > 1]
#    w111 = [i for i in w11 if len(i) > 1]
#    w211 = [i for i in w21 if len(i) > 1]
    #对于同一时间段内的多次观测，星等取中值
    m_mjd, s_mjd = pp(MJD)
    m_w11, std_w11 = pp(w11)
    m_w21, std_w21 = pp(w21)
    
    MAG1.append(m_w11)
    MAG2.append(m_w21)
    MJD1.append(m_mjd)
#    #误差取平均值,降低随机误差(random error)的方法：多次测量 or 取均值
#    m_e11 = [np.median(i) - np.mean(i) for i in w11]
#    m_e21 = [np.median(i) - np.mean(i) for i in w21]
    
    #root mean square
    def RMS(x):
        
        s = 0.0
        #Calculate mean
        m = np.mean(x)
    
        #Calculate sum
        for i in x:
            s += np.square((i - m)) 
        
        #Calculate rms
        rms = math.sqrt((s/(len(x)-1)))
    
        return rms
    #取平方和的均值的平方根
    def ee(x):
        
        s = 0.0
        for i in x:
            s += np.square(i)
        ms = math.sqrt(s/(len(x)))
        return ms
    #以中值作为最似然估计值，求中值的标准偏差
#    def RMS1(x):
#        
#        s = 0.0
#        #Calculate mean
#        m = np.median(x)
#    
#        #Calculate sum
#        for i in x:
#            s += np.square((i - m)) 
#        
#        #Calculate rms
#        rms1 = math.sqrt((s/(len(x)-1)))
#    
#        return rms1    
    #用数据的标准偏差作为数据的误差，如果某个观测时间只有一个数据点，则舍掉不用(#则采用原有的系统误差)
    m_e11 = []
    for i in range(len(w11)):
        if len(w11[i]) > 1:
#            a = RMS(w11[i])*math.sqrt((len(w11[i])-1)/len(w11[i]))
            #SE(median) = 1.2533 * SE(mean) for normal distribution
            a = 1.2533*RMS(w11[i])/math.sqrt(len(w11[i]))
            m_e11.append(a)
#        else:
#            m_e11.append(e11[i])
    ERR1.append(m_e11)
        
    m_e21 = []
    for i in range(len(w21)):
        if len(w21[i]) > 1:
#            a = RMS(w21[i])*math.sqrt((len(w21[i])-1)/len(w21[i]))
            #SE(median) = 1.2533 * SE(mean) for normal distribution
            a = 1.2533*RMS(w21[i])/math.sqrt(len(w21[i]))
            m_e21.append(a)
#        else:
#            m_e21.append(e21[i])
    ERR2.append(m_e21)
    
    w1_list.append(w11)
    w2_list.append(w21)
    mjd_list.append(MJD)
    
#    plt.errorbar(m_mjd,m_w11,yerr=m_e11,label='W1')
#    plt.errorbar(m_mjd,m_w21,yerr=m_e21,label='W2')
#    plt.xlabel('MJD')
#    plt.ylabel('W1&W2(mag)')
#    plt.title(str(RA[0]) + '+' + str(DEC[0]))
#    plt.legend()
##    plt.savefig(str(index[0])+'.pdf')
#    plt.show()
    
    

#MLE
#求取多个源最大似然值<m>,sigma_m,给出99%置信区间
#2.30, 4.61, and 9.21 for probability levels α = 68%, 90%, and 99%, respectively.    

#WISE w1
w1_m_mean_list = []#每个源的w1的MLE值，所有源list
w1_sig_mean_list = []#每个源的w1光变Msigma的MLE值，所有源list

w1_mag_list1 = []#1000步里面，满足99%置信度下，w1可能的数值，对于所有源的list
w1_mag_list2 = []
w1_sigma_list = []

w1_sigma_list_99 = []

RA1 = []
DEC1 = []
MAG_w1 = []
ERR_w1 = []
MJD_w1 = []

w1_list1 = []
mjd_list1 = []
   
n = 0
while n < len(MAG1):
    if len(MAG1[n])>1:
        mag = MAG1[n]
        magerr = ERR1[n]
        
        RA1.append(RA[n])
        DEC1.append(DEC[n])
        MAG_w1.append(MAG1[n])
        ERR_w1.append(ERR1[n])
        MJD_w1.append(MJD1[n])
        w1_list1.append(w1_list[n])
        mjd_list1.append(mjd_list[n])
        
        nn = len(mag)
        L = 1000
        #循环结构求去sima_m,mean_mag
        sigma_m = []
        mag_mean = []
        S = []
        dx = 0.8/L
        i = 1
        
        a = list(range(L))
        #sigma_m = [0.8*(i+1)/L for i in a]
        while i < L:
        #    a = list(range(L))
        #    sigma_m = [0.8*(i+1)/L for i in a]
            sig_m = i*dx
            sigma_m.append(sig_m)
            
            t0 = [sig_m**2 + k**2 for k in magerr]
            
            summag = 0
            j = 0
            while j < len(mag):
                a = mag[j]/t0[j]
                summag = summag + a
                j += 1
            
            sumerr = 0
            j = 0
            while j < len(magerr):
                b = 1/t0[j]
                sumerr = sumerr + b
                j += 1
            
            mag_m = summag/sumerr
            mag_mean.append(mag_m)
            
            #在已知sigma_m and mag_m下给出S函数值
            s1 = nn*math.log(2*math.pi)
            
            mr = [math.log(sig_m**2 + ri**2) for ri in magerr]
            s2 = sum(mr)
            
            s3 = 0
            j = 0
            while j < len(mag):
                s3 = s3 + (mag[j]-mag_m)**2/t0[j]
                j += 1
            
            s = s1 + s2 + s3
            S.append(s)
            i += 1
        

        
    #    if i > 1:
    #        if S[-1] > S[-2]:break
    
    #求出最大似然估计sigma_m and <m>值    
        i = 0
        while i < len(S):
            if S[i+1] > S[i]:break
            i += 1    
        print(i)        
        sig_mean = sigma_m[i]
        m_mean = mag_mean[i]
        s_min = S[i]
        
        w1_m_mean_list.append(m_mean)
        w1_sig_mean_list.append(sig_mean)
#        n += 1
#    else:
#        n += 1
        
        #calculate the m_mean-sigma region of certain confidence level
        #delts=2.3,4.61,9.21 corresponding to 68%, 90%, 99% cinfidence level
        sig_range = []
        mean_range1 = []
        mean_range2 = []
        mean_err_99 = []
        sig_err_99 = []
        
        i = 1
        L = 1000
        dx = 0.8/L
        
        while i < L:
            sig_m = i*dx
            sigma_m.append(sig_m)
            
            t0 = [sig_m**2 + k**2 for k in magerr]
            
            a = [1/k for k in t0]
            A = sum(a)
            
            b = 0
            j = 0
            while j < len(mag):
                b = b + mag[j]/t0[j]
                j += 1
                
            B = -2.0*b
                
            c1 = nn*math.log(2*math.pi) - s_min - 9.21
            
            c21 = [math.log(k) for k in t0]
            c2 = sum(c21)
            
            c3 = 0
            j = 0
            while j < len(mag):
                c3 = c3 + mag[j]**2/t0[j]
                j += 1   
            C = c1 + c2 + c3
            
            DELTA = B**2 - 4.0*A*C
            
            if DELTA > 0:
                x1 = (-1.0*B + math.sqrt(DELTA))/(2.0*A)
                x2 = (-1.0*B - math.sqrt(DELTA))/(2.0*A)
                if x1 > 0:
                    
                    sig_range.append(sig_m)
                    mean_range1.append(x1)
                    mean_range2.append(x2)
             
            i += 1    
                    
#        mean_err1 = max(mean_range1) - m_mean
#        mean_err2 = m_mean - min(mean_range2)
#        sig_err1 = max(sig_range) - sig_mean
#        sig_err2 = sig_mean - min(sig_range)
        
        w1_mag_list1.append(mean_range1)
        w1_mag_list2.append(mean_range2)
        w1_sigma_list.append(sig_range)
        
    #    plt.scatter(m_mean,sig_mean)
    #    plt.scatter(mean_range1,sig_range)
    #    plt.scatter(mean_range2,sig_range)
    #    plt.xlabel('<m>')
    #    plt.ylabel('σ_m')
    
        n += 1
    else:
        n += 1
    
    
#WISE w2
w2_m_mean_list = []
w2_sig_mean_list = []

w2_mag_list1 = []
w2_mag_list2 = []
w2_sigma_list = []

w2_sigma_list_99 = []

RA2 = []
DEC2 = []
MAG_w2 = []
ERR_w2 = []
MJD_w2 = []   

w2_list1 = []

n = 0
while n < len(MAG2):
    if len(MAG2[n])>1:
        mag = MAG2[n]
        magerr = ERR2[n]
        
        RA2.append(RA[n])
        DEC2.append(DEC[n])
        MAG_w2.append(MAG2[n])
        ERR_w2.append(ERR2[n])
        MJD_w2.append(MJD1[n])
        
        w2_list1.append(w2_list[n])
        
        nn = len(mag)
        L = 1000
        #循环结构求去sima_m,mean_mag
        sigma_m = []
        mag_mean = []
        S = []
        dx = 0.8/L
        i = 1
        
        a = list(range(L))
        #sigma_m = [0.8*(i+1)/L for i in a]
        while i < L:
        #    a = list(range(L))
        #    sigma_m = [0.8*(i+1)/L for i in a]
            sig_m = i*dx
            sigma_m.append(sig_m)
            
            t0 = [sig_m**2 + k**2 for k in magerr]
            
            summag = 0
            j = 0
            while j < len(mag):
                a = mag[j]/t0[j]
                summag = summag + a
                j += 1
            
            sumerr = 0
            j = 0
            while j < len(magerr):
                b = 1/t0[j]
                sumerr = sumerr + b
                j += 1
            
            mag_m = summag/sumerr
            mag_mean.append(mag_m)
            
            #在已知sigma_m and mag_m下给出S函数值
            s1 = nn*math.log(2*math.pi)
            
            mr = [math.log(sig_m**2 + ri**2) for ri in magerr]
            s2 = sum(mr)
            
            s3 = 0
            j = 0
            while j < len(mag):
                s3 = s3 + (mag[j]-mag_m)**2/t0[j]
                j += 1
            
            s = s1 + s2 + s3
            S.append(s)
            i += 1
        

        
    #    if i > 1:
    #        if S[-1] > S[-2]:break
    
    #求出最大似然估计sigma_m and <m>值    
        i = 0
        while i < len(S):
            if S[i+1] > S[i]:break
            i += 1    
        print(i)        
        sig_mean = sigma_m[i]
        m_mean = mag_mean[i]
        s_min = S[i]
        
        w2_m_mean_list.append(m_mean)
        w2_sig_mean_list.append(sig_mean)
#        n += 1
#    else:
#        n += 1
        
        #calculate the m_mean-sigma region of certain confidence level
        #delts=2.3,4.61,9.21 corresponding to 68%, 90%, 99% cinfidence level
        sig_range = []
        mean_range1 = []
        mean_range2 = []
        mean_err_99 = []
        sig_err_99 = []
        
        i = 1
        L = 1000
        dx = 0.8/L
        
        while i < L:
            sig_m = i*dx
            sigma_m.append(sig_m)
            
            t0 = [sig_m**2 + k**2 for k in magerr]
            
            a = [1/k for k in t0]
            A = sum(a)
            
            b = 0
            j = 0
            while j < len(mag):
                b = b + mag[j]/t0[j]
                j += 1
                
            B = -2.0*b
                
            c1 = nn*math.log(2*math.pi) - s_min - 9.21
            
            c21 = [math.log(k) for k in t0]
            c2 = sum(c21)
            
            c3 = 0
            j = 0
            while j < len(mag):
                c3 = c3 + mag[j]**2/t0[j]
                j += 1   
            C = c1 + c2 + c3
            
            DELTA = B**2 - 4.0*A*C
            
            if DELTA > 0:
                x1 = (-1.0*B + math.sqrt(DELTA))/(2.0*A)
                x2 = (-1.0*B - math.sqrt(DELTA))/(2.0*A)
                if x1 > 0:
                    
                    sig_range.append(sig_m)
                    mean_range1.append(x1)
                    mean_range2.append(x2)
             
            i += 1    
                    
#        mean_err1 = max(mean_range1) - m_mean
#        mean_err2 = m_mean - min(mean_range2)
#        sig_err1 = max(sig_range) - sig_mean
#        sig_err2 = sig_mean - min(sig_range)
        
        w2_mag_list1.append(mean_range1)
        w2_mag_list2.append(mean_range2)
        w2_sigma_list.append(sig_range)
        
    #    plt.scatter(m_mean,sig_mean)
    #    plt.scatter(mean_range1,sig_range)
    #    plt.scatter(mean_range2,sig_range)
    #    plt.xlabel('<m>')
    #    plt.ylabel('σ_m')
    
        n += 1
    else:
        n += 1




#k = 0
#n = 0
#while k < len(MAG1):
#    s2 = sigma_list[k]
#    s1 = sigma_list[k]
#    s11 = s1[::-1]
#    for i in s11:
#        s2.append(i)
#    s2.append(s2[0])
#    if min(s2)==0.0008:
#        n += 1
#    k += 1
        

#画counter
        
##WISE w2
#w2_m_mean_list = []
#w2_sig_mean_list = []
#
#w2_mag_list1 = []
#w2_mag_list2 = []
#w2_sigma_list = []
#
#w2_sigma_list_99 = []        
        
#画counter     
#j = 0
#w1range = []
#w2range = []
#
#s1range = []
#s2range = []
##w1range = w1_mag_list1 + w1_mag_list2
##w2range = w2_mag_list1 + w2_mag_list2
#while j < len(w1_mag_list1):
#    w1range.append(w1_mag_list1[j]+w1_mag_list2[j][::-1])
#    w2range.append(w2_mag_list1[j]+w2_mag_list2[j][::-1])
#    
#    s1range.append(w1_sigma_list[j]+w1_sigma_list[j][::-1])
#    s2range.append(w2_sigma_list[j]+w2_sigma_list[j][::-1])
#    #w1range.append(w1_mag_list1 + w1_mag_list2)
#    #w2range.append(w2_mag_list1 + w2_mag_list2)
#    j += 1



##统计counter大于0的源    
#j = 0
#n = 0
#rra = []
#ddec = []
#q1 = []
#q2 = []
#while j < len(w1_sigma_list):
#    if min(w1_sigma_list[j])>0.0008 and min(w2_sigma_list[j])>0.0008:
#        n += 1
#        j += 1
###        plt.scatter(w1_mag_list1[j],w1_sigma_list[j],c='black',marker='.')
###        plt.scatter(w1_mag_list1[j],w1_sigma_list[j],c='black',marker='.')
##        plt.scatter(w2_m_mean_list[j],w2_sig_mean_list[j],c='red',marker='+')
##        plt.scatter(w2range[j], s2range[j], c='black',marker='.')
##        plt.show()
#        rra.append(RA1[j])
#        ddec.append(DEC1[j])
#        q1.append(w1_sig_mean_list[j])
#        q2.append(w2_sig_mean_list[j])
#    else:
#        j += 1
#        
#print(n)
        



##画出置信区间
RA3 = []
DEC3 =[]
ss1 = []#w1的光变幅度
ss2 = []#w2的光变幅度
    

k = 0
while k < len(RA1):#MAG1
    if min(w1_sigma_list[k])>0.0008 and min(w2_sigma_list[k])>0.0008:
        
        #sigma value
        ss1.append(w1_sig_mean_list[k])
        ss2.append(w2_sig_mean_list[k])
        
        ###W1
        
        RA3.append(RA1[k])
        DEC3.append(DEC1[k])
        
        
        m2 = w1_mag_list2[k]
        m1 = w1_mag_list1[k]
        m11 = m1[::-1]
        for i in m11:
            m2.append(i)
        m2.append(m2[0])
        
        s2 = w1_sigma_list[k]
        s1 = w1_sigma_list[k]
        s11 = s1[::-1]
        for i in s11:
            s2.append(i)
        s2.append(s2[0])
        
        plt.scatter(w1_m_mean_list[k],w1_sig_mean_list[k],c='red',marker='+')
        plt.scatter(m2,s2,c='black',marker='.')
        plt.xlabel(r'$<m>_{w1}$')
        plt.ylabel(r'$σ_{w1}$(mag)')
        plt.title(str(RA1[k]) + '+' + str(DEC1[k]))
#        plt.savefig(str(k)+'W1'+'.pdf')
        plt.show()
        
        print(w1_sig_mean_list[k])
        ####W2
        
#        RA3.append(RA1[k])
#        DEC3.append(DEC1[k])
        
        m2 = w2_mag_list2[k]
        m1 = w2_mag_list1[k]
        m11 = m1[::-1]
        for i in m11:
            m2.append(i)
        m2.append(m2[0])
        

        s2 = w2_sigma_list[k]
        s1 = w2_sigma_list[k]
        s11 = s1[::-1]
        for i in s11:
            s2.append(i)
        s2.append(s2[0])
        
        plt.scatter(w2_m_mean_list[k],w2_sig_mean_list[k],c='red',marker='+')
        plt.scatter(m2,s2,c='black',marker='.')
        plt.xlabel(r'$<m>_{w2}$')
        plt.ylabel(r'$σ_{w2}$(mag)')
        plt.title(str(RA1[k]) + '+' + str(DEC1[k]))
#        plt.savefig(str(k)+'W2'+'.pdf')
        plt.show()
        
        print(w2_sig_mean_list[k])
        
        plt.errorbar(MJD_w1[k], MAG_w1[k], ERR_w1[k], color='red', linestyle='none', marker='o', label='W1' )
        plt.errorbar(MJD_w2[k], MAG_w2[k], ERR_w2[k], color='blue', linestyle='none', marker='*', label='W2')
        plt.xlabel('MJD')
        plt.ylabel('W1 & W2(mag)')
        plt.legend()
#        plt.savefig(str(k)+'lightcurve'+'.pdf')
        plt.show()
        #save
        
#        print(MAG_w1[k])
#        print(MAG_w2[k])
        print(k)
        
    #    plt.savefig('W1_k_sigma.pdf')
        k += 1
    #    plt.show()
    else:
        k += 1
##        
###Variability objects
#file = open('Reines2013_Var.txt','w')
#i = 0
#while i < len(RA3):
#    a = RA3[i]
#    b = DEC3[i]
#    a1 = '%011.7f' % a
#    b1 = '%011.7f' % b
#    file.write(' ' + a1 + ' ' + b1 + '\n')
#    i += 1
#
#file.close()

#plt.scatter(ss1,ss2)
#plt.xlabel(r'$\sigma_{w1}$')
#plt.ylabel(r'$\sigma_{w2}$')
#plt.savefig('Sartori10_Var.pdf')

    
#dd = 0
#while dd < len(w1_list1[40]):
#    plt.scatter(mjd_list1[40][dd],w1_list1[40][dd],color='black')
#    dd += 1    
    
m = 0
RA4 = []
DEC4 = []
while m < len(RA3):
    if ss1[m]>0.2 and ss2[m]>0.2:
        RA4.append(RA3[m])
        DEC4.append(DEC3[m])
        m += 1
        
    else:
        m += 1
        
#file = open('Reines2013_Var.txt','w')
#i = 0
#while i < len(RA3):
#    a = RA3[i]
#    b = DEC3[i]
#    a1 = '%011.7f' % a
#    b1 = '%011.7f' % b
#    file.write(' ' + a1 + ' ' + b1 + '\n')
#    i += 1
    
    