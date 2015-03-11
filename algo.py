from __future__ import division
import numpy as np
from scipy import special
from progressbar import *
from datetime import datetime
###function needed
def cv_hist_fun(x,k):
	#use cross validation, give the best bin numbers for histogram
	#numpy.ndarray x: p values of non functional loci
	#int k: parameter for cross validation
	#float mbest: the best number of bins for histogram from cross validation
	n = len(x)
	a = 0
	b = 1
	nbins = np.arange(k+1)[1:]    ###number of bins
	h = (b-a)/nbins   ###width of bins
	risk = []
	#pb = txtProgressBar(min=0, max=k, style=3)
	widgets = ['Cross Validation: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
	pbar = ProgressBar(widgets=widgets, maxval=len(nbins)).start()
	for i in nbins:
		###get counts N_j
		br = np.linspace(a,b,nbins[i-1]+1)
		N = np.histogram(x,bins=br,density=False)[0]
		risk.append(np.sum((N**2))/(n**2*h[i-1])  - (2/(h[i-1]*n*(n-1)))*np.sum(N*(N-1)))
		pbar.update(i)
	pbar.finish()
		#sum(N^2)/(n^2*h[i])  - (2/(h[i]*n*(n-1)))*sum(N*(N-1))
		#setTxtProgressBar(pb, i)
	
	risk = np.array(risk)
	hbest = h[risk==risk.min()]
	hbest = hbest[0]  ###in case of tie take first (smallest) one
	mbest = (b-a)/hbest   ###optimal number of bins
	return {'risk':risk,'nbins':nbins,'h':h,'mbest':mbest}


def Get_Density_NHWnon(Input, Density_NHW_non, mbest):
	#calculate the estimated density of p value of non-functional loci
	#numpy.ndarray Input: a vector of numbers ranging from 0 to 1
	#numpy.ndarray Density_NHW_non: estimated density of non-functional p value from histogram
	#float mbest: the best number of bins for histogram from cross validation
	#numpy.ndarray Output: the estimated density at Input
	Index = np.array(map(math.ceil,Input*mbest),dtype=int) - 1
	Output = Density_NHW_non[Index]
	return Output


def Conditional_Exp_NHW(data, theta, Density_NHW_non, mbest):
	#calculate the conditional probability of functional given the GWAS p value
	#numpy.ndarray data: GWAS p value
	#theta: parameter of beta and pi
	#numpy.ndarray Density_NHW_non: estimated density of non-functional p value from histogram
	#float mbest: the best number of bins for histogram from cross validation
	#numpy.ndarray result: conditional probility of being functional given the GWAS p value
	part1 = (data**(theta[1]-1))/special.beta(theta[1], 1)
	part2 = Get_Density_NHWnon(data, Density_NHW_non, mbest)
	result = theta[0]*part1/(theta[0]*part1 + (1-theta[0])*part2)
	return result


def post_GWAS_posterior(data1, data2, thd, cv, ite):
	#calculating the posterior score
	#numpy.ndarray data1: GWAS pvalue
	#numpy.ndarray data2: GenoCanyon score 
	#float thd: threshold for defining functional and non-functional
	#int cv: parameter for cross validation
	#int ite: number of iterations for EM algorithm
	#return posterior score and the trace of parameter Theta in EM
	#numpy.ndarray Posterior: posterior score
	#numpy.ndarray Theta_Trace: record the value of Theta in each iteration
	Pvalue_NHW_func = data1[data2 > thd]  
	Pvalue_NHW_non = data1[data2 <= thd]
	bin_num = cv_hist_fun(Pvalue_NHW_non, cv)['mbest']   ### 949 (under constraint < 1000)
	Breaks = np.linspace(0,1,bin_num+1)
	h_NHW = np.histogram(Pvalue_NHW_non,bins=Breaks,density=True)
	Density_NHW_non = h_NHW[0]
	###   EM   ###
	Theta = np.array([0.01, 0.5])
	Update = np.array([0.,0.])
	Theta_Trace = []
	widgets = ['EM: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
	pbar = ProgressBar(widgets=widgets, maxval=ite).start()
	for j in range(ite):
		T1 = Conditional_Exp_NHW(Pvalue_NHW_func, Theta, Density_NHW_non, bin_num)
		Update[0] = np.mean(T1)
		Update[1] = Update[0]/np.mean(T1*(-np.log(Pvalue_NHW_func)))
		Theta = Update
		Theta_Trace.append(Update)
		pbar.update(j)
	
	pbar.finish()
	###   Calculate the Posterior   ###
	part1 = (data1**(Theta[1]-1))/special.beta(Theta[1], 1)
	part2 = Get_Density_NHWnon(data1,Density_NHW_non,bin_num)
	Prior = Theta[0]*data2
	Posterior = Prior*part1/(Prior*part1 + (1-Prior)*part2)
	return {'Post':Posterior,'Trace':Theta_Trace} 



##############test example##############
###read in data
rs = [] #SNP ID
Chr_NHW = []
Pos_NHW = []
Pvalue_NHW = []
GenoCanyon_10K_NHW = []
with open('/Users/huyiming/Desktop/courses spring 2015/post.GWAS.prioritize/CD_phs000130.pha002847.txt','r') as dat:
	for eachline in dat:
		each = eachline.split('\t')
		if each[3] != '':
			rs.append(each[0])
			Chr_NHW.append(each[2])
			Pos_NHW.append(each[3])
			Pvalue_NHW.append(float(each[1]))

with open('/Users/huyiming/Desktop/courses spring 2015/post.GWAS.prioritize/GenoCanyon_CD_NIDDK.txt','r') as dat1:
	for eachline in dat1:
		each = eachline.split('\n')
		GenoCanyon_10K_NHW.append(float(each[0]))


#Chr_NHW = np.array(Chr_NHW)
#Pos_NHW = np.array(Pos_NHW)
Pvalue_NHW = np.array(Pvalue_NHW)
GenoCanyon_10K_NHW = np.array(GenoCanyon_10K_NHW)

results = post_GWAS_posterior(data1 = Pvalue_NHW, data2 = GenoCanyon_10K_NHW, thd = 0.1, cv = 100, ite = 5000)

res = results['Post']
n = len(Chr_NHW)
output = ['SNP.ID'+'\t'+'Chr.ID'+'\t'+'Chr.Position'+'\t'+'P.value'+'\t'+'Posterior'+'\n']
for i in range(n):
	output.append(rs[i]+'\t'+Chr_NHW[i]+'\t'+Pos_NHW[i]+'\t'+str(Pvalue_NHW[i])+'\t'+str(res[i])+'\n')


with open('/Users/huyiming/Desktop/courses spring 2015/post.GWAS.prioritize/dat.txt','w') as outfile:
	outfile.writelines(output)




 
