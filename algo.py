from __future__ import division
import numpy as np
from scipy import special
from progressbar import *
from decimal import *

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

def Get_Bin(Input, nbins):
	index = np.empty_like(Input, dtype=int)
	np.ceil(Input*nbins, index)
	index = index - 1
	return index


def Get_Density_NHWnon(Input, Density_NHW_non, mbest):
	#calculate the estimated density of p value of non-functional loci
	#numpy.ndarray Input: a vector of numbers ranging from 0 to 1
	#numpy.ndarray Density_NHW_non: estimated density of non-functional p value from histogram
	#float mbest: the best number of bins for histogram from cross validation
	#numpy.ndarray Output: the estimated density at Input
	Index = Get_Bin(Input, mbest)
	#Index = np.array(list(map(math.ceil,Input*mbest)),dtype=int) - 1
	Output = Density_NHW_non[Index]
	return Output


def Conditional_Exp_NHW(data, theta, Density_NHW_non, mbest):
	#calculate the conditional probability of functional given the GWAS p value
	#numpy.ndarray data: GWAS p value
	#theta: parameter of beta and pi
	#numpy.ndarray Density_NHW_non: estimated density of non-functional p value from histogram
	#float mbest: the best number of bins for histogram from cross validation
	#numpy.ndarray result: conditional probility of being functional given the GWAS p value
	#part1 = (data**(theta[1]-1))/special.beta(theta[1], 1)
	part1 = np.power(data, (theta[1]-1))/special.beta(theta[1],1)
	part2 = Get_Density_NHWnon(data, Density_NHW_non, mbest)
	result = theta[0]*part1/(theta[0]*part1 + (1-theta[0])*part2)
	return result


def histogram(data, nbins):
	indeces = Get_Bin(data, nbins)
	hist = np.zeros(nbins, dtype=int)
	for index in indeces:
		hist[index]+=1
	result = np.empty(nbins, dtype=np.float64)
	np.divide(hist*nbins, len(data), result)
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
	Density_NHW_non = histogram(Pvalue_NHW_non, bin_num)

	###   EM   ###
	Theta = np.array([0.01, 0.5], dtype=np.float64)
	Theta_Trace = []
	widgets = ['EM: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
	pbar = ProgressBar(widgets=widgets, maxval=ite).start()
	for j in range(ite):
		T1 = Conditional_Exp_NHW(Pvalue_NHW_func, Theta, Density_NHW_non, bin_num)
		Update = np.empty(2, dtype=np.float64)
		Update[0] = np.mean(T1)
		Update[1] = Update[0]/np.mean(T1*(-np.log(Pvalue_NHW_func)))
		Theta = Update
		Theta_Trace.append(Update)
		pbar.update(j)
	
	pbar.finish()
	###   Calculate the Posterior   ###
	part1 = np.power(data1, (Theta[1]-1))/special.beta(Theta[1], 1)
	part2 = Get_Density_NHWnon(data1,Density_NHW_non,bin_num)
	Prior = Theta[0]*data2
	Posterior = Prior*part1/(Prior*part1 + (1-Prior)*part2)
	return {'Post':Posterior,'Trace':Theta_Trace} 



##############test example##############
###read in data
Chr_NHW = []
Pos_NHW = []
Pvalue_NHW = []
GenoCanyon_10K_NHW = []

with open('/home/dydyd/test/GenoCanyonData/test1/Pvalue','r') as dat:
	for eachline in dat:
		each = eachline.split('\n')
		Pvalue_NHW.append(Decimal(each[0]))

with open('/home/dydyd/test/GenoCanyonData/test1/Pos','r') as dat:
	for eachline in dat:
		each = eachline.split('\n')
		Pos_NHW.append(each[0])

with open('/home/dydyd/test/GenoCanyonData/test1/Chr','r') as dat:
	for eachline in dat:
		each = eachline.split('\n')
		Chr_NHW.append(each[0])

with open('/home/dydyd/test/GenoCanyonData/test1/GenoCanyon','r') as dat:
	for eachline in dat:
		each = eachline.split('\n')
		GenoCanyon_10K_NHW.append(Decimal(each[0]))


#Chr_NHW = np.array(Chr_NHW)
#Pos_NHW = np.array(Pos_NHW)
Pvalue_NHW = np.array(Pvalue_NHW, dtype=np.float64)
GenoCanyon_10K_NHW = np.array(GenoCanyon_10K_NHW, dtype=np.float64)

results = post_GWAS_posterior(data1 = Pvalue_NHW, data2 = GenoCanyon_10K_NHW, thd = 0.1, cv = 100, ite = 5000)

res = results['Post']
n = len(Chr_NHW)
output = ['Chr.ID'+'\t'+'Chr.Position'+'\t'+'P.value'+'\t'+'Posterior'+'\n']
for i in range(n):
	output.append(Chr_NHW[i]+'\t'+Pos_NHW[i]+'\t'+str(Pvalue_NHW[i])+'\t'+str(res[i])+'\n')


with open('result.data','w') as outfile:
	outfile.writelines(output)
# with open('result_theta','w') as outfile:
# 	trace = results["Trace"]
# 	for vector in trace:
# 		outfile.write(str(vector[0])+'\t'+str(vector[1])+'\n')




 
