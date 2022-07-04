import numpy as np
import pandas as pd
from scipy.stats import norm
import itertools
from statistics import stdev as sd
from typing import Any, Dict, List, Optional
from math import factorial
import math
from math import lgamma, floor, ceil

"""
Python port of iNEXT.ind function from iNEXT 

Code: https://github.com/JohnsonHsieh/iNEXT/blob/master/R/iNEXT.r
Paper: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12613

Descriptions: where count is a list of ints == counts/sample (counts is typically equivilant to "reads")
and endpoint is an integer meaning total number of counts or reads to extrapolate to. 
Uses a Baysian approach to interpolate and extrapolate species richness by applying Chao1 matrix calculations
Output is a dict of named lists, each contains n == # species objects for:

['reads', 'method', 'species', 'species_lower_conf',
       'species_upper_conf', 'coverage', 'coverage_lower_conf',
       'coverage_upper_conf']

Usage:
import pynext
import pandas as pd #optional 

count = [1, 5, 10, 17... ] 
endpoint = 500000
out = pynext(count, endpoint)
df  = pd.DataFrame

"""

seq = np.linspace #nod to R

def estimate_community(spec):
    """Use Chao1 estimates to extrapolate relative abundance"""
    sobs = len(spec) #observed species
    n = sum(spec) #sample size
    f1 = sum(i == 1 for i in spec) #singleton
    f2 = sum(i == 2 for i in spec) #doubleton

    #estimation of unseen species via Chao1
    if f2 == 0:
        f0_hat = (n - 1) / n * f1 * (f1 - 1) / 2
    else: 
        f0_hat = (n - 1) / n * f1 **2 / 2 / f2

    if f1 > 0:
        cap_a = n*f0_hat/(n*f0_hat+f1)
    else:
        cap_a = 1
    a = f1/n*cap_a
    b = sum(spec / n * (1 - spec / n) ** n)
  
    #adjusted factor for rare species in the sample
    if f0_hat==0:
        w = 0
    else:
        w = a / b

    #estimation of relative abundance of observed species in the sample
    prob_hat = spec / n * (1 - w * (1 - spec / n) ** n)
    #estimation of relative abundance of unseen species in the sample
    prob_hat_unseen = np.linspace(a/ceil(f0_hat), a/ceil(f0_hat), num=ceil(f0_hat))
    
    #Output: a vector of estimated relative abundance
    seen_unseen = list(prob_hat) + list(prob_hat_unseen)
    return sorted(seen_unseen, reverse = True)

def pymultinom(n, size, prob):
    """port of stats::rmultinom"""
    return np.random.multinomial(size, prob, size = n)

def choose(k, n):
    """
    Port of R choose
    https://www.educative.io/answers/what-is-the-choose-function-in-r
    """
    if k < n:
        return 0
    else: 
        return int(factorial(k)/(factorial(n)*factorial(k-n)))

choose_vec = np.vectorize(choose, excluded=['n'])

def exp(x):
    """port of R exp"""
    return math.e**x

def determine_abundance(x , zero = False ) :
    n = sum ( x ) 
    f1 = sum ( x == 1 )
    f2 = sum ( x == 2 ) 
    f3 = sum ( x == 3 )

    if f2 == 0:
        f1 = max ( f1 - 1 , 0 ) 
        f2 = 1 

    A1 = f1 / n * ( ( n - 1 ) * f1 / ( ( n - 1 ) * f1 + 2 * max ( f2 , 1 ) ) ) 
    A2 = f2 / choose ( n , 2 ) * ( ( n - 2 ) * f2 / ( ( n - 2 ) * f2 + 3 * max ( f3 , 1 ) ) ) ** 2 

    if zero == False :
        x = x [ x > 0 ] 

        def q_solve( q ) :
            e = A1 / sum ( x / n * exp ( - q * x ) ) 
            out = sum ( ( x / n * ( 1 - e * exp ( - q * x ) ) ) ** 2 ) - sum ( choose_vec ( x , 2 ) / choose ( n , 2 ) ) + A2 
            return abs(out)
        
        m = minimize_scalar(q_solve, bounds = (0,1), method = "Bounded")
        if m.success:
            q = m.x
        else:
            q = 1
        e = A1 / sum(x/n * exp(-q*x))
        return [e, q]


def extrapolate_richness(counts, read_breaks):
    """port of iNEXT D0.hat"""
    counts = counts[counts > 0]
    total_reads = sum(counts)
    
    def estimate_sample(sampled_nreads):
        """port of Sub where sampled_nreads == objects in the read_breaks (observed & extrapolated)
        where m == #reads list and k == # reads at sampling point"""
        if sampled_nreads <= total_reads:
            #print(f"k <= n: {k}, {n}")
            #---
            def calculate_gammas(rset):
                """
                port of Fun where rset == # reads in a read_set
                """   
                if rset <= (total_reads - sampled_nreads):
                    lg1 = lgamma(total_reads - rset + 1)
                    lg2 = lgamma(total_reads - sampled_nreads + 1)
                    lg3 = lgamma(total_reads - rset - sampled_nreads + 1)
                    lg4 = lgamma(total_reads + 1)
                    return exp(lg1 + lg2 - lg3 - lg4)
                else:
                    return 0
            #---
            return sum([1 - calculate_gammas(count) for count in counts])
        else: 
            sum_obs = sum(counts > 0)
            f1 = sum(counts == 1)
            f2 = sum(counts == 2)
            if f2 == 0:
                f0_hat = (total_reads-1)/total_reads * f1 * (f1 - 1)/2
            else:
                f0_hat = (total_reads-1)/total_reads * f1**2/2/f2
            alpha = total_reads * f0_hat/(total_reads * f0_hat + f1)
            if f1 == 0:
                return sum_obs
            else:
                return sum_obs + f0_hat *(1-alpha**(sampled_nreads - total_reads))

    return [estimate_sample(read_set) for read_set in read_breaks]


def estimate_coverage(counts, read_set):
    """port of Chat.Ind to estimate coverage"""
    counts = counts[counts>0]
    total_reads = sum(counts)
    f1 = sum(counts == 1)
    f2 = sum(counts == 2)
    #f0_hat = ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2) 
    #estimation of unseen species via Chao1
    if f2 == 0:
        f0_hat = (total_reads - 1) / total_reads * f1 * (f1 - 1) / 2
    else:
        f0_hat = (total_reads - 1) / total_reads * f1 ** 2/ 2 / f2
    if f1 > 0:
        alpha = total_reads*f0_hat/(total_reads*f0_hat+f1)
    else:
        alpha = 1
    
    def calculate_gammas(rset):
        if rset < total_reads:
            sampled_counts = counts[(total_reads-counts)>=rset]
            lg1 = [lgamma(sample) for sample in (total_reads-sampled_counts+1)]
            lg2 = [lgamma(sample) for sample in (total_reads-sampled_counts-rset+1)]
            lg3 = lgamma(total_reads)
            lg4 = lgamma(total_reads-rset)
            weights = sampled_counts / total_reads
            lg_diff = [x1 - x2 - lg3 + lg4 for (x1, x2) in zip(lg1, lg2)]
            lg_exp = [exp(d) for d in lg_diff]
            return 1-sum(weights*lg_exp)
        if rset == total_reads:
            return 1-f1/total_reads*alpha
        if rset > total_reads:
            return 1-f1/total_reads*alpha**(rset-total_reads+1)
    return [calculate_gammas(rset) for rset in read_set]


def pynext(species: list, read_set: list=None, endpoint: int=None, knots=40, se=True, nboot = 200, conf=0.95) -> dict:
    """
    A python port of the main iNEXT function iNEXT.ind from R
    Apply a Baysian estimator to interpolate and extrapolate species/reads based on the Chao1 index
    """
    total_reads = sum(species) #observed
    if not endpoint:
        endpoint = 2*total_reads
    read_set = None #m in iNEXT
    if not read_set:
            if endpoint <= total_reads:
                #Interpolate
                read_set = [floor(sample_set) for sample_set in np.linspace(1, endpoint, num=floor(knots))]
            else: 
                #Extrapolate
                read_set = [floor(sample_set) for sample_set in list(seq(1, sum(species)-1, num=floor(knots/2)-1)) + 
                       [sum(species)] + 
                       list(seq(sum(species)+1, endpoint, num=floor(knots/2)))]
            read_set = [1] + read_set[1:]
    else:
        if max(read_set)>total_reads and len(read_set[read_set==total_reads])==0: 
            read_set = read_set + [total_reads-1] + [total_reads] + [total_reads+1]
            read_set = read_set.sort()           

    richness_model = extrapolate_richness(species, read_set)
    coverage_model = estimate_coverage(species, read_set)

    if se and nboot > 0 and len(species) > 1:
        #community probability matrix
        prob_hat = estimate_community(species)
        #community abundance matrix
        abundance_matrix = pymultinom(nboot, total_reads, prob_hat)

        #boostrap richness +- sd
        boot_counts = map(extrapolate_richness, abundance_matrix, itertools.repeat(read_set, len(abundance_matrix)))
        boot_sd = np.apply_along_axis(sd, 0, list(boot_counts))

        #use quantile normalization to estimate richness confidence intervals
        error = norm.ppf(1-(1-conf)/2)*boot_sd
        left  = richness_model - error
        right = richness_model + error

        #boostrap coverage +- sd
        boot_cov = map(estimate_coverage, abundance_matrix, itertools.repeat(read_set, len(abundance_matrix)))
        boot_cov_sd = np.apply_along_axis(sd, 0, list(boot_cov))
        error_cov = norm.ppf(1-(1-conf)/2)*boot_cov_sd
        left_cov = coverage_model - error_cov
        right_cov = coverage_model + error_cov   

        #interpolation/extrapolation method
        def methhead(m):
            if m > total_reads:
                return "Extrapolated"
            elif m < total_reads:
                return "Interpolated"
            else:
                return "Observed"
        method = [methhead(m) for m in read_set] 

        out = {"reads": read_set, 
               "method": method, 
               "species": richness_model,
               "species_lower_conf": left,
               "species_upper_conf": right,
               "coverage": coverage_model,
               "coverage_lower_conf": left_cov, 
               "coverage_upper_conf": right_cov}
        return out
