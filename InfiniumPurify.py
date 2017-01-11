#!/opt/bin/python

#######   Metadata ####################################################
#  Metadata-Version: 3.1
#  Name: InfiniumPurify
#  Version: 3.0-20160330
#  Summary: Tumor purity estimation for methylation 450K array 
#  Home-page: UNKNOWN
#  Author: Naiqian Zhang & Xiaoqi Zheng
#  Author-email: naiqian@wfu.edu.cn & xqzheng@shnu.edu.cn
#  License: Shanghai normal university, Amory University
#  updates in 20160330:
#  1. Use a universial set of normal samples to get informative DMPs for purity estimation. Now works for most tumor types in TCGA
#  2. Proposed a linear model for differential methylation analysis accounting for tumor purities.
#######################################################################



import sys,re,os,sys,random,math,numpy,scipy,time
from optparse import OptionParser
import logging
from scipy import stats
import rpy2.robjects as robjects
r = robjects.r

#random.seed(1234)

logfhd = open("log","w")
def writelog():
    
    

    logging.basicConfig(level=20,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )

    error   = logging.critical        # function alias
    warn    = logging.warning
    
    return logfhd

def info(a):
    logfhd = writelog()
    logging.info(a)
    logfhd.write(a+"\n")
    logfhd.flush()
    

def get_DMPs(CancerType):
    
    """ load DMPs for a given cancer type. Return two lists of hyper and hypo DMPs"""
    if not os.path.exists("./data/"+CancerType+"_hyper_methylated.sites"):
        print("The current version only works for cancer type: '"+CancerType+"'. Please update iDMC first!")
        sys.exit(1)
    
    hyper_sites=open("./data/"+CancerType+"_hyper_methylated.sites","r")
    hypo_sites=open("./data/"+CancerType+"_hypo_methylated.sites","r")
    hyper=hyper_sites.readlines()
    hypo=hypo_sites.readlines()
    HYPER=[]
    HYPO=[]
    for line in hyper:
        arr = line.strip()
        HYPER.append(arr)
    for line in hypo:
        arr = line.strip()
        HYPO.append(arr)

    hyper_sites.close()
    hypo_sites.close()
    
    return(HYPER,HYPO)

def get_peak(dir_450K,HYPER,HYPO):

    """ Get mode of peak of all DMPs. Input: a 450K array data and hyper/hypo DMPs; Output: Peak mode """
    
    DAT = open(dir_450K,"r")
    DAT.readline()
    DAT.readline()
    value=[]
    for line in DAT:
        arr = line.strip().split("\t")
        cpg = str(arr[0])
        mRatio = arr[1]
        if cpg in HYPER and mRatio!="NA":
            value_n=float(mRatio)
            p=value_n
        elif cpg in HYPO and mRatio!="NA":
            value_n=float(mRatio)
            p=1-value_n
        else: continue
        value.append(p)
    value_float=[float(i) for i in value]
    density = stats.kde.gaussian_kde(value_float)
    x = numpy.arange(0., 1, .01)
    density=density(x)
    M=max(density)
    index=density.tolist().index(M)
    peak=x.tolist()[index]

    DAT.close()
    
    return(peak)



def main():
    usage = "usage: python %prog <-f filename> <-t cancer_type> [...]"
    description = "Get tumor purity from 450K array data.\n For example: python %prog -f 450k_example.txt -c LUAD. -f and -c options are needed!"
    op = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-f","--filename",dest="filename",type="str",
                  help="The file name of 450K array, with directory")
    op.add_option("-c","--CancerType",dest="CancerType",type="str",
                  help="Cancer type, in TCGA abbreviation format, can be chosen from 'BRAC','LUAD','COAD',... (See abbr.txt for detail)")

    (options,args) = op.parse_args()

    
    if not options.filename or not options.CancerType :
        op.print_help()
        sys.exit(1)

    filename = options.filename
    CancerType = options.CancerType

    info("Input filename is: " + filename) 
    info("Cancer Type is: " + CancerType) 

   
    #-------- step 1: Load hyper and hypo DMPs -------#

    info("Loading hyper and hypo iDMCs for: %s" % CancerType)
    (HYPER,HYPO)=get_DMPs(CancerType)
    
    #-------- step 2: Compute mode of peak (purity) -------#

    info("Computating peak mode of iDMCs ...")
    purity = get_peak(filename,HYPER,HYPO)

    info("Estimated purity by InfiniumPurify is: "+ str(purity))

                  

    
            
if __name__ == "__main__":
    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)    
