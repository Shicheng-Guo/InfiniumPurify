#!/opt/bin/python

#######   Metadata ####################################################
#  Metadata-Version: 3.1
#  Name: InfiniumDM
#  Version: 3.0-20160330
#  Summary: Differential methylation analysis accounting for tumor purity 
#  Home-page: UNKNOWN
#  Author: Naiqian Zhang & Xiaoqi Zheng
#  Author-email: naiqian@wfu.edu.cn & xqzheng@shnu.edu.cn
#  License: Shanghai normal university, Amory University
#  updates in 20160915:
#  1. Use a universial set of normal samples to get informative DMPs for purity estimation. Now works for most tumor types in TCGA
#  2. Proposed a linear model for differential methylation analysis accounting for tumor purities.
#  3. Control free DM calling by testing correlation between methylation levels and estimated tumor purities. 
#######################################################################

##################
#  load package
##################

import os,numpy,scipy,random,math,sys,logging
from scipy import stats
from numpy import matrix
import numpy as np
from optparse import OptionParser
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from scipy.stats import t
import rpy2.robjects as robjects
r = robjects.r
from scipy.stats import norm


##################
#  functions
##################

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
    
def get_purity_dict(purity_file):
          # input is the purity file, output is the purity dict.
          PURITY_dict={}
          DAT=open(purity_file,"r")
          DAT.readline()
          for line in DAT:
                    arr=line.strip().split("\t")
                    sample_name=arr[0]
                    purity=arr[1]
                    PURITY_dict[sample_name]=purity
          return(PURITY_dict)


def get_s2(Tumor,Tumor_delt,Normal,Normal_delt):
          ncase=len(Tumor)
          ncntl=len(Normal)
          x1=[1 for i in range(0,len(Tumor)+len(Normal))]
          Tumor_delt.extend(Normal_delt)
          Tumor.extend(Normal)
          y=matrix([Tumor])
          trans_y=y.T
          x=matrix([x1,Tumor_delt])
          trans_x=x.T
          xtrans_x=x*trans_x
          inver_xtrans_x=xtrans_x.I
          hat=inver_xtrans_x*x
          coef=hat*trans_y
          ypred=(trans_x*coef).T
          resi=y-ypred
          resi_case=resi[0,0:ncase]
          resi_case_trans=resi_case.T
          s2case=(resi_case*resi_case_trans)[0,0]/(ncase-2)

          resi_cntl=resi[0,ncase:(ncase+ncntl)]
          resi_cntl_trans=resi_cntl.T
          s2cntl=(resi_cntl*resi_cntl_trans)[0,0]/(ncntl-2)
          
          return(s2case,s2cntl)


def shrinker(vv,Mean):
    tmp = math.log(vv)
    out = math.exp((tmp + Mean)/2.0)
        
    return(out)

def get_statistic_pvalue(Tumor,Tumor_delt,Normal,Normal_delt,S2Case_m,S2Cntl_m):
          ncase=len(Tumor)
          ncntl=len(Normal)
          x1=[1 for i in range(0,len(Tumor)+len(Normal))]
          Tumor_delt.extend(Normal_delt)
          Tumor.extend(Normal)
          y=matrix([Tumor])
          trans_y=y.T
          x=matrix([x1,Tumor_delt])
          trans_x=x.T
          xtrans_x=x*trans_x
          inver_xtrans_x=xtrans_x.I
          hat=inver_xtrans_x*x
          coef=hat*trans_y
          ypred=(trans_x*coef).T
          resi=y-ypred
          resi_case=resi[0,0:ncase]
          resi_case_trans=resi_case.T
          s2case=(resi_case*resi_case_trans)[0,0]/(ncase-2)
          s2case=shrinker(s2case,S2Case_m)
          
          resi_cntl=resi[0,ncase:(ncase+ncntl)]
          resi_cntl_trans=resi_cntl.T
          s2cntl=(resi_cntl*resi_cntl_trans)[0,0]/(ncntl-2)
          s2cntl=shrinker(s2cntl,S2Cntl_m)
          
          hat_1=hat[0:,0:ncase]
          hat_2=hat[0:,ncase:(ncase+ncntl)]
          sum1=((hat_1*(hat_1.T))[1,1])*s2case
          sum2=((hat_2*hat_2.T)[1,1])*s2cntl
          se_new=math.sqrt(sum1+sum2)
          statis=coef[1,0]/se_new
          pvalue=2*t.cdf(-abs(statis),(ncase+ncntl)-2)
          return(statis,pvalue)




def DMCalling(output_file,tumor_file,control_file,PURITY_dict):
          tumor=open(tumor_file,"r")
          tumorname_line=tumor.readline()
          TUMOR_SAMPLE=tumorname_line.strip().split("\t")
          #del TUMOR_SAMPLE[0]
          tumor_list=tumor.readlines()
          normal=open(control_file,"r")
          normal.readline()
          normal_list=normal.readlines()
          cpgsite_sum=len(tumor_list)
          
          S2Case_m = 0
          S2Cntl_m = 0
          step = 0
          for i in range(0,cpgsite_sum):
                    tumor_line=tumor_list[i].strip().split("\t")
                    cpgsite=tumor_line[0]
                    del tumor_line[0]
                    normal_line=normal_list[i].strip().split("\t")
                    del normal_line[0]
                    if tumor_line.count("NA")==len(tumor_line) or normal_line.count("NA")==len(normal_line):
                              continue
                    else:
                              Tumor=[]
                              Tumor_delt=[]
                              for k in range(0,len(TUMOR_SAMPLE)):
                                        sample=TUMOR_SAMPLE[k]
                                        beta=tumor_line[k]
                                        if beta=="NA":
                                                  continue
                                        else:
                                                  transform_beta=math.asin(2*float(beta)-1)
                                                  tumor_delt=float(PURITY_dict[sample])
                                                  Tumor.append(transform_beta)
                                                  Tumor_delt.append(tumor_delt)
                              Normal=[]
                              Normal_delt=[]
                              for arr in normal_line:
                                        if arr=="NA":
                                                  continue
                                        else:
                                                  transform_value=math.asin(2*float(arr)-1)
                                                  Normal.append(transform_value)
                                                  Normal_delt.append(0)

                                                  
                              if len(Tumor)>2 and len(Normal)>2:
                                         
                                        (s2case,s2cntl)=get_s2(Tumor,Tumor_delt,Normal,Normal_delt)
                                        S2Case_m = S2Case_m + math.log(s2case)
                                        S2Cntl_m = S2Cntl_m + math.log(s2cntl)
                                        step = step + 1
                              else:
                                        continue

          S2Case_m = S2Case_m/step
          S2Cntl_m = S2Cntl_m/step
          
          
          out2=open(output_file,"w")
          out2.write("CpgSite\tStatistic\tP-value\n")              
          for i in range(0,cpgsite_sum):
                    tumor_line=tumor_list[i].strip().split("\t")
                    cpgsite=tumor_line[0]
                    del tumor_line[0]
                    normal_line=normal_list[i].strip().split("\t")
                    del normal_line[0]
                    if tumor_line.count("NA")==len(tumor_line) or normal_line.count("NA")==len(normal_line):
                              out2.write(cpgsite+"\tNA\tNA\n")
                    else:
                              Tumor=[]
                              Tumor_delt=[]
                              for k in range(0,len(TUMOR_SAMPLE)):
                                        sample=TUMOR_SAMPLE[k]
                                        beta=tumor_line[k]
                                        if beta=="NA":
                                                  continue
                                        else:
                                                  transform_beta=math.asin(2*float(beta)-1)
                                                  tumor_delt=float(PURITY_dict[sample])
                                                  Tumor.append(transform_beta)
                                                  Tumor_delt.append(tumor_delt)
                              Normal=[]
                              Normal_delt=[]
                              for arr in normal_line:
                                        if arr=="NA":
                                                  continue
                                        else:
                                                  transform_value=math.asin(2*float(arr)-1)
                                                  Normal.append(transform_value)
                                                  Normal_delt.append(0)

                                                  
                              if len(Tumor)>2 and len(Normal)>2:
                                         
                                        (linear_Statistic,linear_Pvalue)=get_statistic_pvalue(Tumor,Tumor_delt,Normal,Normal_delt,S2Case_m,S2Cntl_m)
                                        out2.write(cpgsite+"\t"+str(linear_Statistic)+"\t"+str(linear_Pvalue)+"\n")
                              else:
                                        out2.write(cpgsite+"\tNA\tNA\n")

          out2.close()


## ctrlFree DM calling functions

def DMCalling_ctrlFree(output,tumorFile,PURITY_dict): #(TissueType,DAT_file,TUMOR_dict,TUMOR,PURITY):
    out=open(output,"w")
    out.write("CpgSite\tP-value\n")
    DAT=open(tumorFile,"r")

    TUMOR = DAT.readline().strip().split("\t")

    for line in DAT:
        arr = line.strip().split("\t")
        cpgsite = arr[0]
        
        METHY=[]
        PURITY_NONENA=[]
        for i in range(1,len(arr)):
            sample=TUMOR[i-1]
            methy=arr[i]
            purity=PURITY_dict[sample]
            if methy!="NA":
                METHY.append(float(methy))
                PURITY_NONENA.append(float(purity))
            else:
                continue


        if len(METHY)!=0:
            gradient, intercept, r_value, p_value, std_err = stats.linregress(PURITY_NONENA,METHY)
            beta_mean = np.mean(METHY)
            beta_std = np.std(METHY)

            ## get probability 0408
            beta = gradient
            s2 = 0
            for i in range(len(PURITY_NONENA)):
                s2 += (gradient*PURITY_NONENA[i] + intercept - METHY[i])**2
            s2 = s2/len(PURITY_NONENA)

            xx = 0
            for i in range(len(PURITY_NONENA)):
                xx += (PURITY_NONENA[i] - np.mean(PURITY_NONENA))**2

            s2_beta = s2 / xx

            threshold=0.1
            postprob = norm.cdf(beta - threshold,0,(s2_beta)**(1/2.0)) + 1 - norm.cdf(beta + threshold,0,(s2_beta)**(1/2.0))

        else:

            postprob = "NA"


        out.write(cpgsite+"\t"+str(postprob)+"\n")

    out.close()


          
##################
# main function
##################


def main():
    usage = "usage: python %prog <-t tumorFile> <-p purityFile> [-n normalFile] [...]"
    description = "Differential Methylation analysis accounting for tumor purities. If normal samples are provided, a generalized least square model is executed for DM analysis. If there is no normal control, a control-free DM model is used for testing the regression coeffcient. For example: python %prog -t tumor.txt -n normal.txt -p purity.txt. -t and -p options are needed!"
    op = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-t","--tumorFile",dest="tumorFile",type="str",
                  help="tumor file, rows are CpG sites, columns are tumor samples.")
    op.add_option("-p","--purityFile",dest="purityFile",type="str",
                  help="tumor purities estimated by InfiniumPurify, should have the same number with tumor files.")
    op.add_option("-n","--normalFile",dest="normalFile",type="str",
                  help="normal file, rows are CpG sites, columns are normal samples.")


    (options,args) = op.parse_args()


    if not options.tumorFile or not options.purityFile:
        op.print_help()
        sys.exit(1)


    ## if having normal control file
    if options.normalFile:
        
        tumorFile = options.tumorFile
        purityFile = options.purityFile
        normalFile = options.normalFile
        
        info("Input tumor samples are: " + tumorFile)
        info("Input purities for tumor samples are: " + purityFile)
        info("Input normal samples are: " + normalFile) 

        #-------- step 1: loading purityFile  -------#
        info("Loading purityFile ...")
        PURITY_dict = get_purity_dict(purityFile)

        
        #-------- step 2: DM analysis accounting for tumor purity -------#
        
        info("DM analysis accounting for tumor purity ...")
        DMCalling("DMC.txt",tumorFile,normalFile,PURITY_dict)
        print("Result is saved in: "+os.getcwd()+"/DMC.txt")


    ## if no normal control file        
    else:
        
        tumorFile = options.tumorFile
        purityFile = options.purityFile
        
        info("Input tumor samples are: " + tumorFile)
        info("Input purities for tumor samples are: " + purityFile)

        #-------- step 1: loading purityFile  -------#
        info("Loading purityFile ...")
        PURITY_dict = get_purity_dict(purityFile)

        
        #-------- step 2: DM analysis accounting for tumor purity -------#
        
        info("Control-free DM analysis accounting for tumor purity ...")
        DMCalling_ctrlFree("DMC_ctrlFree.txt",tumorFile,PURITY_dict)
        print("Result is saved in: "+os.getcwd()+"/DMC_ctrlFree.txt")
            

            
if __name__ == "__main__":
    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)    
