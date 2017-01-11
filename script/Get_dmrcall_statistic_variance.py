import os,numpy,scipy,random,math
from scipy import stats
from numpy import matrix
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from scipy.stats import t
import rpy2.robjects as robjects
r = robjects.r
from scipy.stats import norm

##############

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


def get_statistic_pvalue(Tumor,Tumor_delt,Normal,Normal_delt):
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
          hat_1=hat[0:,0:ncase]
          hat_2=hat[0:,ncase:(ncase+ncntl)]
          sum1=((hat_1*(hat_1.T))[1,1])*s2case
          sum2=((hat_2*hat_2.T)[1,1])*s2cntl
          se_new=math.sqrt(sum1+sum2)
          statis=coef[1,0]/se_new
          pvalue=2*t.cdf(-abs(statis),(ncase+ncntl)-2)
          return(statis,pvalue)



def get_statistic_pvalue_file(ranksum_file_name,linearmethod_file,tumor_file,control_file,PURITY_dict):
          tumor=open(tumor_file,"r")
          tumorname_line=tumor.readline()
          TUMOR_SAMPLE=tumorname_line.strip().split(" ")
          del TUMOR_SAMPLE[0]
          tumor_list=tumor.readlines()
          normal=open(control_file,"r")
          normal.readline()
          normal_list=normal.readlines()
          cpgsite_sum=len(tumor_list)
          out1=open(ranksum_file_name,"w")
          out1.write("CpgSite\tStatistic\tP-value\n")
          out2=open(linearmethod_file,"w")
          out2.write("CpgSite\tStatistic\tP-value\n")
          for i in range(0,cpgsite_sum):
                    tumor_line=tumor_list[i].strip().split(" ")
                    cpgsite=tumor_line[0]
                    del tumor_line[0]
                    normal_line=normal_list[i].strip().split(" ")
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
                                                  #Tumor.append(float(beta))
                                                  Tumor_delt.append(tumor_delt)
                              Normal=[]
                              Normal_delt=[]
                              for arr in normal_line:
                                        if arr=="NA":
                                                  continue
                                        else:
                                                  transform_value=math.asin(2*float(arr)-1)
                                                  #Normal.append(float(arr))
                                                  Normal.append(transform_value)
                                                  Normal_delt.append(0)
                              if len(Tumor)>5 and len(Normal)>5:
                                        ranksum = scipy.stats.ranksums(Tumor,Normal)
                                        ranksum_P_value=ranksum[1]
                                        ranksum_Statistic=ranksum[0]
                                        out1.write(cpgsite+"\t"+str(ranksum_Statistic)+"\t"+str(ranksum_P_value)+"\n")
                                        (linear_Statistic,linear_Pvalue)=get_statistic_pvalue(Tumor,Tumor_delt,Normal,Normal_delt)
                                        out2.write(cpgsite+"\t"+str(linear_Statistic)+"\t"+str(linear_Pvalue)+"\n")
                              else:
                                        continue
          out1.close()
          out2.close()


          
#######################################################################
         



def run_pipline(TissueType):
          workdir="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/statistic_data_halfsamples/"
          os.chdir(workdir)
          os.system("mkdir "+TissueType)
          os.chdir(workdir+TissueType)
          tumor_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"+TissueType+"/Tumor_initial.txt"
          control_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"+TissueType+"/Normal_initial.txt"
          purity_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/Infinium_mode/"+TissueType+"/TCGA_"+TissueType+".txt"
          PURITY_dict=get_purity_dict(purity_file)
          get_statistic_pvalue_file("I_Ranksum_Statistic_Pvalue.txt","I_Linear_Statistic_Pvalue_variance.txt",tumor_file,control_file,PURITY_dict)

          
def main():
          file_dir = "/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"
          print "file_dir=:"+file_dir
          TumorType = os.listdir(file_dir)
          for TissueType in TumorType:
                    run_pipline(TissueType)
          



if __name__ == "__main__":
          try:
                    main()
          except KeyboardInterrupt:
                    sys.stderr.write("User interrupt me, see you!\n")
                    sys.exit(0)


