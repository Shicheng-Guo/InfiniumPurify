import os,numpy,scipy,random,math
from scipy import stats
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

import rpy2.robjects as robjects
r = robjects.r


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
                              ranksum = scipy.stats.ranksums(Tumor,Normal)
                              ranksum_P_value=ranksum[1]
                              ranksum_Statistic=ranksum[0]
                              out1.write(cpgsite+"\t"+str(ranksum_Statistic)+"\t"+str(ranksum_P_value)+"\n")
                              
                              Normal.extend(Tumor)
                              Normal_delt.extend(Tumor_delt)
                              robjects.globalenv["x"] = robjects.FloatVector(Normal_delt)
                              robjects.globalenv["y"] = robjects.FloatVector(Normal)
                              lm = r.lm("y ~ x")
                              summary_list=r.summary(lm).rx2('coefficients')
                              line_Statistic=summary_list[5]
                              line_P_value=summary_list[7]
                              out2.write(cpgsite+"\t"+str(line_Statistic)+"\t"+str(line_P_value)+"\n")
          out1.close()
          out2.close()


          
#######################################################################
         



def run_pipline(TissueType):
          workdir="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Statistic_Data/"
          os.chdir(workdir)
          os.system("mkdir "+TissueType)
          os.chdir(workdir+TissueType)

          tumor_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"+TissueType+"/Tumor_initial.txt"
          control_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"+TissueType+"/Normal_initial.txt"
          
          purity_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/Infinium_mode/"+TissueType+"/TCGA_"+TissueType+".txt"
          PURITY_dict=get_purity_dict(purity_file)
          get_statistic_pvalue_file("ranksum_statistic_pvalue.txt","linear_statistic_pvalue.txt",tumor_file,control_file,PURITY_dict)


               
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


