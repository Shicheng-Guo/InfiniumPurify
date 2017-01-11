import os,numpy,scipy,random
from scipy import stats
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr


############################################################################
def get_near_cpgsite_dict(DAT_file):
          DAT=open(DAT_file,"r")
          DAT.readline()
          DAT.readline()
          cpgsite_position_dict={}
          for line in DAT:
                    arr = line.strip().split("\t")
                    cpg = str(arr[0])
                    gene = arr[2]
                    chrom=arr[3]
                    position=int(arr[4])
                    if chrom!="NA":
                              if chrom not in cpgsite_position_dict.keys():
                                        cpgsite_position_dict[chrom]={}
                              else:
                                        cpgsite_position_dict[chrom][cpg]=position
                    else:continue
          near_cpgsite_dict={}
          for chrom in cpgsite_position_dict.keys():
                    d=cpgsite_position_dict[chrom]
                    sorted_cpgsite=sorted(d, key=d.__getitem__)
                    for i in range(0,len(sorted_cpgsite)-1):
                              cpgsite=sorted_cpgsite[i]
                              position=d[cpgsite]
                              next_cpgsite=sorted_cpgsite[i+1]
                              next_position=d[next_cpgsite]
                              if next_position-position< 25:   #######near cutoff#######
                                        near_cpgsite_dict[cpgsite]=next_cpgsite
          return(near_cpgsite_dict)


def get_spatial_correlation(near_cpgsite_dict,statistic_file):
          DAT=open(statistic_file,"r")
          DAT.readline()
          cpgsite_statistic_dict={}
          for line in DAT:
                    arr = line.strip().split("\t")
                    cpg = str(arr[0])
                    statistic = float(arr[1])#
                    cpgsite_statistic_dict[cpg]=statistic
          NEAR_CPGSITE_1=near_cpgsite_dict.keys()
          NEAR_CPGSITE_2=near_cpgsite_dict.values()
          STATISTIC_CPG=cpgsite_statistic_dict.keys()
          NEAR_CPGSITE_1_intersect=list(set(NEAR_CPGSITE_1).intersection(set(STATISTIC_CPG)))
          NEAR_CPGSITE_2_intersect=list(set(NEAR_CPGSITE_2).intersection(set(STATISTIC_CPG)))

          STATISTIC_1=[]
          STATISTIC_2=[]
          for cpgsite in NEAR_CPGSITE_1_intersect:
                    next_cpgsite=near_cpgsite_dict[cpgsite]
                    if next_cpgsite in NEAR_CPGSITE_2_intersect:
                              statistic=cpgsite_statistic_dict[cpgsite]
                              next_statistic=cpgsite_statistic_dict[next_cpgsite]
                              STATISTIC_1.append(statistic)
                              STATISTIC_2.append(next_statistic)
                    else:
                              continue
          correlation=pearsonr(STATISTIC_1,STATISTIC_2)[0]
          return(correlation)


#######################################################################

          



def run_pipline(TissueType):
          m450k_data_dir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/TCGA/"+TissueType+"/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"

          fileName=os.listdir(m450k_data_dir)
          DAT_file=m450k_data_dir+fileName[1]
          near_cpgsite_dict=get_near_cpgsite_dict(DAT_file)
          ranksum_statistic_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Statistic_Data/"+TissueType+"/I_Ranksum_Statistic_Pvalue.txt"
          correlation_ranksum=get_spatial_correlation(near_cpgsite_dict,ranksum_statistic_file)
          #linear_statistic_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Statistic_Data/"+TissueType+"/linear_statistic_pvalue.txt"
          #correlation_linear=get_spatial_correlation(near_cpgsite_dict,linear_statistic_file)
          linear_statistic_file="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Statistic_Data/"+TissueType+"/I_Linear_Statistic_Pvalue_variance.txt"
          correlation_linear=get_spatial_correlation(near_cpgsite_dict,linear_statistic_file)
          result=TissueType+"\t"+str(correlation_ranksum)+"\t"+str(correlation_linear)+"\n"
          print(result)
          return(result)




def main():
          file_dir = "/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"
          print "file_dir=:"+file_dir
          TumorType = os.listdir(file_dir)
          os.chdir("/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/")
          out=open("summary_pearson_variance_I.txt","w")
          out.write("TissueType\tcorrelation_ranksum\tcorrelation_linear\n")
          for TissueType in TumorType:
                    result=run_pipline(TissueType)
                    out.write(result)
          out.close()



if __name__ == "__main__":
          try:
                    main()
          except KeyboardInterrupt:
                    sys.stderr.write("User interrupt me, see you!\n")
                    sys.exit(0)    




