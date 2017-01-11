import os,numpy,scipy,random,math
from scipy import stats


##########
# function
##########

def get_tissue_normal_quantity(data_dir):
          fileName=os.listdir(data_dir)
          NORMAL=[]
          for fn in fileName:
                    sample_name=fn.split('.')[5]
                    index=sample_name.split('-')[3][0:2]
                    if index=="11":
                              NORMAL.append(sample_name)
          return(len(NORMAL))


def get_tumor_normal_dict(data_dir):
          ## input is a 450k data directory, output is a Methylation Dict.
          fileNames = os.listdir(data_dir)
          TUMOR_dict={}
          NORMAL_dict={}
          for fn in fileNames:
                    #print(fn+"...")
                    TCGA_name=fn.split('.')[5]
                    index=int(TCGA_name.split('-')[3][0:2])
                    if index < 10:
                              print (TCGA_name)
                              TUMOR_dict[TCGA_name]={}
                              DAT = open(data_dir+fn,"r")
                              DAT.readline()
                              DAT.readline()
                              for line in DAT:
                                        arr = line.strip().split("\t")
                                        cpg = str(arr[0])
                                        mRatio = arr[1]
                                        TUMOR_dict[TCGA_name][cpg]=mRatio
                    else:
                              print(TCGA_name)
                              NORMAL_dict[TCGA_name]={}
                              DAT = open(data_dir+fn,"r")
                              DAT.readline()
                              DAT.readline()
                              for line in DAT:
                                        arr = line.strip().split("\t")
                                        cpg = str(arr[0])
                                        mRatio = arr[1]
                                        NORMAL_dict[TCGA_name][cpg]=mRatio
          return(TUMOR_dict,NORMAL_dict)
          
"""
def get_tumor_normal_txt(filename,transform_filename,TUMOR_NORMAL_dict,REF_CPG):
          out1=open(filename+".txt","w")
          out1.write("Sample"+"\t"+"\t".join(REF_CPG)+"\n")
          out2=open(transform_filename+".txt","w")
          out2.write("Sample"+"\t"+"\t".join(REF_CPG)+"\n")
          SAMPLE=TUMOR_NORMAL_dict.keys()
          for sample in SAMPLE:
                    sample_betavalue_dict=TUMOR_NORMAL_dict[sample]
                    ALL=[]
                    TRANSFORM_ALL=[]
                    for cpgsite in REF_CPG:
                              bata_value=sample_betavalue_dict[cpgsite]
                              if bata_value=="NA":
                                        value="NA"
                              else:
                                        value=str(bata_value)
                                        transform_value=str(math.asin(2*float(bata_value)-1))
                              ALL.append(value)
                              TRANSFORM_ALL.append(transform_value)
                    line=sample+"\t"+"\t".join(ALL)
                    line2=sample+"\t"+"\t".join(TRANSFORM_ALL)
                    out1.write(line+"\n")
                    out2.write(line2+"\n")
          out1.close()
          out2.close()


"""
def get_tumor_normal_txt(filename,TUMOR_NORMAL_dict,REF_CPG):
          out1=open(filename+".txt","w")
          out1.write("Sample"+"\t"+"\t".join(REF_CPG)+"\n")
          SAMPLE=TUMOR_NORMAL_dict.keys()
          for sample in SAMPLE:
                    sample_betavalue_dict=TUMOR_NORMAL_dict[sample]
                    ALL=[]
                    for cpgsite in REF_CPG:
                              bata_value=str(sample_betavalue_dict[cpgsite])
                              ALL.append(bata_value)
                    line=sample+"\t"+"\t".join(ALL)
                    out1.write(line+"\n")
          out1.close()


def get_transposition(filename):
          os.system("awk '{for(i=1;i<=NF;i++){a[i]=a[i]_FS$i};_FS=FS}END{for(i=1;i<=NF;i++){print a[i]>FILENAME}}' "+filename)


############################################################################

       
def run_pipline(TissueType):
          m450k_data_dir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/TCGA/"+TissueType+"/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
          (TUMOR_dict,NORMAL_dict)=get_tumor_normal_dict(m450k_data_dir)

          REF_CPG=TUMOR_dict[TUMOR_dict.keys()[0]].keys()
          workdir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/Data_Matrix/"
          os.chdir(workdir)
          os.system("mkdir "+TissueType)
          os.chdir(workdir+TissueType)
          get_tumor_normal_txt("Tumor_initial",TUMOR_dict,REF_CPG)
          get_tumor_normal_txt("Normal_initial",NORMAL_dict,REF_CPG)
          
          get_transposition("Tumor_initial.txt")
          get_transposition("Normal_initial.txt")


         
"""          
def main():
          file_dir = "/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/Infinium_mode/"
          print "file_dir=:"+file_dir
          TumorType = os.listdir(file_dir)
          for TissueType in TumorType:
                    m450k_data_dir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/TCGA/"+TissueType+"/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
                    print(TissueType)
                    run_pipline(TissueType)




          
        
def run_pipline(TissueType):
          m450k_data_dir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/TCGA/"+TissueType+"/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
          (TUMOR_dict,NORMAL_dict)=get_tumor_normal_dict(m450k_data_dir)

          REF_CPG=TUMOR_dict[TUMOR_dict.keys()[0]].keys()
          workdir="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Data/"
          os.chdir(workdir)
          os.system("mkdir "+TissueType)
          os.chdir(workdir+TissueType)
          get_tumor_normal_txt("Tumor_initial","Tumor_tansform",TUMOR_dict,REF_CPG)
          get_tumor_normal_txt("Normal_initial","Normal_tansform",NORMAL_dict,REF_CPG)
          
          get_transposition("Tumor_initial.txt")
          get_transposition("Normal_initial.txt")
          get_transposition("Tumor_tansform.txt")
          get_transposition("Normal_tansform.txt")

"""         
          
def main():
          file_dir = "/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/Infinium_mode/"
          print "file_dir=:"+file_dir
          TumorType = os.listdir(file_dir)
          
          for TissueType in TumorType:
                    m450k_data_dir="/mnt/Storage/home/zhengxq/zhengxq/dna_methylation/TCGA/"+TissueType+"/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
                    NORMAL_quantile=get_tissue_normal_quantity(m450k_data_dir)
                    if NORMAL_quantile>20:
                              print(TissueType)
                              run_pipline(TissueType)
                    else:
                              continue
          



if __name__ == "__main__":
          try:
                    main()
          except KeyboardInterrupt:
                    sys.stderr.write("User interrupt me, see you!\n")
                    sys.exit(0)    

