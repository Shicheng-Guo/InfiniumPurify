import os

##############



def run_pipline(TissueType):
          target_dir="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/hao/"
          os.chdir(target_dir)
          os.system("mkdir "+TissueType)
          workdir="/mnt/Storage/home/zhengxq/zhengxq/450Purify/NAR_F/Result/linear_method/Statistic_Data/"
          os.system("cp "+workdir+TissueType+"/Ranksum_Statistic_Pvalue.txt "+target_dir+TissueType)
          os.system("cp "+workdir+TissueType+"/Linear_Statistic_Pvalue_variance.txt "+target_dir+TissueType)
          os.chdir(target_dir+TissueType)
          os.system("mv Linear_Statistic_Pvalue_variance.txt Linear_Statistic_Pvalue.txt")


               
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


