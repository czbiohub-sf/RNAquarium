import os.path
import csv
import argparse
import sys
import linecache
import datetime
import gc
import math
import errno, os, stat, shutil
import pickle

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script parses the .summary files and accumulate into one file')
    parser.add_argument('--dir', default="", type=str, help='path to the folder containing the summarized files', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())
indir = config["dir"].rstrip('/').rstrip('\\')


#####################
##      main       ##
#####################
def main():
    try:
        #check input
        if (config["dir"] is None or config["dir"] == ""):
            log.error(f"please specify --dir")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()

        #load run-accession-to-sample mapping
        with open("/hpc/projects/balla_group/sra_experiments/SRA_metadata/ZF_SRA_accession_list.txt_InfoDict.pkl", 'rb') as handle:
            run_info = pickle.load(handle)
            #print(acc_info)

        #read sample-attr-to-category mapping
        sample_info = dict()
        with open("/hpc/projects/balla_group/sra_experiments/SRA_metadata/attr_by_sample_tissues.txt_MAPPEDWITH_tissue_cat_def.txt.tab", "r") as handle:
            for line in handle:
                fields = line.split("\t")
                Merged_attributes = fields[0].replace(",",";")
                cate = fields[1].replace(",",";")
                sample_ID = fields[2]
                sample_info[sample_ID] = [Merged_attributes,cate]

        #read mapping count data
        counts = dict()
        with open("/hpc/projects/balla_group/sra_experiments/all_zebrafish_RNAseq/read_counts/countsSummaryTable.csv", "r") as handle:
            next(handle)#skip header
            for line in handle:
                fields = line.split(",")
                accession = fields[0]
                mapped_count = fields[1]
                unmapped_count = sum([int(i) for i in fields[2:15]])
                counts[accession] = [mapped_count, unmapped_count]

        outfh = open(f"{indir}.parsed.csv", "w")
        #write header
        outfh.write(f"SRA_run_accession,SRA_sample_accession,tissue category,merged_attr,mapped_count,unmapped_count,readlength,famcvg,fam,score,pctid,depth,aln,glb,length,top,topscore,toplen,topname,seq_divider,seq count,seqcvg_1,seq_1,score_1,pctid_1,depth_1,aln_1,glb_1,length_1,family_1,name_1,seq_divider,seqcvg_2,seq_2,score_2,pctid_2,depth_2,aln_2,glb_2,length_2,family_2,name_2,seq_divider,seqcvg_3,seq_3,score_3,pctid_3,depth_3,aln_3,glb_3,length_3,family_3,name_3,seq_divider,seqcvg_4,seq_4,score_4,pctid_4,depth_4,aln_4,glb_4,length_4,family_4,name_4,,\n")
        file_count = 0
        for file in os.listdir(indir):
            if not file.endswith(".summary"):
                continue

            SRA_ID = file.split(".summary")[0]

            with open(os.path.join(indir,file), "r") as infh:

                match_flag = 0

                family_res = ""
                gene_res = ""
                gene_count = 0

                for line in infh:

                    if line.startswith("readlength"):
                        readlength = line.split("=")[1].split(";")[0]

                    if line.startswith("famcvg"):
                        fields = line.split(";")
                        for key, val in fields.items():
                            fields[key] = val.split("=")[1]
                        famcvg, fam, score, pctid, depth, aln, glb, length, top, topscore, toplen, topname = fields.values()
                        # famcvg = fields[0].split("=")[1]
                        # fam = fields[1].split("=")[1]
                        # score = fields[2].split("=")[1]
                        # pctid = fields[3].split("=")[1]
                        # depth = fields[4].split("=")[1]
                        # aln = fields[5].split("=")[1]
                        # glb = fields[6].split("=")[1]
                        # length = fields[7].split("=")[1]
                        # top = fields[8].split("=")[1]
                        # topscore = fields[9].split("=")[1]
                        # toplen = fields[10].split("=")[1]
                        # topname = fields[11].split("=")[1]
                        match_flag = 1
                        family_res = f"{famcvg},{fam},{score},{pctid},{depth},{aln},{glb},{length},{top},{topscore},{toplen},{topname},"

                    if line.startswith("seqcvg"): #TODO in future there maybe different seqs in each family
                        fields = line.split(";")
                        for key, val in fields.items():
                            fields[key] = val.split("=")[1]
                        Seqcvg, Seq, Score, Pctid, Depth, Aln, Glb, Length, Family, Name = fields.values()
                        # Seqcvg = fields[0].split("=")[1]
                        # Seq = fields[1].split("=")[1]
                        # Score = fields[2].split("=")[1]
                        # Pctid = fields[3].split("=")[1]
                        # Depth = fields[4].split("=")[1]
                        # Aln = fields[5].split("=")[1]
                        # Glb = fields[6].split("=")[1]
                        # Length = fields[7].split("=")[1]
                        # Family = fields[8].split("=")[1]
                        # Name = fields[9].split("=")[1]
                        match_flag = 1
                        gene_count+=1
                        gene_res = gene_res + f"{Seqcvg},{Seq},{Score},{Pctid},{Depth},{Aln},{Glb},{Length},{Family},{Name},|,"



                if match_flag == 1:
                    row = f"{readlength},{family_res},{gene_count},{gene_res}\n"
                else:
                    row = f"{readlength},,,,,,,,,,,,,,,,,,,,,,\n"

                #ANNOTATION
                SAMPLE_ID = run_info[SRA_ID]["Sample"]
                Merged_attr = sample_info[SAMPLE_ID][0].replace(",", ";")
                categories = sample_info[SAMPLE_ID][1].replace(",", ";")

                #counts
                mapped_count=counts[SRA_ID][0]
                unmapped_count=counts[SRA_ID][1]

                row = f"{SRA_ID},{SAMPLE_ID},{categories},{Merged_attr},{mapped_count},{unmapped_count}," + row 


                outfh.write(row)

                match_flag = 0 # reset
                file_count +=1

        outfh.close()

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        print(f"finished in {elapsed_min:.2f} min, processed {file_count} summary files")

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

def handleRemoveReadonly(func, path, exc):
  excvalue = exc[1]
  if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
      os.chmod(path, stat.S_IRWXU| stat.S_IRWXG| stat.S_IRWXO) # 0777
      func(path)
  else:
      raise

if __name__ == "__main__": main()

