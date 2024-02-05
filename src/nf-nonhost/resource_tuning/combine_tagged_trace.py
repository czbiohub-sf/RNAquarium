import glob
from subprocess import call
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import re
from sys import argv


#combine_tagged_trace.py "trace*" zf_54k_read_summary.tsv benchmark_slurm_accounting.txt
matplotlib.use('svg')

trace_files = glob.glob(argv[1])
trace_files.sort()

trace_list = []
for f in trace_files:
    df = pd.read_csv(f, header=0, delimiter="\t")
    trace_list.append(df)
combined = pd.concat(trace_list, axis=0)

jids = combined.native_id.to_string(index=False).split('\n')
get_slurm_acct_script = f"sacct -j {'.batch,'.join(jids)}.batch -o    JobIDRaw,NodeList,State,AveCPU,AveCPUFreq,AveDiskRead,AveDiskWrite,AveRSS,AveVMSize,CPUTimeRAW,ElapsedRaw,MaxDiskRead,MaxDiskWrite,MaxPages,MaxRSS,MaxVMSize,MinCPU,NCPUS,ReqCPUS,ReqMem,TotalCPU --parsable2 > benchmark_slurm_accounting.txt\n"

if len(argv) < 2:
    call(get_slurm_acct_script, shell=True)
    slurm_info = "benchmark_slurm_accounting.txt"
else:
    slurm_info = argv[3]
# with open("get_slurm_acct.sh", "w") as file:
#     file.write("sacct -j ")
#     file.write(f"{'.batch,'.join(jids)}.batch")
#     file.write(" -o    JobIDRaw,NodeList,State,AveCPU,AveCPUFreq,AveDiskRead,AveDiskWrite,AveRSS,AveVMSize,CPUTimeRAW,ElapsedRaw,MaxDiskRead,MaxDiskWrite,MaxPages,MaxRSS,MaxVMSize,MinCPU,NCPUS,ReqCPUS,ReqMem,TotalCPU")
#     file.write(" --parsable2 > benchmark_slurm_accounting.txt\n")


del trace_list
accession_p = re.compile('([SED]RR\\d+)')
process_p = re.compile('(.*) ')
combined.insert(0, 'accession', combined['name'].str.extract(accession_p))
combined.insert(1, 'process', combined['name'].str.extract(process_p))


readstats = pd.read_csv(argv[2], header=0, sep="\t")
readstats.sort_values(by=['reads', 'median'], inplace=True)
readstats['accession'] = readstats['fastq_name'].str.extract(accession_p)
readstats.drop_duplicates(subset='accession', keep='last', inplace=True,
                          ignore_index=True)
full = combined.merge(readstats, on='accession', right_index=False, sort=False)

jobid_p = re.compile('(\\d+)')
nodegroup_p = re.compile('(\\w+-\\w+)-')
slurm_acct = pd.read_csv(argv[3], header=0, sep="|")
slurm_acct['native_id'] = slurm_acct['JobIDRaw'].str.extract(jobid_p)
slurm_acct['native_id'] = pd.to_numeric(slurm_acct['native_id'])
slurm_acct['NodeType'] = slurm_acct['NodeList'].str.extract(nodegroup_p)

full = full.merge(slurm_acct, on='native_id', sort=False)

full = full[full.State == "COMPLETED"]
full.drop_duplicates(subset='name', keep='last', inplace=True)

#columns = None
c = ['accession', 'process', 'native_id', 'name', 'duration', 'realtime', '%cpu',
     'peak_rss', 'peak_vmem', 'rchar', 'wchar', 'median', 'unique', 'reads', 'size']
c = c + ['NodeType', 'NodeList', 'State', 'AveCPU', 'AveCPUFreq', 'AveDiskRead', 'AveDiskWrite', 'AveRSS', 'AveVMSize', 'CPUTimeRAW', 'ElapsedRaw', 'MaxDiskRead', 'MaxDiskWrite', 'MaxPages', 'MaxRSS', 'MaxVMSize', 'MinCPU', 'NCPUS', 'ReqCPUS', 'ReqMem', 'TotalCPU']
full.to_csv("combined_trace_pre_conversion.txt", sep="\t", index=False, columns=c)


def storage_to_bytes(s: str):
    units = {'B': 1, 'KB': 1024, 'MB': 1024*1024, 'GB': 1024*1024*1024,
             'TB': 1024*1024*1024*1024, 'PB': 1024*1024*1024*1024*1024}
    p = re.compile("(\\d+.?(?:\\d+)?) ?(\\w+)")
    m = re.search(p, s)
    if m:
        return pd.to_numeric(m.group(1)) * units[m.group(2)]
    else:
        return np.NaN

def si_to_int(s: str):
    units = {'': 1, 'K': 1024, 'M': 1024*1024, 'G': 1024*1024*1024,
             'T': 1024*1024*1024*1024, 'P': 1024*1024*1024*1024*1024}
    p = re.compile("(\\d+.?(?:\\d+)?) ?(\\w+)")
    m = re.search(p, s)
    if m:
        return pd.to_numeric(m.group(1)) * units[m.group(2)]
    else:
        return np.NaN

def percent_to_float(s: str):
    p = re.compile("(\\d+.?(?:\\d+)?) ?%")
    m = re.search(p, s)
    if m:
        return pd.to_numeric(m.group(1)) / 100
    else:
        return np.NaN
def strip_percent(s: str):
    p = re.compile("(\\d+.?(?:\\d+)?) ?%")
    m = re.search(p, s)
    if m:
        return pd.to_numeric(m.group(1))
    else:
        return np.NaN

def time_to_s(s: str):
    units = {'d': 24*60*60, 'h': 60*60, 'm': 60, 's': 1, '': 1, 'ms': 0.001}
    p = re.compile("([0-9]+\\.?[0-9]*\\w+)")
    p2 = re.compile("([0-9]+\\.?[0-9]*)(\\w+)")
    m = re.findall(p, s)
    sec = 0
    if m:
        for mm in m:
            m2 = re.search(p2, mm)
            if m2:
                sec = sec + pd.to_numeric(m2.group(1)) * units[m2.group(2)]
        return sec
    else:
        return np.NaN

def colontime_to_s(s:str):
    units = {'d': 24*60*60, 'h': 60*60, 'm': 60, 's': 1, '': 1, 'ms': 0.001}
    p = re.compile("(?P<d>\\d+)?-?(?P<h>\\d\\d):(?P<m>\\d\\d):(?P<s>\\d\\d)")
    m = re.search(p, s)
    sec = 0
    if m:
        if m.group('d'):
            sec += pd.to_numeric(m.group('d')) * units['d']
        sec += pd.to_numeric(m.group('h')) * units['h']
        sec += pd.to_numeric(m.group('m')) * units['m']
        sec += pd.to_numeric(m.group('s')) * units['s']
        return sec
    else:
        return np.NaN

    
full['duration'] = full['duration'].map(time_to_s)
full['realtime'] = full['realtime'].map(time_to_s)
full['%cpu'] = full['%cpu'].map(strip_percent)
full['peak_rss'] = full['peak_rss'].map(storage_to_bytes)
full['peak_vmem'] = full['peak_vmem'].map(storage_to_bytes)
full['rchar'] = full['rchar'].map(storage_to_bytes)
full['wchar'] = full['wchar'].map(storage_to_bytes)
full['AveCPU'] = full['AveCPU'].map(colontime_to_s)
full['AveCPUFreq'] = full['AveCPUFreq'].map(si_to_int)
full['AveDiskRead'] = full['AveDiskRead'].map(si_to_int)
full['AveDiskWrite'] = full['AveDiskWrite'].map(si_to_int)
full['AveRSS'] = full['AveRSS'].map(si_to_int)
full['AveVMSize'] = full['AveVMSize'].map(si_to_int)
full['MaxDiskRead'] = full['MaxDiskRead'].map(si_to_int)
full['MaxDiskWrite'] = full['MaxDiskWrite'].map(si_to_int)
full['MaxRSS'] = full['MaxRSS'].map(si_to_int)
full['MaxVMSize'] = full['MaxVMSize'].map(si_to_int)
full['MinCPU'] = full['MinCPU'].map(colontime_to_s)
full['TotalCPU'] = full['TotalCPU'].map(colontime_to_s)

#columns = None
c = ['accession', 'process', 'native_id', 'name', 'duration', 'realtime', '%cpu',
     'peak_rss', 'peak_vmem', 'rchar', 'wchar', 'median', 'unique', 'reads', 'size']
c = c + ['NodeType', 'NodeList', 'State', 'AveCPU', 'AveCPUFreq', 'AveDiskRead', 'AveDiskWrite', 'AveRSS', 'AveVMSize', 'CPUTimeRAW', 'ElapsedRaw', 'MaxDiskRead', 'MaxDiskWrite', 'MaxPages', 'MaxRSS', 'MaxVMSize', 'MinCPU', 'NCPUS', 'ReqCPUS', 'ReqMem', 'TotalCPU']
full.to_csv("combined_trace.txt", sep="\t", index=False, columns=c)