#!/usr/bin/python

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

def substring_after(s, delim):
    return s.partition(delim)[2]

import os
import glob
import pandas as pd
import time
import datetime
import subprocess
from subprocess import STDOUT
from filecmp import cmp
binaries_dict = {}
GPU_binaries_dict = {}

os.chdir("new_binaries")
pathToDir = os.getcwd()
binaries_cpu_new = sorted(glob.glob('GOMC_CPU_*'), key=os.path.getmtime)
binaries_gpu_new = sorted(glob.glob('GOMC_GPU_*'), key=os.path.getmtime)
for binary in binaries_cpu_new:
    binaries_dict[substring_after(binary, "GOMC_")+"_new"] = os.path.join(pathToDir, binary)
for binary in binaries_gpu_new:
    binaries_dict[substring_after(binary, "GOMC_")+"_new"] = os.path.join(pathToDir, binary)

os.chdir("../ref_binaries")
pathToDir = os.getcwd()
binaries_cpu_ref = sorted(glob.glob('GOMC_CPU_*'), key=os.path.getmtime)
binaries_gpu_ref = sorted(glob.glob('GOMC_GPU_*'), key=os.path.getmtime)
for binary in binaries_cpu_ref:
    binaries_dict[substring_after(binary, "GOMC_")+"_ref"] = os.path.join(pathToDir, binary)
for binary in binaries_gpu_ref:
    binaries_dict[substring_after(binary, "GOMC_")+"_ref"] = os.path.join(pathToDir, binary)


os.chdir("../integration")

confFileName = "in.conf"

Log_Template_file = open("IntegrationTest.log", 'w')
Log_Template_file.write("Binaries Dict") 
Log_Template_file.write(str(binaries_dict))
Log_Template_file.write("opened Log")

listOfTests = []

# traverse root directory, and list directories as dirs and files as files
# traverse root directory, and list directories as dirs and files as files
for root, dirs, files in os.walk("."):
    path = root.split(os.sep)
    for file in files:
        if file==confFileName:
            newOrRef = ""
            cpuOrGpu = ""
            if "new_cpu" in path:
                newOrRef = "_new"
                cpuOrGpu = "CPU_"
            elif "ref_cpu" in path:
                newOrRef = "_ref"
                cpuOrGpu = "CPU_"
            elif "new_gpu" in path:
                newOrRef = "_new"
                cpuOrGpu = "GPU_"
            elif "ref_gpu" in path:
                newOrRef = "_ref"
                cpuOrGpu = "GPU_"

            if cpuOrGpu+"NVT"+newOrRef in binaries_dict and 'NVT' in path and 'NVT_GEMC' not in path:
                command = binaries_dict[cpuOrGpu+"NVT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NVT"+newOrRef,cpuOrGpu,newOrRef,"NVT_"+os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"NPT"+newOrRef in binaries_dict and 'NPT' in path and 'NPT_GEMC' not in path:
                command = binaries_dict[cpuOrGpu+"NPT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NPT"+newOrRef,cpuOrGpu,newOrRef,"NPT_"+os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GCMC"+newOrRef in binaries_dict and 'GCMC' in path:
                command = binaries_dict[cpuOrGpu+"GCMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GCMC"+newOrRef,cpuOrGpu,newOrRef,"GCMC_"+os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GEMC"+newOrRef in binaries_dict and 'NVT_GEMC' in path:
                command = binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"NVT_GEMC"+newOrRef,cpuOrGpu,newOrRef,"NVT_GEMC_"+os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GEMC"+newOrRef in binaries_dict and 'NPT_GEMC' in path:
                command = binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"NPT_GEMC"+newOrRef,cpuOrGpu,newOrRef,"NPT_GEMC_"+os.path.basename(root)
                listOfTests.append(command)
Log_Template_file.write("GOMC Call Commands:\n")
Log_Template_file.write(str(listOfTests))
# Create the pandas DataFrame
colNames = ['PathToBinary', 'PathToExample', 'Binary', 'CPU_or_GPU','New_or_Ref','Example']
df = pd.DataFrame(listOfTests, columns = colNames)
df = df.sort_values(by=['Example'])
Log_Template_file.write(str(df))
grouped = df.groupby(['Example'])
all_examples = df['Example'].unique()
CPU_v_GPU_global = True
New_v_Ref_global = True
cross_bool_global = True
CPU_v_GPU_exists_global = False
New_v_Ref_exists_global = False
cross_exists_global = False
for example in all_examples:
    ex_df = grouped.get_group(example)
    for index, row in ex_df.iterrows():
        os.chdir(row['PathToExample'])
        write_log_data = "Changing directory to {}\n".format(row['PathToExample'])
        Log_Template_file.write(str(write_log_data))
        command = (row['PathToBinary'] + " in.conf > out.log")
        write_log_data = "Issuing command: {}\n".format(command)
        Log_Template_file.write(str(write_log_data))
        start = time.time()
        exec_GOMC_run_command = subprocess.Popen(command, shell=True, stderr=STDOUT)
        write_log_data = "Waiting for GOMC Example {} {} to finish.\n".format(row['Binary'],row['Example'])
        Log_Template_file.write(str(write_log_data))
        GOMC_pid_status = os.waitpid(exec_GOMC_run_command.pid, os.WSTOPPED)  # pauses python until box 0 sim done
        end = time.time()
        write_log_data = "Elapsed time: {}.\n".format(datetime.timedelta(seconds=end-start))
        Log_Template_file.write(str(write_log_data))
        write_log_data = "The GOMC Example {} {} has finished.\n".format(row['Binary'],row['Example'])
        Log_Template_file.write(str(write_log_data))
        Log_Template_file.flush()

    # Create a list of the PDB files in this example
    full_path_pdb_files = sorted(glob.glob(os.path.join(row['PathToExample'],'*.pdb')), key=os.path.getmtime)
    just_file_names = []
    for path in full_path_pdb_files:
        just_file_names.append(os.path.basename(path))
    Log_Template_file.write(str(just_file_names))
    cross = ex_df.merge(ex_df, on=['Example'],how='outer')
    Log_Template_file.write(str('---', ex_df['Example'].iloc[0]))
    CPU_v_GPU = True
    New_v_Ref = True
    cross_bool = True
    CPU_v_GPU_exists = False
    New_v_Ref_exists = False
    cross_exists = False
    for pdb_file in just_file_names:
        Log_Template_file.write(str(2 * '---', pdb_file))
        my_tuples = []
        for index, row in cross.iterrows():
            f1 = os.path.join(row['PathToExample_x'],pdb_file)
            f2 = os.path.join(row['PathToExample_y'],pdb_file)
            if ((row['CPU_or_GPU_x'] != row['CPU_or_GPU_y']) and (row['New_or_Ref_x'] == row['New_or_Ref_y'])):
                CPU_v_GPU_exists = True
                CPU_v_GPU_exists_global = True
                result = cmp(f1, f2, shallow=False)
                CPU_v_GPU = CPU_v_GPU and result
                CPU_v_GPU_global = CPU_v_GPU_global and result
            elif ((row['CPU_or_GPU_x'] == row['CPU_or_GPU_y']) and (row['New_or_Ref_x'] != row['New_or_Ref_y'])):
                New_v_Ref_exists = True
                New_v_Ref_exists_global = True
                result = cmp(f1, f2, shallow=False)
                New_v_Ref = New_v_Ref and result
                New_v_Ref_global = New_v_Ref_global and result
            elif ((row['CPU_or_GPU_x'] != row['CPU_or_GPU_y']) and (row['New_or_Ref_x'] != row['New_or_Ref_y'])):
                cross_exists = True
                cross_exists_global = True
                result = cmp(f1, f2, shallow=False)
                cross_bool = cross_bool and result
                cross_bool_global = cross_bool_global and result
    if(CPU_v_GPU_exists):
        if(CPU_v_GPU):
            Log_Template_file.write(str((3 * '---')+"CPU_v_GPU: "+ OKGREEN + "PASS" + ENDC))
        else:
            Log_Template_file.write(str((3 * '---')+"CPU_v_GPU: "+ FAIL + "FAIL" + ENDC))
    if(New_v_Ref_exists):
        if(New_v_Ref):
            Log_Template_file.write(str((3 * '---')+"New vs Ref: "+ OKGREEN + "PASS" + ENDC))
        else:
            Log_Template_file.write(str((3 * '---')+"New vs Ref: "+ FAIL + "FAIL" + ENDC))
    if(cross_exists):
        if(cross_bool):
            Log_Template_file.write(str((3 * '---')+"CPU vs GPU X New vs Ref: "+ OKGREEN + "PASS" + ENDC))
        else:
            Log_Template_file.write(str((3 * '---')+"CPU vs GPU X New vs Ref: "+ FAIL + "FAIL" + ENDC))


if(CPU_v_GPU_exists_global):
    if(CPU_v_GPU_global):
        Log_Template_file.write(str("CPU_v_GPU Global: "+ OKGREEN + "PASS" + ENDC))
    else:
        Log_Template_file.write(str("CPU_v_GPU Global: "+ FAIL + "FAIL" + ENDC))
if(New_v_Ref_exists_global):
    if(New_v_Ref_global):
        Log_Template_file.write(str("New vs Ref Global: "+ OKGREEN + "PASS" + ENDC))
    else:
        Log_Template_file.write(str("New vs Ref Global: "+ FAIL + "FAIL" + ENDC))
if(cross_exists_global):
    if(cross_bool_global):
        Log_Template_file.write(str("CPU vs GPU X New vs Ref Global: "+ OKGREEN + "PASS" + ENDC))
    else:
        Log_Template_file.write(str("CPU vs GPU X New vs Ref Global: "+ FAIL + "FAIL" + ENDC))
