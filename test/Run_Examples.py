#!/usr/bin/python

def substring_after(s, delim):
    return s.partition(delim)[2]

import os
import glob
import sys
import pandas as pd
import time
import datetime
import subprocess
from subprocess import Popen, PIPE, STDOUT
CPU_binaries_dict = {}
GPU_binaries_dict = {}

os.chdir("new_feature_binaries")
pathToDir = os.getcwd()
CPU_binaries_new_feature = sorted(glob.glob('*_CPU_*'), key=os.path.getmtime)
GPU_binaries_new_feature = sorted(glob.glob('*_GPU_*'), key=os.path.getmtime)

for binary in CPU_binaries_new_feature:
    CPU_binaries_dict[substring_after(binary, "GOMC_CPU_")+"_new"] = os.path.join(pathToDir, binary)

for binary in GPU_binaries_new_feature:
    GPU_binaries_dict[substring_after(binary, "GOMC_GPU_")+"_new"] = os.path.join(pathToDir, binary)

os.chdir("../ref_binaries")
pathToDir = os.getcwd()
CPU_binaries_ref = sorted(glob.glob('*_CPU_*'), key=os.path.getmtime)
GPU_binaries_ref = sorted(glob.glob('*_GPU_*'), key=os.path.getmtime)

for binary in CPU_binaries_ref:
    CPU_binaries_dict[substring_after(binary, "GOMC_CPU_")+"_ref"] = os.path.join(pathToDir, binary)

for binary in GPU_binaries_ref:
    GPU_binaries_dict[substring_after(binary, "GOMC_GPU_")+"_ref"] = os.path.join(pathToDir, binary)


print("CPU Binaries Dict", CPU_binaries_dict)
print("GPU Binaries Dict", GPU_binaries_dict)

os.chdir("../integration")

confFileName = "in.conf"

Log_Template_file = open("IntegrationTest.log", 'w')
print("opened Log")

listOfTests = []

# traverse root directory, and list directories as dirs and files as files
# traverse root directory, and list directories as dirs and files as files
for root, dirs, files in os.walk("."):
    path = root.split(os.sep)
#    print((len(path) - 1) * '---', os.path.basename(root))
    for file in files:
#        print(len(path) * '---', file)
        if file==confFileName:
            newOrRef = ""
            if "new_feature_cpu" in path:
                newOrRef = "_new"
            elif "ref_cpu" in path:
                newOrRef = "_ref"

            run = False
            if "NVT"+newOrRef in CPU_binaries_dict and 'NVT' in path and 'GEMC_NVT' not in path:    
                print("Call GOMC")
                command = CPU_binaries_dict["NVT"+newOrRef],os.path.abspath(root),"NVT"+newOrRef,os.path.basename(root)
                print(CPU_binaries_dict["NVT"+newOrRef],os.path.abspath(root),"NVT"+newOrRef,os.path.basename(root))
                listOfTests.append(command)
                run = True
            elif "NPT"+newOrRef in CPU_binaries_dict and 'NPT' in path and 'GEMC_NPT' not in path:    
                print("Call GOMC")
                print(CPU_binaries_dict["NPT"+newOrRef],os.path.abspath(root),"NPT"+newOrRef,os.path.basename(root))
                command = CPU_binaries_dict["NPT"+newOrRef],os.path.abspath(root),"NPT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
                run = True
            elif "GCMC"+newOrRef in CPU_binaries_dict and 'GCMC' in path:
                print("Call GOMC")
                print(CPU_binaries_dict["GCMC"+newOrRef],os.path.abspath(root),"GCMC"+newOrRef,os.path.basename(root))
                command = CPU_binaries_dict["GCMC"+newOrRef],os.path.abspath(root),"GCMC"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
                run = True
            elif "GEMC"+newOrRef in CPU_binaries_dict and 'GEMC_NVT' in path:
                print("Call GOMC")
                print(CPU_binaries_dict["GEMC"+newOrRef],os.path.abspath(root),"GEMC_NVT"+newOrRef,os.path.basename(root))
                command = CPU_binaries_dict["GEMC"+newOrRef],os.path.abspath(root),"GEMC_NVT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
                run = True
            elif "GEMC"+newOrRef in CPU_binaries_dict and 'GEMC_NPT' in path:
                print("Call GOMC")
                print(CPU_binaries_dict["GEMC"+newOrRef],os.path.abspath(root),"GEMC_NPT"+newOrRef,os.path.basename(root))
                command = CPU_binaries_dict["GEMC"+newOrRef],os.path.abspath(root),"GEMC_NPT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
                run = True
            else:
                run = False
            
print(listOfTests)
# Create the pandas DataFrame
df = pd.DataFrame(listOfTests, columns = ['PathToBinary', 'PathToExample', 'Binary', 'Example'])
df = df.sort_values(by=['Example'])
print(df)
for index, row in df.iterrows():
    os.chdir(row['PathToExample'])
    write_log_data = "Changing directory to {}\n".format(row['PathToExample'])
    Log_Template_file.write(str(write_log_data))
    command = (row['PathToBinary'] + " in.conf > out.log")
    write_log_data = "Issuing command: {}\n".format(command)
    Log_Template_file.write(str(write_log_data))
    start = time.time()
    exec_GOMC_run_command = subprocess.Popen(command, shell=True, stderr=STDOUT)
    write_log_data = "Waiting for GOMC Example {} {} to finish.\n".format(row['Binary'],row['Example'])
    print(str(write_log_data))
    Log_Template_file.write(str(write_log_data))
    GOMC_pid_status = os.waitpid(exec_GOMC_run_command.pid, os.WSTOPPED)  # pauses python until box 0 sim done    
    end = time.time()
    write_log_data = "Elapsed time: {}.\n".format(datetime.timedelta(seconds=end-start))
    Log_Template_file.write(str(write_log_data))
    print(str(write_log_data))
    write_log_data = "The GOMC Example {} {} has finished.\n".format(row['Binary'],row['Example'])
    print(str(write_log_data))
    Log_Template_file.write(str(write_log_data))
    Log_Template_file.flush()
