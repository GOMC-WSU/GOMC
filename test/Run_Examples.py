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

print("Binaries Dict", binaries_dict)

os.chdir("../integration")

confFileName = "in.conf"

Log_Template_file = open("IntegrationTest.log", 'w')
print("opened Log")

listOfTests = []

# traverse root directory, and list directories as dirs and files as files
# traverse root directory, and list directories as dirs and files as files
for root, dirs, files in os.walk("."):
    path = root.split(os.sep)
    print((len(path) - 1) * '---', os.path.basename(root))
    for file in files:
        print(len(path) * '---', file)
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

            if cpuOrGpu+"NVT"+newOrRef in binaries_dict and 'NVT' in path and 'GEMC_NVT' not in path:    
                print("Call GOMC")
                command = binaries_dict[cpuOrGpu+"NVT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NVT"+newOrRef,os.path.basename(root)
                print(binaries_dict[cpuOrGpu+"NVT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NVT"+newOrRef,os.path.basename(root))
                listOfTests.append(command)
            elif cpuOrGpu+"NPT"+newOrRef in binaries_dict and 'NPT' in path and 'GEMC_NPT' not in path:    
                print("Call GOMC")
                print(binaries_dict[cpuOrGpu+"NPT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NPT"+newOrRef,os.path.basename(root))
                command = binaries_dict[cpuOrGpu+"NPT"+newOrRef],os.path.abspath(root),cpuOrGpu+"NPT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GCMC"+newOrRef in binaries_dict and 'GCMC' in path:
                print("Call GOMC")
                print(binaries_dict[cpuOrGpu+"GCMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GCMC"+newOrRef,os.path.basename(root))
                command = binaries_dict[cpuOrGpu+"GCMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GCMC"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GEMC"+newOrRef in binaries_dict and 'GEMC_NVT' in path:
                print("Call GOMC")
                print(binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GEMC_NVT"+newOrRef,os.path.basename(root))
                command = binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GEMC_NVT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
            elif cpuOrGpu+"GEMC"+newOrRef in binaries_dict and 'GEMC_NPT' in path:
                print("Call GOMC")
                print(binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GEMC_NPT"+newOrRef,os.path.basename(root))
                command = binaries_dict[cpuOrGpu+"GEMC"+newOrRef],os.path.abspath(root),cpuOrGpu+"GEMC_NPT"+newOrRef,os.path.basename(root)
                listOfTests.append(command)
            else:
                print(cpuOrGpu+"GEMC"+newOrRef)            
print(listOfTests)
# Create the pandas DataFrame
df = pd.DataFrame(listOfTests, columns = ['PathToBinary', 'PathToExample', 'Binary', 'Example'])
df = df.sort_values(by=['Example'])
print(df)
for index, row in df.iterrows():
    print(index)
    print(row)
"""
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
"""
