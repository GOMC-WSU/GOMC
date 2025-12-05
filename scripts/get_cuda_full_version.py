import os

nvcc_version_cmd = 'nvcc -V > output.txt'
os.system(nvcc_version_cmd)

with open('output.txt') as f:
    lines = f.readlines()
    for line in lines:
        if ", release" in line:
            start = line.index(', release') + 10
            end = line.index(',', start)
            result = line[start:end]
            print(result)
            os.remove("output.txt")
            quit()
