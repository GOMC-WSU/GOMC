# Instructions for running integration tests

### Clones GOMC_Examples repo and creates 4 copies 
### CPU_Ref, GPU_Ref, CPU_New, GPU_Ref
### Then is will checkout the most likely branch to merge into
### If you are on development when you run the script,
### the this will be the main branch.
###
### If you are a new feature branch, 
### this will be the development branch
###
### The current branch is compiled and the binaries put into new_binaries 
### The ref branch is compiled and the binaries put into ref_binaries 

$ bash Setup_Examples.sh 

### Then the actual examples are run in serial

$ python Run_Examples.py

### The results are in integration/IntegrationTest.log
### To check for pass/fail in color coded output

$ cat integration/IntegrationTest.log
