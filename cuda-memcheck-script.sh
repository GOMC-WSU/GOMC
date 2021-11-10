#!/bin/bash
#Run script with the executable as the commandline argument
cuda-memcheck --tool memcheck --leak-check full --show-backtrace yes --save memcheck-results%p.log $*
#cuda-memcheck --tool memcheck --leak-check full --show-backtrace yes --save memcheck-results%p.log $*
#cuda-memcheck --tool racecheck --save memcheck-results%p.log $*
#cuda-memcheck --save memcheck-results%p.log $*