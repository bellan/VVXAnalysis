#! /usr/bin/env python

from __future__ import print_function
import sys
import os

def reformat(input_file_name):
    
    input_file = open(input_file_name, 'r')
    output_file = open(output_file_name, 'w')


    counter=0

    output_file.write("process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange( *(\n")
    
    for line in input_file:
        counter=counter+1  
        

        preblock=line.split()[0]
        block=preblock.strip('{"":')
        subs=line.count("[")-1
        #print(line)
         
        for i in range(1,subs+1):
            preindex1=line.split("[")[2]
            index1=preindex1.split(",")[0]
            prepreindex2=line.split("]")[0]
            preindex2=prepreindex2.split(",")[1]
            index2=preindex2.strip(" ")
            #print(index2)
            output_file.write("    '"+block+":"+index1+"-"+block+":"+index2+"',\n") 

    output_file.write("))")
    
    input_file.close()     
    output_file.close()
 

input_file_name = sys.argv[1]
output_file_name = input_file_name.replace(".py","2.py")


reformat(input_file_name)
