#! /usr/bin/env python

from __future__ import print_function
import sys
import os

def filter(input_file_name):
    
    input_file = open(input_file_name, 'r')
    output_file = open(output_file_name, 'w')

    linelist=[]
    in_ev   = 0 # To know if we're inside an event
    in_ev_1 = 0 # The first line after <event> is information so we must skip that as well
    had_num = 0 # Counts the amount of hadrons (excluding protons) in the event
    lep_num = 0 # Counts the amount of leptons in the event
    header  = True

    for line in input_file:

        if line.startswith("#"): continue

        if in_ev_1 == 1:
            in_ev_1 = 0
            in_ev = 1
            continue
    
        if line.startswith("<event>"):
            header = False
            in_ev_1 = 1
            continue
  
        if header: 
            output_file.write(line)
            continue
    
        if in_ev == 1 and line.startswith("</event>"):
            in_ev = 0
            if had_num == 0 and lep_num>3:
                output_file.write("<event>\n")
                for line1 in linelist:     
                    #print(line1)
                    output_file.write(line1)
                output_file.write("</event>\n")
            had_num = 0
            lep_num = 0
            linelist=[]
            #print(linelist)
            continue
            
        if (line.startswith("<") or line.startswith("  <") or line.startswith("   <")):
            continue
    
        if in_ev == 1:
        
            if (abs(int(line.split()[0]))>23 and int(line.split()[0])!=2212):
                had_num = had_num + 1
            elif (abs(int(line.split()[0]))>10 and abs(int(line.split()[0]))<17):
                linelist.append(line)
                lep_num = lep_num + 1
            else:
                linelist.append(line)
             
            #print(line)
    output_file.write("</LesHouchesEvents>")
                      
                      
    input_file.close()     
    temp_file.close()
    output_file.close()
 

input_file_name = sys.argv[1]
temp_file_name = input_file_name.replace(".lhe","temp.lhe")
output_file_name = input_file_name.replace(".lhe","filtered.lhe")


filter(input_file_name)
