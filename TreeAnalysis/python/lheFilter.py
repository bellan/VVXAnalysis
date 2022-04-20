#! /usr/bin/env python

from __future__ import print_function
import sys
import os

def filter(input_file_name):
    
    input_file = open(input_file_name, 'r')
    output_file = open(output_file_name, 'w')

    event_counter = 0        # Counts the amount of event targeted in the ZZ->4l diffractive analysis
    #ll_event_counter = 0     # Counts the amount of events with 4 charged leptons

    linelist=[]
    in_ev   = 0 # To know if we're inside an event
    in_ev_1 = 0 # The first line after <event> is information so we must skip that as well
    header  = True

    had_num = 0 # Counts the amount of hadrons (excluding protons) in the event
    lep_num = 0 # Counts the amount of leptons in the event 
    charg_lep_num = 0 # Counts the amount of charged leptons in the event
    neut_num = 0 # Counts the amount of neutrinos in the event

    
    for line in input_file:

        if line.startswith("#"): continue

        if in_ev_1 == 1:
            in_ev_1 = 0
            in_ev = 1
            linelist.append(line)
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
            if had_num == 0 and charg_lep_num>=4:
                event_counter = event_counter + 1
                #if charg_lep_num >= 4:
                #    ll_event_counter = ll_event_counter + 1
                output_file.write("<event>\n")
                for line1 in linelist: 
                    output_file.write(line1)
                output_file.write("</event>\n")
            had_num = 0
            lep_num = 0
            charg_lep_num = 0
            neut_num = 0
            linelist=[]
            continue
            
        if (line.startswith("<") or line.startswith("  <") or line.startswith("   <")):
            continue
    
        if in_ev == 1:
        
            if (abs(int(line.split()[0]))>23 and int(line.split()[0])!=2212):
                had_num = had_num + 1
            elif (abs(int(line.split()[0]))>=11 and abs(int(line.split()[0]))<=16):
                linelist.append(line)
                lep_num = lep_num + 1
                if ((abs(int(line.split()[0]))==11 or abs(int(line.split()[0]))==13) or abs(int(line.split()[0]))==15):
                    charg_lep_num = charg_lep_num + 1
                else:
                    neut_num= neut_num + 1
            else:
                linelist.append(line)
             
    output_file.write("</LesHouchesEvents>")
    print(event_counter)
    #print(ll_event_counter)
                      
    input_file.close()     
    output_file.close()
 

input_file_name = sys.argv[1]
output_file_name = input_file_name.replace(".lhe","filtered.lhe")


filter(input_file_name)
