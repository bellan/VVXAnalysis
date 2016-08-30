#! /usr/bin/env python
from optparse import OptionParser
import os
import sys
from os import walk
from collections import Counter
import commands




parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-s", "--string1", dest="String1",
                  default="",
                  help="Some string inside a specific job set like the lxplus machine")

parser.add_option("-t", "--string2", dest="String2",
                  default="",
                  help="Some string inside a specific job set like the lxplus machine")


parser.add_option("-d", "--detail", dest="detail",
                  action="store_true",
                  default=False,
                  help="Print details")


(options, args) = parser.parse_args()

String1 = options.String1
String2 = options.String2
detail  = options.detail


AllJobs =  commands.getstatusoutput('bjobs -a')


print AllJobs[1].count('DONE '),"Done"
print AllJobs[1].count('RUN '),"Run"
print AllJobs[1].count('EXIT '),"Exit"
print AllJobs[1].count('PEND '),"Pend"
print AllJobs[1].count('ZOMBI '),"Zombi"

if detail:
    print AllJobs[1]

list1 = AllJobs[1].split()
lxpluslist = []

for s in list1:
    if s not in lxpluslist:
        if "lxplus" in s:
            lxpluslist.append(s)


print "machines used for jobs"
for j in lxpluslist:
    print j


#print type(AllJobs[1]),AllJobs


if String1 != "":
    RunJobs =  commands.getstatusoutput('bjobs -r')
    print RunJobs[1].count(String1),"Run in ",String1
    
    PendJobs =  commands.getstatusoutput('bjobs -p')
    print PendJobs[1].count(String1),"Pending in ",String1

if String2 != "":
    RunJobs =  commands.getstatusoutput('bjobs -r')
    print RunJobs[1].count(String2),"Run in ",String2
    
    PendJobs =  commands.getstatusoutput('bjobs -p')
    print PendJobs[1].count(String2),"Pending in ",String2



