# bash script

#####################################################################
# Automates scram, batch_condor.py and cmsRun of a test sample      #
# This is meant to be run interactively, e.g. `. automated_test.sh` #
#                                                                   #
# Author: A. Mecca (amecca@cern.ch)                                 # 
#####################################################################

events=100
cd $CMSSW_BASE/src ; printf "scram... " ; scram b -j 12 1>/dev/null && echo "Done" || { cd - ; return 1 2>/dev/null || exit 1 ; }
cd $CMSSW_BASE/src/VVXAnalysis/Producers/test || { echo "No VVXAnalysis in CMSSW/src" 1>&2 ; return 1 2>/dev/null || exit 1 ; }
rm -r test_auto || { echo "Failed to rm test_auto" 1>&2 ; return 1 2>/dev/null || exit 1 ; }
batch_Condor.py ../python/samples_PKU_2016UL_MC.csv -i ../python/analyzer_VVjj.py -o test_auto | grep -Pv "N:\s+\d*[1-9]\s" || { return 1 2>/dev/null || exit 1 ; }
cd test_auto/ZGToLLG-v17_Chunk0/ || { echo "cd into Chunk0 failed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }
line_n=$(( $(grep -n "process.maxEvents" run_cfg.py | grep -oP "([^:]+)(?=:)") + 1 ))
sed -i "${line_n}s/int32(.*)/int32($events)/" run_cfg.py >/dev/null || { echo "sed failed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }
[ $(grep -A1 "process\.maxEvents" run_cfg.py | grep -oP "(?<=int32\()[^)*]+(?=\))") -ne $events ] && { echo "Wrong sed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }
cmsRun run_cfg.py &> cmsrun.log &
less +F cmsrun.log
cd ../..
