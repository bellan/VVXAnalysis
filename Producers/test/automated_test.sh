# bash script

#####################################################################
# Automates scram, batch_condor.py and cmsRun of a test sample      #
# This is meant to be run interactively, e.g. `. automated_test.sh` #
#                                                                   #
# Author: A. Mecca (amecca@cern.ch)                                 # 
#####################################################################

analyzer=../python/analyzer_VVjj.py
specific_file=""  # es. file:/eos/home-a/amecca/samples/.../somefile.root
csv_file=../python/samples_PKU_2016UL_MC.csv
sample_ID=ZGToLLG-v17
events=100


# First, recompile CMSSW and remove any results from last execution
cd $CMSSW_BASE/src ; printf "scram... " ; scram b -j 12 1>/dev/null && echo "Done" || { cd - ; return 1 2>/dev/null || exit 1 ; }
cd $CMSSW_BASE/src/VVXAnalysis/Producers/test || { echo "No VVXAnalysis in CMSSW/src" 1>&2 ; return 1 2>/dev/null || exit 1 ; }
rm -r test_auto || { echo "Failed to rm test_auto" 1>&2 ; return 1 2>/dev/null || exit 1 ; }

# Create the config using the specified template ("analyzer") for all the uncommented lines in csv_file
batch_Condor.py ${csv_file} -i $analyzer -o test_auto | grep -Pv "N:\s+\d*[1-9]\s" || { return 1 2>/dev/null || exit 1 ; }

# Take the first chunk, which should always exist
cd test_auto/${sample_ID}_Chunk0/ || { echo "cd into Chunk0 failed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }

# process.maxEvents is split in more than one line: figure out at which line number does it start...
line_n=$(( $(grep -n "process.maxEvents" run_cfg.py | grep -oP "([^:]+)(?=:)") + 1 ))
# ... then edit the next line...
sed -i "${line_n}s/int32(.*)/int32($events)/" run_cfg.py >/dev/null || { echo "sed maxEvents failed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }
# and check that everything went well
[ $(grep -A1 "process\.maxEvents" run_cfg.py | grep -oP "(?<=int32\()[^)*]+(?=\))") -ne $events ] && { echo "Wrong sed maxEvents." 1>&2 ; return 1 2>/dev/null || exit 1 ; }

# Optional: insert a custom filename
[ -n $specific_file ] && sed -i -r "s/(fileNames = cms\.untracked\.vstring).+$/\1('${specific_file}')/" run_cfg.py 1>/dev/null

# Run in background, redirect stdout and stderr to a file and read from it (so that the log is permanent and you can ^C without killing cmsRun)
cmsRun run_cfg.py &> cmsrun.log &
less +F cmsrun.log
cd ../..
