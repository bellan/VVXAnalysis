#!/bin/bash

################################################################
#  Loops through subdirectories to find jobs-containing ones   #
#  (they contain a condor.sub file), then checks if there are  #
#  still unfinished jobs. Otherwise, it hadds all the chunks   #
#  in AAAOK using haddChunks.py                                #
#  It also synchronizes the hadd-ed file to a folder in eos    #
#                                                              #
#  Author: A. Mecca  (alberto.mecca@cern.ch)                   #
################################################################

campaign=$(git describe --tags --abbrev=0)
initial=$(printf "%.1s" $USER)
eosdir="/eos/user/$initial/$USER/samples/$campaign"

source _findJobDirs.sh ; jobdirs=$(findJobDirs $@)  # All subfolders of the arguments containing a file named condor.sub

dirsExist () {
    for d in $@ ; do [ -d $d ] || return 1 ; done
}

for dir in $jobdirs ; do
    (
	cd $dir
	printf "%s \t--> " $dir
        # Were all the samples in this jobdir completed?
	if ls | grep -q Chunk ; then
	    # ... no
	    j_done=$( [ -d AAAOK ] && ls AAAOK/ | grep Chunk | wc -l || echo "0" )
	    j_todo=$(ls        | grep Chunk | wc -l)
	    printf "Still not done (%d/%d)\n" $j_done $(($j_done + $j_todo))
	    exit 0
	fi
	# ... yes

        # There may be more than one sample in a jobdir
        samples=$(ls AAAOK/ | grep -oP ".+(?=_Chunk)" | uniq)
        echo $samples | sed "s/\n/ /g" | xargs printf "%s "

	# Do we store rootfiles in AAAOK or is there a symlink at trasferPath?
	has_transfer_path=false
	for sample in $samples ; do
	    if [ -L transferPath ] ; then
		has_transfer_path=true
	    elif has_transfer_path ; then
		echo
		echo "WARNING: some samples in $dir have a transferPath symlink, but $sample does not!" 1>&2
		exit 1
	    fi
	done

	# Hadd either with a transfer path (repeat for each sample) or with the classic AAAOK directory (once for all the samples in this jobdir)
	jobdir="$(basename "${dir}")"
	if $has_transfer_path ; then
	    for sample in $samples ; do
		target_hadd="transferPath/${jobdir}/"
		# Was hadd already called?
		# if [ -e "${target_hadd}/${sample}/ZZ4lAnalysis.root" ] ; then
		if [ -d "${target_hadd}/${sample}" ] ; then
		    echo "(Hadd-ed)"
		else
		    haddChunks.py "${target_hadd}"
		fi
	    done
	else
	    # Was hadd already called?
	    if dirsExist $(echo "$samples" | sed "s:^:AAAOK/:g") ; then
		echo "(Hadd-ed)"
	    else
		haddChunks.py AAAOK
	    fi
	fi

	# Now transfer the rootfiles to eosdir
        for sample in $samples ; do
	    if ${has_transfer_path} ; then
		sourceFile="transferPath/${jobdir}/${sample}/ZZ4lAnalysis.root"
	    else
		sourceFile="AAAOK/${sample}/ZZ4lAnalysis.root"
	    fi
	    stripDir=$(echo $dir | sed -r "s:/[^/]+$::g")
	    destFile="$eosdir/$stripDir/${sample}.root"

	    # Check that the source file exists
	    if ! [ -e "$sourceFile" ] ; then
		if [ -e "$destFile" ] ; then
		    exit 0  # No source, but there is the destination. Probably the source was cleaned up
		else
		    printf "ERROR: Tried to hadd but no rootfile has been created at %s\n" "$sourceFile" 1>&2
		    exit 1
		fi
	    fi

	    # Transfer output
            mkdir -p "$eosdir/$stripDir"
	    if ! [ -e $destFile ] ; then
		if [ $(stat -c "%d" "$eosdir") -eq $(stat -c "%d" "$sourceFile") ] ; then
		    # Same filesystem, a simple mv will suffice
		    mv "$sourceFile" "$destFile" && echo "Copied to $destFile"
		else
		    rsync -au --progress "$sourceFile" "$destFile" && echo "Copied to $destFile"
		fi
	    fi
        done
    )
done
