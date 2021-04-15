#!/bin/tcsh
# Search for missing root files
# parameter to be specified: production path

#set echo 

#set basename=$1

foreach chunk ( *Chunk* )
 set fail="false"

 # Check that root file is existing and not empty
 set filename=${chunk}/$1
 if ( ! -e $filename ) then
   echo "Missing root file in " ${chunk}
   set fail="true"
 else if ( -z $filename ) then
#   echo "Empty file: " $filename
   set fail="true"
 endif

 if ( ! -e ${chunk}/exitStatus.txt) then
   echo "Missing txt file in " ${chunk}
   set fail="true"
 endif

#echo $chunk `ls ${chunk} | grep job | awk '{print $9}'`


# if ( $fail == "true" ) then
#   if ( -e ${filename}.recovered ) echo "   found: " ${filename}.recovered
#   if ( -e ${filename}.corrupted ) echo "   found: " ${filename}.corrupted
# else
#   set exitstatus=`grep "2012 with exit status" ${chunk}/*.txt | awk '{print $NF}'`
#   if ( $exitstatus != 0 ) then

 # Check job exit status
 set exitStatus = 0
 if ( -es ${chunk}/exitStatus.txt ) then
   set exitStatus=`cat ${chunk}/exitStatus.txt`
   set fail="true"     
 endif

 # Archive succesful jobs, or report failure
 if ( $fail == "false" ) then
  mkdir -p AAAOK
  mv $chunk AAAOK/
 else 
   set description=""
   switch($exitStatus)  # From https://twiki.cern.ch/twiki/bin/view/CMSPublic/StandardExitCodes
     case 0:   ; set description="(unknown)" ; breaksw
     case 84:  ; set description="(missing input file)" ; breaksw
     case 85:  ; set description="(failed to open local and fallback files)" ; breaksw
     case 134: ; set description="(Crashed)" ; breaksw
     case 152: ; set description="(Exceeded CPU time)" ; breaksw
     case 153: ; set description="(Exceeded File size limit)" ; breaksw
   endsw
   echo $chunk ": failed, exit status = " $exitStatus $description
 endif
end
# cd -
#end
