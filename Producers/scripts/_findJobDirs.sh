findJobDirs() {
    [ $# -ge 1 ] && top=$@ || top=.
    echo $(find $top -maxdepth 4 -name "*Chunk*" -prune -o -name condor.sub | grep -oP ".+(?=/condor.sub$)" | sed "s|^\./||g" )
}
