#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks=8
#SBATCH --partition=debug

# function to change Par_file
change_parfile() {
    local param=$1
    local value=$2 
    local file="DATA/Par_file"

    local oldstr=`grep "^$param " $file`
    local newstr="$param     =       $value"

    sed  "s?$oldstr?$newstr?g" $file  > tmp
    mv tmp $file 
}


# run force
change_parfile USE_FORCE_POINT_SOURCE .true.
change_parfile USE_CMT_AND_FORCESOLUTION .false.
bash run.sh 
mv OUTPUT_FILES OUTPUT_FILES.FORCE

# run cmt
change_parfile USE_FORCE_POINT_SOURCE .false. 
change_parfile USE_CMT_AND_FORCESOLUTION .false.
bash run.sh 
mv OUTPUT_FILES OUTPUT_FILES.CMT

change_parfile USE_CMT_AND_FORCESOLUTION .true.
bash run.sh 

