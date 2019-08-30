#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=sfs.pairwise
#################  
#a file for job output, you can check job progress
#SBATCH --output=sfs.pairwise.%j.out
#################
# a file for errors from the job
#SBATCH --error=sfs.pairwise.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=6:00:00
#################
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=6
#################
#SBATCH --mem=12G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=BEGIN,END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=xxx
#################
#now run normal batch commands
##################
#echo commands to stdout

set -x

n=`ls *.saf.idx | wc -l`

angsd="~/programs/angsd"

for (( i = 1; i < n; i++ ))
do
    for (( j = i+1; j <= n; j++ ))
    do
        idxi=`ls -1 *.saf.idx | sed -n "$i"p`
        idxj=`ls -1 *.saf.idx | sed -n "$j"p`
	popi=`echo $idxi | cut -f1 -d'.'`
	popj=`echo $idxj | cut -f1 -d'.'`
	
	# calculate 2dsfs from saf files
	"$angsd"/misc/realSFS "$idxi" "$idxj" > fst/"$popi"."$popj".ml
	# prepare the fst for easy analysis
	"$angsd"/misc/realSFS fst index "$idxi" "$idxj" -sfs "$popi"."$popj".ml -fstout "$popi"."$popj"
	# global estimate of fst
	"$angsd"/misc/realSFS fst stats "$popi"."$popj".fst.idx > fst/"$popi"."$popj".fst.global.txt
	# sliding window: 50kb window, 25kb slide
	"$angsd"/misc/realSFS fst stats2 "$popi"."$popj".fst.idx -win 50000 -step 25000 -whichFST 1 > fst/"$popi"."$popj".fst.50kbin.txt
    done
done

## After this, combine the fst.global.txt files for isolation by distance analysis in R or whatever
# $ cat `ls *.global.txt` > Fst.global.txt
# $ echo `ls *.global.txt | cut -f1,2 -d'.'` > Pops.pairwise.txt
# $ paste Pops.pairwise.txt Fst.global.txt > Pairwise.fst.txt
## In Pairwise.fst.txt, col1 is pairwise pop names, col2 is fst unweighted, col3 is fst weighted (see Reynolds et al. 1983 for weighted vs unweighted)
