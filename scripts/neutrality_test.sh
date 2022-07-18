#!/bin/bash

#SBATCH -A mdiop
#SBATCH -J FuLi
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem=15000
#SBATCH -w fmem2
#SBATCH --array=1-7
#SBATCH --error="FuLi.error"
#SBATCH --output="FuLi_%A_%a.out" 
#SBATCH --mail-type=BEGIN,FAIL,END          
#SBATCH --mail-user=mdiop@mrc.gm

ml java/10.0.2
ml bcftools/1.9
ml R/3.5.0
ml vcftools/0.1.15
ml tabix/0.2.6

function usage 
{
	printf "\n\t\t\t\t\tUsage: %s: -v <value> -g <value> -o <value> \n" $0  
	 echo -e "\n\n\t\t\t\t\t\t+-------------------------------------------------+"
     echo -e "\t\t\t\t\t\t|                   Essential parameters              |"
	 echo -e "\t\t\t\t\t\t|                -v <value>: Path to VCF Input File   |"
	 echo -e "\t\t\t\t\t\t|                -g <value>: Path to gtf file         |"
	 echo -e "\t\t\t\t\t\t|                -o <value>: output directory         |"
     echo -e "\t\t\t\t\t\t+-----------------------------------------------------+"
}

function validate 
{
	if [[ -d $v ]]
	then
		echo ""
	else
		echo -e $v "\n----- The file doesn't exist -----\n"
		usage
        exit
	fi
	
	if [[ -f $g ]]
	then
		echo ""
	else
		echo -e $g "\n----- The file doesn't exist -----\n"
		usage
        exit
	fi

	if [[ -d $o ]]
	then
		echo ""
	else
		echo -e "\n----- Incorrect value for -o option : Not a directory -----\n"
		usage
        exit
	fi
}

while getopts v:g:o: name
do
    case $name in
		v) v="$OPTARG";;
        g) g="$OPTARG";;
		o) o="$OPTARG";;
		h) usage;;
        ?) usage;;
    esac
done
shift $(($OPTIND - 1))
mkdir $o
validate

workingDir=`pwd`

seed=1
files=`ls $v/*gz`
file1=`echo $files | cut -d" " -f $SLURM_ARRAY_TASK_ID`

chmod +x neutrality_test.R
${workingDir}/neutrality_test.R $v $g $o

sleep 1
