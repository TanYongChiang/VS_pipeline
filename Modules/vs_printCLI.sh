#!/bin/bash

set +a

function CLI_print_line ()
{
# 1- job name to print : count
# 2- current : count
    if [ $2 -lt $1 ]
    then
        CLI_status="\e[33mIN QUEUE..\e[0m"
    elif [ $2 == $1 ]
    then
        CLI_status="\e[5mRUNNING.. \e[0m"
    elif [[ $2 -gt $1 ]]
    then
        CLI_status="\e[92mFINISHED  \e[0m"
    fi
    
    case $1 in
        "1") local CLI_job="Preparing receptor PDBQT";;
        "2") local CLI_job="Generating DPF";;
        "3") local CLI_job="AutoSite: Pocket prediction";;
        "4") local CLI_job="Generating GPF";;
        "5") local CLI_job="Preparing ligand PDBQT";;
        "6") local CLI_job="General warmup";; # VS STARTS HERE
        "7") local CLI_job="Grid processing";;
        "8") local CLI_job="Virtual Screening";;
        "9") local CLI_job="Generating Summaries";;
    esac

    printf "\t|     | \e[2m>> %-31s\e[0m | $CLI_status  |     |\n" "$CLI_job"
}

while getopts c: flag
do
    case "${flag}" in
        c) current_CLI=${OPTARG}
    esac
done

clear
printf "\n"
printf "\t+--------------------------------------------------------------+\n"
printf "\t|                                                              |\n"
printf "\t|    \e[1m\e[96mP I P E \e[30m\e[106m L I N E \e[0m                                         |\n"
printf "\t|    ---------------                                           |\n"
printf "\t|    \e[2mDeveloped by Tan Yong Chiang, Sunway University\e[0m           |\n"
printf "\t|                                                              |\n"
printf "\t|     +------------------------------------+-------------|     |\n"
printf "\t|     | \e[1m\e[96mJOBS\e[0m                               | \e[1m\e[96mSTATUS\e[0m      |     |\n"
printf "\t|     +------------------------------------+-------------|     |\n"
printf "\t|     | \e[1mPREPARATION\e[0m                                      |     |\n"

if [[ $prepare_receptor == "y" || $prepare_receptor == "Y" ]]
then
    CLI_print_line "1" "$current_CLI"
fi

if [[ $AD_ALGO =~ ^[1-3]$ ]]
then
    CLI_print_line "2" "$current_CLI"
fi

if [[ $AUTO_SITE -eq 2 ]]
then
    CLI_print_line "3" "$current_CLI"
fi

if [[ $AUTO_SITE -eq 3 ]]
then
    CLI_print_line "4" "$current_CLI"
fi

if [ $(extension_present "pdbqt" "$ligand_dir") == false ]
then
    CLI_print_line "5" "$current_CLI"
fi

printf "\t|     |                                                  |     |\n"
printf "\t|     | \e[1mVIRTUAL SCREENING\e[0m                                |     |\n"
CLI_print_line "6" "$current_CLI"

if [[ $ALGO == "AD" ]]
then
    CLI_print_line "7" "$current_CLI"
fi

CLI_print_line "8" "$current_CLI"

CLI_print_line "9" "$current_CLI"


printf "\t|     +------------------------------------+-------------+     |\n"
printf "\t|                                                              |\n"
printf "\t+--------------------------------------------------------------+\n"
printf "\n\n"