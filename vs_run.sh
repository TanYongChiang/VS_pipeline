#!/bin/bash
set -e
set +a

SCRIPTPATH="$( cd -- "$(dirname "$0")" > /dev/null 2>&1 ; pwd -P )"

printf "\n\t\t\t\e[96m\e[1mVirtual Screening Pipeline \e[0m\n\n"
printf "\t\e[2mDeveloper: Tan Yong Chiang\n"
printf "\tAffiliation: Sunway University\n"
printf "\tMail: alexyc991@gmail.com\e[0m\n\n"

export SCRIPTPATH

bash $SCRIPTPATH/Modules/vs_prepare.sh