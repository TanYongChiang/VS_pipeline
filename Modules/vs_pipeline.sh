#!/bin/bash
set -e
set +a

    #region add trap
log() { printf '%s\n' "$*"; }
error() { log "ERROR: $*" >&2; }
fatal() { error "$@"; exit 1; }

# appends a command to a trap
#
# - 1st arg:  code to add
# - remaining args:  names of traps to modify
#
trap_add() {
    trap_add_cmd=$1; shift || fatal "${FUNCNAME} usage error"
    for trap_add_name in "$@"; do
        trap -- "$(
            # helper fn to get existing trap command from output
            # of trap -p
            extract_trap_cmd() { printf '%s\n' "$3"; }
            # print existing trap command with newline
            eval "extract_trap_cmd $(trap -p "${trap_add_name}")"
            # print the new trap command
            printf '%s\n' "${trap_add_cmd}"
        )" "${trap_add_name}" \
            || fatal "unable to add to trap ${trap_add_name}"
    done
}
# set the trace attribute for the above function.  this is
# required to modify DEBUG or RETURN traps because functions don't
# inherit them unless the trace attribute is set
declare -f -t trap_add
    # endregion add trap
        # usage: trap_add 'command' EXIT

# region others 
# upon exit run this function
function remove_tmp_file_parallel {
    # remove temp_file for suppress warning only
    rm ${temp_file}
}
trap remove_tmp_file_parallel EXIT

# create temporary file to suppress warning in parallel -j0 only
temp_file=$(mktemp)

# Set prot_ext
## protein files pdbqt gpf dpf or pdbqt txt
if [[ $ALGO == "AD" ]]
then
    prot_ext=".pdbqt .gpf .dpf"
elif [[ $ALGO == "VINA" ]]
then
    prot_ext=".pdbqt .txt"
fi
# endregion

# region functions 
function check_exist () 
{
    if [ -e "$1" ]; then
        :
    else
        print_err_color "$1 does not exist"
        exit 1
    fi
}

function vina_check_cpu () 
{
    # check if number of cpu has been declared in config.txt
    # if present, then delete the cpu no. declaration
    # 1- receptor name
    # 2- protein_dir
    cpu_present=$(cat $2/$1.txt | awk '/(cpu =|cpu=)/ {print}')
    # if output is null
    if [ -z $cpu_present ]
    then 
        :
    else 
        cat $2/$1.txt | awk '!/(cpu =|cpu=)/ {print}' > $2/$1.txt
    fi
}

function prepareGPF () 
{
# dir, receptor name, lig name, MGL_dir
    cd $1
    $4/bin/pythonsh $4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l $3.pdbqt -r $2.pdbqt -i $2.gpf -o grid.gpf > /dev/null 2>&1
}

function prepareDPF ()
{
# dir, receptor name, lig name, MGL_dir
    cd $1
    if [[ $(grep -w "do_local_only" $2.dpf | wc -l) -ge 1 ]]
    then
        param="-L"
    elif [[ $(grep -w "simanneal" $2.dpf | wc -l) -ge 1 ]]
    then
        param="-S"
    elif [[ $(grep -w "do_global_only" $2.dpf | wc -l) -ge 1 ]]
    then
        printf "\nga_run $(grep "do_global_only" $2.dpf | awk '{print $2}')\n" >> $2.dpf
    fi
    $4/bin/pythonsh $4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -l $3.pdbqt -r $2.pdbqt $param -i $2.dpf -o dock.dpf > /dev/null 2>&1
    if [[ $(grep -w "do_local_only" $2.dpf | wc -l) -ge 1 ]] || [[ $(grep -w "simanneal" $2.dpf | wc -l) -ge 1 ]]
    then
        grep -v "unbound_model" dock.dpf > dock.temp
        cat dock.temp > dock.dpf
        rm -f dock.temp
    fi
}

function autogrid () 
{
# receptor name, lig name, working dir
    cd $3/docking/$1/$2/
    autogrid4 -p grid.gpf -l grid.glg > /dev/null 2>&1
}

function autodock () 
{
# receptor name, lig name, working dir
    cd $3/docking/$1/$2/
    autodock4 -p dock.dpf -l $1__$2.dlg > /dev/null 2>&1
}

function ad_vina ()
{
    #---------------------------
    # 1- receptor name, 
    # 2- lig name, 
    # 3- working dir, 
    #---------------------------
    cd $3/docking/$1/$2/
    vina --receptor $1.pdbqt --ligand $2.pdbqt --config $1.txt --out out.pdbqt --log log.pdbqt > /dev/null 2>&1
}

function vina_summaries ()
{
    # 1- receptor name, 
    # 2- lig name, 
    # 3- working dir,
    # 4- result_dir
    cd $3/docking/$1/$2/
    result_byline=$(awk 'BEGIN {FS=" "; OFS=","}; /REMARK VINA RESULT/ {print "'$1'","'$2'",$4}' out.pdbqt | head -n 1)
    printf "$result_byline" | tr -d [:space:] >> $4/summary_virtualscreening.csv
    printf "\n" >> $4/summary_virtualscreening.csv
}

function vs_genenergymatrix ()
{
# 1- working_dir
# 2- result_dir
# 3- ALGO
    printf "Binding Energies" > $2/vs_energymatrix.csv
    
    # x axis
    while read receptor_byline || [ -n "$line" ]
    do
        printf "," >> $2/vs_energymatrix.csv
        printf "$receptor_byline" >> $2/vs_energymatrix.csv
    done < $1/proteins.csv
    printf "\n" >> $2/vs_energymatrix.csv

    while read ligand_byline || [ -n "$line" ]
    do
        printf "$ligand_byline" >> $2/vs_energymatrix.csv
        while read receptor_byline || [ -n "$line" ]
        do
            printf "," >> $2/vs_energymatrix.csv
            cd $1/docking/$receptor_byline/$ligand_byline/
            if [[ $3 == "AD" ]]
            then
                awk 'BEGIN {FS=","; OFS=","}; NR==2 {print $5}' summary | tr -d [:space:] >> $2/vs_energymatrix.csv
            elif [[ $3 == "VINA" ]]
            then
                local result_byline=$(awk 'BEGIN {FS=" "; OFS=","}; /REMARK VINA RESULT/ {print $4}' out.pdbqt | head -n 1)
                printf "$result_byline" | tr -d [:space:] >> $2/vs_energymatrix.csv
            fi
        done < $1/proteins.csv
        printf "\n" >> $2/vs_energymatrix.csv
    done < $1/ligands.csv
}

function vs_summaries ()
{
    # 1- receptor name, 
    # 2- lig name, 
    # 3- working dir, 
    # 4- MGL_dir, 
    # 5- result_dir
    cd $3/docking/$1/$2/
    $4/bin/pythonsh $4/MGLToolsPckgs/AutoDockTools/Utilities24/summarize_results4.py -d $3/docking/$1/$2/ -b -k -r $1.pdbqt -o summary > /dev/null 2>&1
    awk 'BEGIN {FS=","; OFS=","}; NR==2 {print "'$1'","'$2'",$5,$7}' summary | tr -d [:space:] >> $5/summary_virtualscreening.csv
    printf "\n" >> $5/summary_virtualscreening.csv
}

function extension_present ()
{
# 1- suffix name
# 2- directory
# return true or false
    cd $2
    local filecount=$(find . -maxdepth 1 -type f -name "*.$1" | wc -l)
    if [[ $filecount -ge 1 ]]
    then
        echo true
    else
        echo false
    fi
}

function print_err_color ()
{
    printf "\n\n\e[1m\e[91mERROR: \e[39m$1\e[0m\n\n"
}

function gen_best_structure ()
{
# 1- protein name
# 2- ligand name
# 3- working_dir
# 4- result_dir
    local dlg_file="$3/docking/$1/$2/$1__$2.dlg"
    mkdir -p $4/best_structures/$1/$2/
    cat $3/docking/$1/$2/$1.pdbqt | awk '/^(ATOM|TER)/ {print}' > $4/best_structures/$1/$2/DOCKED_COMPLEX.pdb
    trigger=0
    while read line || [ -n "$line" ]
    do
        if [[ "$line" =~ ^MODEL ]]
        then
            trigger=$(( trigger+1 ))
        fi
        if [ "$trigger" == 1 ]
        then
            if [[ "$line" =~ ^ATOM ]] || [[ "$line" =~ ^TER ]]
            then
                echo "$line" >> $4/best_structures/$1/$2/DOCKED_LIGAND.pdb
                echo "$line" >> $4/best_structures/$1/$2/DOCKED_COMPLEX.pdb
                if [[ "$line" =~ ^TER ]]
                then
                    trigger=$(( trigger+1 ))
                fi
            fi
        fi
    done < ${dlg_file}
}
# endregion functions

# region preparations 

    # region for CLI
export prepare_receptor
export AD_ALGO
export AUTO_SITE
export -f extension_present
export ligand_dir
export ALGO
    # endregion for CLI

# create list of proteins and ligands in csv

bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 6

parallel -j 0 echo {/.} ::: $(ls $protein_dir/*.pdbqt) > $working_dir/proteins.csv
parallel -j 0 echo {/.} ::: $(ls $ligand_dir/*.pdbqt) > $working_dir/ligands.csv
trap_add "rm -f $working_dir/proteins.csv" EXIT
trap_add "rm -f $working_dir/ligands.csv" EXIT

echo "Checking: proper preparation of protein files"
export -f check_exist
export -f vina_check_cpu
parallel -j 0 check_exist "$protein_dir/{1}{2}" :::: $working_dir/proteins.csv ::: $prot_ext

if [[ $ALGO == "VINA" ]]
then
    parallel -j 0 vina_check_cpu "{1}" "$protein_dir" :::: $working_dir/proteins.csv
fi

echo "Virtual Screening Preparation: Allocating input files"
# read csv, the allocate files *working_dir/*protein/*ligand/*files*
cat ${temp_file} | parallel -j 0 mkdir -p $working_dir/docking/{1}/{2} :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
cat ${temp_file} | parallel -j 0 cp -a $protein_dir/{1}{3} $working_dir/docking/{1}/{2} :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv ::: $prot_ext
cat ${temp_file} | parallel -j 0 cp -a $ligand_dir/{2}.pdbqt $working_dir/docking/{1}/{2} :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv

if [[ $ALGO == "AD" ]]
then
    echo "Virtual Screening Preparation: Configuring GPF "
    export -f prepareGPF
    parallel -j $number_jobs --eta prepareGPF "$working_dir/docking/{1}/{2}" "{1}" "{2}" "$MGL_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
    echo "Virtual Screening Preparation: Configuring DPF"    
    for protein_files in $(ls $protein_dir/*.pdbqt)
    do 
        protein_file=$(printf '%s\n' "${protein_files%.pdbqt}")
        if [[ $(grep -w "do_local_only" $protein_file.dpf | wc -l) -ge 1 ]] || [[ $(grep -w "simanneal" $protein_file.dpf | wc -l) -ge 1 ]]
        then
            if [[ $(grep -w "unbound_model" $protein_file.dpf | wc -l) -ge 1 ]] && [[ $(grep -w "unbound_model" $protein_file.dpf | awk '{print $2}') -ne "bound" ]]
            then
                printf "\n\e[2m\e[31mWARNING: unbound_model other than bound is not supported for non-GA algorithms, unbound_model bound is automatically specified instead. \e[0m\n"
                break
            fi
        fi
        if [[ $(grep -w "do_global_only" $protein_file.dpf | wc -l) -ge 1 ]]
        then
            printf "\n\e[2m\e[31mWARNING: Sole genetic algorithm is not currently supported by AutoDock4 virtual screening, converted into Lamarckian Genetic Algorithm with same number of runs. \e[0m\n"
            break
        fi
    done
    export -f prepareDPF
    parallel -j $number_jobs --eta prepareDPF "$working_dir/docking/{1}/{2}" "{1}" "{2}" "$MGL_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
fi

# endregion

# region virtual screening

if [[ $ALGO == "AD" ]]
then
    # Run autodock4
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 7
    echo "Virtual Screening: grid maps processing... "
    export -f autogrid
    parallel -j $number_jobs --eta autogrid "{1}" "{2}" "$working_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv

    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 8
    echo "Virtual Screening: virtual screening... "
    export -f autodock
    parallel -j $number_jobs --eta autodock "{1}" "{2}" "$working_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv

elif [[ $ALGO == "VINA" ]]
then
    # Run vina
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 8
    echo "Virtual Screening: virtual screening... "
    export -f ad_vina
    parallel -j $number_jobs --eta ad_vina "{1}" "{2}" "$working_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv

fi

# endregion

# region analysis 
bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 9
if [[ $ALGO == "AD" ]]
then
    echo "Virtual Screening: generating summaries... "
    echo "Receptor,Ligand,Binding Energy (kcal/mol),Intermolecular hydrogen bonds" > $result_dir/summary_virtualscreening.csv
    export -f vs_summaries
    cat ${temp_file} | parallel -j 1 --eta vs_summaries "{1}" "{2}" "$working_dir" "$MGL_dir" "$result_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
elif [[ $ALGO == "VINA" ]]
then
    echo "Virtual Screening: generating summaries... "
    echo "Receptor,Ligand,Binding Energy (kcal/mol)" > $result_dir/summary_virtualscreening.csv
    export -f vina_summaries
    cat ${temp_file} | parallel -j 1 --eta vina_summaries "{1}" "{2}" "$working_dir" "$result_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
fi
vs_genenergymatrix "$working_dir" "$result_dir" "$ALGO"

if [[ $ALGO == "AD" ]]
then
    export -f gen_best_structure
    cat ${temp_file} | parallel -j 0 gen_best_structure "{1}" "{2}" "$working_dir" "$result_dir" :::: $working_dir/proteins.csv :::: $working_dir/ligands.csv
fi

bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 10
printf "DONE! Check the RESULT folder in your working directory!\n\n"
# endregion

# region unset functions 
unset -f check_exist
unset -f vina_check_cpu
unset -f prepareGPF
unset -f prepareDPF
unset -f autogrid
unset -f autodock
unset -f ad_vina
unset -f vs_summaries
unset -f vina_summaries
unset prepare_receptor
unset AD_ALGO
unset AUTO_SITE
unset -f extension_present
unset ligand_dir
unset ALGO
# endregion