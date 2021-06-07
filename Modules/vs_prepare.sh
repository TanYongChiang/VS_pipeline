#!/bin/bash

set +a
set -e

# region functions 
function check_cmd ()
{
    if ! command -v $1 &> /dev/null
    then
        print_err_color "$1 is not installed or callable"
        exit 1
    fi
}

function check_exist () 
{
    if [ -e "$1" ]; then
        :
    else
        print_err_color "$1 does not exist"
        exit 1
    fi
}

function check_spaces ()
{
    string="$@"
    pattern=" "
    if [[ $string =~ $pattern ]]
    then
        print_err_color "Directory $string contain spaces, please delete the spaces or replace them with _"
        exit 1
    else
        :
    fi
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

function prepare_receptor4_py ()
{
# 1- receptor name without .pdb
# 2- MGL_dir
# 3- protein_dir
# convert to pdbqs then pdbqt because kollman inaccessible for pdbqt but pdbqs can
    $2/bin/pythonsh $2/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor.py -r $3/$1.pdb -o $3/$1.pdbqs > /dev/null 2>&1
    $2/bin/pythonsh $2/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqs_to_pdbqt.py -s $3/$1 -C -o $3/$1.pdbqt > /dev/null 2>&1
    rm -f $3/$1.pdbqs
}

function generate_dpf ()
{
# 1- AD_ALGO
# 2- num_runs
# 3- protein_dir
# 4- protein name
    if [[ $1 -eq 1 ]]
    then
        printf "set_ga\n" > $3/$4.dpf
        printf "ga_run $2\n" >> $3/$4.dpf
    elif [[ $1 -eq 2 ]]
    then
        printf "runs $2\n" > $3/$4.dpf
        printf "simanneal\n" >> $3/$4.dpf
    elif [[ $1 -eq 3 ]]
    then
        printf "do_local_only $2\n" > $3/$4.dpf
    fi
}

function autosite_pred ()
{
# 1- MGL2_dir
# 2- protein_dir
# 3- protein name
# 4- size of box (A)
# 5- ALGO
# 6- num_runs
# 7- AUTO_SITE_SHAPE
    $1/bin/pythonsh $1/MGLToolsPckgs/AutoSite/bin/AS.py -r $2/$3.pdbqt > /dev/null 2>&1
    unset start_col
    unset end_col
    unset dim_min
    unset dim_max
    unset dim_midpt
    
    declare -A start_col
    declare -A end_col
    declare -A dim_min
    declare -A dim_max
    declare -A dim_midpt
    for dim in $(printf "x y z")
    do  
        case $dim in
            "x")
                start_col["$dim"]="31"
                end_col["$dim"]="38";;
            "y") 
                start_col["$dim"]="39"
                end_col["$dim"]="46";;
            "z") 
                start_col["$dim"]="47"
                end_col["$dim"]="54";;
        esac
        coords=$(cat $2/$3_cl_1.txt | awk 'BEGIN {FS=""};{for(i='"${start_col["$dim"]}"';i<='"${end_col["$dim"]}"';i++) printf $i; print ""}')
        # dim_avrg["$dim"]=$(printf "$coords" | awk 'BEGIN {sum=0};{sum+=$1};END {printf "%f\n", sum/NR}')
        dim_min["$dim"]=$(printf "$coords" | awk 'NR==1 {min=$1};{if ($1<min) min=$1};END {printf "%f\n", min}')
        dim_max["$dim"]=$(printf "$coords" | awk 'NR==1 {max=$1};{if ($1>max) max=$1};END {printf "%f\n", max}')
        dim_midpt["$dim"]=$(printf "${dim_min["$dim"]} ${dim_max["$dim"]}" | awk '{print ($1+$2)/2}')
    done
    if [[ $5 == "AD" ]]
    then
        printf "npts " > $2/$3.gpf
        for i in {1..3..1}
        do
            printf "$(awk 'BEGIN {printf "%.0f ", '"$4/0.375"'}') " >> $2/$3.gpf
        done
        printf "\nspacing 0.375\n" >> $2/$3.gpf
        printf "gridcenter " >> $2/$3.gpf
        for dim in $(printf "x y z")
        do
            printf "%.4f " ${dim_midpt["$dim"]} >> $2/$3.gpf
        done
    elif [[ $5 == "VINA" ]]
    then
        printf "" > $2/$3.txt
        for dim in $(printf "x y z")
        do
            printf "center_$dim = ${dim_midpt["$dim"]}\n" >> $2/$3.txt
        done
        printf "\n" >> $2/$3.txt
        for dim in $(printf "x y z")
        do
            printf "size_$dim = $4\n" >> $2/$3.txt
        done
        printf "\n" >> $2/$3.txt
        printf "num_modes = $6" >> $2/$3.txt
    fi
    rm -f $(ls $2/*.pdb | grep "$3\_fp\_")
    cl_filter=$(ls -d $2/*.pdb | grep "$3\_cl\_")
    rm -f $(printf "$cl_filter" | grep -v "$3\_cl\_1\.pdb")
}

function pdbqt_xyz ()
{
# 1- protein_dir
# 2- resiconf_file
# 3- ALGO
# 4- expand xyz in Angstrom (gpf_length)
# 5- num_runs
    while read line || [ -n "$line" ]
    do
        prot_name="$(printf -- "$line" | awk -F "," '{print $1}')"
        chn_res_s="$(printf -- "$line" | awk -F "," '{print $2}')"
        get_ATOM="$(cat $1/$prot_name.pdbqt | grep ^ATOM)"
        printf -- "" > $1/$prot_name.pdbqtxyz

        # detect type of res input, chn_res or just res
        chn_res_val=$(echo "$chn_res_s" | awk -F " " '{print $1}')
        col2_val=$(printf -- "$chn_res_val" | awk -F "_" '{print $2}')

        if ! [[ -n "$col2_val" ]]
        then
            for chn_res in $chn_res_s
            do
                res=$(printf -- "$chn_res")
                filterchn_res_s=$(printf -- "$get_ATOM" | awk 'BEGIN {FS=""};{for(i=23;i<=26;i++) printf $i; print ""}')
                filterres_line_s=$(printf -- "$filterchn_res_s" | awk '{if($1 == "'"$res"'") print FNR}')
                for filterres_line in $filterres_line_s
                do
                    printf -- "$get_ATOM" | sed "${filterres_line}q;d" >> $1/$prot_name.pdbqtxyz
                done
            done
        else
            for chn_res in $chn_res_s
            do
                chn=$(printf -- "$chn_res" | awk -F "_" '{print $1}')
                res=$(printf -- "$chn_res" | awk -F "_" '{print $2}')
                filterchn=$(printf -- "$get_ATOM" | awk 'BEGIN {FS=""};{if($22 == "'"$chn"'") print}')
                
                filterchn_res_s=$(printf -- "$filterchn" | awk 'BEGIN {FS=""};{for(i=23;i<=26;i++) printf $i; print ""}')
                filterres_line_s=$(printf -- "$filterchn_res_s" | awk '{if($1 == "'"$res"'") print FNR}')
                for filterres_line in $filterres_line_s
                do
                    printf -- "$filterchn" | sed "${filterres_line}q;d" >> $1/$prot_name.pdbqtxyz
                done
            done
        fi

        # generate start end midpt
        unset start_col
        unset end_col
        unset dim_min
        unset dim_max
        unset dim_midpt
        unset dim_npts

        declare -A start_col
        declare -A end_col
        declare -A dim_min
        declare -A dim_max
        declare -A dim_midpt
        declare -A dim_npts
        for dim in $(printf "x y z")
        do  
            case $dim in
                "x")
                    start_col["$dim"]="31"
                    end_col["$dim"]="38";;
                "y") 
                    start_col["$dim"]="39"
                    end_col["$dim"]="46";;
                "z") 
                    start_col["$dim"]="47"
                    end_col["$dim"]="54";;
            esac
            coords=$(cat $1/$prot_name.pdbqtxyz | awk 'BEGIN {FS=""};{for(i='"${start_col["$dim"]}"';i<='"${end_col["$dim"]}"';i++) printf $i; print ""}')
            dim_min["$dim"]="$(printf -- "$coords" | awk 'NR==1 {min=$1};{if ($1<min) min=$1};END {printf "%f\n", min}')"
            dim_max["$dim"]="$(printf -- "$coords" | awk 'NR==1 {max=$1};{if ($1>max) max=$1};END {printf "%f\n", max}')"
            dim_midpt["$dim"]="$(printf -- "${dim_min["$dim"]} ${dim_max["$dim"]}" | awk '{print ($1+$2)/2}')"
            dim_npts["$dim"]="$(printf -- "${dim_max["$dim"]} ${dim_midpt["$dim"]} $4" | awk '{print (($1-$2)*2)+$3}')"
        done

        if [[ $3 == "AD" ]]
        then
            printf "npts " > $1/$prot_name.gpf
            for dim in $(printf "x y z")
            do
                printf -- "$(awk 'BEGIN {printf "%.0f ", '"${dim_npts["$dim"]}/0.375"'}') " >> $1/$prot_name.gpf
            done
            printf "\nspacing 0.375\n" >> $1/$prot_name.gpf
            printf "gridcenter " >> $1/$prot_name.gpf
            for dim in $(printf "x y z")
            do
                printf -- "%.4f " "${dim_midpt["$dim"]}" >> $1/$prot_name.gpf
            done
        elif [[ $3 == "VINA" ]]
        then
            printf "" > $1/$prot_name.txt
            for dim in $(printf "x y z")
            do
                printf "center_$dim = ${dim_midpt["$dim"]}\n" >> $1/$prot_name.txt
            done
            printf "\n" >> $1/$prot_name.txt
            for dim in $(printf "x y z")
            do
                printf "size_$dim = ${dim_npts["$dim"]}\n" >> $1/$prot_name.txt
            done
            printf "\n" >> $1/$prot_name.txt
            printf "num_modes = $5" >> $1/$prot_name.txt
        fi
        # rm -f $1/$prot_name.pdbqtxyz
    done < $2
}

function obabel_addh ()
{
# 1- bab_addh
# 2- ligand full location + suffix

    if [[ "$1" == 1 ]]
    then
        local bab_addh='-p 7.4'
    elif [[ "$1" == 2 ]]
    then
        local bab_addh='-h'
    fi
    obabel "$2" -O "$2" $bab_addh &> /dev/null
}

function obabel_gen3d ()
{
# 1- bab_gen3d_speed
# 2- ligand full location with suffix

    case "$1" in
        "1") local bab_gen3d_speed='--fastest';;
        "2") local bab_gen3d_speed='--fast';;
        "3") local bab_gen3d_speed='--medium';;
        "4") local bab_gen3d_speed='--better';;
        "5") local bab_gen3d_speed='--best';;
    esac
    local lig_loc_nosuffix=$(echo "$2" | cut -f 1 -d '.')

    obabel "$2" -O "$lig_loc_nosuffix.pdb" --gen3d $bab_gen3d_speed &> /dev/null
}

function obabel_conformer ()
{
# 1- bab_conformer_ff
# 2- bab_conformer_score
# 3- bab_conformer_no
# 4- ligand full location + suffix

    case "$1" in
        "1") local bab_conformer_ff='--systematic';;
        "2") local bab_conformer_ff='--random';;
        "3") local bab_conformer_ff='--weighted';;
        "4") local bab_conformer_ff='--ff MMFF94';;
        "6") return 0;;
    esac

    if [[ "$1" == 5 ]]
    then
        case "$2" in
            "1") local bab_conformer_ff='--score rmsd';;
            "2") local bab_conformer_ff='--score energy';;
        esac
    fi
    obabel "$4" -O "$4" --conformer -nconf $3 $bab_conformer_ff &> /dev/null
}

function obabel_minimization ()
{
# 1- bab_minim_sd_nsteps
# 2- bab_minim_sd_conv
# 3- bab_minim_cg_nsteps
# 4- bab_minim_cg_conv
# 5- bab_minim_ff
# 6- ligand full location + suffix

    case "$5" in
        "1") local bab_minim_ff='-ff MMFF94';;
        "2") local bab_minim_ff='-ff MMFF94s';;
        "3") local bab_minim_ff='-ff Ghemical';;
        "4") local bab_minim_ff='-ff Gaff';;
        "5") local bab_minim_ff='-ff Uff';;
    esac

    local lig_loc_nosuffix=$(echo "$6" | cut -f 1 -d '.')

    if [[ $1 -gt 0 ]]
    then
        obminimize -sd -n $1 -c $2 $bab_minim_ff $6 > $lig_loc_nosuffix.bblmintmp 2> /dev/null
        cat $lig_loc_nosuffix.bblmintmp > $6
    fi
    if [[ $3 -gt 0 ]]
    then
        obminimize -cg -n $3 -c $4 $bab_minim_ff $6 > $lig_loc_nosuffix.bblmintmp 2> /dev/null
        cat $lig_loc_nosuffix.bblmintmp > $6
    fi
    rm -f $lig_loc_nosuffix.bblmintmp
}

function obabel_mol2_pdb ()
{
# 1- ligand full location + .mol2
    local lig_loc_nosuffix=$(echo "$1" | cut -f 1 -d '.')
    obabel "$1" -O "$lig_loc_nosuffix.pdb" &> /dev/null
}

function prepare_ligand4_py ()
{
# 1- MGL_dir
# 2- ligand full location + suffix
    local lig_loc_nosuffix=$(echo "$2" | cut -f 1 -d '.')
    $1/bin/pythonsh $1/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l $2 -o $lig_loc_nosuffix.pdbqt > /dev/null 2>&1
}

function print_err_color ()
{
    printf "\n\n\e[1m\e[91mERROR: \e[39m$1\e[0m\n\n"
}

function print_ask_head ()
{
    printf "\n\n\e[95m\e[92m$1\n\n\e[0m"
}

function print_ask_body ()
{
    printf "\e[22m\e[32m$1\n\e[0m"
}

function ligand_check_atm ()
{
# 1- ligand_dir
# 2- ligand.pdbqt full path + suffix
    local ligand_atmtype=$(cat $2 | awk -F "" '/^(HETATM|ATOM)/ {for(i=78;i<=79;i++) printf $i; print ""}')
    for atm_line in $ligand_atmtype
    do
        if ! [[ "$atm_line" =~ ^(H|HD|HS|C|A|N|NA|NS|O|OA|OS|F|Mg|MG|P|SA|S|Cl|CL|Ca|CA|Mn|MN|Fe|FE|Zn|ZN|Br|BR|I)$ ]]
        then
            mkdir -p $1/Incompatible
            mv $2 $1/Incompatible/
        fi
    done
}
# endregion functions

# region IMPORTANT CHECKING 

    # region checking essential environment 
printf "\nChecking: essential environments\n"

# check autodock and autogrid command
check_cmd "parallel"
if [[ $ALGO == "AD" ]]
then
    for command in autodock4 autogrid4
    do
        check_cmd "$command"
    done
elif [[ $ALGO == "VINA" ]]
then
    check_cmd "vina"
fi
check_cmd "obabel"
    # endregion checking essential environment

# endregion CHECKING

# region GET important directories 
print_ask_head "Enter directory containing receptor files: "
while true
do
    read protein_dir
    if ! [ -z $protein_dir ]
    then
        if [ -e "$protein_dir" ]
        then
            break
        fi
    fi
    print_err_color "Directory $protein_dir not found:"
done

print_ask_head "Enter directory containing ligand files: "
while true
do
    read ligand_dir
    if ! [ -z $ligand_dir ]
    then
        if [ -e "$ligand_dir" ]
        then
            break
        fi
    fi
    print_err_color "Directory $ligand_dir not found:"
done

print_ask_head "Enter MGLTools install directory: (e.g. /xx/xx/MGLTools-1.5.6)"
while true
do
    read MGL_dir
    if ! [ -z $MGL_dir ]
    then
        if [ -e "$MGL_dir" ]
        then
            break
        fi
    fi
    print_err_color "Directory $MGL_dir not found:"
done
# check MGL directories
export -f check_exist
parallel -j 0 check_exist "{1}" ::: $MGL_dir $MGL_dir/bin/pythonsh $MGL_dir/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py $MGL_dir/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py

print_ask_head "Enter working directory:"
while true
do
    read real_working_dir
    if ! [ -z $real_working_dir ]
    then
        if [ -e "$real_working_dir" ]
        then
            if [ -n "$(ls -A $real_working_dir)" ]; then
                print_err_color "Directory $real_working_dir not empty:"
            else
                break
            fi
        fi
    fi
    print_err_color "Directory $real_working_dir not found:"
done

print_ask_head "Virtual Screening Algorithm:"
print_ask_body " 1: AutoDock 4"
print_ask_body " 2: AutoDock Vina"
printf "\n"
while true
do
    read ALGO
    if ! [ -z $ALGO ]
    then
        if [ "$ALGO" -eq 1 ]
        then
            ALGO="AD"
            break
        elif [ "$ALGO" -eq 2 ]
        then
            ALGO="VINA"
            break
        fi
    fi
    print_err_color "Please enter 1 or 2:"
done

mkdir -p $real_working_dir/WORKING
working_dir="$real_working_dir/WORKING"
mkdir -p $real_working_dir/RESULT
result_dir="$real_working_dir/RESULT"

echo -e "\nChecking: input directories\n"
export -f check_spaces
parallel -j 0 check_spaces "{}" ::: "$ligand_dir" "$protein_dir" "$working_dir" "$result_dir" "$MGL_dir"

if [ -n "$(ls -A $working_dir)" ]; then
    print_err_color "Directory $working_dir not empty"
    exit 1
fi
if [ -n "$(ls -A $result_dir)" ]; then
    print_err_color "Directory $result_dir not empty"
    exit 1
fi

# endregion GET important directories

print_ask_head "Enter number of jobs to be run in parallel: "
while true
do
    read number_jobs
    if ! [ -z $number_jobs ]
    then
        if [[ $number_jobs =~ ^[0-9]+$ ]]
        then
            break
        fi
    fi
    print_err_color "Please enter an integer: "
done

# region RECEPTOR 

cd $protein_dir
printf "Receptor: checking receptor directory\n"

    # region protein pdbqt absent
if [ $(extension_present "pdbqt" "$protein_dir") == false ]
then
    printf "Receptor: pdbqt files undetected\n"
    if [ $(extension_present "pdb" "$protein_dir") == true ]
    then
        printf "Receptor: pdb files detected\n"
        print_ask_head "Prepare receptor through prepare_receptor.py? (Y/N)"
        while true
        do
            read prepare_receptor
            if [[ "$prepare_receptor" =~ ^(Y|y|N|n)$ ]]
            then
                break
            fi
            print_err_color "Please enter Y or N:"
        done
        if [[ $prepare_receptor == "n" || $prepare_receptor == "N" ]]
        then
            printf "\n\n\e[93mPrepare the receptor .pdbqt files manually.\e[0m\n\n"
            exit 1
        fi
    elif [ $(extension_present "pdb" "$protein_dir") == false ]
    then
        print_err_color "Receptor structure files not detected in specified directory."
        exit 1
    fi
fi
    # endregion protein pdbqt present 

    # region check parameters 
if [[ $ALGO == "AD" ]]
then
        # region check dpf (AD only) 
    printf "\nParameters: checking dpf files\n"
    if [ $(extension_present "dpf" "$protein_dir") == true ]
    then
        export check_exist
        parallel -j 0 check_exist "{1.}.dpf" ::: $(ls $protein_dir/*.pdbqt)
    elif [ $(extension_present "dpf" "$protein_dir") == false ]
    then
        printf "Parameters: dpf files undetected, preparing dpf files\n"
        print_ask_head "Select AutoDock algorithm:"
        print_ask_body " 1: Lamarckian Genetic Algorithm (LGA)"
        print_ask_body " 2: Simulated Annealing (SA)"
        print_ask_body " 3: Local Search (LS)"
        print_ask_body " 4: Exit and manually prepare dpf files\n"
        while true
        do
            read AD_ALGO
            if ! [ -z $AD_ALGO ]
            then
                if [[ $AD_ALGO =~ ^[1-4]$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter a number between 1 to 4:"
        done
        if [[ $AD_ALGO -eq 4 ]]
        then
            printf "\n\n\e[93mPrepare the receptor .pdbqt files manually.\e[0m\n\n"
            exit 1
        fi
        print_ask_head "Enter number of runs: "
        while true
        do
            read num_runs
            if ! [ -z $num_runs ]
            then
                if [[ $num_runs =~ ^[0-9]+$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an integer: "
        done
    fi
        # endregion check dpf
            # return AD_ALGO, 1- LGA, 2- GA, 3- SA, 4- LS, 5- EXIT
            # return num_runs
        # region check gpf (AD only) 
    printf "\nParameters: checking gpf files\n"
    if [ $(extension_present "gpf" "$protein_dir") == true ]
    then
        export check_exist
        parallel -j 0 check_exist "{1.}.gpf" ::: $(ls $protein_dir/*.pdb)
    elif [ $(extension_present "gpf" "$protein_dir") == false ]
    then
        printf "Parameters: gpf files undetected, preparing gpf files\n"
        print_ask_head "Generate gpf (search grid):"
        print_ask_body " 1: Manually"
        print_ask_body " 2: Pocket prediction using AutoSite (developing..)"
        print_ask_body " 3: Specify residues\n"
        while true
        do
            read AUTO_SITE
            if ! [ -z $AUTO_SITE ]
            then
                if [[ $AUTO_SITE -eq 1 ]]
                then
                    printf "\n\n\e[93mKindly prepare the .gpf files manually.\e[0m\n\n"
                    exit 1
                elif [[ $AUTO_SITE -eq 2 ]] || [[ $AUTO_SITE -eq 3 ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter a number between 1 to 3:"
        done
    fi
        # endregion check gpf
            # return AUTO_SITE, 2- yes, 3- specify residue
            # return gpf_length if AUTO_SITE=2
elif [[ $ALGO == "VINA" ]]
then
        # region check config 
    printf "\nParameters: checking config files\n"
    if [ $(extension_present "txt" "$protein_dir") == true ]
    then
        export check_exist
        parallel -j 0 check_exist "{1.}.txt" ::: $(ls $protein_dir/*.pdbqt)
    elif [ $(extension_present "txt" "$protein_dir") == false ]
    then
        printf "Parameters: config files undetected, preparing config files\n"
        print_ask_head "Specify search grid:"
        print_ask_body " 1: Manually"
        print_ask_body " 2: Pocket prediction using AutoSite (developing..)"
        print_ask_body " 3: Specify residues\n"
        while true
        do
            read AUTO_SITE
            if ! [ -z $AUTO_SITE ]
            then
                if [[ $AUTO_SITE -eq 1 ]]
                then
                    printf "\n\n\e[93mKindly prepare the .gpf files manually.\e[0m\n\n"
                    exit 1
                elif [[ $AUTO_SITE -eq 2 ]] || [[ $AUTO_SITE -eq 3 ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter a number between 1 to 3:"
        done
        print_ask_head "Enter number of runs: "
        while true
        do
            read num_runs
            if ! [ -z $num_runs ]
            then
                if [[ $num_runs =~ ^[0-9]+$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an integer: "
        done
    fi
        # endregion check config
            # return num_runs
            # return AUTO_SITE, 2- yes
            # return gpf_length if AUTO_SITE=2
fi
        # region param autosite
if [[ $AUTO_SITE -eq 2 ]]
then
            # region get MGL2_dir
    print_ask_head "Enter directory of MGLTools2:"
    while true
    do
        read MGL2_dir
        if ! [ -z $MGL2_dir ]
        then
            if [ -e "$MGL2_dir" ]
            then
                break
            fi
        fi
        print_err_color "Directory $MGL2_dir not found:"
    done
            # endregion get MGL2_dir
    print_ask_head "AutoSite: search grid geometry:"
    print_ask_body " 1: Cube "
    print_ask_body " 2: Auto positioning and sizing \n"
    while true
    do
        read AUTO_SITE_SHAPE
        if ! [ -z $AUTO_SITE_SHAPE ]
        then
            if [[ $AUTO_SITE_SHAPE -eq 1 ]] || [[ $AUTO_SITE_SHAPE -eq 2 ]]
            then
                break
            fi
        fi
        print_err_color "Please enter a number between 1 to 2:"
    done
    if [[ $AUTO_SITE_SHAPE -eq 1 ]]
    then
        print_ask_head "Side length of cubic search space (angstrom):"
        while true
        do
            read gpf_length
            if ! [ -z $gpf_length ]
            then
                if [[ $gpf_length =~ ^[0-9]+([.][0-9]+)?$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an positive numerical value:"
        done
    elif [[ $AUTO_SITE_SHAPE -eq 2 ]]
    then
        :



    fi
elif [[ $AUTO_SITE -eq 3 ]]
then
    print_ask_head "Enter location of configuration file containing specified residues:"
    print_ask_body "(Receptor directory is $protein_dir)\n"
    while true
    do
        read resiconf_file
        if ! [ -z $resiconf_file ]
        then
            if [ -e "$resiconf_file" ]
            then
                break
            fi
        fi
        print_err_color "File $resiconf_file not found:"
    done

    print_ask_head "Expand search grid by: (in Angstrom) "
    while true
    do
        read gpf_length
        if ! [ -z $gpf_length ]
        then
            if [[ $gpf_length =~ ^[0-9]+([.][0-9]+)?$ ]]
            then
                break
            fi
        fi
        print_err_color "Please enter an positive numerical value:"
    done
fi
        # endregion param autosite
    # endregion check parameters 
# endregion RECEPTOR 

# region LIGAND 
cd $ligand_dir
    # region ligand pdbqt absent
printf "\nLigand: checking ligand directory\n"
if [ $(extension_present "pdbqt" "$ligand_dir") == false ]
then
    printf "Ligand: .pdbqt files undetected\n"

    # if input is smiles
    if [ $(extension_present "smi" "$ligand_dir") == true ]
    then
        printf "Ligand preparation: .smi files detected\n"
    fi
    # if input is sdf
    if [ $(extension_present "sdf" "$ligand_dir") == true ]
    then
        printf "Ligand preparation: .sdf files detected\n"
    fi
    # if input is mol2
    if [ $(extension_present "mol2" "$ligand_dir") == true ]
    then
        printf "Ligand preparation: .mol2 files detected\n"
    fi
    # if input is pdb
    if [ $(extension_present "pdb" "$ligand_dir") == true ]
    then
        printf "Ligand preparation: .pdb files detected\n"
    fi
    # ask 3D gen
    if [ $(extension_present "smi" "$ligand_dir") == true ] || [ $(extension_present "sdf" "$ligand_dir") == true ]
    then
        print_ask_head "Ligand preparation: Enter 3D structure conversion speed: "
        print_ask_body " 1: fastest"
        print_ask_body " 2: fast"
        print_ask_body " 3: medium (default)"
        print_ask_body " 4: better"
        print_ask_body " 5: best\n"
        while true
        do
            read bab_gen3d_speed
            if ! [ -z $bab_gen3d_speed ]
            then
                if [[ $bab_gen3d_speed =~ ^[1-5]$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an integer between 1 to 5:"
        done
    fi
    print_ask_head "Ligand preparation: Protonation: " 
    print_ask_body " 1: at physiological pH (pH 7.4)"
    print_ask_body " 2: explicit"
    print_ask_body " 3: do not add hydrogens\n"
    while true
    do
        read bab_addh
        if ! [ -z $bab_addh ]
        then
            if [[ $bab_addh =~ ^[1-3]$ ]]
            then
                break
            fi
        fi
        print_err_color "Please enter an integer between 1 to 3:"
    done
    print_ask_head "Ligand preparation: Conformer search method: " 
    print_ask_body " 1: systematic"
    print_ask_body " 2: random"
    print_ask_body " 3: weighted rotor search"
    print_ask_body " 4: MMFF94 force field"
    print_ask_body " 5: genetic algorithm"
    print_ask_body " 6: do not generate conformer\n"
    while true
    do
        read bab_conformer_ff
        if ! [ -z $bab_conformer_ff ]
        then
            if [[ $bab_conformer_ff =~ ^[1-6]$ ]]
            then
                break
            fi
        fi
        print_err_color "Please enter an integer between 1 to 6:"
    done
    if [[ $bab_conformer_ff == 5 ]]
    then
        print_ask_head "Ligand preparation: Conformer scoring function: " 
        print_ask_body " 1: rmsd"
        print_ask_body " 2: energy"
        while true
        do
            read bab_conformer_score
            if ! [ -z $bab_conformer_score ]
            then
                if [[ $bab_conformer_score =~ ^[1-2]$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter 1 or 2:"
        done
    fi
    if [[ $bab_conformer_ff =~ ^[1-5]$ ]]
    then
        print_ask_head "Ligand preparation: Number of conformers: "
        while true
        do
            read bab_conformer_no
            if ! [ -z $bab_conformer_no ]
            then
                if [[ $bab_conformer_no =~ ^[0-9]+$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an integer:"
        done
    fi
    print_ask_head "Ligand preparation: Energy minimization - Steepest descent: "
    print_ask_body "Ligand preparation: Steepest descent - Enter maximum number of steps: "
    print_ask_body "(to skip this step, enter 0) \n"
    while true
    do
        read bab_minim_sd_nsteps
        if ! [ -z $bab_minim_sd_nsteps ]
        then
            if [[ $bab_minim_sd_nsteps =~ ^[0-9]+$ ]]
            then
                break
            fi
        fi
        print_err_color "Please enter an integer:"
    done
    if ! [[ $bab_minim_sd_nsteps == 0 ]]
    then
        print_ask_head "Ligand preparation: Steepest descent - Enter convergence criteria: \n"
        while true
        do
            read bab_minim_sd_conv
            if ! [ -z $bab_minim_sd_conv ]
            then
                break
            fi
            print_err_color "Please enter an integer:"
        done
    fi
    print_ask_head "Ligand preparation: Energy minimization - Conjugate gradient: "
    print_ask_body "Ligand preparation: Conjugate gradient - Enter maximum number of steps: "
    print_ask_body "(to skip this step, enter 0) \n"
    while true
    do
        read bab_minim_cg_nsteps
        if ! [ -z $bab_minim_cg_nsteps ]
        then
            if [[ $bab_minim_cg_nsteps =~ ^[0-9]+$ ]]
            then
                break
            fi
        fi
        print_err_color "Please enter an integer:"
    done
    if ! [[ $bab_minim_cg_nsteps == 0 ]]
    then
        print_ask_head "Ligand preparation: Conjugate gradient - Enter convergence criteria: "
        while true
        do
            read bab_minim_cg_conv
            if ! [ -z $bab_minim_cg_conv ]
            then
                break
            fi
            print_err_color "Please enter an integer:"
        done
    fi
    if [[ $bab_minim_sd_nsteps != 0 ]] || [[ $bab_minim_cg_nsteps != 0 ]]
    then
        print_ask_head "Ligand preparation: Energy minimization force field: " 
        print_ask_body " 1: MMFF94"
        print_ask_body " 2: MMFF94s"
        print_ask_body " 3: Ghemical"
        print_ask_body " 4: Gaff"
        print_ask_body " 5: Uff"
        print_ask_body " 6: None\n"
        while true
        do
            read bab_minim_ff
            if ! [ -z $bab_minim_ff ]
            then
                if [[ $bab_minim_ff =~ ^[1-6]$ ]]
                then
                    break
                fi
            fi
            print_err_color "Please enter an integer between 1 to 6:"
        done
    fi
fi
    # endregion ligand pdbqt present
# endregion LIGAND

# region run

    # region export for CLI
export prepare_receptor
export AD_ALGO
export AUTO_SITE
export -f extension_present
export ligand_dir
export ALGO
    # endregion export for CLI
    # region prepare pdbqt 

if [[ $prepare_receptor == "y" || $prepare_receptor == "Y" ]]
then
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 1
    printf "Receptor: Preparing receptor pdbqt files from pdb\n"
    export -f prepare_receptor4_py
    parallel -j 0 --eta prepare_receptor4_py "{1/.}" "{2}" "$protein_dir" ::: $(ls $protein_dir/*.pdb) ::: $MGL_dir
fi
    # endregion prepare pdbqt

    # region generate dpf 

if [[ $AD_ALGO =~ ^[1-3]$ ]]
then
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 2
    printf "Parameters: generating dpf files\n"
    export -f generate_dpf
    parallel -j 0 --eta generate_dpf "$AD_ALGO" "{2}" "$protein_dir" "{1/.}" ::: $(ls $protein_dir/*.pdbqt) ::: "$num_runs"
fi
    # endregion generate dpf

    # region generate gpf 
if [[ $AUTO_SITE -eq 2 ]]
then
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 3
    printf "Pocket Prediction: autosite running\n"
    export -f autosite_pred
    parallel -j 0 --eta autosite_pred "$MGL2_dir" "$protein_dir" "{1/.}" "{2}" "$ALGO" "{3}" "{4}" ::: $(ls $protein_dir/*.pdbqt) ::: "$gpf_length" ::: "$num_runs" ::: "$AUTO_SITE_SHAPE"
elif [[ $AUTO_SITE -eq 3 ]]
then
    bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 4
    printf "Parameters: generating gpf files"
    pdbqt_xyz "$protein_dir" "$resiconf_file" "$ALGO" "$gpf_length" "$num_runs"
fi
    # endregion generate gpf

    # region prepareligand
bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
if [ $(extension_present "pdbqt" "$ligand_dir") == false ]
then
    export -f obabel_addh
    export -f obabel_gen3d
    export -f obabel_conformer
    export -f obabel_minimization
    export -f obabel_mol2_pdb
    export -f prepare_ligand4_py

    # if input is smiles
    if [[ $(ls $ligand_dir/*.smi | wc -l) -ge 1 ]]
    then
        if [[ $(ls $ligand_dir/*.smi | wc -l) == 1 ]]
        then
            # split smi
            smi_file="$(ls $ligand_dir/*.smi)"
            mkdir -p $ligand_dir/LIGAND_SMI
            tot_smi=$(cat $smi_file | wc -l)
            n=1
            while read line || [ -n "$line" ]
            do
                smi_suffix=$(printf "%0*d" "${#tot_smi}" "$n")
                printf "$line" | awk '{print $1}' > $ligand_dir/LIGAND_SMI/LG_$smi_suffix.smi
                n=$((n+1))
            done < ${smi_file}

        elif [[ $(ls $ligand_dir/*.smi | wc -l) -gt 1 ]]
        then
            mkdir -p $ligand_dir/LIGAND_SMI
            mv $ligand_dir/*.smi $ligand_dir/LIGAND_SMI
        fi
    fi
    # if input is sdf
    if [ $(extension_present "sdf" "$ligand_dir") == true ]
    then
        mkdir -p $ligand_dir/LIGAND_SDF
        mv $ligand_dir/*.sdf $ligand_dir/LIGAND_SDF
    fi
    # if input is mol2
    if [ $(extension_present "mol2" "$ligand_dir") == true ]
    then
        mkdir -p $ligand_dir/LIGAND_MOL2
        mv $ligand_dir/*.mol2 $ligand_dir/LIGAND_MOL2
    fi
    # if input is pdb
    if [ $(extension_present "pdb" "$ligand_dir") == true ]
    then
        mkdir -p $ligand_dir/LIGAND_PDB
        mv $ligand_dir/*.pdb $ligand_dir/LIGAND_PDB
    fi

    if [ -e "$ligand_dir/LIGAND_SMI" ]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: 3D structures generation\n"
        if [[ $(ls $ligand_dir/LIGAND_SMI/*.smi | wc -l) -ge 1 ]]
        then
            parallel -j 0 --eta obabel_gen3d "$bab_gen3d_speed" "{1}" ::: $(ls $ligand_dir/LIGAND_SMI/*.smi)
        fi
        mkdir -p $ligand_dir/LIGAND_PDB
        mv $ligand_dir/LIGAND_SMI/*.pdb $ligand_dir/LIGAND_PDB
    elif [ -e "$ligand_dir/LIGAND_SDF" ]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: 3D structures generation\n"
        parallel -j 0 --eta obabel_gen3d "$bab_gen3d_speed" "{1}" ::: $(ls $ligand_dir/LIGAND_SDF/*.sdf)
        mkdir -p $ligand_dir/LIGAND_PDB
        mv $ligand_dir/LIGAND_SDF/*.pdb $ligand_dir/LIGAND_PDB
    elif [ -e "$ligand_dir/LIGAND_MOL2" ]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: conversion from mol2 into pdb\n"
        parallel -j 0 --eta obabel_mol2_pdb "{1}" ::: $(ls $ligand_dir/LIGAND_MOL2/*.mol2)
        mkdir -p $ligand_dir/LIGAND_PDB
        mv $ligand_dir/LIGAND_MOL2/*.pdb $ligand_dir/LIGAND_PDB
    fi

    if ! [[ $bab_addh == 3 ]]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: protonation "
        if [[ $bab_addh == 1 ]]
        then
            printf "at pH of 7.4 \n"
        elif [[ $bab_addh == 2 ]]
        then
            printf "explicitly \n"
        fi
        parallel -j 0 --eta obabel_addh "$bab_addh" "{1}" ::: $(ls $ligand_dir/LIGAND_PDB/*.pdb)
    fi

    if [[ $bab_conformer_ff =~ ^[1-5]$ ]]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: conformer generation \n"
        parallel -j 0 --eta obabel_conformer "$bab_conformer_ff" ::: "$bab_conformer_score" ::: "$bab_conformer_no" ::: $(ls $ligand_dir/LIGAND_PDB/*.pdb)
    fi

    if [[ $bab_minim_sd_nsteps -gt 0 ]] || [[ $bab_minim_cg_nsteps -gt 0 ]]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: energy minimization \n"
        parallel -j 0 --eta obabel_minimization "$bab_minim_sd_nsteps" ::: "$bab_minim_sd_conv" ::: "$bab_minim_cg_nsteps" ::: "$bab_minim_cg_conv" ::: "$bab_minim_ff" ::: $(ls $ligand_dir/LIGAND_PDB/*.pdb)
    fi

    if [ -e "$ligand_dir/LIGAND_PDB" ]
    then
        bash $SCRIPTPATH/Modules/vs_printCLI.sh -c 5
        printf "\nLigand preparation: pdbqt generation \n"
        parallel -j 0 --eta prepare_ligand4_py "$MGL_dir" "{1}" ::: $(ls $ligand_dir/LIGAND_PDB/*.pdb)
        mv $ligand_dir/LIGAND_PDB/*.pdbqt $ligand_dir
    fi
fi

printf "\n\nLigand preparation: pdbqt checking \n"
export -f ligand_check_atm
parallel -j 0 ligand_check_atm "$ligand_dir" "{1}" ::: $(ls $ligand_dir/*.pdbqt)

    # endregion prepare ligand
# endregion run

# region unset functions
unset -f check_exist
unset -f check_spaces
unset -f prepare_receptor4_py
unset -f generate_dpf
unset -f autosite_pred
unset -f obabel_addh
unset -f obabel_gen3d
unset -f obabel_conformer
unset -f obabel_minimization
unset -f obabel_mol2_pdb
unset -f prepare_ligand4_py
unset -f ligand_check_atm
# endregion unset functions

# your turn vs_pipeline
export ligand_dir
export result_dir
export protein_dir
export working_dir
export number_jobs
export MGL_dir
export ALGO
export SCRIPTPATH
bash $SCRIPTPATH/Modules/vs_pipeline.sh