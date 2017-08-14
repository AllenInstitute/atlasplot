#!/usr/bin/env bash
# creates a list of the function and which file they're in for an R package
# it should be run in the same directory as the package root directory
#
# eg.
# $ bash create_function_list.sh atlasplot

r_pkg="$1/R"  # first parameter is the package name

# ensure R directory exists; exit if not
if [[ ! -d "$r_pkg"  ]]
then
    echo "The specified file is not an R package or there is no `R` directory"
    exit 1
fi

# create file; ensure it is empty to avoid appending to old func list
func_list="$1/function_list.md"
touch "$func_list"
echo "# Functions By File" > "$func_list"
for i in $( ls "$r_pkg" )
do
    # make sure that we're not added a file for sysdata.rda
    if [[  "${i#*.}" = "R" ]] || [[ "${i#*.}" = "r" ]]
    then
        echo "#### $r_pkg/$i" >> "$func_list"
        funcs=$( grep "<-\s*function" "$r_pkg/$i" | sed "s/<-.*$//g" | sed '/^#/d' | sed 's/^/* /g' )
        for f in "$funcs"
        do
            echo "$f" >> "$func_list"
        done
        printf "\n" >> "$func_list"
    fi
done
