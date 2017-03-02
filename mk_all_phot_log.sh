#!/bin/bash   
#  commentary 
## test module

if [ $# -lt 1 ]
then
   echo "                                                                   "
   echo -e "\033[1m mk_all_phot_log.sh\033[0m                               "
   echo " Program  prepares a log for further photometric                   "
   echo " analysis. Program checks status of all files listed in individual "
   echo " logs and modifies them in order to exclude nonexistent files.     "
   echo "                                                                   "
   echo " You should use this program in a directory that contains at least "
   echo " one data directory, i.e.: Berkeley59_15_03_03_'PATTERN',...       "
   echo " where 'PATTERN' is set inside your 'master.conf' file.            "
   echo " All prepared data should come from pipeline 'ccdphot'.            "
   echo " Program looks for a given pattern in the directory structure.     "
   echo "                                                                   "
   echo " Requirements: bash, ccdphot (Z. Kolaczkowski, IA UWr)             "
   echo "                                                                   "
   echo " Usage:                                                            "
   echo " ./mk_all_phot_log.sh <object_name> [-n] [-p PATTERN] [-h]         "
   echo " Options: -n        : no pattern - search in all folders           "
   echo "          -p        : set searching PATTERN (default: 'RESULTS')   "
   echo "          -h        : see help                                     "
   echo "                                                                   "
   echo " Version: 2017.03.02, (PM: http://github.com/astromiki)            "
   echo "                                                                   "
   exit
fi

# importing existing configuration
source get_conf.sh
get_variables

# setting basic variables
object_name="$1"
all_phot_log_filename=${object_name}"_all_phot.log"

# determining rest of the arguments
if [[ $# -gt 1 ]]
then
   if [[ $2 == "-n" ]]
   then
      pattern=""
   elif [[ $2 == "-p" ]]
   then
      if [[ $# != 3 ]]
      then
         echo "# ERROR! No pattern entered!"
         exit 1
      else
         pattern="$3"
      fi
   elif [[ $2 == "-h" ]]
   then
      $(basename "$0")
      exit
   else
      echo "# ERROR! Wrong option entered!"
      exit 1
   fi
else
   pattern=$PATTERN
fi
## echo "OBJECT NAME: "${object_name} "PATTERN: "${pattern} "OUTPUT: "${all_phot_log_filename} && exit 0

# starting the process
STARTTIME=$(date +%s)

# looking for directories matching set pattern
if [[ ${pattern} == "" ]]
then
   echo "./" > .dirs
else
   ls -d *${pattern}*/ > .dirs
fi
nr_dirs=`wc -l .dirs | awk '{print $1}'`
echo ""
echo " $ Number of directories : "$nr_dirs

# declaring temporary variables
declared=0
found=0
unique_log=0
unique_dir=0

echo -n -e " > Wait..."

# looking for files with complete photometry (all_phot)
i=1
while [[ $i -le $nr_dirs ]] ; do
   dir=`cat .dirs | sed $i'q;d' | awk '{print $1}'` 
   phot_log=`ls ${dir}*phot.log`
   cat ${phot_log} | sed s/.fits/.all_phot/g | awk -v d=$dir 'NR > 2 {print d$1}' > .subs_in_log-${i}
   sub_lines=`wc -l .subs_in_log-${i} | awk '{print $1}'`
   declared=$[declared + sub_lines]
   cat ${phot_log} | sed s/.fits/.all_phot/g | sed 's/^ *//' | awk -v d=$dir 'NR > 2 {print d $0}' > .temp
   ls ${dir}*all_phot > .subs_in_dir-${i}
   sub_files=`wc -l .subs_in_dir-${i} | awk '{print $1}'`
   found=$[found + sub_files]
   comm -12 .subs_in_log-${i} .subs_in_dir-${i} > .common
   comm -23 .subs_in_log-${i} .subs_in_dir-${i} > .log-unique-${i}
   if [ $(wc -l .log-unique-${i} | awk '{print $1}') -ne 0 ]
   then
      unique_log=$[unique_log + $(wc -l .log-unique-${i} | awk '{print $1}')]
   fi
   comm -23 .subs_in_dir-${i} .subs_in_log-${i} > .dir-unique-${i}
   if [ $(wc -l .dir-unique-${i} | awk '{print $1}') -ne 0 ]
   then
      unique_dir=$[unique_dir + $(wc -l .dir-unique-${i} | awk '{print $1}')]
   fi
   nr_comm=`wc -l .common | awk '{print $1}'`
   j=1
   while [[ $j -le $nr_comm ]] ; do
      common=`cat .common | sed $j'q;d' | awk '{print $1}'`
      cat .temp | grep ${common} >> .phot_log-${i}
      j=$[j + 1]
   done
   rm .subs_in_log-${i} .subs_in_dir-${i} .common .temp
   i=$[i + 1]
done

# creating list of files existing only in logs and only in directiories
cat .phot_log-* | sed s/-sub//g > ${all_phot_log_filename}
cat ${all_phot_log_filename} | sort -n -k 14 > a && mv a ${all_phot_log_filename}
cat .log-unique* > LOG-unique
cat .dir-unique* > DIR-unique

# communicating with user
rm .dirs .phot_log-* .log-unique-* .dir-unique-*
ENDTIME=$(date +%s)
echo " "
echo " $ Frames declared       : "${declared}
echo " $ Frames found          : "${found}
if [ $unique_log -ne 0 ]
then
   echo " # Warning! See file     : LOG-unique!"
else
   rm LOG-unique
fi
if [ $unique_dir -ne 0 ]
then
   echo " # Warning! See file     : DIR-unique!"
else
   rm DIR-unique
fi
echo " $ Log written to file   : "${all_phot_log_filename}"!"
echo " $ Job DONE in           : "$(($ENDTIME - $STARTTIME))"s!"
echo ""
exit 0
