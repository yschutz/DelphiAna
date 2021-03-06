#! /bin/bash
help()
{
	echo "usage ana.sh  h    = help"
	echo "             ch ps folder  = make control histograms and compare with hbook histograms for selected particles"
	echo "             cd ps folder  = make control histograms and nmakes relative difference with hbook histograms for selected particles"
	echo "             co ps folder  = make correlations of selected particles" 
	echo "                ps = a(all), c(charged), h(hadrons), ch(charged hadrons)"
	echo "             sh folder     = make single histogram"
	exit
}
if [ "$#" -eq 0 ]
then 
	help
	exit
else
	CPATH=`pwd`
	cd $CPATH/library	
	path=$CPATH/../data/$3
	eventfile=$path/event.root
	echo `pwd` $eventfile
	if [ ! -f $eventfile ]
	then 
		echo $eventfile not found
		exit
	fi
	export DATAPATH=$path
	export OPT1=`echo $1| awk '{print toupper($0)}'`
	if [ "$1" = "h" ]
	then
		help
	elif [ "$1" = "ch" -o "$1" = "cd" -o "$1" = "co" ]; then
		if [ "$#" -eq 1 ]; then 
			echo "!!! Missing argument"
			help
		elif [ "$2" = "a" -o "$2" = "c" -o "$2" = "h" -o "$2" = "ch" ]; then	
			export OPT2=`echo $2| awk '{print toupper($0)}'`
		fi
	elif [ "$1" = "sh" ]; then 
		echo "Choose among: "
		echo "             1->Thrust"
		echo "             2->Multiplicity"
        read -p "Enter parameter: " parameter
		export PARAM=$parameter
		echo $PARAM
	else
		echo "!!! ana.sh $1: unknown option"
		help
		exit
	fi	
	if [ "$1" = "co" ]; then 
		read -p "Enter [0,low], [high, infinite] multipicity and mixing level (RET = 20,30,10)" str
		IFS=,
		ary=($str)
		lowm=${ary[0]}
		highm=${ary[1]}
		mix=${ary[2]}
		if [ "$lowm" = "" ]; then 
			lowm=20		
		fi
		if [ "$highm" = "" ]; then 
			highm=30
		fi
		if [ "$mix" = "" ]; then 
			mix=10
		fi
		export HIGHM=$highm
		export LOWM=$lowm
		export MIX=$mix
	fi
fi
LIBDIR=$CPATH/library
if [ -d "$LIBDIR" ]; then
	cd $CPATH/library
	make
    retVal=$?
 	if [ $retVal -ne 0 ]; then
    	echo "Compilation Error!"
		cd ..
		exit $retVal
	fi
else
	echo "ERROR: $LIBDIR not found!"
	exit 1
fi 
echo `pwd`
root -l $CPATH/DAna.cxx++