#!/bin/bash
SCRIPTLOC="/scratch/st-mjju-1/tffnytse/Codev1/matPdoct/"

LOC="/scratch/st-mjju-1/tffnytse/"
PROJECT='QUEUE'
##PROJECT='RP-EN'

DISP_ROI1=15
DISP_ROI2=500

##ADAPT=0


echo Project location $LOC$PROJECT, running with dispROI $DISP_ROI1 - $DISP_ROI2 $'\n'

cd $LOC$PROJECT

FOLDER=($(ls -d */))

echo $'\n' "List of folders: ${FOLDER[@]}" $'\n'

################################### Loop through all the folders ######################################
for F in "${FOLDER[@]}"
do
	echo Current folder: $F
	cd $F	

	################## Check for sub-directories ################
	# EYE=($(ls -d */))
	shopt -s nullglob
	EYE=(*/)
	shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
	
	if (( ${#EYE[@]} == 0 )); then
		echo No subdirectories found
		DATALOC="${LOC}${PROJECT}/${F}"

	else
		################# Loop through all sub-directories #################
		for E in "${EYE[@]}"
		do
			DATALOC="${LOC}${PROJECT}/${F}${E}"

			# Check if there is even a unp file in this busdirectory
			RETURN=$(pwd)
			cd $DATALOC
			shopt -s nullglob
			UNP=(*.unp)
			shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
			if (( ${#UNP[@]} == 0)); then # There are  no .unp files at this level, skip to next iteration
				echo "          No .unp files at $DATALOC"
				cd $RETURN
				continue
			fi
			cd $RETURN
			
			# run the job if there are .unp files
			echo "          $DATALOC"
			JOB="${PROJECT}-${F::-1}-${E::-1}"
			echo $JOB
			cd $SCRIPTLOC
			qsub -v DATALOC="${DATALOC}",dROI1=${DISP_ROI1},dROI2=${DISP_ROI2},JOB="${JOB}" run_salmon_array.pbs
			cd $RETURN
			sleep 2

		done
	fi
	echo Exiting location $LOC$PROJECT/$F $'\n'$'\n'
	cd ..
done
