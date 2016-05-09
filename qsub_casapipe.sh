#!/bin/bash

# set script directories
shopt -s expand_aliases
alias CASA=/csr/mwa/casapy/casapy-41.0.24668-001-64b/casapy
MITPIPE_DIR=/csr/mwa/pipeline
MWATOOLS_DIR=/csr/mwa/python/mwa_git/mwatools_setup/bin
MITBIN_DIR=/csr/mwa/bin

# get command line arguments
program=`basename $0`
config=$1
obsfile=$2

# print usage
if [ $# -ne 2 ]; then
    echo "Usage:" ${program} "<config.txt> <obsfile.list>"
    echo "    Example config file: /csr/mwa/pipeline_dev/example_config.txt"
    exit
fi

# verify arguments
if [ ! -f ${config} -o ! -f ${obsfile} ]; then
    echo "One or more of your inputs do not exist. Quit script."
    exit
fi

# parse config file
this_dir=`pwd | sed -e 's/\/nfs//g' | awk '{print "/nfs" $1}'`
source ${this_dir}/${config}

# create scratch space in ram
if [ $SHM = "True" ]; then
    scratch=/dev/shm/casapipe
else
    scratch=${this_dir}
fi

if [ ! -d "${scratch}" ]; then
    mkdir -v ${scratch}
fi

# iterate over each dataset
while read -r obsid
do

    # create directory for obsid
    if [ ! -d "${obsid}" ]; then
	mkdir -v ${obsid}
    fi

#    # move uvfits to obsid directory
#    if [ -f "${obsid}.uvfits" ]; then
#	mv -v ${obsid}.uvfits ${obsid}/${obsid}.uvfits
#    fi

    # process in scratch space
    cd ${scratch}; pwd
    if [ ! -d "${obsid}" ]; then
	mkdir -v ${obsid}
    fi

    # copy over uvfits
    if [ $PREPDATA = "True" ]; then
	cp -v ${this_dir}/${obsid}.uvfits ${scratch}/${obsid}
    fi

    # copy over metafits
    uvpath=`readlink ${this_dir}/${obsid}.uvfits`
    metadir=$(dirname ${uvpath})
    cp -v ${metadir}/${obsid}.metafits ${scratch}/${obsid}

    # copy config file
    timestamp=`date +%Y%m%d%H%M%S`
    obsconfig=${config}_${timestamp}
    cp -v ${this_dir}/${config} ${scratch}/${obsid}/${obsconfig}

    # create logfile directory
    cd ${obsid}; pwd
    if [ ! -d logfiles ]; then
	mkdir -v logfiles
    fi

    # import data
    if [ $PREPDATA = "True" ]; then

	echo "Import" ${obsid}".uvfits into CASA"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "prep" ${obsid} > logfiles/dataprep_${timestamp}.log
	mv -v casapy*.log logfiles/dataprep_${timestamp}.casalog
	rm ${scratch}/${obsid}/${obsid}.uvfits  # remove uvfits after ms is imported

	echo "Fix measurement set"
	${MITBIN_DIR}/fixmwams ${obsid}.ms ${obsid}.metafits
	${MWATOOLS_DIR}/fixms.sh ${obsid}.ms ${MWATOOLS_DIR}

	echo "Plot beam map"
	len=${#obsid}
	if [ $len -eq 10 ]; then
	    time python ${MWATOOLS_DIR}/get_observation_info.py -g ${obsid} -i
	else
	    time python ${MWATOOLS_DIR}/get_observation_info.py -f ${obsid} -i
	fi
	mv *MHz.png beam_${obsid}.png
    fi

    # change phase center
    if [ $CHGCENTRE = "True" ]; then

	echo Change phase center to $PHASERA and $PHASEDEC
	${MITBIN_DIR}/chgcentre ${obsid}.ms $PHASERA $PHASEDEC

    fi
    
    # run flag statistics
    if [ $ACAVG = "True" ]; then

	# UNTESTED
	echo "Run flag statistics"
	time casa_acflag ${obsid}
	cp ${MITPIPE_DIR}/show_ac.macro .
	ssh grs1915 "cd ${this_dir}; chmod +x make_sm_plots.go; ./make_sm_plots.go"

    fi

    # flag additional data
    if [ $FLAGGING = "True" ]; then

	echo "Flag additional data"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "flag" ${obsid} > logfiles/flagdata_${timestamp}.log
	mv -v casapy*.log logfiles/flagdata_${timestamp}.casalog

    fi

    # run delaycal
    if [ $DELAYCAL = "True" ]; then
	
	echo "Run delaycal"
	time CASA --nologger -c ${MITPIPE_DIR}/bash_delaycal.py --cal=${CALFILE} --src=${CALSRC} --cat=${CATALOG} --nsrcs=${NSRCS} --fov=${FOV} --practice --overwrite --refant=${REFANT} ${obsid}.ms > logfiles/delaycal_${timestamp}.log
	mv -v casapy*.log logfiles/delaycal_${timestamp}.casalog
	rm -r *.cl.fits *.cl.im

    fi

    # run imagecal
    if [ $IMAGECAL = "True" ]; then

	echo "Run image calibration"
	echo "Use image " ${IMCALFILE}.fits
	cp ${MITPIPE_DIR}/${IMCALFILE}.fits .
	python ${MITPIPE_DIR}/update_header.py -o ${obsid} -i ${IMCALFILE}.fits
	#delays=`python ${MITPIPE_DIR}/get_beam_delays.py ${obsid}`
	delays=`python ${MITPIPE_DIR}/get_beam_delays_metafits.py ${obsid}.metafits`
	python ${MWATOOLS_DIR}/make_beam.py -f ${IMCALFILE}.fits -d ${delays} -v
	python ${MWATOOLS_DIR}/imarith.py -i ${IMCALFILE}_beamXX.fits,${IMCALFILE}_beamYY.fits -o ${IMCALFILE}_beamI.fits -e "(im0+im1)/2" -v
	python ${MWATOOLS_DIR}/imarith.py -i ${IMCALFILE}.fits,${IMCALFILE}_beamI.fits -o ${IMCALFILE}_cal.fits -e "im0*im1" -v	
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "imagecal" ${obsid} > logfiles/imagecal_${timestamp}.log
	mv -v casapy*.log logfiles/imagecal_${timestamp}.casalog
	rm -rv ${IMCALFILE}.fits ${IMCALFILE}_beamXX.fits ${IMCALFILE}_beamYY.fits imcal.model

    fi

    # generate calibration solutions
    if [ $GENCAL = "True" ]; then

	echo "Generate calibration solutions"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "gencal" ${obsid} > logfiles/gencal_${timestamp}.log
	mv -v casapy*.log logfiles/gencal_${timestamp}.casalog

    fi

    # apply calibration solutions
    if [ $APPLYCAL = "True" ]; then

	echo "Apply calibration solutions"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "applycal" ${obsid} > logfiles/applycal_${timestamp}.log
	mv -v casapy*.log logfiles/applycal_${timestamp}.casalog

    fi

    # apply FHD calibration solutions
    if [ $FHDCAL = "True" ]; then

	echo "Apply FHD .npz calibration solutions"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "fhdcal" ${obsid} > logfiles/fhdcal_${timestamp}.log
	mv -v casapy*.log logfiles/fhdcal_${timestamp}.casalog

    fi

    # apply self calibration
    if [ $SELFCAL = "True" ]; then

	echo "Apply self calibration with" $NLOOPS "loops"
	if [ $WSCLEAN = "True" ]; then
	    for (( loop = 1; loop <= $NLOOPS; loop++ ))
	    do
		echo "Use wsclean loop" $loop
		time wsclean -name selfcal_xx_$loop -pol xx -size $WIMSIZE $WIMSIZE -scale $WCELL -niter $WNITERSELF -threshold $WTHRESHOLDSELF -mgain $WMGAIN -weight $WWEIGHTING -datacolumn CORRECTED_DATA -absmem $WSMEM -j $NCORES ${obsid}.ms
		time wsclean -name selfcal_yy_$loop -pol yy -size $WIMSIZE $WIMSIZE -scale $WCELL -niter $WNITERSELF -threshold $WTHRESHOLDSELF -mgain $WMGAIN -weight $WWEIGHTING -datacolumn CORRECTED_DATA -absmem $WSMEM -j $NCORES ${obsid}.ms
		time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py -l $loop ${obsconfig} "selfcal_wsclean" ${obsid} > logfiles/selfcal_${timestamp}.log
		mv -v casapy*.log logfiles/selfcal_${timestamp}.casalog
	    done
	else
	    echo "Use CASA clean"
	    time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "selfcal" ${obsid} > logfiles/selfcal_${timestamp}.log
            mv -v casapy*.log logfiles/selfcal_${timestamp}.casalog
	fi

    fi

    # export calibrated uvfits
    if [ $EXPORTUVFITS = "True" ]; then

	echo "Export calibrated uvfits"
	time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "exportuvfits" ${obsid} > logfiles/export_${timestamp}.log
	mv -v casapy*.log logfiles/export_${timestamp}.casalog

    fi

    # image
    if [ $CLEAN = "True" ]; then

	echo "Image"
	if [ $WSCLEAN = "True" ]; then

	    if [ $ZENITH = "True" ]; then
		echo "Rephase to zenith"
		time chgcentre -minw ${obsid}.ms
	    fi

	    echo "Use wsclean"
            time wsclean -name ${obsid} -size $WIMSIZE $WIMSIZE -scale $WCELL -niter $WNITER -threshold $WTHRESHOLD -gain $WGAIN -mgain $WMGAIN -gridmode $WGRIDMODE -pol xx,yy -weight $WWEIGHTING -datacolumn CORRECTED_DATA -absmem $WSMEM -j $NCORES -channelrange $STARTCHAN $ENDCHAN ${obsid}.ms
            mv -v ${obsid}-XX-image.fits clean_${IMAGEMODE}_${obsid}_XX.fits
            mv -v ${obsid}-YY-image.fits clean_${IMAGEMODE}_${obsid}_YY.fits
	    rm -v *-tmp.fits  # wsclean things?
	else
	    echo "Use CASA clean"
            time CASA --nologger -c ${MITPIPE_DIR}/casapipe.py ${obsconfig} "clean" ${obsid} > logfiles/imaging_${timestamp}.log
            mv -v casapy*.log logfiles/imaging_${timestamp}.casalog
	fi

    fi

    # correct for primary beam
    if [ $PBCOR = "True" ]; then

	echo "Correct for primary beam"

	python ${MITPIPE_DIR}/obsdate_expfix.py -i clean_${IMAGEMODE}_${obsid}_XX.fits -m ${obsid}.metafits -v -e 
	python ${MITPIPE_DIR}/obsdate_expfix.py -i clean_${IMAGEMODE}_${obsid}_YY.fits -m ${obsid}.metafits -v -e 

	if [ $ANALYTIC = "True" ]; then
	    echo "Use analytic beam model"
	    time python ${MWATOOLS_DIR}/make_beam.py --analytic -f clean_${IMAGEMODE}_${obsid}_XX.fits -m ${obsid}.metafits -v
	    bmmodel=analytic
	else
	    echo "Use curtin numerical beam model"
	    time python ${MWATOOLS_DIR}/make_beam.py -f clean_${IMAGEMODE}_${obsid}_XX.fits -m ${obsid}.metafits -v
	    bmmodel=curtin
	fi
	
	# set up file names because we don't want to overwrite different beam models
	bmxx=clean_${IMAGEMODE}_${obsid}_XX_beam_${bmmodel}.fits
        bmyy=clean_${IMAGEMODE}_${obsid}_YY_beam_${bmmodel}.fits
	iim=pbcor_clean_${IMAGEMODE}_${obsid}_I_${bmmodel}.fits
	qim=pbcor_clean_${IMAGEMODE}_${obsid}_Q_${bmmodel}.fits
	
	# make I/Q images
	mv -v clean_${IMAGEMODE}_${obsid}_XX_beamXX.fits ${bmxx}
	mv -v clean_${IMAGEMODE}_${obsid}_XX_beamYY.fits ${bmyy}
	time python ${MWATOOLS_DIR}/imarith.py -i clean_${IMAGEMODE}_${obsid}_XX.fits,clean_${IMAGEMODE}_${obsid}_YY.fits,${bmxx},${bmyy} -o ${iim} -e "(im0*im2+im1*im3)/(im2*im2+im3*im3)" -v
	time python ${MWATOOLS_DIR}/imarith.py -i clean_${IMAGEMODE}_${obsid}_XX.fits,clean_${IMAGEMODE}_${obsid}_YY.fits,${bmxx},${bmyy} -o ${qim} -e "(im0*im2-im1*im3)/(im2*im2+im3*im3)" -v
	
    fi

    # run source extractor
    if [ $SEXT = "True" ]; then

	echo "Run source extractor aegean"
	ls *I*.fits | sed -e 's/.fits//g' > images.list
	while read -r image
	do

	    echo "Save image RMS"
	    time python ${MITBIN_DIR}/aegean.py --cores=$NCORES --save_background ${image}.fits
	    echo "Mask large RMS values"
	    time python ${MITPIPE_DIR}/mask_large_rms.py --maxrms $MAXRMS
	    echo "Find sources"
	    time python ${MITBIN_DIR}/aegean.py --cores=$NCORES --nonegative --outfile=${image}.aeg_cat --rmsin=aegean-rms-masked.fits --bkgin=aegean-background.fits ${image}.fits
            echo "Plot RMS"
            time python ${MITPIPE_DIR}/imstat.py -i ${image}.fits -b clean_${IMAGEMODE}_${obsid}_XX_beam.fits

	done < images.list

    fi

    # clean up
    rm *.last ipython*.log
    if [ $ARCHIVECLEAN = "True" ]; then
	rm -r *.ms*
    fi

    # return to parent directory
    if [ $SHM = "True" ]; then
	time cp -r ${scratch}/${obsid}/* ${this_dir}/${obsid}
	rm -rf ${scratch}/${obsid}
    fi
    cd ${this_dir}; pwd


done < "${obsfile}"

if [ $SHM = "True" ]; then
    rmdir ${scratch}  # clean up
fi
