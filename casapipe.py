# 2013/10/29

####################################################################################################

# import statements
import sys, optparse
import numpy as np

###################################################################################################

# create command line arguments
parser = optparse.OptionParser()
parser.set_usage("Usage: casapy [--nologger] -c casapipe.py [options] <config.txt> <mode> <obsid>")
parser.add_option("-l", "--loop", default=1, type="int",
                  help="Number of selfcal loops. [default: %default]")

# parse command line arguments
casa_id = sys.argv.index('-c') + 2  # remove CASA options
(opts, args) = parser.parse_args(sys.argv[casa_id:])
obsid, mode, configfile = args[-1], args[-2], args[-3]


####################################################################################################

def read_config(filename):
    """
    Read the pipeline config file <filename>.
    Return the config parameters parsed into a dictionary.

    """

    # read pipeline parameters
    print "=" * 50
    print "Read configuration file " + filename
    config = {}
    with open(filename, "r") as cfile:
        for line in cfile:
            if line[0] != "#":
                params = line.split("#")[0].split("=")  # split off comments
                config[params[0]] = eval(params[1])
    return config


####################################################################################################

def import_data(obsid):
    """
    Import <obsid>.uvfits into <obsid>.ms.
    Output an observation summary file <obsid>.obs.

    """

    # parse filenames
    myuvfits = obsid + ".uvfits"
    myvis = obsid + ".ms"
    mylistfile = obsid + ".obs"

    # import data into CASA
    print "=" * 50
    print "Import " + myuvfits + " into " + myvis
    importuvfits(fitsfile=myuvfits, vis=myvis)
    listobs(vis=myvis, listfile=mylistfile, verbose=True, selectdata=False)
    flagdata(vis=myvis, mode="summary")


####################################################################################################

def prep_rfi_files(obsid):
    """
    Prepare the RFI analysis files for <obsid>.ms
    All output files have <obsid> in their filenames.

    """

    ## begin code that prepares files for RFI analysis
    ## modified from Sam Simmons' code

    # parse filenames
    myvis = obsid + ".ms"
    acvis = obsid + "_ac.ms"

    # split out ms that only contains autocorrelations
    print "=" * 50
    print "Prepare " + myvis + " for RFI analysis"
    split(vis=myvis, outputvis=acvis, datacolumn="data", antenna="*&&&", keepflags=True)

    # make column vectors with times and antenna numbers for each entry
    print "Build xxfull/yyfull files"
    tb.open(acvis)
    ant = tb.getcol(columnname="ANTENNA1")
    length = len(ant)
    ant = ant.reshape(length, 1)
    dumbtime = tb.getcol(columnname="TIME")
    dumbtime = dumbtime.reshape(length, 1)

    # initialize arrays for all data pols
    findatyy = np.zeros((length, 768))
    findatxx = np.zeros((length, 768))
    findatxy = np.zeros((length, 768))
    findatyx = np.zeros((length, 768))

    # get data from table row by row 
    # (1 row is the data for 1 ac and one time, ac number cycles fastest)
    for z in range(length):

        # data is a 4x768 where 4 rows are 4 pols
        acdata = tb.getcell(columnname="DATA", rownr=z)
        xxmagsquare = np.absolute(acdata[0,:])
        yymagsquare = np.absolute(acdata[3,:])
        xxmagsquare1 = xxmagsquare.reshape(1, 768)
        yymagsquare1 = yymagsquare.reshape(1, 768)

        # each findat is length(32*times)x768 and contains data from all antennas 
        findatyy[z] = yymagsquare1
        findatxx[z] = xxmagsquare1

    # concatenate times and antennas
    # write out length x780 array 
    # where columns are time, antenna number, chan1, chan2,...,chan768. 
    outxx = np.concatenate((dumbtime, ant, findatxx), axis=1)
    outyy = np.concatenate((dumbtime, ant, findatyy), axis=1)

    # write out data for diagnostics
    np.savetxt("yyfull_"+obsid+".txt",outyy)
    np.savetxt("xxfull_"+obsid+".txt",outxx)

    # close table
    tb.close()


####################################################################################################

def ftimage(obsid):
    """
    FT a calibration image.

    """

    myvis = obsid + ".ms"
    imcalfile = config["IMCALFILE"]+"_cal.fits"
    importfits(fitsimage=imcalfile, imagename="imcal.model")
    im.open(myvis)
    im.defineimage()
    im.setoptions(ftmachine="wproject", wprojplanes=config["WPROJPLANES"])
    im.ft(model="imcal.model")
    im.close()


####################################################################################################

def generate_calibration(obsid, bcalname, selfcal=False, wsclean=False):
    """
    Generate <obsid>.bcal from <obsid>.cl.

    """

    # parse filenames
    myvis = obsid + ".ms"
    if selfcal or wsclean:
        mymodel = bcalname + "_pos.model"
        mybcal = bcalname + ".bcal"
    else:
        mycomplist = obsid + "_" +str(config["NSRCS"]) + ".cl"
        mybcal = obsid + "_" + str(config["NSRCS"]) + ".bcal"

    # clear previous calibration
    if not wsclean:
        print "=" * 50
        print "Clear previous calibration for " + myvis
        clearcal(myvis)

    # create model visibilities
    if selfcal and not wsclean:
        print "Create model visibilities using model image " + mymodel
        im.open(myvis)
        im.defineimage()
        im.setoptions(ftmachine="wproject", wprojplanes=config["WPROJPLANES"])
        im.ft(model=mymodel)
        im.close()
    if not selfcal and not wsclean:
        print "Create model visibilities using component list " + mycomplist
        ft(vis=myvis, complist=mycomplist)

    # calculate calibration solutions
    print "Calculate bandpass solutions.\nCalibration table: " + mybcal
    print "Bandpass parameters: "
    print "\tminsnr = " + str(config["MINSNR"])
    print "\trefant = " + str(config["REFANT"])
    print "\tminblperant = " + str(config["MINBLPERANT"])
    print "\tbandtype = " + str(config["BANDTYPE"])
    print "\tsolint = " + str(config["SOLINT"])
    print "\tsolnorm = " + str(config["SOLNORM"])
    print "\tuvrange = " + str(config["UVRANGE"])
    bandpass(vis=myvis, caltable=mybcal, minsnr=config["MINSNR"], refant=config["REFANT"], minblperant=config["MINBLPERANT"], bandtype=config["BANDTYPE"], solint=config["SOLINT"], solnorm=config["SOLNORM"], selectdata=True, uvrange=config["UVRANGE"], append=False)


####################################################################################################

def apply_calibration(obsid, selfcal=False, wsclean=False, bcalname=None):
    """
    Apply .bcal to <obsid>.ms.
    bandpass parameters specified in config. 

    """

    # parse filenames
    myvis = obsid + ".ms"
    if config["CALTABLE"] == "default":
        if selfcal or wsclean:
            mybcal = bcalname + ".bcal"
        else:
            mybcal = obsid + "_" + str(config["NSRCS"]) + ".bcal"
    else:
        mybcal = config["CALTABLE"] + ".bcal"

    # clear previous calibration
    if not wsclean:
        print "=" * 50
        print "Clear previous calibration for " + myvis
        clearcal(myvis)

    # apply calibration
    print "=" * 50
    print "Apply calibration table " + mybcal
    print "\tuvrange = " + str(config["UVRANGE"])
    applycal(vis=myvis, gaintable=mybcal, selectdata=True, uvrange=config["UVRANGE"], flagbackup=True)


####################################################################################################

def fhd_cal(obsid):

    # get ms info
    myvis = obsid + ".ms"
    ms.open(myvis)
    metadata = ms.metadata()
    nchan = metadata.nchan(0)  # number of frequency channels in ms
    ms.done()

    # parse fhd calibration filename
    fhdcaldir = config["FHDCALDIR"]
    npzfile = fhdcaldir + "/" + obsid + "_FHD_cal.npz"

    # load fhd calibration file
    print "Load " + npzfile
    fhdcal = np.load(npzfile)
    gains = fhdcal['G'] # shape (2,384,128)
    invgains = 1./gains
    invgains[np.isinf(invgains)]=0
    calflags = fhdcal['mask'] # shape (2,384,128)
    fhdchan = gains.shape[1]  # number of frequency channels in fhd calsoln
    offset = int(nchan / fhdchan)
    print "Measurement set has " + str(nchan) + " fine channels"
    print "FHD calibration has " + str(fhdchan) + " fine channels"
    print "Will apply the same FHD calibration solutions to every " + str(offset) + " consecutive channels"
    
    # clear previous calibration
    print 'Run clearcal to add CORRECTED_DATA column to measurement set'
    clearcal(vis=myvis)
    tb.open(myvis, nomodify=False)
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    
    print 'Apply FHD calibration and flags to measurement set'
    for ch in range(fhdchan):
        ch1 = ch * offset
        ch2 = (ch+1) * offset
	print 'Apply FHD calibration for channels ' + str(ch1) + " to " + str(ch2)
        while ch1 < ch2:
            data = tb.getcolslice('DATA',[0,ch1],[3,ch1],[1,1]) # shape (4, 1, 455168)

            data[0,0,:] = data[0,0,:] * np.transpose(invgains[0,ch,ant1]*np.conj(invgains[0,ch,ant2])) # XX
            data[3,0,:] = data[3,0,:] * np.transpose(invgains[1,ch,ant1]*np.conj(invgains[1,ch,ant2])) # YY
            tb.putcolslice('CORRECTED_DATA',data,[0,ch1],[3,ch1],[1,1])
            
            flags = tb.getcolslice('FLAG',[0,ch1],[3,ch1],[1,1])
            flags[0,0,:] = np.logical_or(np.logical_or(flags[0,0,:],np.transpose(calflags[0,ch,ant1])),np.transpose(calflags[0,ch,ant2]))
            flags[3,0,:] = np.logical_or(np.logical_or(flags[3,0,:],np.transpose(calflags[1,ch,ant1])),np.transpose(calflags[1,ch,ant2]))
            tb.putcolslice('FLAG',flags,[0,ch1],[3,ch1],[1,1])
            ch1 += 1
        
    tb.flush()
    tb.close()


####################################################################################################

def selfcal(obsid, wsclean=False):
    """
    Enter selfcal loop. 

    """
    
    # parse filenames
    myvis = obsid + ".ms"

    # apply selfcal
    if wsclean:
        print "=" * 50
        print "WSCLEAN SELFCAL"
        generate_calibration(obsid, "selfcal_"+str(opts.loop), wsclean=wsclean)
        apply_calibration(obsid, wsclean=wsclean, bcalname="selfcal_"+str(opts.loop))
    else:
        print "=" * 50
        print "CASA CLEAN SELFCAL"
        for i in xrange(config["NLOOPS"]):
            tmp_im = "selfcal_" + str(i)
            clean_image(obsid, selfcal=True, imname=tmp_im)
            rm_negative_pix(tmp_im+".model", tmp_im+"_pos.model")
            generate_calibration(obsid, tmp_im, selfcal=True)
            apply_calibration(obsid, selfcal=True, bcalname=tmp_im)


def rm_negative_pix(oldmodelfile, newmodelfile):
    """
    Remove negative pixels from the clean model. 

    """

    tb.open(oldmodelfile)
    tb.copy(newtablename=newmodelfile)
    tb.close()
    
    tb.open(newmodelfile, nomodify=False)
    d = tb.getcol("map")
    d[d<0] = 0
    tb.putcol("map", d)
    tb.flush()
    tb.close()
    

####################################################################################################

def export_calibrated_uvfits(obsid):
    """
    Export calibrated <obsid>.ms to <obsid>_cal.uvfits.

    """

    # parse filenames
    myvis = obsid + ".ms"
    myuvfits = obsid + "_cal.uvfits"

    # export to uvfits
    print "=" * 50
    print "Export " + myvis + " to " + myuvfits
    exportuvfits(vis=myvis, fitsfile=myuvfits, datacolumn="corrected", multisource=False)


####################################################################################################

def clean_image(obsid, selfcal=False, imname=None):
    """
    Clean <obsid>.ms.
    clean parameters specified in config.

    """
    # parse filenames
    myvis=obsid+".ms"

    # echo parameter inputs
    print "="*50
    print "Make images for " + myvis
    print "\t- number of iterations: " + str(config["NITER"])
    print "\t- threshold: " + str(config["THRESHOLD"])
    print "\t- gridmode: " + str(config["GRIDMODE"])
    print "\t- wprojplanes: " + str(config["WPROJPLANES"])
    print "\t- gain: " + str(config["GAIN"])
    print "\t- psfmode: " + str(config["PSFMODE"])
    print "\t- imagermode: " + str(config["IMAGERMODE"])
    print "\t- cyclefactor: " + str(config["CYCLEFACTOR"])
    print "\t- imagesize: " + str(config["IMSIZE"])
    print "\t- cellsize: " + str(config["CELL"])
    print "\t- weighting: " + str(config["WEIGHTING"])
    print "\t- uvrange: " + str(config["UVRANGE"])

    # make clean image
    imagelist = []
    if selfcal:
        clean_image = imname
        pols = ["XXYY"]
    else:
        clean_image = "clean_" + config["IMAGEMODE"] + "_" +obsid
        if config["IMAGESUFFIX"] != "none":
            clean_image = clean_image + "_" + config["IMAGESUFFIX"]
        pols = ["XX", "YY"]
    for mystokes in pols:
        if selfcal:
            myimagename = clean_image
        else:
            myimagename = clean_image + "_" + mystokes
        imagelist.append(myimagename)
        print "\nImage for Stokes " + mystokes + " named " + myimagename
        if config["IMAGEMODE"] == "fullband":
            mymode = "mfs"
            mynchan = -1
            mystart = 0
            print "\t == FULLBAND == "
            print "\t- mode: " + str(mymode)
        elif config["IMAGEMODE"] == "imgcube":
            mymode="channel"
            mynchan = 768
            mystart = 0
            print "\t == IMAGE CUBE =="
            print "\t- mode: " + str(mymode)
            print "\t- nchan: " + str(mynchan)
            print "\t- start: " + str(mystart)
        clean(vis=myvis, imagename=myimagename, niter=config["NITER"], mode=mymode, nchan=mynchan, start=mystart, gridmode=config["GRIDMODE"], wprojplanes=config["WPROJPLANES"], gain=config["GAIN"], threshold=config["THRESHOLD"], psfmode=config["PSFMODE"], imagermode=config["IMAGERMODE"], cyclefactor=config["CYCLEFACTOR"], imsize=[config["IMSIZE"], config["IMSIZE"]], cell=config["CELL"], stokes=mystokes, weighting=config["WEIGHTING"], selectdata=True, uvrange=config["UVRANGE"], interactive=False)

    return imagelist


####################################################################################################

def export_images_to_fits(imagelist):
    """
    Export <imagename>.image to <imagename>.fits.

    """

    for image in imagelist:
        print "=" * 50
        print "Export " + image + ".image to fits"
        exportfits(imagename=image+".image", fitsimage=image+".fits", velocity=False, stokeslast=False, overwrite=True)


####################################################################################################

# begin script

# read config file
config = read_config(configfile)

if mode == "prep":
    import_data(obsid)
    if config["ACAVG"]:
        prep_rfi_files(obsid)

if mode == "flag":
    if config["FLAGFILE"] == "default":
        flagfile = obsid + ".flg"
    else:
        flagfile = config["FLAGFILE"]
    print "=" * 50
    print "Apply flag file " + flagfile
    myvis = obsid + ".ms"
    execfile(flagfile)
    flagdata(vis=myvis, mode="summary")

if mode == "imagecal": 
    ftimage(obsid)

if mode == "gencal":
    generate_calibration(obsid, None)

if mode == "applycal":
    apply_calibration(obsid)

if mode == "fhdcal":
    fhd_cal(obsid)

if mode == "selfcal_wsclean":
    selfcal(obsid, wsclean=True)

if mode == "selfcal":
    selfcal(obsid, wsclean=False)

if mode == "exportuvfits":
    export_calibrated_uvfits(obsid)

if mode == "clean":
    imagelist = clean_image(obsid)
    export_images_to_fits(imagelist)  # no pb correction
