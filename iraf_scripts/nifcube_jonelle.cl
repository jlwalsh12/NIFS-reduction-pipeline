# Cpyright(c) 2006-2013 Association of Universities for Research in Astronomy, Inc.
#
# NIFCUBE -- Transform long slit images to rectified and wavelength
# calibrated datacubes.  This requires the 2D dispersion and distortion
# functions be previously fit by NSFITCOORDS.
#
# Version: Sept 23, 2005 FV  First version
#
# Complete revision history found in CVS
#
# JLW edits to nifcube.cl (12/29/17):
#
# nifcube is a wrapper to the more general gemcube. here, nifcube has
# been modified so that more information could be passed in and
# returned, such as the extension that is being worked on, the primary
# header information that the final cube should have, and the output
# weight image from gemcube. in addition, the input images have a
# header keyword that provides the location/file name of a mask. this
# mask is used to exclude bad edges around the geometrically rectified
# slices. the call to gemcube (where all the heavy lifting is done)
# remains the same, except here the output weight image is asked to be
# returned, and the mask for each slice is provided.

procedure nifcube_jonelle (inimages, inext, outphu) 

char    inimages    {prompt = "Input NIFS spectra"}
char    inext       {prompt = "Extension Type (SCI/VAR/DQ)"}
char    outphu      {prompt = "Fits file name; will copy primary header to output cube"}
char    outcubes    {"", prompt = "Output cubes"}
char    outprefix   {"c", prompt = "Prefix for output cubes"}
char    outweight   {"", prompt = "Output weight cube"}
char    reference   {"", prompt = "Reference output cube"}
int     waxis       {3, prompt = "Axis for wavelength"}
int     saxis       {2, prompt = "Axis for coordinate along slices"}
int     taxis       {1, prompt = "Axis for coordinate across slices"}
real    wmin        {INDEF, prompt = "Minimum wavelength (A)"}
real    wmax        {INDEF, prompt = "Maximum wavelength (A)"}
char    wscale      {"+INDEF", prompt = "Wavelength scale (A/pix)"}
char    sscale      {"0.05", prompt = "Spatial scale (arcsec/pix)"}
char    drizscale   {"1.", prompt = "Drizzle scale factors"}
real    blank       {0., prompt = "Value when there is no data"}
char    logfile     {"", prompt = "Logfile"}                        # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                       # OLDP-2
int     status      {0, prompt = "Exit status (0=Good)"}            # OLDP-4
struct  *scanfile1  {"", prompt = "Internal use"}                   # OLDP-4

begin
    char    l_inimages = ""
    char    l_inext = ""
    char    l_outphu = ""
    char    l_outcubes = ""
    char    l_outprefix = ""
    char    l_outweight = ""
    char    l_reference = ""
    char    l_wscale = ""
    char    l_sscale = ""
    struct  l_drizscale = ""
    char    l_logfile = ""
    int     l_waxis, l_saxis, l_taxis
    real    l_wmin, l_wmax, l_blank
    bool    l_verbose

    char    l_ext = ""      # Obtained from nsheaders

    bool    intdbg
    int     nin, nout, junk
    int     axis1, axis2, axis3
    real    cmin1, cmin2, cmin3
    real    cmax1, cmax2, cmax3
    char    cd1_1, cd2_2, cd3_3
    char    badhdr
    char    imgin, imgout, line
    char    tmpin, tmpout, tmpfile, tmpcube
    struct  l_date, sline
        
    # Define temporary files
    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout") 
    tmpfile = mktemp ("tmpfile") 

    # Cache parameter files.
    cache ("gemextn", "ckinput", "gemdate")
        
    # Initialize
    status = 1
    intdbg = no
    axis1 = INDEF; axis2 = INDEF; axis3 = INDEF
    cmin1 = INDEF; cmin2 = INDEF; cmin3 = INDEF
    cmax1 = INDEF; cmax2 = INDEF; cmax3 = INDEF
    cd1_1 = "INDEF"; cd2_2 = "INDEF"; cd3_3 = "INDEF"

    junk = fscan (  inimages, l_inimages)
    junk = fscan (  inext, l_inext)
    junk = fscan (  outphu, l_outphu)
    junk = fscan (  outcubes, l_outcubes)
    junk = fscan (  outprefix, l_outprefix)
    junk = fscan (  outweight, l_outweight)
    junk = fscan (  reference, l_reference)
    l_waxis = waxis
    l_saxis = saxis
    l_taxis = taxis
    l_wmin = wmin
    l_wmax = wmax
    junk = fscan (  wscale, l_wscale)
    junk = fscan (  sscale, l_sscale)
    junk = fscan (  drizscale, l_drizscale)
    l_blank = blank
    junk = fscan (  logfile, l_logfile)
    l_verbose = verbose

    badhdr = ""
    if (l_inext == "SCI") {
       junk = fscan (nsheaders.sci_ext, l_ext)
       if ("" == l_ext) badhdr = badhdr + " sci_ext"
       }
    if (l_inext == "VAR") {
       junk = fscan (nsheaders.var_ext, l_ext)
       if ("" == l_ext) badhdr = badhdr + " var_ext"
       }
    if (l_inext == "DQ") {
       junk = fscan (nsheaders.dq_ext, l_ext)
       if ("" == l_ext) badhdr = badhdr + " dq_ext"
       }

    if (l_logfile == "") {
        junk = fscan (nifs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "nifs.log"
            printlog ("WARNING - NIFCUBE:  Both logfile and nifs.logfile",
                l_logfile, verbose+) 
            printlog ("                        are empty.  Using " //l_logfile,
                l_logfile, verbose+)
        }
    }

    printlog ("---------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)
    date | scan (l_date)
    printlog ("NIFCUBE -- " // l_date, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    if ("" != badhdr) {
        printlog ("ERROR - NIFCUBE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+)
        goto clean
    }

    # The input and output list handling is a little tricky because
    # we might want to combine all extensions into one cube or combine
    # the extensions from each input mef into a separate cube.  In
    # either case the input to TRANSCUBE will be a list of extensions
    # and the output is a single cube for that list.

    # Expand and verify input.  First check all input and then
    # make list of mef files.
    ckinput (input=l_inimages, output="", prefix="", outlist=tmpin,
        name="NIFCUBE_JONELLE", dependflag="NSFITCOO", procflag="",
        sci_ext=l_ext, logfile=l_logfile, verbose=l_verbose,
        vverbose=intdbg)
    if (ckinput.status == 1)
        goto clean
    nin = ckinput.nfiles

    # Expand and verify output.
    if (intdbg) print ("expansion and verification of output")
    if (l_outcubes == "" && l_outprefix == "") {
        printlog ("ERROR - NIFCUBE: Output not specified",
            l_logfile, verbose+)
        goto clean
    }
    gemextn (l_outcubes, check="absent", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension,kernel", replace="",
        outfile=tmpout, logfile="", glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NIFCUBE: Problems with output.",
            l_logfile, verbose+)
        goto clean
    }
    if (gemextn.count == 0) {
        printf ("%%^%%%s%%@%s\n", l_outprefix, tmpin) | scan (line)
        gemextn (line, check="absent", process="none", index="", extname="",
            extversion="", ikparams="", omit="kernel,exten", replace="",
            outfile=tmpout, logfile="", glogpars="", verbose=l_verbose)
        if (gemextn.count == 0 || gemextn.fail_count != 0) {
            printlog ("ERROR - NIFCUBE: No or incorrectly formatted output \
                files\n", l_logfile, verbose+)
            goto clean
        }
    }
    nout = gemextn.count

    # Check number of input and output images
    if (nout > 1 && nin != nout) {
        printlog ("ERROR - NIFCUBE: Input and output files don't match.",
            l_logfile, verbose+)
        goto clean
    }

    # Log output files.
    printlog ("Using output files:", l_logfile, l_verbose)
    if (l_verbose) type (tmpout)
        type (tmpout, >> l_logfile)

    # Make file list.
    if (intdbg) print ("constructing file list")
    if (nin == nout) {
        rename (tmpout, tmpfile, field="all")
        joinlines (tmpfile, tmpin, output=tmpout, delim=" ", missing="-", \
            maxchar=161, shortest+, verbose-)
        delete (tmpfile, verify-)
    }
    delete (tmpin, verify-)

    # Main loop for making cubes.
    scanfile1 = tmpout
    while (fscan (scanfile1, imgout, imgin) != EOF) {
        if (intdbg) print ("creating " // imgout)
        if (nin == nout)
            gemextn (imgin, check="", process="expand", index="",
                extname=l_ext, extversion="1-", ikparams="", omit="",
                replace="", outfile=tmpin, logfile="", glogpars="",
                verbose=l_verbose)
        else
            gemextn (l_inimages, check="", process="expand", index="",
                extname=l_ext, extversion="1-", ikparams="", omit="",
                replace="", outfile=tmpin, logfile="", glogpars="",
                verbose=l_verbose)

        if (intdbg) print ("gemcube ...")
        if (verbose)
            line = "STDOUT,"//l_logfile
        else
            line = l_logfile

        # Set output axis parameters.
        if (l_waxis == 1) {
            axis1 = 1
            cmin1 = l_wmin
            cmax1 = l_wmax
            cd1_1 = l_wscale
        } else if (l_waxis == 2) {
            axis2 = 1
            cmin2 = l_wmin
            cmax2 = l_wmax
            cd2_2 = l_wscale
        } else if (l_waxis == 3) {
            axis3 = 1
            cmin3 = l_wmin
            cmax3 = l_wmax
            cd3_3 = l_wscale
        }
        if (l_saxis == 1) {
            axis1 = 2
            cd1_1 = l_sscale
        } else if (l_saxis == 2) {
            axis2 = 2
            cd2_2 = l_sscale
        } else if (l_saxis == 3) {
            axis3 = 2
            cd3_3 = l_sscale
        }
        if (l_taxis == 1) {
            axis1 = 3
            cd1_1 = l_sscale
        } else if (l_taxis == 2) {
            axis2 = 3
            cd2_2 = l_sscale
        } else if (l_taxis == 3) {
            axis3 = 3
            cd3_3 = l_sscale
        }

        tmpcube = mktemp ("tmpcube")
	gemcube.geofunc = "gfwcs"
        gemcube ("@"//tmpin, tmpcube, masks="", weights=l_outweight, logfiles=line, \
            bpm="!GEMMASK", scale="", wt="", wcsreference=l_reference, \
            wttype="drizzle", drizscale=l_drizscale, blank=l_blank, \
            memalloc=200., \
            geofunc.axis1=axis1, geofunc.cmin1=cmin1, geofunc.cmax1=cmax1, \
            geofunc.cd1_1=cd1_1, \
            geofunc.axis2=axis2, geofunc.cmin2=cmin2, geofunc.cmax2=cmax2, \
            geofunc.cd2_2=cd2_2, \
            geofunc.axis3=axis3, geofunc.cmin3=cmin3, geofunc.cmax3=cmax3, \
            geofunc.cd3_3=cd3_3, \
            geofunc.nonspatial=l_waxis, geofunc.square=l_saxis//" "//l_taxis, \
            geofunc.rotate=no)
        gemhedit (tmpcube, "dispaxis", l_waxis, "", delete-)

	delete (tmpin, verify-, >& "dev$null")

        # Create a MEF file
        printlog (tmpcube//" -> "//imgout, l_logfile, l_verbose)
	fxcopy (l_outphu//"[0]", imgout//".fits", >& "dev$null")
        fxinsert (tmpcube//".fits", imgout//".fits[1]", groups="0", ver-)

        # Fix headers
        gemhedit (imgout//"[1]", "EXTVER", 1, "Extension version", delete-)
        gemhedit (imgout//"[1]", "EXTNAME", l_ext, "Extension name",
            delete-)

        # From Kathleen...
        # Note to James and Insoek about the PHU headers:
        #   I guess there's a whole bunch of WCS stuff in the first extension
        #     that needs to be copied to the PHU.
        #   There are probably more header fixing to be done.  That's for
        #     you to find out.

        gemdate ()
        gemhedit (imgout//"[0]", "NIFCUBE", gemdate.outdate,
            "UT Time stamp for NIFCUBE", delete-)
        gemhedit (imgout//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)


        imdelete (tmpcube, verify-, >& "dev$null")
    }

    # Completed successfully
    status = 0

clean:
    if (status == 0) {
        printlog ("", l_logfile, no)
        printlog ("NIFCUBE  Exit status good", l_logfile, l_verbose)
        printlog ("--------------------------------------------------------\
            -----------------------", l_logfile, l_verbose)
    }
    scanfile1 = ""
    delete (tmpin, ver-, >& "dev$null")
    delete (tmpout, ver-, >& "dev$null")
    delete (tmpfile, ver-, >& "dev$null")
end
