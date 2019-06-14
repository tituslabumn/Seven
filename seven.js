importPackage(Packages.ij); 
importPackage(Packages.java.io); 
importPackage(Packages.java.awt); 
importPackage(Packages.ij.measure); 
load(IJ.getDirectory("ImageJ")+"complex.js");
 
//// seven.js 
//// a library to analyze filopodia using tip markers 
//// 
//// ImageJ API documentation: 
//// http://rsb.info.nih.gov/ij/developer/api/ 
//// 
//// Java API documentation: 
//// http://docs.oracle.com/javase/7/docs/api/overview-summary.html 
//// 
//// JavaScript documentation for the Rhino interpreter: 
//// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference 
 
// Useful constants 
var newline = "\n"; 
var sep = File.separator; 
var WATER_AUTO = -2; 
var LOCAL_AUTO = -1; 
var GLOBAL_AUTO = 0; 
var digits = 6; // number of significant digits for logging
var DEBUG = false; 
var once = true;
var dt = 0.793; // default timestep in seconds if not available from logfile
var boxwidth_um = 0.8; // height in microns of linescan 
var minarea_um2 = 1; // minimum area in um^2; prevents glitches caused by very small ROIs 
var interpolation = 0.5; // scaling factor <= 1 for interpolation of the cell perimeter ROI
var MEASUREMENTS = Measurements.AREA + Measurements.MEAN + 
	Measurements.MEDIAN + Measurements.MIN_MAX + Measurements.MODE + 
	Measurements.SHAPE_DESCRIPTORS + Measurements.PERIMETER +
	Measurements.SKEWNESS + Measurements.KURTOSIS;
 
// Assay-specific options: check these before running! 
var format = "tif"; // file extension 
var maskformat = "zip"; // file extension for mask file, if different 
var resultsformat = "csv"; // file extension for results files
var linescantag = "newzeroscan"; 
var anaversion = "7"; 
var frame = 1; // frame of the raw data to analyze (starting with 1)
var minfp = 1; // set to 1 unless you know what this does
var banded = true; // Whether to look at cell perimeter band 
var thresholdcode = LOCAL_AUTO; // see list above 
var skiptips = false; // skip filopod tip search 
var semi = true; // set to true to wait for manual corrections - recommended with local thresholds 
var editlinescans = false; // manual editing of manual ROIs after they are loaded 
var onlyreviewlinescans = false; // manual review of linescans but don't overwrite
 
// Threshold images to determine cell boundaries 
// generates two images: the binary mask and the mask applied to the original image 
function ThresholdCells(anadir, acqname, thresholds, prefix) { 
	var acqlen = acqname.length; 
	var sublist = anadir.list(); 
	var imagefile = null; 
	var subname = ""; 
 
	// Preliminary step - copy any images in subdirectories to the root 
	for (var i=0; i < sublist.length; i++) { 
		IJ.log(sublist[i]); 
		var subdir = new File(anadir.getCanonicalPath(), sublist[i]); 
		if (subdir.isDirectory() && acqlen > 0) { 
			var imagelist = subdir.list(); 
			for (var j=0; j < imagelist.length; j++) { 
				var imagename = imagelist[j]; 
				var imagefile = new File(subdir.getCanonicalPath(), imagename); 
				if (Packages.java.lang.Integer.parseInt(imagename.length()) > acqlen) { 
					if (imagename.substring(0,acqlen) == acqname) { 
						CopyFile(imagefile, anadir, imagename); 
					} 
				} 
			} 
		} 
	} 
 
	// Apply thresholds to create a mask: 
	// loop over the list of thresholds, then over the images 
	// Masks are stored in subdirectories named with the form prefix_a_i 
	for (var a = 0; a < thresholds.length; a++) { 
		var th = thresholds[a]; 
		for (var i = 0; i < sublist.length; i++) { 
			var imagefile = new File(anadir.getCanonicalPath(), sublist[i]); 
			var imagename = sublist[i]; 
			
			// Raw images are identified by numeric ID in the original file name 
			// Because the ID is handled as a string, first count the characters (digits) 
			var digits = 0; 
			for (var j = acqlen; j < imagename.length(); j++) { 
				var testchar = imagename.substring(j,j+1); 
 
				if (testchar == "_" ||  testchar == ".") { 
					break; 
				} 
				digits++; 
			} 
 
			// Read the image ID from the file name and create the analysis subdirectory 
			if (Packages.java.lang.Integer.parseInt(imagename.length()) > acqlen+digits && !imagefile.isDirectory()) { 
				if (imagename.substring(0,acqlen) == acqname && getExt(imagename) != resultsformat) { 
					// number subdirectories using the number following the acqname prefix 
					subdir_index = imagename.substring(acqlen,acqlen+digits); 
					if (filterInt(subdir_index) > 0) { 
						subname = prefix+IJ.d2s(th,0)+"_"+subdir_index; 
					} else { 
						subname = prefix+IJ.d2s(th,0)+"_"+ 
							IJ.d2s(i+1,0);	// in case image file names are not numbered 
					} 
					var subdir = new File(root.getCanonicalPath(), subname); 
					if (!subdir.isDirectory()) { 
						subdir.mkdir(); 
					} 
					// Subdirectory contains a hint to the original file name 
					hint = subdir.getCanonicalPath() + sep + "hint.txt"; 
					// Save the hint as a text file
					saveText(hint, imagename, false); 
 
					// If the mask file already exists, skip processing 
					var f = new File(subdir.getCanonicalPath(), 
						"Capture-mask-"+anaversion+"."+maskformat); 
					if (!f.exists()) { 
 
						// Process image 
						var mask = IJ.openImage(imagefile); 
 						if (mask.getNSlices() >= frame)
 							mask.setSlice(frame); 
						var img = mask.duplicate(); 
 						if (img.getNSlices() >= frame)
 							img.setSlice(frame); 
						// Read the pixel size, converting cm->um 
					    var cal = mask.getCalibration(); 
					    cm2um(cal); 
					    var dx = cal.getX(1); 
				    	// Cell boundaries are determined by median filter followed by thresholding. 
				    	// Length scale for the filter is approx. 2 um (originally 10 pixels @ 63X mag.) 
				    	var mask_width_um = 1;
					    var radius_med = Math.ceil(mask_width_um/dx); 
						IJ.run(mask, "Median...", "radius="+IJ.d2s(radius_med,0)+" slice stack"); 
						// Input image is 16-bit, but some methods need 8-bit image, 
						// so the image lookup table must be normalized before processing 
						IJ.run(mask, "Enhance Contrast", "saturated=0.35"); 
 
						// Apply threshold to create the mask 
						switch(th) { 
							case GLOBAL_AUTO:	// Global threshold (calculated using entire XY image) 
								IJ.run(mask, "8-bit", ""); 
								IJ.setAutoThreshold(mask, "Moments stack"); 
								IJ.run(mask, "Convert to Mask", "method=Moments background=Light stack"); 
								break; 
 
							case WATER_AUTO: 
							case LOCAL_AUTO:	// Local threshold (Bernsen algorithm) 
								IJ.run(mask, "8-bit", ""); 
								// Length scale for local thresholds = 2x the mask radius 
								var bernsen_width_um = 4; 
								var radius_th = Math.ceil(bernsen_width_um/dx); 
								IJ.run(mask, "Auto Local Threshold", "method=Bernsen radius="+ 
									IJ.d2s(radius_th,0)+" parameter_1=0 parameter_2=0 white stack"); 
								IJ.run(mask, "Convert to Mask", "method=Moments background=Light stack"); 
								if (th == WATER_AUTO) { 
									// Use watershed to separate neighboring cells 
									IJ.run(mask, "Watershed", "stack"); 
								} 
								invertImage(mask); 
								break; 
 
							default:		// Global threshold, using a manually specified value 
								if (th < 0) 
									IJ.showMessage("Invalid image threshold!"); 
								IJ.setThreshold(mask, 0, th); 
								IJ.run(mask, "8-bit", ""); 
								IJ.run(mask, "Convert to Mask", "method=Moments background=Light stack"); 
							} 
 
						//Return to the original image and normalize 
						IJ.run(img, "Enhance Contrast", "saturated=0.35"); 

						// Let the user manually correct the mask file 
						img.show(); 
						if (img.getNSlices() >= frame) // reset image position due to potential side effects
							img.setSlice(frame);
						mask.show(); 
						if (mask.getNSlices() >= frame) // reset image position due to potential side effects
							mask.setSlice(frame);

						if (semi) { 
							new Packages.ij.gui.WaitForUserDialog( 
								"Manual Corrections", "Please press OK when done.").show(); 
						} 
						saveImage(mask, maskformat, anadir, subname+"/Capture-mask", Packages.java.lang.Integer.parseInt(anaversion)); 
 
						// Apply the mask to the image (not required, but useful for visual inspection) 
						ic = new Packages.ij.plugin.ImageCalculator(); 
						img2 = ic.run("Min stack create", mask, img); 
						if (img2.getNSlices() >= frame)
							img2.setSlice(frame);

						invertImage(img2); 
						saveImage(img2, format, anadir, subname+"/Capture-threshold", Packages.java.lang.Integer.parseInt(anaversion)); 
 
						//Clean up 
						mask.changes = false; 
						img.changes = false; 
						img2.changes = false; 
						if (!DEBUG) { 
							mask.close(); 
							img.close(); 
							img2.close(); 
						} 
					} 
				} 
			} 
		} 
	} 
} 
 
// Analyze cells using the mask to calculate ROIs, cell area and mean intensity 
function AnalyzeCells(img1, anadir, depth, resultname, intensities, areas, cell_xs, cell_ys) { 
	if (img1.getNSlices() >= frame)
		img1.setSlice(frame);

	var img3 = img1.duplicate();  		
	if (img3.getNSlices() >= frame)
 		img3.setSlice(frame);

	var h = img1.getHeight(); 
	var w = img1.getWidth(); 
	// set image scale cm->um 
    var cal = img1.getCalibration(); 
	cm2um(cal); 
	var stats = img1.getStatistics(MEASUREMENTS); 
	var fullarea = stats.area; 
	var bgval = stats.dmode;
	 
	//Open mask file 
    var maskfile = anadir + "Capture-mask-"+anaversion+"." + maskformat; 
    var mask = openIf(maskfile, maskformat); 
	setMask(mask, depth, false); 
	 
    // get ROI of cell body 
    var rs = new RoiSet(); 
	var rsname = anadir+"analysis-RoiSet"+anaversion+".zip"; 
	var rsfile = new File(rsname); 
	// if rsfile does not exist, can still use the rsmanager window 
	if (rsfile.exists()) 
		rs.runCommand("Open", rsname); 
 
    if (rs.getName(1) != "Cells") { 
        IJ.showMessage("Error", "ROI Manager not populated"); 
        exit; 
    } 
    // get xy coordinates to identify cells with wand tool 
    var rois = rs.getRoisAsArray(); // returns Java array of ij.gui.Roi[]
	var xs = rois[1].getPolygon().xpoints; // returns Java array of int[]
    var ys = rois[1].getPolygon().ypoints; // returns Java array of int[]
 
    // invert image 
    invertImage(img1); 
    var ic = new Packages.ij.plugin.ImageCalculator(); 
    var img2 = ic.run("Max create stack", img1, mask); 
	if (img2.getNSlices() >= frame)
		img2.setSlice(frame);
 
    // invert image 
    invertImage(img2); 
    saveImage(img2, format, anadir, "Capture-mask-"+resultname+anaversion, 0); 
 
	// return to original image for quantitation 
	var rs = new RoiSet(); 
 
    // set image scale cm->um 
    var maskcal = mask.getCalibration(); 
	cm2um(maskcal); 
	var maskstats = mask.getStatistics(MEASUREMENTS); 
	// this is where the cell outlines are stored as ROIs
    for (var i = 0; i < xs.length; i++) { 
        IJ.doWand(mask, xs[i], ys[i], 0, "4-connected"); 
		StoreRoi(mask, mask.getRoi(), resultname+" "+IJ.d2s(i,0), rs); 
    } 
 
    // save rois 
    rs.runCommand("Select All"); 
    rs.runCommand("Save", anadir + "cell-"+resultname+"-RoiSet"+anaversion+".zip"); 
    var rois = rs.getRoisAsArray(); 

    // set image scale cm->um 
    var cal = img3.getCalibration(); 
	cm2um(cal); 
	// Measure mean intensity in original unmodified image      
	for (var i = 0; i < rois.length; i++) { 
		img3.setRoi(rois[i]); 
        rs.select(img3, i); 
        var cellstats = img3.getStatistics(MEASUREMENTS); 
        // copy cell intensity and xy position to global data arrays
        intensities[i] = cellstats.mean - bgval; // mean background-corrected cell intensity
        if (xs[i] > 0 && ys[i] > 0) {
        	cell_xs[i] = xs[i];
        	cell_ys[i] = ys[i];
        }
 
        // Measure cell area and copy to global data array
        areas[i] = 0; 
        if (cellstats.area < fullarea) // Incorrectly drawn ROIs will cover the full frame 
	        areas[i] = cellstats.area; 
        IJ.log(IJ.d2s(intensities[i], 0)); 
    } 
 
	if (!DEBUG) { 
		img1.close(); 
		img2.close(); 
		img3.close(); 
		mask.close(); 
	} 
} 
 
 
// Analysis of filopod line scans 
function AnalyzeScans(img1, imagefile, anadir, boxwidth_um) { 
	if (img1.getNSlices() >= frame)
		img1.setSlice(frame);

	var h = img1.getHeight(); 
	var w = img1.getWidth(); 
	// set image scale cm->um 
    var cal = img1.getCalibration(); 
	cm2um(cal); 
	var stats = img1.getStatistics(MEASUREMENTS); 
	var fullarea = stats.area; 
	var timestep = dt; 
	var kymodir = new File(anadir+sep+"Kymograph"); 
	var paramtab = new Packages.ij.measure.ResultsTable(); 
	paramtab.setPrecision(digits); 
 
    // get ROI of cell body 
    var rs = new RoiSet(); 
    var roifile = new File(anadir+"scan.zip"); 
    if (!roifile.exists() || !roifile.isFile()) { 
    	roifile = new File(anadir+"scan.roi"); 
    } 
    if (roifile.exists() && roifile.isFile()) { 
		if (!kymodir.isDirectory()) { 
			kymodir.mkdir(); 
		} 
		timestep = getTimingFromLogfile(img1, imagefile, timestep); 
		LoadRois(rs, roifile, ""); 
 
		if (editlinescans || onlyreviewlinescans) { 
		new Packages.ij.gui.WaitForUserDialog( 
				"Manual Corrections to ROI", "Please press OK when done.").show(); 
			if (editlinescans) { 
				rs.runCommand("Select All"); 
				rs.runCommand("Save",roifile.getCanonicalPath()); 
			} 
		} 
 
	    // Read pixel size, convert cm->um 
	    var cal = img1.getCalibration(); 
		cm2um(cal); 
	    var dx = cal.getX(1); 
	
		kymograph(img1, rs, timestep, boxwidth_um, paramtab, kymodir); 
	     
	    // Measure intensity in the masked image 
	    var masked = new File(anadir+sep+"Capture-mask-body"+anaversion+"-0."+format); 
	    if (masked.isFile()) { 
	    	var img2 = IJ.openImage(masked); 
 			if (img2.getNSlices() >= frame)
 				img2.setSlice(frame); 
		    var rois = rs.getRoisAsArray(); 
		    for (var i = 0; i < rois.length; i++) { 
				rois[i].setStrokeWidth(1); 
				img2.setRoi(rois[i]); 
		        rs.select(img2, i); 
				ScanLine(img2, anadir, "linescan", dx, i); 
		    } 
		    img2.changes = false; 
		    if (!DEBUG) { img2.close(); } 
	    } else { 
	    	IJ.showMessage("Could not file masked file to match "+imagefile); 
	    } 
    } 
 
	// Save the results table as a delimited text file
	paramtab.save(kymodir.getPath()+sep+"kymograph_results."+resultsformat); 
	img1.changes = false; 
	if (!DEBUG) { img1.close(); } 
} 
 
function ScanLine(img, dir, name, pixelwidth, index) { 
	var profileplot = new Packages.ij.gui.ProfilePlot(img); 
	var linescan = profileplot.getProfile(); 
	//var text = "Distance (um)\tIntensity (A.U.)"+newline; 
	var text = "";
	if (linescan != null) {
		for (var j = 0; j < linescan.length; j++) { 
			text += IJ.d2s(j*pixelwidth, digits)+"\t"+IJ.d2s(linescan[j],0)+newline; 
			//IJ.log(IJ.d2s(j*pixelwidth, digits)+"\t"+IJ.d2s(linescan[j],0)); 
		} 
		var textfile = dir+name+"-"+IJ.d2s(index,0)+"-"+anaversion+".txt"; 
		//var text = IJ.getLog(); 

		// Save the linescan as text file
		saveText(textfile, text, false); 
	}
} 
 
// Analysis for counting filopodia tips and registering to cells 
function AnalyzeTips(img1, imagefile, anadir, imagetab, boxwidth_um, firstpass) { 
	if (img1.getNSlices() >= frame)
		img1.setSlice(frame);

	// local array for filopodia
	var fparray = new Array();

	// get raw image stats	
	var stats = img1.getStatistics(MEASUREMENTS); 
	var bgval = stats.dmode; // use instead of minval - this is black level of whole image
	var minSNR = 3; 
	var noise = noiseThreshold(stats, minSNR); 
	
	// pixel size, convert cm->um 
    var cal = img1.getCalibration(); 
 	cm2um(cal); 
    var dx = cal.getX(1); 
	var boxheight = Math.ceil(boxwidth_um/dx); // height in pixels of linescan 
    var pixelsize = IJ.d2s(dx, digits); 
 
	//Open mask file 
	var maskfile = anadir + "Capture-mask-"+anaversion+"."+maskformat; 
	var mask = openIf(maskfile, maskformat); 
	setMask(mask, 0, true); 

	// get image dimensions
	var h = img1.getHeight(); 
	var w = img1.getWidth(); 
	var nFrames = img1.getNSlices(); // ImageJ parses images as XYZ by default (treat as XYT)
	var timestep = getTimingFromLogfile(img1, imagefile, dt); 
	
	// invert mask 			
	var ic = new Packages.ij.plugin.ImageCalculator(); 
	IJ.run(img1, "Subtract...", "value=1 stack"); // subtract 1 from raw image to avoid white pixels 
	var img2 = ic.run("Max create stack", img1, mask); 
	if (img2.getNSlices() >= frame)
		img2.setSlice(frame);

	var ip2 = img2.getProcessor(); 
	// Save an intermediate image with cells masked in white; this gets used in "linescanonly" mode 
	saveImage(img2, format, anadir, "Capture-mask-body"+anaversion, 0); 
 
	// Clean up from masking 
	if (!DEBUG) { mask.close(); } 
 
	IJ.run(img2, "Find Maxima...", 
		"noise="+IJ.d2s(noise,0)+" output=[Point Selection]"); 
	var rRaw = img2.getRoi(); 
	var rois = new RoiSet(); 
	StoreRoi(img2, rRaw, "Raw", rois); 
	var nCells = 0; 
	var nRawTips = 0; 
 
	var tipx = rRaw.getPolygon().xpoints; 
	var tipy = rRaw.getPolygon().ypoints; 
	var nPoints = rRaw.getNCoordinates(); 
	var maxval = Math.pow(2, 16)-1; 
 
	// Count cells 
	var pCells = new Polygon(); 
	var pRawTips = new Polygon(); 
	var pTips = new Polygon(); 
 
	for (var i=0; i<nPoints; i++) { 
		val = getValue(img2, tipx[i], tipy[i]); 
		if (val == maxval) { 
			nCells++; 
			pCells.addPoint(tipx[i],tipy[i]); 
		} else { 
			nRawTips++; 
			pRawTips.addPoint(tipx[i],tipy[i]); 
		} 
	} 

	// Copy cell position to global array
	var cell_xpos = new Array(nCells);
	var cell_ypos = new Array(nCells);
	for (var i=0; i<nCells; i++) {
		cell_xpos[i] = pCells.xpoints[i];
		cell_ypos[i] = pCells.ypoints[i];
	}

	if (pCells.npoints > 0) { // only write if populated 
		var rCells = new Packages.ij.gui.PointRoi(pCells); 
		rCells.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
		StoreRoi(img2, rCells, "Cells", rois); 
	} 
	if (pRawTips.npoints > 0) { // only write if populated 
		var rRawTips = new Packages.ij.gui.PointRoi(pRawTips); 
		rRawTips.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
		StoreRoi(img2, rRawTips, "RawTips", rois); 
	} 
 
	if (!skiptips) { 
		//Calculate radial distribution 
		var tipdata = "";
		var sampletab = new Packages.ij.measure.ResultsTable(); 
		var mpp = 0.21164;
		sampletab.setPrecision(digits);
		img2.setColor(Color.WHITE); 
		// Scan in radial pattern (approx 5 um diameter; original radius 12 px at 63X) 
		var scan_width = 5; 
		var radius_scan = Packages.java.lang.Integer.parseInt(Math.floor(scan_width/(dx*2))); 
		var fp = 0; 
		tipx = pRawTips.xpoints; 
		tipy = pRawTips.ypoints; 
		
		// Black out tips to improve radial search (approx 1 um square; originally 3x3 at 63X) 
		var blackout_width = 1; 
		var radius_black = Math.floor(blackout_width/(dx*2)); 
		for (var i=0; i < nRawTips; i++) { 
			var x = tipx[i]; 
			var y = tipy[i]; 
			var xlo = x - radius_black; 
			var xhi = x + radius_black; 
			var ylo = y - radius_black; 
			var yhi = y + radius_black; 
			if (xlo < 0) { xlo = 0; } 
			if (ylo < 0) { ylo = 0; } 
			if (xhi > (w-1)) { xhi = w-1; } 
			if (yhi > (h-1)) { yhi = h-1; } 
			for (var u=xlo; u<=xhi; u++) { 
				for (var v=ylo; v<=yhi; v++) { 
					ip2.putPixel(u, v, stats.mean); // mean value of image instead of bgval or black = 0
				} 
			} 
		} 
		 
		// Scan tips to see if they are connected to cells 
		var pTips = new Polygon(); 
		for (var a=0; a<nRawTips; a++) { 
			if (tipx[a] >= radius_scan && tipx[a] < (w-radius_scan) && 
				tipy[a] >= radius_scan && tipy[a] < (h-radius_scan)) {
				var radius_scanroi = new Packages.ij.gui.Line( 
					tipx[a], tipy[a], Packages.java.lang.Integer.parseInt(tipx[a]+radius_scan), tipy[a]); 
				radius_scanroi.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
				img2.setRoi(radius_scanroi);
				var thetastep = 3; 
				var nSteps = Math.ceil(360/thetastep); 
				img2.hide(); 
				IJ.run(img2, "Set Scale...", "distance=1 known="+pixelsize+" pixel=1 unit=um"); 
				IJ.run(img2, "Radial Reslice", "angle=360 degrees_per_slice="+ 
					IJ.d2s(thetastep, 0)+" direction=Clockwise"); 
				var radial_3d = IJ.getImage(); 
				radial_3d.setRoi(new Packages.ij.gui.Line(0, 0, Packages.java.lang.Integer.parseInt(radius_scan-1), 0)); 
				IJ.run(radial_3d, "Reslice [/]...", "output="+pixelsize+" slice_count=1 avoid");
				var radial_2d = IJ.getImage(); 
				IJ.run(radial_2d, "Rotate 90 Degrees Right", ""); 
				radial_2d.setRoi(new Rectangle(0, 0, nSteps, Packages.java.lang.Integer.parseInt(radius_scan-1))); 
				var radial_dist = new Array(nSteps); 
	 
				for (var u = 0; u < nSteps; u++) { 
					radial_dist[u] = 0; 
					for (var v = 0; v < radius_scan; v++) { 
						radial_dist[u] += getValue(radial_2d, u, v); 
					} 
					radial_dist[u] /= radius_scan; 
				} 
				if (!DEBUG) { 
					radial_3d.close(); 
					radial_2d.close(); 
				} 
	 
				var skew = skewness(radial_dist); 
				 
				if (skew > 1) { 
					pTips.addPoint(tipx[a],tipy[a]); 
					fp++; // increment filo count
				} 
			} 
		} 
		if (pTips.npoints > 0) { // only write if populated 
			var rTips = new Packages.ij.gui.PointRoi(pTips); 
			rTips.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
			StoreRoi(img2, rTips, "Tips", rois); 
		} 
	 
		// Register FP to cells 
		var maxdistance_um = 20; // maximum distance from cell center in microns 
		var maxpixeldist = maxdistance_um/dx;  
		var regfp = fp; 
		var regcross = 0;
		var cells_with_fp = nCells; 
		var fp_per_cell = new Array(nCells); // create local array for bookkeeping
		// initialize to zero - probably not necessary?
		for (var q=0; q<nCells; q++)
			fp_per_cell[q] = 0; 

		// assign global array dimensions - probably not necessary (JavaScript very forgiving)
		cell_per_fp = new Array(fp);     	// reassign dimensions of this global array
		intensity_per_fp = new Array(fp);   // reassign dimensions of this global array
		xpos_per_fp = new Array(fp);     	// reassign dimensions of this global array
		ypos_per_fp = new Array(fp);    	// reassign dimensions of this global array
		
		for (var u=0; u<fp; u++) { 
			var score = 1e6; 
			var reg = -1; 
			for (var v=0; v<nCells; v++) { 
				var disp_x = pCells.xpoints[v] - pTips.xpoints[u]; 
				var disp_y = pCells.ypoints[v] - pTips.ypoints[u]; 
				var disp = Math.sqrt(Math.pow(disp_x,2)+Math.pow(disp_y,2)); 
	 
				if (score > disp && maxpixeldist > disp) { 
					score = Math.sqrt(Math.pow(disp_x,2)+Math.pow(disp_y,2)); 
					reg = v; 
				} 
			} 
			if (reg >= 0) { 
				// Select the line connecting tip to cell center-of-area 
				var tiproi = new Packages.ij.gui.Line(pTips.xpoints[u],pTips.ypoints[u], 
						pCells.xpoints[reg],pCells.ypoints[reg]); 
				tiproi.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
				img2.setRoi(tiproi); 

				// Update the count of registered filopodia
				fp_per_cell[reg]++; 
				// update global data array for intensity using the unmasked img1 
				intensity_per_fp[u] = getValue(img1, pTips.xpoints[u], pTips.ypoints[u]) 
					- bgval; // background-corrected peak intensity at tip 
				// update global data arrays for cell registration and xy position 
				cell_per_fp[u] = reg;
				xpos_per_fp[u] = pTips.xpoints[u];
				ypos_per_fp[u] = pTips.ypoints[u];
			} else { 
				// the tip cannot be registered to a cell, i.e., it is a false positive detection event.
				intensity_per_fp[u] = 0; 
			    cell_per_fp[u] = -1; 
				xpos_per_fp[u] = -1;
				ypos_per_fp[u] = -1;
				regfp--; 
			} 
			//IJ.showMessage("index "+IJ.d2s(u,0)+", cell_per_fp = "+IJ.d2s(cell_per_fp[u],0));
		} 
			 
		// local array for filopodia
	 	fparray = new Array(regfp);
	 	var v = 0;
		for (var u = 0; u < fp; u++) {
			if (cell_per_fp[u] >= 0) {
				fparray[v] = new Filopod(xpos_per_fp[u]*dx, ypos_per_fp[u]*dx, 
					pCells.xpoints[cell_per_fp[u]]*dx, pCells.ypoints[cell_per_fp[u]]*dx,
					cell_per_fp[u], intensity_per_fp[u]);
				if (fparray[v] == null)
					IJ.error("Seven.js", "could not initialize filopod "+IJ.d2s(v,0));
				v++;
			}
		}
		if (v != regfp)
			IJ.error("Seven.js", "Failed to initialize Filopod array");
		
	    // need to set line color here for contrast with white mask in img2
		img2.setColor(Color.GRAY); 

		// Second pass to save zeroscans of filopod length 
		for (var u=0; u<fp; u++) { 
			index = cell_per_fp[u]; 
			// IJ.log("Cell index count =" + index + "Number of filopod attached =" + fp_per_cell[index]);
			//IJ.log("Filopodia u count =" + u);
			if (index >= 0 && fp_per_cell[index] >= minfp) { 
				var scanroi = new Packages.ij.gui.Line(pTips.xpoints[u],pTips.ypoints[u], 
						pCells.xpoints[index],pCells.ypoints[index]); 
				//var angle = scanroi.getAngle();
				//IJ.log("Angle value = "+ angle);
				scanroi.setPosition(1, frame, 1); // specify which frame (slice) to analyze
				scanroi.setStrokeWidth(1); 
				img2.setRoi(scanroi); 
	 
				// Calculate line scan using the masked img2 
				ScanLine(img2, anadir, linescantag, dx, u);
	 
				// Mark with a line on the masked img2 
				ip2.draw(scanroi); 
			} 
		} 
	 
		// Trim the count of cells 
		for (var z=0;z<nCells;z++) { 
			if (fp_per_cell[z] == 0) { cells_with_fp--; } 
		}
	}
 
	// Save ROIs on first pass 
	rois.runCommand("Select All"); 
	rois.runCommand("Save",anadir+"analysis-RoiSet"+anaversion+".zip"); 
 
	// Save ROIs on first pass, on second pass save log as well 
	if (!firstpass) { 
		// table for reporting results	 
		var celltab = new Packages.ij.measure.ResultsTable(); 
		celltab.setPrecision(digits); 
	    var spacingtab = new Packages.ij.measure.ResultsTable(); 
		spacingtab.setPrecision(digits); 
	    var crossingtab = new Packages.ij.measure.ResultsTable(); 
		crossingtab.setPrecision(digits); 
		 
	    // Spacing analysis of cell perimeter band 
	    // Load ROIs for cells 
	    var rs = new RoiSet(); 
	    var roifile = new File(anadir+"cell-body-RoiSet"+anaversion+".zip"); 
	    if (!roifile.exists()) { 
	    	roifile = new File(anadir+"cell-body-RoiSet"+anaversion+"+.roi"); 
	    } 
		
		// begin perimeter band analysis
	    if (banded && roifile.exists()) { 
			rs.runCommand("Open",roifile.getCanonicalPath()); 
			var rois = rs.getRoisAsArray(); 
			if (rs.getCount() != cell_xpos.length)
				IJ.error("seven.js", "ROIs != nCells: "+IJ.d2s(rs.getCount(),0)+" != "+
					IJ.d2s(cell_xpos.length,0));

			// iterate over cells 
			for (var i = 0; i < rs.getCount(); i++) { 
				// begin reporting 
				celltab.incrementCounter(); 
				celltab.addValue("Cell ID+1", i+1); 
				celltab.addValue("Cell Intensity", cell_body[i]); 
				celltab.addValue("Cell Area (um^2)", cell_area_body[i]); 
			    //if (fp_per_cell[i] > 0) { // optionally only look at cells with fp 
				spacingtab.incrementCounter(); 
				 
				// analyze cells by ROI 
				img1.setRoi(rois[i]); 
				rs.select(img1, i); 
 
				if (cell_area_body[i] > minarea_um2) { 
					// convert magic wand ROI to line ROI 
					//IJ.run(img1, "Area to Line", ""); // this leaves a gap at end 
					//var interp_polygon = img1.getRoi().getPolygon();
					var interp_polygon = img1.getRoi().getInterpolatedPolygon(interpolation, true);
					var outline_points = interp_polygon.xpoints.length;
					var contour = new Array(outline_points);
					var startroi = new Packages.ij.gui.Arrow(interp_polygon.xpoints[0],  
						interp_polygon.ypoints[0],interp_polygon.xpoints[0]+1,  
						interp_polygon.ypoints[0]+1);
					startroi.setPosition(1, frame, 1); // specify which frame (slice) to analyze 
					img2.setRoi(startroi); 
					// Draw a 1px line on the masked img2 
					ip2.draw(startroi); 
					
					var narrow = new Packages.ij.plugin.Straightener(); 
					var outline_roi = outliner(interp_polygon, outline_points, frame);
					var band = new ImagePlus();
					
					// initialize contour for lookup of u-coordinate from outline point index
					for (var j = 0; j<outline_points; j++) {
						if (j>1) {
							var this_outline_roi = outliner(interp_polygon, j, frame);
							img1.setRoi(this_outline_roi);
							bandp = narrow.straighten(img1, this_outline_roi, 1);
							//if (j==2)
							//	IJ.showMessage(IJ.d2s(this_outline_roi.getLength(),3)+"/"+IJ.d2s(outline_points,3));
							contour[j] = this_outline_roi.getLength();
						} else {
							contour[j] = 0;
						}
					}

					// store crossing xy point and outline point index in filopod data
					for (var j = 0; j<outline_points; j++) {
						var outline_x = interp_polygon.xpoints[j]*dx;
						var outline_y = interp_polygon.ypoints[j]*dx;

						// iterate over filopodia
						for (var k = 0; k < fparray.length; k++) {
							//IJ.showMessage("Index "+IJ.d2s(i,0)+" ? "+IJ.d2s(fparray[k].outline_index,0));
							if (i == fparray[k].cell_index && colinear(
								fparray[k].x, fparray[k].y, 
								fparray[k].cell_x, fparray[k].cell_y,
								outline_x, outline_y)) {
									var old_cross_x = fparray[k].cross_x;
									var old_cross_y = fparray[k].cross_y;
									var old_extension = fparray[k].extension();

									fparray[k].cross_x = outline_x;
									fparray[k].cross_y = outline_y;
									//if (fparray[k].outline_index >= 0)
									//	IJ.showMessage(IJ.d2s(fparray[k].extension(),digits)+" vs "+
									//		IJ.d2s(old_extension,digits)+" at k="+IJ.d2s(k,0));
									
									if (fparray[k].outline_index < 0) {
										fparray[k].outline_index = j; // set the index for the first time
										regcross++; // only count once per fp
									} else if (fparray[k].extension() < old_extension) {
										fparray[k].outline_index = j; // save the new index
									} else {
										// revert to previously detected crossing point and index
										fparray[k].cross_x = old_cross_x;
										fparray[k].cross_y = old_cross_y;	
									}

									if (!fparray[k].colinear())
										IJ.error("Seven.js", "Failed to register filopod #"+IJ.d2s(k,0));
							}
								
						}
					}

					// check for filopodia without crossing assignments
					var noncrossing = 0;
					for (var k = 0; k < fparray.length; k++) {
						if  (i == fparray[k].cell_index && fparray[k].outline_index <0) {
						//	IJ.showMessage("Cell ("+IJ.d2s(fparray[k].cell_x,2)+", "+
						//		IJ.d2s(fparray[k].cell_y,2)+"); Tip ("+IJ.d2s(fparray[k].x,2)+
						//		", "+IJ.d2s(fparray[k].y,2)+")");
							noncrossing++;
						}
					}
					if (noncrossing > 0)
						IJ.showMessage("missing assignments: "+IJ.d2s(noncrossing,0));	
					
					// Generate banded image from original 
					img1.setRoi(outline_roi);
					bandp = narrow.straighten(img1, outline_roi, boxheight);
					band = new ImagePlus("Perimeter band", bandp);
					var bandperimeter = outline_roi.getLength(); // maximum perimeter distance in um
					//var bandperimeter = band.getWidth()*dx; // alternate definition
				
	 				// Set intensity threshold
					var bandnoise = noise; // defined above based on whole image 
					var bandx = new Array(); // points of interest on the band (1-D position in um)
					var filo_spacing = true;
					var fps = 0;

					if (filo_spacing) {
						// set up for the main action found after the end of this block
						bandx = new Array(fp);
						for (var k = 0; k < fparray.length; k++) {
							// parametric crossing position in pixels
							if (i == fparray[k].cell_index) {
								bandx[fps] = contour[fparray[k].outline_index];
								fps++;
							}
						}
					} else {	
						// Find maxima 
						var bandoptions = "noise="+IJ.d2s(bandnoise,0)+" output=[Point Selection] exclude";
						IJ.run(band, "Find Maxima...", bandoptions);
						var bandpoints = band.getRoi(); 
						if (bandpoints != null && bandpoints.getPolygon() != null)
							bandx = bandpoints.getPolygon().xpoints; 
						
						if (bandx != null && bandx.length > 0) {
							var radstep = 2*Math.PI/bandperimeter*interpolation*dx; // parametric distance of one pixel in radians
							var skewx = skewness(bandx); 
							var neighborx = neighbor(bandx, bandperimeter, 1, true);	
							var leftneighborx = neighbor(bandx, bandperimeter, 1, false);	
							var neighbormeanx = 0;  
							var meanx = bandperimeter/bandx.length; // perimeter in um divided by # detection events
							
							for (var j = 0;j<neighborx.length; j++) { 
								crossingtab.incrementCounter(); 
	
								// report data to table
								crossingtab.addValue("Spacing (Delta um)", neighborx[j]);
								crossingtab.addValue("Left Spacing (Delta um)", leftneighborx[j]);
								//crossingtab.addValue("Crossing angle (deg)", cross_angle.deg);
	
								// accumulate values for averaging
								neighbormeanx += neighborx[j]; 
							} 
							// calculate average
							neighbormeanx /= neighborx.length; 
						}
					}

					// main filo analysis
					if (bandx != null && bandx.length > 0) {
						// start of neighbor analysis
						var radstep = 2*Math.PI/bandperimeter*interpolation*dx; // parametric distance of one pixel in radians
						
						if (DEBUG) {
							// report the mapping from raw pixel (xy, before interpolation) to contour (parametric)
							for (var k = 0; k < fparray.length; k++) {
								if (i == fparray[k].cell_index) {
									// Load the parametric crossing position in um from the fparray
									//temparray[fps] = contour[fparray[k].outline_index]; fps++;
									// ----> already done in previous block
										IJ.showMessage("mapping x (real px) -> u (contour px): \n"+
											IJ.d2s(fparray[k].outline_index*interpolation,0)+
											" -> "+IJ.d2s(contour[fparray[k].outline_index]/dx,0));
								}
							}
						}
							
						// Copy to java array for compatibility with neighbor()
						if (IJ.isJava17()) {
							var FloatArray = Java.type("float[]");
							var fp_u = new FloatArray(fps);
						} else {
							var fp_u = new Packages.java.lang.reflect.Array.newInstance(java.lang.Float, fps); 
						}
						for (var k = 0; k < fps; k++)
							fp_u[k] = bandx[k];

						var skewx = null;
						var neighborx = null;
						var leftneighborx = null;
						var neighbormeanx = 0;  
						var meanx = 0;
						if (fps > 0) {
							skewx = skewness(fp_u); 
							neighborx = neighbor(fp_u, bandperimeter, 1, true);	
							leftneighborx = neighbor(fp_u, bandperimeter, 1, false);
							meanx = bandperimeter/fp_u.length; // perimeter in um divided by # detection events
							//IJ.showMessage(IJ.d2s(neighborx.length,0)+" neighbors found in this cell");
							
							for (var k = 0;k<neighborx.length; k++)
								neighbormeanx += neighborx[k]; 
							// calculate average
							neighbormeanx /= neighborx.length; 
						}
						
						// reset fp counter
						var fps = 0;
						for (var k = 0; k < fparray.length; k++) {	
							if (i == fparray[k].cell_index) {		
								var complex_pos = new Complex(
									Math.cos(fparray[k].outline_index*radstep), 
									Math.sin(fparray[k].outline_index*radstep));
								// parametric crossing angle in range [0..2PI] radians
								var cross_angle = new Angle(complex_pos.Arg(), bandperimeter, 1); 
								
								fparray[k].cell_perimeter = rois[fparray[k].cell_index].getLength();
								fparray[k].cross_u = cross_angle.dist;
								fparray[k].cross_rad = cross_angle.rad;
								fparray[k].cross_deg = cross_angle.deg;
								if (neighborx != null) {
									fparray[k].neighbor_near_u = neighborx[fps];
									fparray[k].neighbor_left_u = leftneighborx[fps];
									fps++;							
								}

								if (DEBUG && i == 0)
									IJ.showMessage("Compare cross_u = "+IJ.d2s(cross_angle.dist,3)+
										" vs contour coordinate = "+IJ.d2s(contour[fparray[k].outline_index],3));
								}
						}
						
						// Cui .. Rice, J Chem Phys, 2002 
						// expected probability density = 2n exp(-n * 2R) where n is number density 
						// mean = integrated density = integral of -exp(-2nR) 0->inf = 1/2n 
						var expected = meanx/2;  // distance in um
						 
						// generate exponential random data 
						var simpoints = 1000; 
						var expdata = new Array(simpoints); 
						var sumdata = 0; 
						for (var j = 0; j<simpoints; j++) { 
							var expvar = 0; 
							// y = -exp(-2/meanx*R)  
							// ln (-y) = 2/meanx*R 
							//         ----> R = meanx/2*ln (-y) 
							// let y = [-1, 0]; -y = [0, 1] = random variable 
							do { 
								expvar = Math.random(); 
								expdata[j] = -expected*Math.log(expvar); // distance in um
							} 
							while (expdata[j] > bandperimeter/2); // enforce circularity 
								sumdata += expdata[j]; 
						} 
						// can calculate expected value from simulation -  
						// however this is not reliable at low density (n < 5 or so) 
						var simexpected = sumdata / simpoints; 
						var skew_neighbor = null;
						if (neighborx != null)
							skew_neighbor = skewness(neighborx); 
						var skew_sim = skewness(expdata); 
	 
						// Reporting 
						spacingtab.addValue("Number of Particles", bandx.length); 
						spacingtab.addValue("Band perimeter (um)", bandperimeter); 
						spacingtab.addValue("Band perimeter/N (um)", meanx); 
						spacingtab.addValue("Neighbor spacing (um)", neighbormeanx); 
						if (neighborx != null)
							spacingtab.addValue("Neighbor skewness", skew_neighbor); 
						spacingtab.addValue("Neighbor spacing - expected (um)", IJ.d2s(simexpected,digits)); 
						spacingtab.addValue("Neighbor skewness - expected", skew_sim);  
					} 
					
					// Analyze banded image intensity 
					IJ.run(band, "Select All", ""); 
					var bs = band.getStatistics(MEASUREMENTS); 
					cell_area_band[i] = bs.area; 
					//cell_band[i] = bs.mean-bgval; // calculated below from peak cortex intensity
					ScanLine(band, anadir, "spacing", dx, i); 

					var bandprofile = new Packages.ij.gui.ProfilePlot(band);
					var borderpixels = bandprofile.getProfile(); 
					var nPixels = borderpixels.length;
					var gfpval = 200; // ? standardize intensity weights - need to determine suitable value
					var borderweights = new Array(nPixels); 

					// use intensity weights
					for (var j = 0; j<nPixels; j++) { 
						// perform background subtraction &
						// normalize intensity weights
						borderweights[j] = (borderpixels[j] - bgval) / gfpval;
					}
					
					// smoothing
					var USESMOOTH = false;	
					if (USESMOOTH) {
						var smoothwidth_um = 1; // width of the smoothing window in micrometers
						var smoothwidth = Math.ceil(smoothwidth_um/dx); // width in pixels
						borderpixels = smooth(borderpixels, smoothwidth);
					}

					// get min and max of pixel intensities
					var minval = Math.pow(2,16)-1;
					var maxval = 0;
					var xmin = 0; // position of min intensity in pixels
					var xmax = 0; // position of max intensity in pixels
					for (var j = 0; j<nPixels; j++) {
						if (borderpixels[j] < minval) {
							minval = borderpixels[j];
							xmin = j;
						} else if (borderpixels[j] > maxval) {
							maxval = borderpixels[j];
							xmax = j;
						}
					}
					cell_band[i] = maxval;

					// rescale intensity values to conservatively estimate sample size
					var sumweights = 0;
					var mx1 = new Complex(0, 0);
					var mx2 = new Complex(0, 0);
					var kappa = 4; // width parameter of von Mises distribution = 1/dispersion
					var bessel = 11.30192; // besseli(4,0) - 95% density from -1.08 to +1.08 rad = 124 deg
					var radstep = 2*Math.PI/nPixels;
					
					for (var j = 0; j<nPixels; j++) {
						var cj = new Complex(Math.cos(j*radstep), Math.sin(j*radstep));
						var cj2 = cj.RealPow(2);
						
						// Simulate von Mises distribution
						//borderweights[j] = Math.exp(kappa * Math.cos(j*radstep)) / (2*Math.PI*bessel);
						
						mx1 = mx1.Add(cj.Scale(borderweights[j]));
						mx2 = mx2.Add(cj2.Scale(borderweights[j]));
						sumweights += borderweights[j];
					}
					mx1 = mx1.Scale(1/sumweights); // first moment vector
					mx2 = mx2.Scale(1/sumweights); // second moment vector
					var mxangle = new Angle(mx1.Arg(), nPixels, dx); // mean angle in range [-PI..PI] radians
					//cell_band[i] = borderpixels[mxangle.pixel]; // alt def = intensity at mean angle

					// Confidence limits of the mean
					//var chi2CI = 1; // 68% CI; chi2(1, 1) ~ 0.68
					//var chi2CI = 2.706; // 90% CI; chi2inv(.90, 1) ~ 2.706
					var chi2CI = 3.841; // 95% CI
					var lim = Math.acos(
						Math.sqrt(2*sumweights*(2*Math.pow(sumweights*mx1.Abs(),2)-sumweights*chi2CI)/
						(4*sumweights-chi2CI))/
						(sumweights*mx1.Abs()));
					var mxanglelower = new Angle(mx1.Arg() - lim, nPixels, dx);
					var mxangleupper = new Angle(mx1.Arg() + lim, nPixels, dx);

					// Circular & angular variance
					var mxanglevar = 2*(1 - mx1.Abs()); // range [0..2]
					var mxangledev = 0;
					var mxvar = 0; // range [0..Inf]
					var mxdev = 0;
					
					if (mx1.Abs() > 0) {
						mxvar = -2*Math.log(mx1.Abs());
					}
					if (mxvar >= 0 && mxanglevar >= 0) {
						mxdev = Math.sqrt(mxvar);
						mxangledev = Math.sqrt(mxanglevar);
					}

					// Circular dispersion
					var mxdisp = 0;
					if (mx1.Abs() > 0) {
						mxdisp = (1 - mx2.Abs())/(2*mx1.Abs());
					}

					// Pewsey skewness statistic
					var pewseyskew = 0;
					for (var j = 0; j<nPixels; j++) {
						pewseyskew += borderpixels[j] * 
							Math.sin(2*ComplexDist(j*radstep, mxangle.rad));
					}
					pewseyskew /= sumweights;
					
					// Cicular skewness - see Statistical analysis of circular data, Fisher, p. 34
					var mxskew = mx2.real * Math.sin(ComplexDist(mx2.imag, 2*mx1.Arg())) 
						/ Math.pow((1 - mx1.Abs()), 1.5);						

					// Rayleigh test, Zar 2010 eq. 27.4
					//var radicand = 1 + 4*(sumweights + Math.pow(sumweights,2) - Math.pow(mx1.Abs(),2));
					//var rayleighprob = 1;
					//if (radicand >= 0)
					//	rayleighprob = Math.exp(Math.sqrt(radicand)	- (1 + 2*sumweights));

					// Hodges-Ajne test, Zar 2010, eq. 27.8
					var m = nPixels;
					for (var j = 0; j<nPixels; j++) {
						var halfcircle = Math.ceil(nPixels/2);
						var temp = 0;
						for (var k = 0; k<halfcircle; k++) {
							var ind = j+k;
							if (ind >= nPixels)
								ind -= nPixels;
							temp += borderpixels[ind];
						}
						if (temp < m)
							m = temp;
					}
					var hodgesA = Math.PI *  Math.sqrt(sumweights) / (2 * (sumweights - 2*m));
					var hodgesprob = Math.sqrt(2*Math.PI) / hodgesA * 
						Math.exp(-Math.pow(Math.PI, 2) / (8 * Math.pow(hodgesA, 2)));
						
					// Intensity-weighted mean position and angle on cell perimeter
					celltab.addValue("Cell Perimeter (um)", nPixels*dx);
					celltab.addValue("Lower Position (um)", mxanglelower.dist);
					celltab.addValue("Mean Position (um)", mxangle.dist);
					celltab.addValue("Upper Position (um)", mxangleupper.dist);
					celltab.addValue("Mean Angle (deg)", mxangle.deg);
					celltab.addValue("Width (deg)", lim*2*180/Math.PI); // on 2 sides of C.I. defined by +/- lim
					celltab.addValue("Peak Position (um)", xmax*dx);
					celltab.addValue("Dispersion", mxdisp);
					//celltab.addValue("Angular Variance (rad^2)", mxanglevar);
					//celltab.addValue("Angular Deviation (deg)", mxangledev*180/Math.PI);
					//celltab.addValue("Circular Variance (rad^2)", mxvar);
					//celltab.addValue("Circular standard deviation (deg)", mxdev*180/Math.PI);
					//celltab.addValue("Circular Skewness", mxskew);
					//celltab.addValue("Pewsey Skewness", pewseyskew);
					//celltab.addValue("Rayleigh probability", rayleighprob);
					//celltab.addValue("Hodges-Ajne probability", hodgesprob);
					
			    	// Report band intensity 
					celltab.addValue("Cortex Intensity", cell_band[i]); 
					//celltab.addValue("Band Area (um^2)", cell_area_band[i]); 
					celltab.addValue("Band Variance", Math.pow(bs.stdDev, 2));
					celltab.addValue("Band Skewness", bs.skewness);
					celltab.addValue("Band Kurtosis", bs.kurtosis);
					 
					// clean up 
					saveImage(band, format, anadir, "spacing-"+i, anaversion); 
					band.changes = false; 
					if (!DEBUG) { band.close(); } 
				} else {
					IJ.showMessage("Cell area < minimum area: cell ID "+IJ.d2s(i,0));
				}
			} 
			//standardize(cell_band, cell_body, cell_area_band, cell_area_body);  
		} 
		
		// Clean up from first pass 
		img1.changes = false; 
		if (!DEBUG) { img1.close(); } 
 
		if (!skiptips) { 
			if (regcross != regfp)
				IJ.error("Seven.js", "Number of crossing points ("+
					IJ.d2s(regcross,0)+") != Number of registered filopodia ("+
					IJ.d2s(regfp,0)+")"); 

			IJ.log("Cell count = "+IJ.d2s(nCells,0)); 
			IJ.log("Raw FP count = "+IJ.d2s(nRawTips,0)); 
			IJ.log("FP count = "+IJ.d2s(fp,0)); 
			IJ.log("Registered FP count = "+IJ.d2s(regfp,0)); 
			IJ.log("Cells with FP = "+IJ.d2s(cells_with_fp,0) 
				+", fraction = "+IJ.d2s(cells_with_fp/nCells,2) 
				+" - summary follows:"); 
	 
			imagetab.addValue("Directory",imagefile.getParent()); 
			imagetab.addValue("Filename",imagefile.getName()); 
			imagetab.addValue("Number of cells", nCells); 
			imagetab.addValue("Raw number of filopodia", nRawTips); 
			imagetab.addValue("Number of tips before registration", fp); 
			imagetab.addValue("Number of filopodia registered", regfp); 
			imagetab.addValue("Number of cells with registered filopodia", cells_with_fp); 
			imagetab.addValue("Fraction of cells with registered filopodia", 
				IJ.d2s(cells_with_fp/nCells,3)); 
			imagetab.addValue("Frame time (s)", IJ.d2s(timestep, digits)); 
			imagetab.addValue("Elapsed time (s)", IJ.d2s(timestep*nFrames, digits)); 
			//imagetab.addValue("Cortex:Cell Ratio", IJ.d2s(cell_band/cell_body, digits)); 

			// Log per-cell information
			//IJ.log("\tCell ID\tNumber of Tips\tMean Cell Intensity\tMean Tip Intensity");	 
			
			// build results table
			for (var u=0; u<nCells; u++) { 
				cell_tipavg[u] = 0; 
				for (var v= 0; v < fp; v++) { 
					if (cell_per_fp[v] == u) { 
						cell_tipavg[u] += intensity_per_fp[v]; 
					} 
				} 
				if (fp_per_cell[u] > 0) { 
					cell_tipavg[u] /= fp_per_cell[u]; 
				} 
				IJ.log("Cell ID = "+IJ.d2s(u,0)+" Number of Tip = " 
					+ IJ.d2s(fp_per_cell[u],0)+" Mean Cell Intensity = " 
					+ IJ.d2s(cell_body[u],0)+"Mean Tip Intensity = " 
					+ IJ.d2s(cell_tipavg[u],0)); 
					 
	            celltab.setValue("Cell X (um)", u, cell_xpos[u]*dx); 
	            celltab.setValue("Cell Y (um)", u, cell_ypos[u]*dx); 
				celltab.setValue("Number of filopodia", u, fp_per_cell[u]); 
				celltab.setValue("Mean filopod tip intensity", u, cell_tipavg[u]); 
				celltab.setValue("Tip:Cell Ratio", u, cell_tipavg[u]/cell_body[u]); 
				celltab.setValue("Tip:Cortex Ratio", u, cell_tipavg[u]/cell_band[u]); 
				celltab.setValue("Cortex:Cell Ratio", u, cell_band[u]/cell_body[u]); 
			} 
			 
			IJ.log("Filopod information (summary follows)"); 
			var reg = 0; 
			var filotab = new Packages.ij.measure.ResultsTable(); 
			filotab.setPrecision(digits);

			// log per-tip information
		    //IJ.log("\tTip ID\tTip Intensity\tCell Intensity\tCell ID\tTip:Body Ratio\tTip:Band Ratio");

			// Build results table
			for (var v = 0; v < fparray.length; v++) { 
			    filotab.incrementCounter(); 
			    var cellid = fparray[v].cell_index; 
			    if (cellid >= 0) { 
			    	reg++; 
				    var tipbody = fparray[v].intensity / cell_body[cellid]; 
				    var tipband = fparray[v].intensity / cell_band[cellid]; 
				 
				    IJ.log("\t" + IJ.d2s(reg, 0) + "\t" 
				    	+ IJ.d2s(fparray[v].intensity, 0) + "\t" 
			            + IJ.d2s(cell_body[cellid],0) + "\t" 
			            + IJ.d2s(cellid, 0) + "\t" 
			            + IJ.d2s(tipbody, digits) + "\t" 
			            + IJ.d2s(tipband, digits)); 
			             
		            filotab.addValue("Filopod ID+1", v+1); 
		            filotab.addValue("Cell ID+1", cellid+1); 
		            filotab.addValue("Tip X (um)", fparray[v].x); 
		            filotab.addValue("Tip Y (um)", fparray[v].y); 
		            filotab.addValue("Cross X (um)", fparray[v].cross_x); 
		            filotab.addValue("Cross Y (um)", fparray[v].cross_y); 
		            filotab.addValue("Extension Length (um)", fparray[v].extension()); 
		            filotab.addValue("Parametric Distance (um))", fparray[v].cross_u);
		            filotab.addValue("Parametric Angle (deg))", fparray[v].cross_deg);
		            filotab.addValue("Cell X (um)", fparray[v].cell_x); 
		            filotab.addValue("Cell Y (um)", fparray[v].cell_y); 
		            filotab.addValue("Perimeter (um)", fparray[v].cell_perimeter); 
		            filotab.addValue("Nearest neighbor (Delta um)", fparray[v].neighbor_near_u); 
		            filotab.addValue("Nearest neighbor (Delta deg)", fparray[v].neighbor_near_deg()); 
		            filotab.addValue("Left neighbor (Delta um)", fparray[v].neighbor_left_u); 
		            filotab.addValue("Left neighbor (Delta deg)", fparray[v].neighbor_left_deg()); 
		            filotab.addValue("Tip intensity", fparray[v].intensity); 
		            filotab.addValue("Cell intensity", cell_body[cellid]); 
		            filotab.addValue("Tip:Cell Ratio", tipbody); 
		            filotab.addValue("Tip:Band Ratio", tipband); 
			    } 
			}
 
			// Save the results table as a delimited text file 
			celltab.save(anadir+"analysis-results-"+anaversion+"."+resultsformat); 
			filotab.save(anadir+"filopod-results-"+anaversion+"."+resultsformat); 
		} 
 
		// Save marked img2 
		saveImage(img2, format, anadir, "Analyzed-"+anaversion, nCells); 
 
		// Save the message log as a text file
		var textfile = anadir+"analysis-log-"+anaversion+".txt"; 
		var text = IJ.getLog(); 
		saveText(textfile, text, false); 
 
		// Save the results table as a delimited text file 
		spacingtab.save(anadir+"spacing-results-"+anaversion+"."+resultsformat); 
		 
		// Clean up from second pass 
		img2.changes = false; 
		if (!DEBUG) { img2.close(); } 
 
		// Clean up stray images - copies of img1 that don't close due to an unidentified issue
		var strayimg = WindowManager.getCurrentImage(); 
		if (strayimg != null && !DEBUG) { strayimg.close(); } 
	} 
} 

function seven_multi(root, acqname) { 
	var dirlist = root.list(); 
	var acqlen = acqname.length; 
	var prefix = "sliced_frame"+IJ.pad(frame,4)+"_"; 
	var linescanonly = false; //Skip full analysis and process linescans instead 
	// To compare multiple thresholding methods, specify an array: 
	//var thresholds = [-1,0,2000, 1400, 800]; 
	// Otherwise use the global default value: 
	var thresholds = [thresholdcode]; 
	var imagetab = new Packages.ij.measure.ResultsTable(); 
	imagetab.setPrecision(digits); 
 
	// Create threshold files if not already present 
	ThresholdCells(root, acqname, thresholds, prefix); 
 
	// Repeat the listing to discover newly created directories 
	dirlist = root.list(); 
 
	// Process thresholded images 
	for (var th = 0; th < thresholds.length; th++) { 
		for (var i = 0; i<dirlist.length; i++) { 
			var subdirname = dirlist[i]; 
			var subdir = new File(root.getCanonicalPath(), subdirname); 
			if (subdir.isDirectory()) { 
				// Check if the directory matches the threshold configuration 
				// Painfully written because java strings and JavaScript "Strings" are mixed 
				var process = false; 
				var thstring= IJ.d2s(thresholds[th],0); 
				var minlength = prefix.length + thstring.length(); 
				if (Packages.java.lang.Integer.parseInt(subdirname.length()) > minlength) { 
					var foundthreshold = subdirname.substring(prefix.length, minlength); 
					if (Packages.java.lang.Integer.parseInt(foundthreshold) == thresholds[th]) { 
						process = true; 
					} 
				} 
				if (process) { 
					var anadir = subdir.getCanonicalPath()+sep; 
 
					var sublist = subdir.list(); 
					var imagefile = null; 
 
					for (var j=0; j<sublist.length; j++) { 
						var temp = sublist[j]; 
						if (linescanonly && Packages.java.lang.Integer.parseInt(temp.length()) > 17) { 
							if (temp.substring(0,17) == "Capture-mask-body") { 
								imagefile = new File(subdir.getCanonicalPath(), sublist[j]); 
							} 
						} else { 
							var hint = readPaths(subdir.getCanonicalPath()+sep+"hint.txt"); 
							if (hint.length == 1) { 
								var imagename = hint[0]; 
								imagefile = new File(root.getCanonicalPath(), imagename); 
							} 
						} 
					} 
 
					if (imagefile != null) { 
 
						if (linescanonly) { 
							seven_scans(imagefile, anadir); 
						} else { 
							seven_run(imagefile, anadir, imagetab); 
						} 
					} 
				} 
			} 
		} 
	} 
	// Save the results table as a delimited text file 
	imagetab.save(root+sep+"experiment-results-"+anaversion+"."+resultsformat); 
} 
 
function seven_scans(imagefile, anadir) { 
	ClearLog(); 
 	var img = IJ.openImage(img);
	if (img.getNSlices() >= frame)
		img.setSlice(frame); 
		
	AnalyzeScans(img, imagefile, anadir, boxwidth_um); 
} 
 
function seven_run(imagefile, anadir, imagetab) { 
	if (imagefile.isFile()) {
		ClearLog(); 
		IJ.log(imagefile); 
		imagetab.incrementCounter(); 
		 
		// Name the image img0 because it will be duplicated for each function call 
		var img0 = IJ.openImage(imagefile);
 		if (img0.getNSlices() >= frame)
 			img0.setSlice(frame);
	 
		// Analyze tips (first pass) 
		AnalyzeTips(img0.duplicate(), imagefile, anadir, imagetab, boxwidth_um, true); 
	 
		// analyze cell body intensity 
		// passing of global arrays by reference probably not necessary & somewhat confusing?
		AnalyzeCells(img0.duplicate(), anadir, 0, "body", cell_body, cell_area_body, cell_xpos, cell_ypos); 
		AnalyzeScans(img0.duplicate(), imagefile, anadir, boxwidth_um); 
	 
		// Analyze tips (second pass) 
		AnalyzeTips(img0.duplicate(), imagefile, anadir, imagetab, boxwidth_um, false); 
			 
		// Clean up 
		img0.changes = false; 
		if (!DEBUG) { img0.close(); } 
	}
} 

function Filopod(x, y, cell_x, cell_y, cell_index, intensity) {
	// object to store filopod data 
	// default values assigned at initialization
	this.x = x;
	this.y = y;
	this.cell_x = cell_x;
	this.cell_y = cell_y;
	this.cell_index = cell_index;
	this.intensity = intensity;
	if (x<0 || y<0 || cell_x<0 || cell_y<0 || cell_index<0 || intensity<0 )
		IJ.error("Seven.js", "invalid Filopod() initialized with negative values");

	// default values for properties with true values assigned later
	this.cell_perimeter = 0;
	this.outline_index = -1;
	this.cross_x = -1;
	this.cross_y = -1;
	this.cross_u = -1;
	this.cross_rad = NaN;
	this.cross_deg = NaN;
	this.neighbor_near_u = -1;
	this.neighbor_left_u = -1;

	// calculate filopod length/extension distance
	// (cannot call it length because length is a reserved word)
	this.extension = function () {
		var disp_x = this.x - this.cross_x;
		var disp_y = this.y - this.cross_y;
		return Math.sqrt(Math.pow(disp_x, 2) + Math.pow(disp_y, 2));
	}

	// calculate parametric angles in radians	
	this.neighbor_near_rad = function () {
		if (this.neighbor_near_u >= 0 && this.cell_perimeter > 0)
			return this.neighbor_near_u/this.cell_perimeter*2*Math.PI;
		else
			return NaN;
	}
	
	this.neighbor_left_rad = function () {
		if (this.neighbor_left_u >= 0 && this.cell_perimeter > 0)
			return this.neighbor_left_u/this.cell_perimeter*2*Math.PI;
		else
			return NaN;
	}
	
	// calculate parametric angles in degrees (evaluate at call time)
	this.neighbor_near_deg = function () {
		return this.neighbor_near_rad() * 180 / Math.PI;
	}
	this.neighbor_left_deg = function () {
		return this.neighbor_left_rad() * 180 / Math.PI;
	}
	
	// consistency check to verify colinearity
	this.colinear = function() {
		var tolerance = 0.4;
		var theta = Math.PI/4; // 45 degrees
		var ratio = colinear2(this.x, this.y, this.cell_x, this.cell_y, this.cross_x, this.cross_y,
			Math.cos(theta), Math.sin(theta));
		//var ratio2 = colinear2(this.cell_x, this.cell_y, this.x, this.y, this.cross_x, this.cross_y,
		//	Math.cos(theta), Math.sin(theta));
		//var ratio3 = colinear2(this.cross_x, this.cross_y, this.cell_x, this.cell_y, this.x, this.y,
		//	Math.cos(theta), Math.sin(theta));
		
		//if (Math.abs(ratio-1) < tolerance || Math.abs(ratio2-1) < tolerance || Math.abs(ratio3-1) < tolerance)
		if (Math.abs(ratio-1) < tolerance)
			return true;
		else
			return false;
	}

}

function colinear(x0, y0, x1, y1, x2, y2) {
	var tolerance = 0.4;
	var theta = Math.PI/4; // 45 degrees
	var ratio = colinear2(x0, y0, x1, y1, x2, y2, Math.cos(theta), Math.sin(theta));
	//var ratio2 = colinear2(x1, y1, x0, y0, x2, y2, costh, sinth);
	//var ratio3 = colinear2(x2, y2, x1, y1, x0, y0, costh, sinth);
	
	//if (Math.abs(ratio-1) < tolerance || Math.abs(ratio2-1) < tolerance || Math.abs(ratio3-1) < tolerance)
	if (Math.abs(ratio-1) < tolerance)
		return true;
	else
		return false;
}

function colinear2(x0, y0, x1, y1, x2, y2, costh, sinth) {
	var errval = 10;
	var xr0 = rotx(x0, y0, costh, sinth);
	var yr0 = roty(x0, y0, costh, sinth);
	var xr1 = rotx(x1, y1, costh, sinth);
	var yr1 = roty(x1, y1, costh, sinth);
	var xr2 = rotx(x2, y2, costh, sinth);
	var yr2 = roty(x2, y2, costh, sinth);

	if ((x0 == x1 && x0 == x2) || (y0 == y1 && y0 == y2) ||
		(xr0 == xr1 && xr0 == xr2) || (yr0 == yr1 && yr0 == yr2))
		return 1; // the values are colinear

	var denom = ((y1-y0)*(x2-x0));
	var denom_rot = ((yr1-yr0)*(xr2-xr0));
	var numerator = (denom == 0) ? ((yr2-yr0)*(xr1-xr0)) : ((y2-y0)*(x1-x0));
	
	if (denom != 0)
		return numerator / denom; // exit and return the ratio using untransformed coordinates

	if (denom_rot != 0)
			return numerator / denom_rot; // exit and return the ratio using transformed coordinates

	return errval; // exit and return the error value >> 1
}

function rotx(x, y, costheta, sintheta) {
	return x * costheta + y * sintheta;
}

function roty(x, y, costheta, sintheta) {
	return -x * sintheta + y * costheta;
}

function lookup(cumulative, input) {
	var output = null;
	var i = 0;
	
	while (output == null && i < cumulative.length)
		if (cumulative[i] >= input)
			output = i;

	return output;
}

function outliner(polygon, npoints, frame) {				
	var outline = new Polygon(); 
	if (npoints > polygon.xpoints.length)
		IJ.error("Seven.js", "Invalid argument to outliner(): npoints too large");
	
	for (var j = 0; j<npoints; j++)
		outline.addPoint(polygon.xpoints[j], polygon.ypoints[j]);
 
	// By default the band ROI has "C" shape. Now correct this to an "O" shape:
	if (npoints == polygon.xpoints.length)
		outline.addPoint(polygon.xpoints[0], polygon.ypoints[0]); 
	var outline_roi = new Packages.ij.gui.PolygonRoi(outline,  
		Packages.ij.gui.Roi.FREELINE); 
	outline_roi.setPosition(1, frame, 1); // specify which frame (slice) to analyze

	return outline_roi;
} 

function ClearLog() { 
	if (IJ.getLog() != null) { 
		IJ.selectWindow("Log"); 
		IJ.run("Close"); 
	} 
} 
 
// Save a copy of a file using streams 
function CopyFile(datafile, outputdir, outputname) { 
	var copyfile = new File(outputdir+sep+outputname); 
	if (!copyfile.exists()) { 
		copyfile.createNewFile(); 
		var instream = FileInputStream(datafile); 
		var outstream = FileOutputStream(copyfile); 
		var buf = new Packages.java.lang.reflect.Array.newInstance(java.lang.Byte.TYPE, 1024); 
		var bytesRead = 0; 
		while ((bytesRead = instream.read(buf)) > 0) { 
	            outstream.write(buf, 0, bytesRead); 
	    } 
	    instream.close(); 
	    outstream.close(); 
	} 
} 
 
function StoreRoi (img, roi, name, roiset) { 
	roiset.add(img,roi,-1); 
	roiset.select(roiset.getCount()-1); 
	roiset.runCommand("Rename", name); 
} 
 
function readBetween(logpath, path, startstring, endstring) { 
	var datafile = new File(path); 
	var thisLine = null; 
 
     // open file as input stream 
     var fr = new FileReader(path); 
     var br = new BufferedReader(fr); 
     var startread = false; 
     var endread = false; 
     var text = ""; 
 
     // bypass mode (read entire stream) 
     if (startstring == null || endstring == null) { 
     	startread = true; 
     } 
 
	// read from stream 
     while ((thisLine = br.readLine()) != null) { 
     	// check for end 
     	if (endstring != null && endstring.length <= thisLine.length() && 
     		thisLine.substring(0,endstring.length) == endstring) { 
     		endread = true; 
     	} 
 
     	// read a line and write it out 
     	if (startread && !endread) { 
     		text = thisLine+newline; 
 
			if (text.length > 0) { 
				IJ.log(text); 
				if (logpath != null) { 
					saveText(logpath, text, true); 
				} 
			} 
     	} 
 
     	// check for start 
     	if (startstring != null && startstring.length <= thisLine.length() && 
     		thisLine.substring(0,startstring.length) == startstring) { 
     		startread = true; 
     	} 
     } 
} 
 
// read a numeric value from a text file 
function readValue(textfile, logline, logstring) { 
	// example of how to parse files for numeric values 
	var val = 0; 
	if (textfile != null && logline >= 0 && logstring != null) { 
		var data = readPaths(textfile); 
		if (data.length > logline) { 
			var row = data[logline]; 
			if (row.length() >= logstring.length && 
				row.substring(0, logstring.length) == logstring) { 
				// read the rest of the line following the logstring 
				var temp = Packages.java.lang.Float.parseFloat(row.substring(logstring.length,row.length())); 
				if (temp != null && temp > 0) { 
					val = temp; 
				} 
			} 
		} 
	} 
	return val; 
} 
 
// read text file as array, one string per line 
function readPaths(list) { 
     var fr = new FileReader(list); 
     var br = new BufferedReader(fr); 
	 var i = 0; 
	 var thisLine = null; 
	 var dirs = []; 
 
	 // read stream into array 
	 while((thisLine = br.readLine()) != null) { 
	 	dirs[i] = thisLine; 
	 	i++; 
	 } 
 
	// return a string array 
	 return dirs; 
} 
 
// Strict number parsing function from dev.mozilla.org 
filterInt = function (value) { 
  if(/^(\-|\+)?([0-9]+|Infinity)$/.test(value)) 
    return Number(value); 
  return NaN; 
} 
 
// get pixel intensity as 24-bit value 
function getValue(img, x, y) { 
	return img.getProcessor().getPixel(x, y);
} 
 
function dotnote(a, b) { 
	// returns a string of the form a.b 
	a = Math.abs(a); 
	b = Math.abs(b); 
	var c = Math.ceil(Math.log(b)/Math.log(10)); 
	var d = a + b/Math.pow(10, c); 
	var note = IJ.d2s(d, c); 
 
	return note; 
} 
 
// Standardize by area = outer-inner 
function standardize(outermean, innermean, outerarea, innerarea) { 
    for (var i = 0; i < outermean.length; i++) { 
		// check for bad data 
		if (outerarea[i] > 0 && innerarea[i] > 0) { 
			// standardize using the following formula: 
			// outermean = (outerint - innerint)/(outerarea - innerarea) 
			// = (outermean*outerarea - innermean*innerarea)/(outerarea - innerarea) 
			outermean[i] = (outermean[i]*outerarea[i] - innermean[i]*innerarea[i]) 
				/ (outerarea[i]-innerarea[i]); 
			outerarea[i] -= innerarea[i]; 
		} else { 
			outermean[i] = 0; 
		} 
    } 
} 
 
function invertImage(img) { 
	// This operation is not reversible - it resets the black level! 
    IJ.run(img, "Invert", "stack"); 
} 
 
function setMask(img, depth, invert) { 
	// manipulate a binary mask for 16-bit  
	// image operations 
    // CAUTION - binary images only! 
 
	var direction = (depth > 0) ? "Dilate" : "Erode"; 
 
	for (var i = 0; i < Math.abs(depth) ; i++) 
	    IJ.run(img, direction, "stack"); 
 
	// Scale 8-bit mask up to 16-bit 
	if (img.getType() == ImagePlus.GRAY8) { 
	    IJ.run(img, "16-bit", ""); 
	    IJ.run(img, "Multiply...", "value=1000 stack"); 
	} 
	 
	// Optionally invert the mask 
	if (invert) 
	    invertImage(img); 
} 
 
function saveText(path, text, append) { 
	var wr = new FileWriter(path, append); 
	var pr = new PrintWriter(wr); 
	pr.print(text); 
	pr.close(); 
} 
 
function openIf(file, format) { 
	var opener = new Packages.ij.io.Opener(); 
	var image = null; 
	if (format == "zip") { 
		image = opener.openZip(file); 
	} else { 
		image = opener.openImage(file); 
	} 
	if (image != null && image.getNSlices() >= frame)
		image.setSlice(frame); 
	return image;
} 
 
function saveIf(img, format, path) { 
	// save image one time only 
	var f = new File(path); 
	if (!f.exists()) 
		IJ.saveAs(img, format, path); 
} 
 
function saveImage(img, format, dir, name, n) { 
	// save image with a custom file name 
	var s = n; 
	if (typeof n == "number") 
		s = IJ.d2s(Math.abs(n), 0); 
	var path = dir+sep+name+"-"+s+"."+format; 
 
	if (DEBUG) { IJ.showMessage(img.getTitle()); } 
	saveIf(img, format, path); 
	img.show(); 
} 
 
function getFileList(dir) { 
    var f = new File(dir); 
    if (!f.exists() || !f.isDirectory()) 
        return null; 
    list = f.list(); 
    if (list==null) 
        return null; 
    var f2 = new File(dir, list[0]); 
    var hidden = 0; 
    for (var i=0; i<list.length; i++) { 
        if (list[i].startsWith(".") || list[i].equals("Thumbs.db")) { 
            list[i] = null; 
            hidden++; 
        } else { 
            f2 = new File(dir, list[i]); 
            if (f2.isDirectory()) 
                list[i] = list[i] + sep; 
        } 
    } 
    var n = list.length-hidden; 
    if (n<=0) 
        return null; 
    if (hidden>0) { 
        list2 = new Array(n); 
        var j = 0; 
        for (var i=0; i<list.length; i++) { 
            if (list[i]!=null) {
                list2[j] = list[i]; 
                j++;
            }
        } 
        list = list2; 
    } 
    array = new Array(n); 
    for (var i=0; i<n; i++) 
        array[i] = new Array(list[i]); 
    return array; 
} 
 
function RoiSet() { 
	// Open the ROI manager and reset it 
	importClass(Packages.ij.plugin.frame.RoiManager); 
	var r = RoiManager.getInstance(); 
 
	if (r != null) { 
		r.close(); 
	} 
	r = new RoiManager(); 
	return r; 
} 
 
function getExt(filename) { 
	var result = null; 
	var start = filename.lastIndexOf(".")+1; 
	var end = filename.length(); 
	if (end > (start+3)) { 
		result = filename.substring(start, start+3); 
	} else { 
		result = filename.substring(start, end); 
	} 
	return result; 
} 
 
function noiseThreshold(stats, minSNR) { 
	// standard noise floor = 20 gray levels in a mix/maxed image 
	var noise = (stats.max-stats.min) * 20/255; 
	var baseline = 0; 
	if (stats.dmode > 0) { 
		baseline = stats.min; // historically used, works when dmode > 0, min > 0 
	} else { 
		baseline = stats.median; // In case black level was subtracted (min = 0) 
	} 
	if (noise < minSNR*baseline) { 
		IJ.log("Noisy image, resetting noise floor... Noise level = "+ 
			IJ.d2s(noise,0)+", median = "+IJ.d2s(stats.median,0)); 
		noise = minSNR*baseline; 
	} 
	return noise; 
} 
 
function neighbor(vararray, maxdist, scale, minimize) { 
	// this function copies an array of x values, ensures they are sorted and then
	// returns an array with the nearest neighbor distance (delta x) for each position
	// assuming a circular x coordinate.
	//
	// Input values
	// vararray : a Java array of 1-D positions in pixels
	// maxdist	: perimeter distance in pixels
	// scale	: width of a pixel in um
	// minimize	: true = return smallest distance, false = return left hand distance
	//
	// Output value is a distance in um

	// ensure the array is sorted and non-null
	var n = 0;
	for (var i = 0; i<vararray.length; i++)
		if (vararray[i] != null)
			n++;
			
	var order = new Array(n);
	var sorted = new Array(n); 
	var temparray = new Array(n); 
	n = 0;
	for (var i = 0; i<vararray.length; i++) {
		if (vararray[i] != null) {
			temparray[n] = Packages.java.lang.Float.parseFloat(vararray[i]);
			n++;
		}
	}
	
	for (var i = 0; i<n; i++) { 
		var min_val = Infinity;
		for (var j = 0; j<n; j++) {
			if (temparray[j] < min_val) {
				min_val = Packages.java.lang.Float.parseFloat(temparray[j]);
				sorted[i] = min_val;
				order[i] = j;
			}
		}
		if (order[i] != null)
			temparray[order[i]] = Infinity; // set to find the next smallest value next round
		//IJ.showMessage(IJ.d2s(order[i],0)+": "+IJ.d2s(sorted[i],3));
	}
	
	// shift the array left and right, ensuring non-negative values
	var jsarray = new Array(sorted.length); 
	var result = new Array(sorted.length); 
	var left = new Array(sorted.length); 
	var right = new Array(sorted.length); 
	for (var i = 0; i<sorted.length; i++) { 
		jsarray[i] = (sorted[i] > 0) ? sorted[i] : 0; 
		left[i] = jsarray[i];
		right[i] = jsarray[i];	
		
		if (jsarray[i] > maxdist)
			IJ.showMessage("bad maxdist parameter, "+IJ.d2s(jsarray[i],digits)+" > "+IJ.d2s(maxdist,digits));		
	} 
 
	// shift the arrays to compare left and right neighbors
	var first = left.shift(); 
	left.push(first); 
	var last = right.pop(); 
	right.unshift(last); 
	 
	// search for nearest neighbor 
	for (var i = 0; i<jsarray.length; i++) { 
			 
		// Implement circular array 
		var dleft = 0; 
		var dright = 0; 
	 
		// calculate immediate neighbors absolute distance 
		var circ = 0; // implement circularity 
		var dist = 0; // distance in um
		circ = (jsarray[i] > left[i]) ? maxdist : 0; 
		dleft = (left[i] - jsarray[i] + circ) * scale; 
		circ = (jsarray[i] < right[i]) ? maxdist : 0; 
		dright = (jsarray[i] - right[i] + circ) * scale; 
		if (minimize || dright == 0) // if the distance is zero then there is no right or left neighbor
 			dist = (dleft < dright) ? dleft : dright; 
 		else
			dist = dleft;
 
		// enforce distance > 0 
		if (jsarray[i] < 0) { 
			var minimize_string = minimize ? "True" : "False";
			IJ.showMessage("Neighbor distance failed, input nonpositive\n Minimize is "+minimize_string+"\n"+
				IJ.d2s(jsarray[i],digits)+", "+IJ.d2s(left[i],digits)+", "+IJ.d2s(maxdist,digits)); 
		} 

		// copy result to array
		// result[i] = dist; // old behavior - output in sorted order
		result[order[i]] = dist; // new behavior - output in original order		
	} 
	return result; 
} 
 
function skewness(vararray) { 
	var skew = 0; 
	var sumsq = 0; 
	var sum = 0; 
	var n = vararray.length; 
 
	if (n > 2) { 
		// Calculate standard deviation 
		for (var j=0; j<n; j++) { 
				sum += vararray[j]; 
				sumsq += Math.pow(vararray[j], 2); 
		} 
		var mean = sum/n; 
		var sdev = Math.sqrt(sumsq/n - Math.pow(mean, 2)); 
	 
		// Calculate skewness 
		for (var j=0; j<n; j++) { 
				skew += Math.pow((vararray[j] - mean)/sdev, 3); 
		} 
		skew *= n/((n-1)*(n-2)); 
	} else { 
		skew = 0; 
	} 
	return skew; 
} 
 
function kymograph(img0, rs, dt, boxwidth_um, paramtab, resultsdir) { 
	var rois = rs.getRoisAsArray(); 
	var cal = img0.getCalibration(); 
	cm2um(cal); // set micrometer scale 
	var dx = cal.getX(1); 
	var boxheight = Math.ceil(boxwidth_um/dx); // height in pixels of linescan 
 
	// set dimensions for the resized image 
	var hw = Packages.java.lang.Integer.parseInt(Math.floor(1.5*img0.getWidth())); 
	var hh = Packages.java.lang.Integer.parseInt(Math.floor(1.5*img0.getHeight())); 
	var woffset = Packages.java.lang.Integer.parseInt(Math.floor(0.25*img0.getWidth())); 
	var hoffset = Packages.java.lang.Integer.parseInt(Math.floor(0.25*img0.getHeight())); 
 
	// Resize image before loading ROIs - note the offset is added to ROIs later 
	resize(img0, hw, hh, woffset, hoffset); 
 
	for (var j=0; j < rois.length; j++) { 
		var img = img0.duplicate(); 
		if (img.getNSlices() >= frame)
 			img.setSlice(frame);
 
		img.show(); 
		img.setRoi(rois[j]); 
		rs.select(j); 
		var n = img.getSlice(); 
		var outname = img0.getTitle()+"_"+IJ.d2s(j,0)+"_"+IJ.d2s(n,0); 
		rs.runCommand("Rename", outname); 
		var line = img.getRoi(); 
		// align the filopod with x-axis and transform ROI coordinates to the final image 
		if (line.isLine() && line.getPolygon().xpoints.length == 2) { // straight line only 
			var outputtab = new Packages.ij.measure.ResultsTable(); 
			// compute angle of filopod 
			var theta_deg = line.getAngle(); 
 
			// Rotate the image 
			img.show(); 
			IJ.run(img, "Rotate... ", 
				"interpolation=Bicubic stack angle="+ 
				IJ.d2s(theta_deg, digits));			// Rotate image 
 
			// Draw bounding rectangle for 3D-kymograph in the rotated image stack 
			getRectangle(img, line, hw/2, hh/2, woffset, hoffset, boxheight, theta_deg); 
 
			// Reslice image to generate 3D-kymograph, project into 2D-kymograph 
			IJ.run(img, "Reslice [/]...", "output=1.000 slice_count=1 avoid"); // Make kymograph 
			var reslice = IJ.getImage(); 
 
			// Project max intensity kymograph across the boxheight direction 
			IJ.run(reslice, "Z Project...", "projection=[Max Intensity]"); 
			var kymograph = IJ.getImage(); 
			IJ.run(kymograph, "Rotate 90 Degrees Left", ""); 
			saveImage(kymograph, format, resultsdir, outname, j); 
			var kymograph2 = kymograph.duplicate(); 
 
			processKymograph(kymograph2); 
 
			paramtab.incrementCounter(); 
			paramtab.addValue("Filename", outname); 
			analyzeKymograph(kymograph2, outputtab, paramtab, dt, dx, digits); 

			// Save the kymograph as a delimited text file
			outputtab.save(resultsdir+sep+outname+"."+resultsformat); 
			drawRoi(kymograph, kymograph2.getRoi()); 
			saveImage(kymograph, format, resultsdir, outname+"_FIT", j); 
			saveImage(kymograph2, maskformat, resultsdir, outname, j); 
 
			reslice.changes = false; 
			kymograph.changes = false; 
			kymograph2.changes = false; 
			if (!DEBUG) { 
				reslice.close(); 
				kymograph.close(); 
				kymograph2.close(); 
			} 
		} 
		img.changes = false; 
		if (!DEBUG) { img.close(); } 
	} 
} 
 
function getRectangle(img, line, hw, hh, woffset, hoffset, rectheight, theta_deg) { 
	var theta = theta_deg*Math.PI/180.0; 
	//IJ.showMessage(IJ.d2s(line.x1,0)+","+IJ.d2s(line.y1,0)+";"+IJ.d2s(line.x2,0)+","+IJ.d2s(line.y2,0)); 
	 
	if (line.isLine() && line.getPolygon().xpoints.length == 2) { // straight line only 
		centerRoi(line, hw, hh, woffset, hoffset, 1); 
		rotateRoi(line, -theta); 
		centerRoi(line, hw, hh, 0, 0, -1); 
 
		// Calculate bounding rectangle 
		var xleft = line.x1; 
		var rectwidth = line.x2 - line.x1; 
		var yup = line.y1 - Math.floor(rectheight/2); 
 
		// Check geometry is correct 
		if (rectwidth <= 0 || rectheight <= 0) { 
			IJ.showMessage("ROI is malformed"); 
			exit; 
		} 
		//if (xleft < 0 || xleft+rectwidth >= hw || 
		//	yup < 0 || yup+rectheight >= hh) { 
		//	IJ.showMessage("ROI out of bounds"); 
		//	exit; 
		//} 
 
		img.setRoi(xleft, yup, rectwidth, rectheight); // Set bounding rectangle 
	} 
} 
 
function rotateRoi(line,theta) { 
	var x1 = line.x1 * Math.cos(theta) + line.y1 * -1 * Math.sin(theta); 
	var x2 = line.x2 * Math.cos(theta) + line.y2 * -1 * Math.sin(theta); 
	var y1 = line.x1 * Math.sin(theta) + line.y1 * Math.cos(theta); 
	var y2 = line.x2 * Math.sin(theta) + line.y2 * Math.cos(theta); 
 
	// Fix values 
	line.x1 = x1; 
	line.x2 = x2; 
	line.y1 = y1; 
	line.y2 = y2; 
} 
 
function centerRoi(line,hw,hh,woffset,hoffset,sense) { 
	var x1 = woffset + line.x1 - hw * sense; 
	var x2 = woffset + line.x2 - hw * sense; 
	var y1 = hh - (hoffset + line.y1); 
	var y2 = hh - (hoffset + line.y2); 
 
	// Fix values 
	line.x1 = x1; 
	line.x2 = x2; 
	line.y1 = y1; 
	line.y2 = y2; 
} 
 
function processKymograph(img) { 
	// Process the kymograph to extract the tip location 
	var imp = img.getProcessor(); 
	imp.findEdges(); 
	img.setProcessor(imp); 
	IJ.run(img, "Enhance Contrast", "saturated=0.35"); 
	IJ.run(img, "8-bit", ""); 
	IJ.setAutoThreshold(img, "Moments stack"); 
	IJ.run(img, "Convert to Mask", "method=Moments background=Light stack"); 
} 
 
function analyzeKymograph(img, outputtab, paramtab, dt, dx, digits) { 
	var duration = img.getWidth(); 
	var distance = img.getHeight(); 
	var maskval = 0; 
	var ts = new Array(duration); 
	var xs = new Array(duration); 
	outputtab.setPrecision(digits); 
 
	for (var k=0; k<duration; k++) { 
		var x = 0; 
		var count = 0; 
		// Scan image for masked pixels 
		for (var m=0; m<distance; m++) { 
			if (getValue(img, k, m) == maskval) { 
				x += m; 
				count++; 
			} 
		} 
		// Compute first moment of distance (= position of tip) 
		ts[k] = dt*k; 
		if (count > 0) { 
			xs[k] = dx*Math.round(x/count); 
		} else { 
			xs[k] = 0; 
		} 
 
		// Write to table 
		outputtab.incrementCounter(); 
		outputtab.addValue("Time (s)", ts[k]); 
		outputtab.addValue("Distance (um)", xs[k]); 
	} 
 
	var minsize = 5; // minimum segment length in pixels 
	var minslope = 0; // only return positive slopes 
	var minrsq = .98; 
	var params = new Array(); 
 
	params = fitKymograph(ts, xs, minsize, minslope, minrsq); 
 
	paramtab.addValue("Slope (um/s)", params[1]); 
	paramtab.addValue("RSquared", params[3]); 
 
	// Calculate ROI for fitted line in image coordinates 
	var x1 = params[4]; 
	var x2 = x1 + params[5]; 
	var y1 = (params[0]+params[1]*dt*x1)/dx; 
	var y2 = (params[0]+params[1]*dt*x2)/dx; 
	line = new Packages.ij.gui.Line(x1, y1, x2, y2); 
	img.setRoi(line); 
} 
 
function fitKymograph(x, y, minsize, minslope, minrsq) { 
	var LINREG = Packages.ij.measure.CurveFitter.STRAIGHT_LINE; 
	var maxsize = minsize; 
	// create empty array that stores the positions of the fitted segments 
	var startpoints = new Array(x.length); 
	for (var j = 0; j<x.length; j++) { 
		startpoints[j] = 0; 
	} 
 
	// Greedy search - find the longest segment with R2 > minsize 
	for (var i = minsize; i<x.length; i++) { 
		for (var j = 0; j<x.length; j++) { 
			var cf = new Packages.ij.measure.CurveFitter( 
				x.slice(j-i,j), y.slice(j-i,j)); 
			cf.doFit(LINREG); 
 
			if (cf.getRSquared() > minrsq && cf.getParams()[1] > minslope) { 
				// length of segment recorded at left index position 
				startpoints[j-i] = i; 
				maxsize = i; 
			} 
		} 
	} 
 
	// In case of ties for longest segment, use the leftmost segment 
	var linestart = 0; 
	for (var j = x.length-1; j>0; j--) { 
		if (startpoints[j] == maxsize) { 
			linestart = j; 
		} 
	} 
	var linestop = linestart+maxsize; 
	if (linestop > x.length) { 
		linestop = x.length; 
	} 
 
	// Run the fit a final time with the chosen segment 
	var cf = new Packages.ij.measure.CurveFitter( 
		x.slice(linestart,linestop), y.slice(linestart,linestop)); 
	cf.doFit(LINREG); 
 
	if (DEBUG) { 
		IJ.showMessage("Start position = "+IJ.d2s(x[linestart],digits)); 
		IJ.showMessage("Length = "+IJ.d2s(maxsize,0)); 
		IJ.showMessage(cf.getResultString()); 
	} 
 
	// return the results of the fit y = a + bx 
	// params[0] = a, the y-intercept in scaled units 
	// params[1] = b, the slope in scaled units 
	// params[2] = SSE 
	var params = cf.getParams(); 
	// CurveFitter returns a Java native array, but a JS array is more convenient 
	var result = new Array(params.length); 
	for (var k = 0; k<params.length; k++) { 
		result[k] = params[k]; 
	} 
	result.push(cf.getRSquared()); //result[3] 
	result.push(linestart); // result[4], start of segment in pixels 
	result.push(maxsize); // result[5], length of segment in pixels 
	return result; 
} 
 
function resize(img, hw, hh, woffset, hoffset) { 
	var resizer = new Packages.ij.plugin.CanvasResizer(); 
	var ims = img.getStack(); 
	IJ.setBackgroundColor(0, 0, 0); // fill expanded image with black pixels 
	ims = resizer.expandStack(ims, hw, hh, woffset, hoffset); 
	img.setStack(ims); 
} 
 
function cm2um(cal) { 
	if (cal.getXUnit() == "cm" && cal.getYUnit() == "cm") { 
		cal.setXUnit("um"); 
		cal.setYUnit("um"); 
		cal.pixelWidth = 1e4*cal.getX(1); 
		cal.pixelHeight = 1e4*cal.getY(1); 
	} 
} 
 
function LoadRois(rs, roifile, prompt) { 
	if (roifile == null) { 
		var rc = new OpenDialog("Choose ROIs to use with image: "+prompt); 
		roifile = new File(rc.getPath()); 
	} 
	if (roifile.isFile() && (getExt(roifile.getName()) == "roi" || 
		getExt(roifile.getName()) == "zip")) { 
		rs.runCommand("Open", roifile); 
	} 
} 
 
function drawRoi(img, roi) { 
	IJ.run(img, "Enhance Contrast", "saturated=0.35"); 
	IJ.run(img, "RGB Color", ""); 
	img.setColor(Packages.java.awt.Color.YELLOW); 
	var ip = img.getProcessor(); 
	ip.draw(roi); 
	img.setProcessor(ip); 
} 

function getTimingFromLogfile(img1, imagefile, timestep) { 
	// SlideBook-specific commands 
	// Attempt to open the SlideBook log file 
	var parent = new File(imagefile.getParent()); 
	var logfilename = imagefile.getName().substring(0,imagefile.getName().lastIndexOf("."))+".log"; 
	logfilename = logfilename.replace("T00_", "T0_"); // workaround for a bug in Slidebook file names 
	var metadatafile = new File(parent.getParent(), "logs"+sep+logfilename); 
 
	if (metadatafile.exists()) { 
		// Read nSlices, nFrames from the log file 
		var logslices = readValue(metadatafile, 2, "Z Planes: "); 
		var logframes = readValue(metadatafile, 3, "Time Points: "); 
		// Assume nChannels = 1 for exported TIFFs 
		//var logchans = readValue(metadatafile, 4, "Channels: "); 
 
		// Swap slices and frames (workaround for a SlideBook limitation) 
		if (logslices*logframes == img1.getNSlices()*img1.getNFrames()) 
			img1.setDimensions(img1.getNChannels(), logslices, logframes); 
 
		// Read the timestep from the log file 
		timestep = readValue(metadatafile, 7, "Average Timelapse Interval: "); 
		// convert ms->s 
		timestep = timestep/1e3; 
	    if (DEBUG) { IJ.showMessage("timestep is now: "+IJ.d2s(timestep, digits)); } 
		// Add time scale back to image 
	    var cal = img1.getCalibration(); 
	    cal.setTimeUnit("s"); 
		cal.frameInterval = timestep; 
	    img1.setCalibration(cal); 
	} else { 
		IJ.log("No logfile available for "+imagefile.getName()+ 
			newline+"Using default timestep "+IJ.d2s(timestep, digits)); 
	} 
	return timestep; 
}

function smooth(intensities, smoothwidth) {
	var n = intensities.length;
	var smoothedpixels = new Array(n);
	if (n > 2) {
		// perform circular smoothing
		for (var j = 0; j<n; j++) {
			var temp = 0;
			for (var k = 1; k<=smoothwidth; k++) {
				var ind = j-k+Math.ceil(k/2);
				if (ind < 0)
					ind += n;
				if (ind >= n)
					ind -= n;
				temp += intensities[ind];
			}
			smoothedpixels[j] = temp/smoothwidth;
		}
		return smoothedpixels;
	} else {
		IJ.showMessage("Array is too small for circular smoothing");
		return null;
	}
}

function Angle(rad, steps, scale) {
	// This object reports angular values in radians or degrees and pathlength values in pixels or um
	//
	// rad		: 	angle in radians
	// steps	:	total perimeter distance in pixels
	// scale	:	width of a pixel in um
	
	this.rad = (rad < 0) ? rad+2*Math.PI : rad; // shift to range [0..2PI]rad;
	this.deg = this.rad * 180/Math.PI;

	this.pixel = this.rad * steps/(2*Math.PI);
	
	if (scale != 1)	
		this.pixel = Math.floor(this.rad); // enforce integer values when using pixel units
			
	this.dist = this.pixel * scale;
	
	// Convert to a string ...
	this.toString = function() {
		return IJ.d2s(this.deg, digits) + " deg.";
	}
}
