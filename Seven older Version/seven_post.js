importPackage(Packages.ij);
importPackage(Packages.ij.io);
importPackage(Packages.java.io);
load(IJ.getDirectory("ImageJ")+"seven.js");

var reanalyze = true; // If false, summarize data but don't reanalyze
var linescantag = "newzeroscan"; // Prefix for length scans
var anaversion = "6";
var sep = File.separator;

//var dc = new Packages.ij.io.DirectoryChooser("Choose the experiment analysis directory:");
//var root = new Packages.java.io.File(dc.getDirectory());
//var dirlist = root.list();

var od = OpenDialog("Choose the list of input directories:");

var list = od.getPath();
var datadir = new File(list+"_Results-"+anaversion);
datadir.mkdir();
var tipdir = new File(datadir+sep+"TipAnalysis-"+anaversion);
tipdir.mkdir();
var lengthdir = new File(datadir+sep+"LengthAnalysis_"+linescantag+"-"+anaversion);
lengthdir.mkdir();
var spacedir = new File(datadir+sep+"SpacingAnalysis-"+anaversion);
spacedir.mkdir();
var exptlist = readPaths(list);
exptlist.sort();
var mastercellpath = list+"-"+anaversion+"-fp-per-cell.txt";
var masterfilopath = list+"-"+anaversion+"-filos.txt";
var masterkymopath = list+"-"+anaversion+"-kymograph.csv";
var masterimagepath = list+"-"+anaversion+"-images.csv";
saveText(mastercellpath, "", false);
saveText(masterfilopath, "", false);
saveText(masterkymopath, "", false);
saveText(masterimagepath, "", false);
IJ.showMessage("Number of expts = "+IJ.d2s(exptlist.length,0));

for (var u = 0; u < exptlist.length; u++) {
	var root = new File(exptlist[u]);
	var cell_per_fp = new Array(); // global array for registering filopodia to cells
	var cell_body = new Array(); // global array for cell body intensity
	var cell_tipavg = new Array(); // global array for cell average tip intensity
	var cell_band = new Array(); // global array for cell band intensity
	var area_body = new Array(); // global array for area of cell body
	var area_band = new Array(); // global array for area of band 
	var intensity_per_fp = new Array(); // global array for filopod tip intensity	
	var acqname = "Capture "; // Default prefix for input files
	//var acqname = root.getName()+"."; // Default prefix for input files

	// Hook to run all the basic analysis before collating the results
	if (reanalyze) {
		seven_multi(root, acqname);
	}

	// Create summary files for all expts
	var dirlist = root.list();
	dirlist.sort();
	var summarypath = root.getCanonicalPath()+sep+"summary-"+anaversion+".txt";
	var cellpath = root.getCanonicalPath()+sep+"summary-fp-per-cell-"+anaversion+".txt";
	var filopath = root.getCanonicalPath()+sep+"summary-filos-"+anaversion+".txt";
	var kymopath = root.getCanonicalPath()+sep+"summary-kymograph-"+anaversion+".csv";
	var imagepath = root.getCanonicalPath()+sep+"experiment-results-"+anaversion+".csv";
	saveText(summarypath, "", false);
	saveText(cellpath, "", false);
	saveText(filopath, "", false);
	saveText(kymopath, "", false);
	IJ.log(IJ.d2s(dirlist.length,0));

	// Create per-experiment results directories
	var tipexdir = new File(tipdir+sep+root.getName()+"_"+IJ.d2s(u,0));
	tipexdir.mkdir();
	var lengthexdir = new File(lengthdir+sep+root.getName()+"_"+IJ.d2s(u,0));
	lengthexdir.mkdir();
	var spaceexdir = new File(spacedir+sep+root.getName()+"_"+IJ.d2s(u,0));
	spaceexdir.mkdir();
	
	for (var i = 0; i<dirlist.length; i++) {
		var subdir = new File(root.getCanonicalPath(), dirlist[i]);
		var anadir = subdir.getCanonicalPath()+sep;
		
		if (subdir.isDirectory()) {
			var sublist = subdir.list();
			sublist.sort();
			var logfile = null;
			var imagetable = null;
			var celltable = null;
			var spacingtable = null;
			var filotable = null;
			
			for (var j=0; j<sublist.length; j++) {
				var tempfile = new File(subdir.getCanonicalPath(), sublist[j]);
				var temp = sublist[j];
				var lastchar = temp.length();
				var versionchar = lastchar-anaversion.length-4; // length of anaversion+".txt"

				// Copy files from source subdir to results dir
				if (parseInt(temp.length()) > 8) {
					if (temp == "analysis-log-"+anaversion+".txt") {
						logfile = tempfile;
					} else if (temp == "analysis-results-"+anaversion+".csv") {
						celltable = tempfile;
					} else if (temp == "spacing-results-"+anaversion+".csv") {
						spacingtable = tempfile;
					} else if (temp == "filopod-results-"+anaversion+".csv") {
						filotable = tempfile;
					} else if (temp.substring(versionchar, lastchar) == anaversion+".txt") {
						if (temp.substring(0,8) == "linescan") {
							CopyFile(tempfile, tipexdir, subdir.getName()+"-"+temp);
						} else if (temp.substring(0,8) == linescantag.substring(0,8)) {
							CopyFile(tempfile, lengthexdir, subdir.getName()+"-"+temp);
						} else if (temp.substring(0,7) == "spacing") {
							CopyFile(tempfile, spaceexdir, subdir.getName()+"-"+temp);
						}
					}
				}
			}
			var tempfile = new File(subdir.getCanonicalPath()+sep+"Kymograph", "kymograph_results.csv");
			var skipheader = (i == 0) ? null : "1";
			if (tempfile.isFile()) {
				readBetween(kymopath, tempfile, skipheader, null);
			}
	
			if (logfile != null) {
				seven_post(logfile.getCanonicalPath(), anadir);			
			}
		}
	}

	var skipheader = (u == 0) ? null : "1";
	readBetween(masterimagepath, imagepath, null, null);
	readBetween(masterkymopath, kymopath, null, null);
	readBetween(mastercellpath, cellpath, null, null);
	readBetween(masterfilopath, filopath, null, null);
}

function seven_post(datafile, anadir) {
		IJ.log(datafile);
	if (IJ.getLog() != null) {
		IJ.selectWindow("Log");
		IJ.run("Close");
	}

	saveText(summarypath, datafile+newline, true);
	readBetween(summarypath, datafile, "Cell", "Filopod");
	readBetween(cellpath, datafile, "Cells", "Filopod");
	readBetween(filopath, datafile, "Filopod", "End");
	IJ.log(summarypath);
}
