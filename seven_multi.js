importPackage(Packages.ij);
importPackage(Packages.ij.io);
importPackage(Packages.java.io);
load(IJ.getDirectory("ImageJ")+"seven.js");

var dc = new DirectoryChooser("Choose the experiment analysis directory:");
var root = new File(dc.getDirectory());

	var cell_per_fp = new Array(); // global array for registering filopodia to cells
	var cell_body = new Array(); // global array for cell body intensity
	var cell_tipavg = new Array(); // global array for cell average tip intensity
	var cell_band = new Array(); // global array for cell band intensity
	var area_body = new Array(); // global array for area of cell body
	var area_band = new Array(); // global array for area of band 
	var intensity_per_fp = new Array(); // global array for filopod tip intensity
	var fp_len = new Array(); // global array for filopod length in micrometers
//	var acqname = "Capture "; // KJP convention - Default prefix for input files
	var acqname = root.getName()+"."; // ALA convention - Default prefix for input files

//    IJ.showMessage("acqname", acqname); 
seven_multi(root, acqname);