importPackage(Packages.ij);
importPackage(Packages.ij.io);
importPackage(Packages.java.io);
load(IJ.getDirectory("ImageJ")+"seven.js");

var dc = new DirectoryChooser("Choose the experiment analysis directory:");
var dirpath = dc.getDirectory();
var root = (dirpath != null) ? new File(dirpath) : null;

	// fp data - probably not necessary as global (only referenced in AnalyzeTips())
	var cell_per_fp = new Array(); // global array for registering filopodia to cells
	var intensity_per_fp = new Array(); // global array for filopod tip intensity
	var xpos_per_fp = new Array(); // global array for filopod raw X coordinate (in pixels)
	var ypos_per_fp = new Array(); // global array for filopod raw Y coordinate (in pixels)

	// cell data - assigned in AnalyzeCells() and referenced again in AnalyzeTips()
	var cell_body = new Array(); // global array for cell body intensity
	var cell_xpos = new Array(); // global array for cell raw X coordinate (in pixels)
	var cell_ypos = new Array(); // global array for cell raw Y coordinate (in pixels)
	var cell_tipavg = new Array(); // global array for cell average tip intensity
	var cell_band = new Array(); // global array for cell band intensity
	var cell_area_body = new Array(); // global array for area of cell body
	var cell_area_band = new Array(); // global array for area of band 
//	var acqname = "Capture "; // KJP convention - Default prefix for input files
	var acqname = root.getName()+" - "; // ALA convention - Default prefix for input files

//    IJ.showMessage("acqname", acqname); 

if (root != null)
	seven_multi(root, acqname);
else
	IJ.error("Seven.js", "Invalid directory path");