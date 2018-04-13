   dir = getDirectory("Choose a Directory ");
   orgdir = dir;
   maskedval = "65535";
   count = 0;
   countFiles(dir);
   n = 0;
   macrodir = getDirectory("imagej");
   runMacro(macrodir+"replace_files.ijm", dir);
   processFiles(dir);
   print(count+" files processed");
   
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              countFiles(""+dir+list[i]);
          else
              count++;
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/")) {
              processFiles(""+dir+list[i]);
              // print(list[i]);
          }
          else {
             //showProgress(n++, count);
             path = dir+list[i];
             processFile(path);
            // print(list[i]);
          }
      }
  }

  function processFile(path) {
  	     scanpath = newArray();
       if (matches(path, ".*new.*")) {
       	   scanpath = Array.concat(scanpath, path);
      } 
       for (i=0; i<scanpath.length;i++) {
       	    tipspath = scanpath[i];
       	   // print(tipspath);
       		findtipslength(tipspath);
       }
  }
function findtipslength(tipspath) {
    filestring=File.openAsString(tipspath);
	dirpath = File.directory; 
	rows=split(filestring, "\n");
	 x = newArray();
	 y = newArray();
	//print(rows.length); 
	for(i=0; i<rows.length; i++){ 
	columns=split(rows[i],"\t");
	x = Array.concat(x, columns[0]);
	y = Array.concat(y, columns[1]);
	} 
	//Array.print(x);
	//Array.print(y);
	index = 0;
	for (j=0; j<y.length;j++) {	
	  if (y[j] == maskedval )	{
	  	break;
	  }
	  index++;
	}
	savelog = dir +"tips_results.txt";
	savelog2 = orgdir + "overall_tips_results.txt";
	File.append(x[index-1], savelog);
	File.append(x[index-1], savelog2);
} 