   dir = getDirectory("Choose a Directory ");
   orgdir = dir;
   macrodir = getDirectory("imagej");
   runMacro(macrodir+"replace_angles.ijm", dir);
   count = 0;
   countFiles(dir);
   n = 0;
   processFiles(dir);
   print(count+" files processed");
   angleValue = newArray();
   
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
       if (matches(path, ".*sample-.*")) {
       	   scanpath = Array.concat(scanpath, path);
      } 
       for (i=0; i<scanpath.length;i++) {
       	    tipspath = scanpath[i];
       	   // print(tipspath);
       		findtipsAngle(tipspath);
       }
  }
function findtipsAngle(tipspath) {
    filestring=File.openAsString(tipspath);
	dirpath = File.directory; 
	val = newArray();
	rows=split(filestring, "\n");
	 cell_x = newArray();
	 cell_y = newArray();
	 tip_x = newArray();
	 tip_y = newArray();
	 cell_rad = newArray();
//	Array.print(rows);
	for(i=0; i<rows.length; i++){ 
	columns=split(rows[i],"\t");
	cell_x = Array.concat(cell_x, columns[0]);
	cell_y = Array.concat(cell_y, columns[1]);
	tip_x = Array.concat (tip_x, columns[2]);
	tip_y = Array.concat (tip_y, columns[3]);
	cell_rad = Array.concat (cell_rad, columns[4]);
	} 
	
//	print(tip_x[1]);
//	print(tip_y[1]);
//	print(cell_x[1]);
//	print(cell_y[1]);
//	print (tip_x[2]);
//	print (tip_y[2]);
//	xCoord = newArray(parseInt(tip_x[1]),parseInt(cell_x[1]),parseInt(tip_x[2]));
//	yCoord = newArray(parseInt(tip_y[1]),parseInt(cell_y[1]),parseInt(tip_y[2]));
//	calAngle(xCoord, yCoord);

	for (j=0; j<cell_x.length;j++) {
     xCoord = newArray;
     yCoord = newArray;
	 xCoord = newArray(parseInt(tip_x[0]),parseInt(cell_x[0]),parseInt(tip_x[j]));
	 yCoord = newArray(parseInt(tip_y[0]),parseInt(cell_y[0]),parseInt(tip_y[j]));
	 val = Array.concat (val, d2s(calAngle(xCoord, yCoord),2));
	 length_res = Array.concat(length_res, d2s(calArc(xCoord, yCoord, cell_rad[0]),2));
	 //val = calAngle(xCoord, yCoord, cell_rad[0]);
	}
   // Array.print(val);
   // Array.print(length_res);
   //  File.append(val, tipspath);
	Anglelog = dir + "overall_tips_Angle_results.txt";

	for(i=0; i<rows.length; i++){ 
	//	File.append(cell_x[i] + "\t" + cell_y[i] +"\t" + tip_x[i] +"\t" + tip_y[i] +"\t" + val[i], Anglelog);
		 File.append(cell_x[i] + "\t" + cell_y[i] +"\t" + tip_x[i] +"\t" + tip_y[i] +"\t" + val[i] + "\t" + length_res[i+1], Anglelog);
		} 
	
} 

function calAngle (xCoord, yCoord) {
	//result_array = newArray();
	x1=xCoord[0]; y1=yCoord[0]; 
	x2=xCoord[1]; y2=yCoord[1]; 
	x3=xCoord[2]; y3=yCoord[2]; 

	vx1 = x1-x2; 
	vy1 = y1-y2; 
	vx2 = x3-x2; 
	vy2 = y3-y2; 

	scalarProduct=(vx1*vx2 + vy1*vy2); 
	lengthProduct =sqrt((pow(vx1, 2)+pow(vy1, 2))) * sqrt((pow(vx2, 2)+pow(vy2, 2))); 
	costheta = scalarProduct/lengthProduct ; 

	thetadegrees = acos(costheta)*180/PI; 
	return thetadegrees;
	}

	function calArc(xCoord, yCoord, cell_rad){
    thetadegrees = calAngle(xCoord, yCoord);
    rad = cell_rad;
    	length_arc = 2*PI*rad*(thetadegrees/360);
	    return length_arc;
	    }