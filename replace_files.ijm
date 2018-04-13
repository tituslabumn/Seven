 // dir = getDirectory("imagej");
   dir = getArgument();
   count = 0;
   countFiles(dir);
   n = 0;
   processFiles(dir);
   //print(count+" files processed");
   
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
          if (endsWith(list[i], "/"))
              processFiles(""+dir+list[i]);
          else {
            // showProgress(n++, count);
             path = dir+list[i];
             deleteFiles(path);
          }
      }
  }

  function deleteFiles(path) {
  	  delpath = newArray();
      if (matches(path, ".*tips_results.*")) {
       	   delpath = Array.concat(delpath, path);
      } 
       for (i=0; i<delpath.length;i++) {
       	    tipspath = delpath[i];
       	  //print(tipspath);
       	    deletepath(tipspath);
       }
  }

  function deletepath(tipspath) {
  	  if (File.exists(tipspath)) {
           File.delete(tipspath); 	  	
  	}
  }
