# Seven

## About 
Seven is image analysis plugin developed by Karl Petersen and its published in this research paper [here](http://www.pnas.org/content/pnas/113/50/E8059.full.pdf). This image analysis script utilizes imageJ application called [FIJI](https://fiji.sc/). The script allows you to exports the stats of the image (Amoeba cell image) such as number of cells, cell area, cell intensity, filopodia per cell, filopodia length, etc.

## How to use Seven
1. Seven utilizes Fiji/ImageJ application. If you already don't have the latest version, Please download it from [here](https://fiji.sc/).

2. Install the updated version of seven from Github repository and move all the seven files to the Fiji.app root folder.

3. Move the tiff files from an experiment to a folder with the name of the file string beginning i.e. 95.1 - 95.25 into a folder 95 or so.

4.	From Fiji/Image J Select :

        a.  Plugins>Macros>Edit>seven_multi
        b.  Run seven_multi file.
        c.  Select the folder with tiff files which you want to analyze. 

5.	Seven Runs through the folder and creates a ZIP mask for each image.

6. The seven will ask you to Manually edit the image file incase the tiff image is not clear.

7. Manually erase any erroneous cells by comparing with the original tiff file, if image is not clearly inverted.

8. For each image in the experiment folder, seven will create a sliced folder ( for ex. sliced_1-1 or so). These folder will contained all the stats exported from an image which will be helpful in conducting analysis.



&copy; Regents of the University of Minnesota

