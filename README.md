# Tioga
## Create trap files 
2018: SegmentTransects_2018 and 2018TiogaElkCoordinatesUTMSex.csv  
*Note: 2018TrapFile.R was previously used.*  

2019: SegmentTransects_2019 and 2019TiogaTracksCombined.csv  
*Note: 2019TrapFile.R was previously used.*  
This is the messiest code of all the scripts...  

## Create capture history files 
2018: SegmentTransects_2018 and 2018McKElkLocations.csv  
*Note: 2018TiogaCaptHistFileCode.R was previously used.*  
  
2019: SegmentTransects_2019 and 2019TiogaElkCoordinatesUTMSex.txt  
*Note: 2019TiogaCaptHistFileCode.R was previously used.*  

## Create mask file  
TiogaMaskScript.R  
I have all the required Rasters and shapefiles if needed.  

## Files to Run Models  
  Located in Nelson_secrfiles
  * 2018TiogaTrapFile_300Msegments_secr.txt (created with SegmentTransects_2018.R)  
  * 2019TiogaTrapFile_300Msegments_secr.txt (created with SegmentTransects_2019.R)  
  * combinedcaptfile_300Msegments.txt (created by combining 2018 TiogaCH_300Msegments.txt and 2019TiogaCH_300Msegments.txt created with SegmentTransects_2018 and SegmentTransects_2019, respectively)  
  * 2018.2019.TiogaTelemetry.300_300Msegments.txt (300 random locations from telemetered individuals)  
  * TiogaMask.comb.km_8000buffer500res.rds (created with TiogaMaskScript.R)
  

