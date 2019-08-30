# mms_shockstat

## Brief description
Some tools to do plot overviews and do statistics on bow shock crossings by MMS.

## Compatibility
Requires Matlab 2016b or newer. Requires the [irfu-matlab](https://github.com/irfu/irfu-matlab) package. Is only compatible with Mac and Linux. 

## Running in Matlab
Run the script
```matlab
ShockStatMMS
```

All other routines are called from this script by prompts.

## Lists
There are several lists compiled from [MMS SITL reports](http://www.ssl.berkeley.edu/~moka/eva/sitl_report.html) in this repo. 

shock_list.txt - Main file containing many shock crossings  
notshock_list.txt - Events that are not shock crossings (solar wind CSs, HFAs, IP shocks...)  
shock2_list.txt - Hand-picked events with stable OMNI data  
shock_list_accepted.txt - Events where the OMNI data is relatively stable  


## Detailed description

### Operations 
#### Operations on SITL txt lists
- Calculate shock parameters: select a txt list and calculate various shock parameters for the events. Relies on reading a large amount of MMS data and is therefore slow. The function relies on a SITL text list and a time interval list in csv format. The function  typically writes a mat file containing saved data.
- Set time intervals for up- and downstream: Visually set up- and downstream intervals for every event. Writes a csv file with time intervals.
- Get line numbers of events with stable OMNI data: Automatic selection of events with stable solar wind conditions. Writes new SITL txt list. 
        
#### Operations on MAT lists
- Load shock parameters: loads mat file with certain constraints.
- Plot various shock parameters: plots a large amount of figures
doing various statistics.
- Construct bow shock model: Fits a conic section to all spacecraft positions.
 
#### Plots from SITL txt lists
All functions save figures as pngs:

- MMS vs OMNI plots: comparison between MMS and OMNI data.
- Ion burst overview plots
- Electron burst overview plots
- Ion composition overview plots: uses HPCA data.
  
#### Other
- Write shock parameters to CSV file: takes a mat file and writes certain parameters to a csv file
 
 
### Files
The programs mainly deals with three different types of data files:
- SITL files in txt format: Contains lines copied from SITL reports, preferrably where segments are merged into one line. Each line constitutes one shock crossing.
- Time interval files in cvs format: Contains time intervals specifying one upstrean and one downstream time interval for each shock crossing. This file can only be constructed by visually selecting up- and downstream (option in ShockStatMMS).
- Files containing calculated shock parameters in mat format: Calculated using the other two types of files. Data analysis is best done on these files.

