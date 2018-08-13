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
