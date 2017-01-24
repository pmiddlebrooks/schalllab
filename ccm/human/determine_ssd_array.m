function SSDArrayScreenFlips = determine_ssd_array(meanSSD)

meanSSD = round(meanSSD);

SSDArrayScreenFlips = meanSSD - 6 : 4 : meanSSD + 14;