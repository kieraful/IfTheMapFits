# ENGO500_Capstone

# Executive Summary

In indoor mapping applications, the determination of position utilizes three different sensor systems; GNSS, INS and LiDAR. 
The problem with having multiple sensors in the system is that they are all collecting different kinds of data with respect to 
their own coordinate frames and reference center. When determining the location of the user, we are getting location data from 
three different sensors in three different reference frames. To solve this problem, we must perform a boresight calibration 
to scale, translate and orient all three sensor systems to be in the same reference frame. After collecting simultaneous 
GNSS, INS and LiDAR data, and obtaining calibration parameters to orient all three systems to the same frame, we can test the 
stability of the parameters. The computed calibration parameters will be tested in different environments, as well as with 
varying temporal conditions to obtain quality measures of our result. 
