#################################################
#                                               #
# Control file for COMCOT program (v1.7)        #
#                                               #
#################################################
#--+-----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
#===============================================:===============================
# General Parameters for Simulation             : Value Field                  |
#===============================================:===============================
#Job Description: NZ30sec bathymetry, Spherical Coordinates for code testing
 Total run time (Wall clock, seconds)           :  14400.0
 Time interval to Save Data    ( unit: sec )    :  600.0
 Output Zmax & TS (0-Max Z;1-Timeseries;2-Both) :     0
 Start Type (0-Cold start; 1-Hot start)         :     0
 Resuming Time If hot start (Seconds)           :   2400.00
 Specify Min WaterDepth offshore  (meter)       :    0.00
 Initial Cond. (0:FLT,1:File,2:WM,3:LS,4:FLT+LS):     0
 Specify BC  (0-Open;1-Sponge;2-Wall;3-FACTS)   :     0
 Specify Input Z filename (for BC=3, FACTS)     : mw94_n22_nz_ha.xyt
 Specify Input U filename (for BC=3, FACTS)     : mw94_n22_nz_ua.xyt
 Specify Input V filename (for BC=3, FACTS)     : mw94_n22_nz_va.xyt

#===============================================:===============================
# Parameters for Fault Model (Segment 01)       :Values                        |
#===============================================:===============================
 No. of FLT Planes (With fault_multi.ctl if >1) :   236
 Fault Rupture Time (seconds)                   :   16.0
 Faulting Option (0: Model; 1- Data;)           :   0
 Focal Depth                             (meter):   16990.000
 Length of source area                   (meter):   25000.000
 Width of source area                    (meter):   17750.000
 Dislocation of fault plate              (meter):    6.390
 Strike direction (theta)               (degree):   18.700
 Dip  angle       (delta)               (degree):    17.5000
 Slip angle       (lamda)               (degree):    122.99480
 Origin of Comp. Domain (Layer 01) (Lat, degree):   -43.240000
 Origin of Comp. Domain (Layer 01) (Lon, degree):   280.000000
 Epicenter: Latitude                    (degree):   -36.3477
 Epicenter: Longitude                   (degree):   286.4160
 File Name of Deformation Data                  : segment_parameter.dat
 Data Format Option (0-COMCOT; 1-MOST; 2-XYZ)   :     2

#===============================================:===============================
#  Parameters for Wave Maker                    :Values                        |
#===============================================:===============================
 Wave type  ( 1:Solit, 2:given, 3:focusing )    :     1
 FileName of Customized Input (for Type=2)      : fse.dat
 Incident direction( 1:top,2:bt,3:lf,4:rt,5:ob ):     2
 Characteristic Wave Amplitude        (meter)   :     0.500
 Typical Water depth                  (meter)   :  2000.000 
 
#===============================================:===============================
#  Parameters for Submarine LS/Transient Motion :ValUes                        |
#===============================================:===============================
 X Coord. of Left/West Edge of Landlide Area    :  177.00
 X Coord. of Right/East Edge of Landlide Area   :  179.00
 Y Coord. of Bottom/South Edge of Landlide Area :  -41.00
 Y Coord. of Top/North Edge of Landlide Area    :  -39.00
 File Name of landslide Data                    : landslide_test.dat
 Data Format Option (0-Old; 1-XYT; 2-Function)  :     2
 
#===============================================:===============================
# Configurations for all grids                  :Values                        |
#===============================================:===============================
# Parameters for 1st-level grid -- layer 01     :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     0
 Coordinate System    (0:spherical, 1:cartesian):     0
 Governing Equations  (0:linear,    1:nonlinear):     0
 Grid Size  (dx, sph:minute, cart:meter)        :     2.16
 Time step                            ( second ):     1.0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 X_start                                        :  280.000000
 X_end                                          :  289.000000
 Y_Start                                        :  -43.240000
 Y_end                                          :  -31.000000
 File Name of Bathymetry Data                   : global_salida.xyz
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    01
 Grid Level                                     :     1
 Parent Grid's ID Number                        :    -1

#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 02    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     0
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        : 280.911600
 X_end                                          : 281.442600
 Y_start                                        : -33.895800
 Y_end                                          : -33.370200
 FileName of Water depth data                   : jfer0018_conisla.txt
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    02
 Grid Level                                     :     2
 Parent Grid's ID Number                        :     1

#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 03    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :  280.911600
 X_end                                          :  281.442600
 Y_start                                        :  -33.895800
 Y_end                                          :  -33.370200
 FileName of Water depth data                   : jfer0018_conisla.txt
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    03
 Grid Level                                     :     3
 Parent Grid's ID Number                        :    02 
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 04    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     0
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     0
 GridSize Ratio of Parent layer to current layer:     6
 X_start                                        :    281.162530 
 X_end                                          :    281.190980
 Y_start                                        :    -33.644460
 Y_end                                          :    -33.621600
 FileName of Water depth data                   : juanfernandez0003.txt
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     2
 Grid Identification Number                     :    04
 Grid Level                                     :     4
 Parent Grid's ID number                        :    03 
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 05    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    61.
 X_end                                          :    80.
 Y_start                                        :    61.
 Y_end                                          :    80.
 FileName of Water depth data                   : layer24.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    05
 Grid Level                                     :     3
 Parent Grid's ID number                        :    02 
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 06    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   115.
 X_end                                          :   233.
 Y_start                                        :   407.
 Y_end                                          :   573.
 FileName of Water depth data                   : layer31.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    06
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 07    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   140.
 X_end                                          :   233.
 Y_start                                        :   143.
 Y_end                                          :   310.
 FileName of Water depth data                   : layer32.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    07
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01 

#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 08    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Use Bottom friction ?(only cart,nonlin,0:y,1:n):     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   274.
 X_end                                          :   329.
 Y_start                                        :   143.
 Y_end                                          :   235.
 FileName of Water depth data                   : layer33.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    08
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 09    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    41.
 X_end                                          :    60.
 Y_start                                        :    41.
 Y_end                                          :    60.
 FileName of Water depth data                   : layer34.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    09
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 10    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   129.
 X_end                                          :   247.
 Y_start                                        :   471.
 Y_end                                          :   588.
 FileName of Water depth data                   : layer41.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    10
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 11    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :   858.
 X_end                                          :   968.
 Y_start                                        :   246.
 Y_end                                          :   388.
 FileName of Water depth data                   : layer42.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    11
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 12    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     1
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :  2154.
 X_end                                          :  2258.
 Y_start                                        :   671.
 Y_end                                          :   825.
 FileName of Water depth data                   : layer43.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    12
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01
 
#===============================================:=============================== 
#  Parameters for Sub-level grid -- layer 13    :Values                        |
#===============================================:===============================
 Run This Layer ?       (0:Yes,       1:No     ):     1
 Coordinate           (0:spherical, 1:cartesian):     1
 Governing Eqn.       (0:linear,    1:nonlinear):     0
 Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1
 Manning's Roughness Coef. (For fric.option=0)  :     0.013
 Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1
 GridSize Ratio of Parent layer to current layer:     5
 X_start                                        :    41.
 X_end                                          :    60.
 Y_start                                        :    41.
 Y_end                                          :    60.
 FileName of Water depth data                   : layer44.dep
 Data Format Option (0-OLD;1-MOST;2-XYZ;3-ETOPO):     0 
 Grid Identification Number                     :    13
 Grid Level                                     :     2
 Parent Grid's ID number                        :    01