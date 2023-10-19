
# This is earlyloc's parameter file

# Basic Earthworm setup:
#
MyModuleId         MOD_EARLYLOC # module id for this instance of template
InRingName         PICK_RING    # shared memory ring for input picking information
OutRingName        HYPO_RING    # shared memory ring for output report string

LogFile            1            # 0 to turn off disk log file; 1 to turn it on
                                # to log to module log but not stderr/stdout
HeartBeatInterval  15           # seconds between heartbeats

# Output related settings:
#
#
ReportPath          /home/EARLYLOC_REPORT
OutputPostfix       postfix
ReportTermNum       50         # The last report should be less than this number.
OutputPickRing      1          # 0 to turn off output pick result to OutRing; 1 to turn it on.
OutputRejectPick    1          # 0 to turn off output rejected pick to report file and OutRing; 1 to turn it on.

# Pickings related settings:
#
#
IgnoreWeightP       3           # P-phase picks with weight larger than & equal to this number will be ignored.
IgnoreWeightS       0           # S-phase picks with weight larger than & equal to this number will be ignored.
PickAliveTime       60.0        # Survival time of each picking, it is second between the picking receiving time and current time
PickFetchStrategy   step        # Fetching strategy of picks, 'greedy' (default) which will fetch all picks until there is not any more;
                                # 'step' mode which will process hypo after every fetching of valid pick.
#
# Hypo & Clustering related settings:
#
TriggerPicks        6
HypoAliveTime       60.0        # Survival time of each hypo pool, it is second between the last hypoing time and current time
#
ClusterTimeDiff     5.0        # The same phase arrival time between each clustered station
ClusterDist         40.0       # Distances between each clustered station
ClusterPicks        2

# Wave velocity layer model:
#
#
# P wave velocity model:
#
PWaveBoundary        40.0                # boundary of shallow and deep layers
ShallowPWaveVel      5.10298             # initial velocity in shallow layer
ShallowPWaveGrad     0.06659             # gradient velocity in shallow layer
DeepPWaveVel         7.80479             # initial velocity in deep layer
DeepPWaveGrad        0.00457             # gradient velocity in deep layer
# S wave velocity model:
#
SWaveBoundary        50.0                # boundary of shallow and deep layers
ShallowSWaveVel      2.9105              # initial velocity in shallow layer
ShallowSWaveGrad     0.0365              # gradient velocity in shallow layer
DeepSWaveVel         4.5374              # initial velocity in deep layer
DeepSWaveGrad        0.0023              # gradient velocity in deep layer

# 3D Wave velocity model (Optional):
#
#
3DVelocityModelFile     /home/3D_VELOCITY_MODEL

# MySQL server information:
#
# If you setup the follow parameter especially SQLHost, this program will fetch
# list from MySQL server or you can just comment all of them, then it will turn
# off this function.
#
#SQLHost         127.0.0.1         # The maximum length is 36 words
#SQLPort         3306              # Port number between 1 to 65536
#SQLDatabase     EEW               # The maximum length is 36 words

# Login information example
#
# SQLUser       test
# SQLPassword   123456
#@LoginInfo_sql                    # Please keep the security of the SQL login information

# List the stations lists that will grab from MySQL server
#
# Even when you using MySQL server to fetch station information.
#
#SQLStationTable    PalertList
#SQLStationTable    SecondList
#SQLStationTable    ThirdList

# List the message logos to grab from transport ring
#              Installation       Module        Message
GetEventsFrom  INST_WILDCARD    MOD_WILDCARD    TYPE_EARLY_PICK

# Local SNL list:
#
# The local list for stations that will receive. By the way, the priority of local list
# is higher than the one from remote data. And the layout should be like these example:
#
#         Station   Network   Location    Latitude      Longitude       Elevation(m)
# UseSNL   TEST       TW         --       23.050514     121.215483      1250.0             # example
