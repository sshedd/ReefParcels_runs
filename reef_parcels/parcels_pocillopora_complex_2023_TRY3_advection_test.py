#playing around with Ocean Parcels
#before running, make sure parcels is activated using 'conda activate py3_parcels'

import os
import numpy as np
from datetime import datetime, timedelta
import random
import glob
import datetime
import netCDF4
from parcels import (FieldSet, Field, NestedField, ParticleSet, JITParticle, StatusCode,
    AdvectionRK4_3D, StatusCode, Variable, ParcelsRandom)
from lanaau_build_spawn_2023 import getSpawningSites, getSpawningEvents
from lanaau_functions import (HandleBoundaries, RWDiffusion_Lanaau, deleteParticle, submergeParticle,
    SettleFish2, InitDepth)
from lanaau_functions_complex import Behavior, Behavior_vertonly, GrowFishComplex, Buoyancy, KillFish, AdvectionRK4_3D_Lanaau
from LarvalParticleComplex import define_Larval_Class

def KeepInOcean(particle, fieldset, time):
    # If a particle would go through the surface, neutralize the vertical step
    if particle.state == StatusCode.ErrorThroughSurface:
        # In v3 kernels, adjust *increments* (particle_ddepth), not the absolute depth
        particle_ddepth = 0.0
        particle.state = StatusCode.Success

def CheckOutOfBounds(particle, fieldset, time):
    # If interpolation failed because the particle left the domain, delete it
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()

def setFields_nonest(dt):

    dim = {'lon':'Longitude_u', 'lat':'Latitude_t', 'depth':'Depth_w'} #u, t works for 1km; t,t works for 4km
    ts = np.expand_dims(np.arange(0, 86400, 86400), axis = 1) #set up timestamps (corresponds to day of model)
        
    #combine everything into nested fieldset
    fieldset_nested = FieldSet.from_mitgcm(["/drives/hdd10t/users/evan/emily_parcels/current_products/0000000720_cropped.cdf"],
    variables = {'U':'u', 'V':'v', 'W':'w'}, dimensions=dim, deferred_load=True)

    #add fields for coastlines
    fieldset_nested.add_constant("dt", dt)
    northCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/northcoasts_1km.nc", 
            variable=("F_N", "northcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/southcoasts_1km.nc", 
            variable=("F_S", "southcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/eastcoasts_1km.nc", 
            variable=("F_E", "eastcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/westcoasts_1km.nc", 
            variable=("F_W", "westcoast_1km"), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    
    for coast in [northCoasts_1km, southCoasts_1km, eastCoasts_1km, westCoasts_1km]:
        fieldset_nested.add_field(coast)
        
    #set up field to check interpolation
    F = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2021_final_1km/1km_mask.nc", 
            variable=('F','F1'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    fieldset_nested.add_field(F)
    
    #set up bathymetry field
    bathy = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2021_final_1km/1km_bathy.nc", 
            variable=('seafloor','depth1km'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    return(fieldset_nested)

def setFields_from_fieldset(dt):
    #trying nested interpolation fields from fieldsets
    #extra step but reading in just fieldsets (setFields) doesn't seem to work
    dim_u = {'lon':'Longitude_u', 'lat':'Latitude_t', 'depth':'Depth_t'}
    dim_v = {'lon':'Longitude_t', 'lat':'Latitude_v', 'depth':'Depth_t'}
    dim_w = {'lon':'Longitude_t', 'lat':'Latitude_t', 'depth':'Depth_w'}
    ts = np.expand_dims(np.arange(0, 86400, 86400), axis = 1) #set up timestamps (corresponds to day of model in seconds)

    fieldset_1km = FieldSet.from_mitgcm(["/drives/hdd10t/users/evan/emily_parcels/current_products/0000000720_cropped.cdf"],
    variables = {'U1':'u', 'V1':'v', 'W1':'w'}, dimensions = {'U1':dim_u, 'V1':dim_v, 'W1':dim_w}, 
    timestamps=ts, deferred_load=True)
    
    fieldset_4km = FieldSet.from_mitgcm(["/drives/hdd10t/users/evan/emily_parcels/current_products/0000000720_cropped.cdf"],
    variables = {'U2':'u', 'V2':'v', 'W2':'w'}, dimensions = {'U2':dim_u, 'V2':dim_v, 'W2':dim_w}, 
    timestamps=ts, deferred_load=True)
    
    U_nested = NestedField('U', [fieldset_1km.U1, fieldset_4km.U2])
    V_nested = NestedField('V', [fieldset_1km.V1, fieldset_4km.V2])
    W_nested = NestedField('W', [fieldset_1km.W1, fieldset_4km.W2])
    
    #combine everything into nested fieldset
    fieldset_nested = fieldset_4km
    fieldset_nested.add_field(U_nested)
    fieldset_nested.add_field(V_nested)
    fieldset_nested.add_field(W_nested)
    
    #add separate fine fields for silliness near coasts
    fieldset_nested.add_field(fieldset_1km.U1)
    fieldset_nested.add_field(fieldset_1km.V1)
    fieldset_nested.add_field(fieldset_1km.W1)
    
    #set up field to check interpolation
    F1 = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2021_final_1km/1km_mask.nc", 
            variable=('F1','F1'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F2 = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2021_final_4km/4km_mask.nc", 
            variable=('F2','F2'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    F = NestedField('F', [F1, F2])
    fieldset_nested.add_field(F)
    
    #set up bathymetry field
    B1 = Field.from_netcdf("/home/ssheddmu/working/emily_parcels/spawn_files/2021_final_1km/1km_bathy.nc", 
            variable='depth1km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    B2 = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2021_final_4km/4km_bathy.nc", 
            variable='depth4km', dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    bathy = NestedField('seafloor', [B1, B2])
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    #add fields for coastlines
    fieldset_nested.add_constant("dt", dt)
    northCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/northcoasts_1km.nc", 
            variable="northcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/southcoasts_1km.nc", 
            variable="southcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/eastcoasts_1km.nc", 
            variable="eastcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_1km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/westcoasts_1km.nc", 
            variable="westcoast_1km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    
    northCoasts_4km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/northcoasts_4km.nc", 
            variable="northcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    southCoasts_4km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/southcoasts_4km.nc", 
            variable="southcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    eastCoasts_4km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/eastcoasts_4km.nc", 
            variable="eastcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    westCoasts_4km = Field.from_netcdf("/home/sshedd/working/emily_parcels/coast_files/westcoasts_4km.nc", 
            variable="westcoast_4km", dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)

    fieldset_nested.add_field(NestedField('F_N', [northCoasts_1km, northCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_S', [southCoasts_1km, southCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_E', [eastCoasts_1km, eastCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_W', [westCoasts_1km, westCoasts_4km]))

    return(fieldset_nested)
    
#def get_unique_time_files(folder, pattern="*.cdf", exclude_substring=None):
#    import os
#    import glob
#    import xarray as xr

#    files = sorted(glob.glob(os.path.join(folder, pattern)))
#    selected = []
#    seen_times = set()

#    for fn in files:
#        base = os.path.basename(fn)
#        if exclude_substring is not None and exclude_substring in base:
#            continue

#        ds = xr.open_dataset(fn)
#        t = ds["Time"].values[0]

#        if t not in seen_times:
#            selected.append(fn)
#            seen_times.add(t)
            
          # sort by internal time, not filename
#    sorted_times = sorted(time_to_file.keys())
#    sorted_files = [time_to_file[t] for t in sorted_times]

#    return selected
    
def setFields(dt):
    import os
    import glob

    # current files on Scylla
    fine_dir = "/home/sshedd/working/emily_parcels/current_files/1km_MITgcm"
    coarse_dir = "/home/sshedd/working/emily_parcels/current_files/4km_MITgcm"

    fine_files = sorted(glob.glob(os.path.join(fine_dir, "*.cdf")))
    coarse_files = sorted(glob.glob(os.path.join(coarse_dir, "*.cdf")))
    #fine_files = get_unique_time_files(fine_dir, "*.cdf", exclude_substring="_zeros")
    #coarse_files = get_unique_time_files(coarse_dir, "*.cdf")
   
    print(f"Found {len(fine_files)} 1 km current files")
    print(f"Found {len(coarse_files)} 4 km current files")
    print("Example 1 km file:", fine_files[:1])
    print("Example 4 km file:", coarse_files[:1])

    if not fine_files:
        raise FileNotFoundError(f"No .cdf files found in {fine_dir}")
    if not coarse_files:
        raise FileNotFoundError(f"No .cdf files found in {coarse_dir}")

    dim_u = {'lon': 'Longitude_u', 'lat': 'Latitude_t', 'depth': 'Depth_w', 'time':'Time'}
    dim_v = {'lon': 'Longitude_t', 'lat': 'Latitude_v', 'depth': 'Depth_w', 'time':'Time'}
    dim_w = {'lon': 'Longitude_t', 'lat': 'Latitude_t', 'depth': 'Depth_w', 'time':'Time'}

    U_fine = Field.from_netcdf(
        fine_files,
        variable=('U1', 'u'),
        dimensions=dim_u,
        deferred_load=True,
        gridindexingtype="mitgcm",
        fieldtype="U"
    )
    V_fine = Field.from_netcdf(
        fine_files,
        variable=('V1', 'v'),
        dimensions=dim_v,
        deferred_load=True,
        gridindexingtype="mitgcm",
        fieldtype="V"
    )
    W_fine = Field.from_netcdf(
        fine_files,
        variable=('W1', 'w'),
        dimensions=dim_w,
        deferred_load=True,
        gridindexingtype="mitgcm"
    )

    U_coarse = Field.from_netcdf(
        coarse_files,
        variable=('U2', 'u'),
        dimensions=dim_u,
        deferred_load=True,
        gridindexingtype="mitgcm",
        fieldtype="U"
    )
    V_coarse = Field.from_netcdf(
        coarse_files,
        variable=('V2', 'v'),
        dimensions=dim_v,
        deferred_load=True,
        gridindexingtype="mitgcm",
        fieldtype="V"
    )
    W_coarse = Field.from_netcdf(
        coarse_files,
        variable=('W2', 'w'),
        dimensions=dim_w,
        deferred_load=True,
        gridindexingtype="mitgcm"
    )

    U_nested = NestedField('U', [U_fine, U_coarse])
    V_nested = NestedField('V', [V_fine, V_coarse])
    W_nested = NestedField('W', [W_fine, W_coarse])

    fieldset_nested = FieldSet(U_nested, V_nested)
    fieldset_nested.add_field(W_nested)

    fieldset_nested.add_field(U_fine)
    fieldset_nested.add_field(V_fine)
    fieldset_nested.add_field(W_fine)

    F1 = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/spawn_files/2021_final_1km/1km_mask.nc",
        variable=('F1', 'F1'),
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    F2 = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/spawn_files/2021_final_4km/4km_mask.nc",
        variable=('F2', 'F2'),
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    F = NestedField('F', [F1, F2])
    fieldset_nested.add_field(F)

    B1 = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/spawn_files/2021_final_1km/1km_bathy.nc",
        variable='depth1km',
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    B2 = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/spawn_files/2021_final_4km/4km_bathy.nc",
        variable='depth4km',
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    bathy = NestedField('seafloor', [B1, B2])
    fieldset_nested.add_field(bathy)
    fieldset_nested.seafloor.interp_method = "nearest"

    fieldset_nested.add_constant("dt", dt)

    northCoasts_1km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/northcoasts_1km.nc",
        variable="northcoast_1km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    southCoasts_1km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/southcoasts_1km.nc",
        variable="southcoast_1km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    eastCoasts_1km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/eastcoasts_1km.nc",
        variable="eastcoast_1km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    westCoasts_1km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/westcoasts_1km.nc",
        variable="westcoast_1km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )

    northCoasts_4km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/northcoasts_4km.nc",
        variable="northcoast_4km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    southCoasts_4km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/southcoasts_4km.nc",
        variable="southcoast_4km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    eastCoasts_4km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/eastcoasts_4km.nc",
        variable="eastcoast_4km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )
    westCoasts_4km = Field.from_netcdf(
        "/home/sshedd/working/emily_parcels/coast_files/westcoasts_4km.nc",
        variable="westcoast_4km",
        dimensions={'lon': 'lon', 'lat': 'lat'},
        allow_time_extrapolation=True
    )

    fieldset_nested.add_field(NestedField('F_N', [northCoasts_1km, northCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_S', [southCoasts_1km, southCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_E', [eastCoasts_1km, eastCoasts_4km]))
    fieldset_nested.add_field(NestedField('F_W', [westCoasts_1km, westCoasts_4km]))

    return fieldset_nested
  
def main():
    
    #set random seed
    #666, 808
    random.seed(808)
    np.random.seed(808)
    
    #set up timestep - at 12 minutes we don't skip over any cells at max velocity
    dt = 60 #in minutes
    simdays = 1

    ################ VARIABLES FOR CREATING PARTICLE SETS AND SPAWNING EVENTS
    
    n_eggs = 1 #eggs per timestep per site
    typeSpawner = 1 #broadcast:0, brood:1
    monthSpawner = True #year-round (False) or seasonal reproduction (True)
    spawnMonths = [1]
    timedSpawner = True #spawning during a certain time of day (True) or not (False)
    timeSpawnBegin = datetime.timedelta(hours = 23) #in hours
    timeSpawnEnd = datetime.timedelta(hours = 24) #in hours
    peakSpawner = False #peak spawning during a certain time of year
    peakSpawnMonths = "NA"
    moonSpawner = True #spawning only at certain moon phases (True) or not (False)
    spawnPhase = ["New moon"]

    fieldset = setFields(dt) #get fieldset                                                                                      # Add constants so kernels can access them
    fieldset.add_constant("dt", dt)  # time step, if not already added
    fieldset.add_constant("typeSpawner", typeSpawner)  # 0=broadcast, 1=brood
    
    ################ FIELDSET CONSTANTS
    
    #spawning
    maxAdultDepth = 27
    fieldset.add_constant("maxAdultDepth", maxAdultDepth) #max adult depth in meters
    
    #egg
    hatchTime = 0.83
    
    firstBuoy = -0.001
    fieldset.add_constant("hatchTime", hatchTime) #hatch time in days

    #pelagic larval duration in days
    minPLD = 1
    maxPLD = 100
    fieldset.add_constant("minPLD", minPLD + hatchTime)
    fieldset.add_constant("maxPLD", maxPLD + hatchTime)
    
    #for linear growth
    fieldset.add_constant("growth_m", 0) #slope for growth function
    fieldset.add_constant("growth_c", 0) #y-intercept for growth function
    
    #for random-walk swimming
    fieldset.add_constant("k", 3) #kappa
    fieldset.add_constant("B", 4) #maximum sensing distance in km
    
    #mortality
    tsMortality = (1.0/maxPLD) / 120 #assumes that by max PLD, all larvae are dead
    fieldset.add_constant("eggMort", tsMortality)
    fieldset.add_constant("larvaeMort", tsMortality)
    
    #development
    #stage 0-1 larvae
    fieldset.add_constant("zeroBuoyancy", 0.0001) #in m/s
    fieldset.add_constant("zeroHorizontal", 1) #0 for no behavior, 1 for yes behavior
    fieldset.add_constant("zeroVertical", 5) #0 for no behavior, 1 for positive phototaxis, 2 for negative phototaxis, 3 for neg. during full moon, 4 for depth only, 5 for constant down, 6 for constant up
    fieldset.add_constant("zeroGrowthType", 2) #0 for age-dependent speed, 1 for length-dependent speed, 2 for absolute speed
    fieldset.add_constant("zeroSpeedMin", 0.002) 
    fieldset.add_constant("zeroSpeedMax", 0.002)
    fieldset.add_constant("zeroSpeed_b", 0) #if function and not absolute
    fieldset.add_constant("zeroMinDepth", -100) #in meters from surface
    fieldset.add_constant("zeroMaxDepth", -100) #in meters from surface
    fieldset.add_constant("devOneTime", -100) #when we advance to next stage (in days)
     
    #stage 1-2 larvae
    fieldset.add_constant("oneBuoyancy", -100) #in m/s
    fieldset.add_constant("oneHorizontal", -100) #0 for no behavior, 1 for yes behavior
    fieldset.add_constant("oneVertical", -100) #0 for no behavior, 1 for positive phototaxis, 2 for negative phototaxis, 3 for neg. during full moon, 4 for depth only
    fieldset.add_constant("oneGrowthType", -100) #0 for age-dependent speed, 1 for length-dependent speed, 2 for absolute speed
    fieldset.add_constant("oneSpeedMin", -100) 
    fieldset.add_constant("oneSpeedMax", -100)
    fieldset.add_constant("oneSpeed_b", -100) #if function and not absolute
    fieldset.add_constant("oneMinDepth", -100) #in meters from surface
    fieldset.add_constant("oneMaxDepth", -100) #in meters from surface
    fieldset.add_constant("devTwoTime", -100) #when we advance to next stage (in days)
    
    #stage 2-3 larvae
    fieldset.add_constant("twoBuoyancy", -100) #in m/s
    fieldset.add_constant("twoHorizontal", -100) #0 for no behavior, 1 for yes behavior
    fieldset.add_constant("twoVertical", -100) #0 for no behavior, 1 for positive phototaxis, 2 for negative phototaxis, 3 for neg. during full moon, 4 for depth only
    fieldset.add_constant("twoGrowthType", -100) #0 for age-dependent speed, 1 for length-dependent speed, 2 for absolute speed
    fieldset.add_constant("twoSpeedMin", -100) 
    fieldset.add_constant("twoSpeedMax", -100)
    fieldset.add_constant("twoSpeed_b", -100) #if function and not absolute
    fieldset.add_constant("twoMinDepth", -100) #in meters from surface
    fieldset.add_constant("twoMaxDepth", -100) #in meters from surface
    fieldset.add_constant("devThreeTime", -100)
    
    #stage 3-4 larvae
    fieldset.add_constant("threeBuoyancy", -100) #in m/s
    fieldset.add_constant("threeHorizontal", -100) #0 for no behavior, 1 for yes behavior
    fieldset.add_constant("threeVertical", -100) #0 for no behavior, 1 for positive phototaxis, 2 for negative phototaxis, 3 for neg. during full moon, 4 for depth only
    fieldset.add_constant("threeGrowthType", -100) #0 for age-dependent speed, 1 for length-dependent speed, 2 for absolute speed
    fieldset.add_constant("threeSpeedMin", -100) 
    fieldset.add_constant("threeSpeedMax", -100)
    fieldset.add_constant("threeSpeed_b", -100) #if function and not absolute
    fieldset.add_constant("threeMinDepth", -100) #in meters from surface
    fieldset.add_constant("threeMaxDepth", -100) #in meters from surface
    fieldset.add_constant("devFourTime", -100)
    
    #stage 4-5 larvae
    fieldset.add_constant("fourBuoyancy", -100) #in m/s
    fieldset.add_constant("fourHorizontal", -100) #0 for no behavior, 1 for yes behavior
    fieldset.add_constant("fourVertical", -100) #0 for no behavior, 1 for positive phototaxis, 2 for negative phototaxis, 3 for neg. during full moon, 4 for depth only
    fieldset.add_constant("fourGrowthType", -100) #0 for age-dependent speed, 1 for length-dependent speed, 2 for absolute speed
    fieldset.add_constant("fourSpeedMin", -100) 
    fieldset.add_constant("fourSpeedMax", -100)
    fieldset.add_constant("fourSpeed_b", -100) #if function and not absolute
    fieldset.add_constant("fourMinDepth", -100) #in meters from surface
    fieldset.add_constant("fourMaxDepth", -100) #in meters from surface
    fieldset.add_constant("devFiveTime", -100)
    
    #full moons
    #fieldset.add_constant("fullMoons", 0)
    
    #field limits
    fieldset.add_constant("lonmin", min(fieldset.U[1].grid.lon))
    fieldset.add_constant("lonmax", max(fieldset.U[1].grid.lon))
    fieldset.add_constant("latmin", min(fieldset.V[1].grid.lat))
    fieldset.add_constant("latmax", max(fieldset.V[1].grid.lat))
    
    ################ ADDING FIELDS
    
    #add settlement boundary to fieldset
    settle = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2023_spawnsettle/p_meandrina_settle.nc", 
            variable=('settle','settle'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
            
    #add fields for random walk to fieldset
    distance = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2023_spawnsettle/p_mea_distance_final_2023.nc", 
            variable=('distance','distance'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)
    angle = Field.from_netcdf("/home/sshedd/working/emily_parcels/spawn_files/2023_spawnsettle/p_mea_angle_final_2023.nc", 
            variable=('angle','angle'), dimensions={'lon':'lon', 'lat':'lat'}, allow_time_extrapolation=True)

    fieldset.add_field(settle)
    fieldset.settle.interp_method = "nearest"
    fieldset.add_field(distance)
    fieldset.distance.interp_method = "linear"
    fieldset.add_field(angle)
    fieldset.angle.interp_method = "nearest"

    #save spawning files
    spawn_file = "/home/sshedd/working/emily_parcels/spawn_files/2023_spawnsettle/p_meandrina_spawn.nc"
    depth_file = "/home/sshedd/working/emily_parcels/spawn_files/2023_habitat/Depth.nc"
    
    #set up spawning sites
    #lat, lon, dep = getSpawningSites(spawn_file, 
     #                                depth_file,
      #                               maxAdultDepth,
       #                              n_eggs)

    # ---- CUSTOM SPAWN SITES FROM CSV (bypass getSpawningSites) ----
    import pandas as pd
    csv_path = "/home/sshedd/working/emily_parcels/spawn_files/my_spawn_sites.csv.2"

    sites_df = pd.read_csv(csv_path)
    cols_lower = {c.lower(): c for c in sites_df.columns}
    if 'lat' not in cols_lower or 'lon' not in cols_lower:
        raise ValueError(f"CSV {csv_path} must have columns 'lat' and 'lon' (depth optional). Found columns: {list(sites_df.columns)}")

    lat = sites_df[cols_lower['lat']].astype(float).tolist()
    lon = sites_df[cols_lower['lon']].astype(float).tolist()
    if 'depth' in cols_lower:
        dep = sites_df[cols_lower['depth']].astype(float).fillna(5.0).tolist()
    else:
        dep = [5.0] * len(lat)

    # Convert negative longitudes to 0–360°E if your grid is 0–360 (MITgcm in Hawaiʻi commonly is)
    if any(L < 0 for L in lon):
        lon = [L + 360.0 if L < 0 else L for L in lon]

    # OPTIONAL: downsample for quick test
    # N = 10; lat, lon, dep = lat[:N], lon[:N], dep[:N]

    assert len(lat) == len(lon) == len(dep), "lat, lon, dep must be same length"
    print(f"Loaded {len(lat)} custom spawn sites from {csv_path}")
    # ---- END CUSTOM SPAWN SITES ----



    #set up spawning events
    #starttime, events = getSpawningEvents(
     #   monthSpawner,
      #  moonSpawner,
      #  timedSpawner,
      #  dt,
      #  simdays
    #)
    # --- set up spawning events (using lanaau_build_spawn_2023 signature) ---
    # Ensure types are what the helper expects
    _peakMonths = [] if (isinstance(peakSpawnMonths, str) and peakSpawnMonths.upper() == "NA") else peakSpawnMonths
    if isinstance(_peakMonths, int):
        _peakMonths = [_peakMonths]
    if isinstance(spawnPhase, str):
        spawnPhase = [spawnPhase]

    starttime, events = getSpawningEvents(
        monthSpawner,         # True/False: seasonal vs year-round
        spawnMonths,          # e.g., [4, 5] for April–May
        timedSpawner,         # True/False: restrict time-of-day window
        timeSpawnBegin,       # datetime.timedelta(hours=...)
        timeSpawnEnd,         # datetime.timedelta(hours=...)
        peakSpawner,          # True/False: has peak season?
        _peakMonths,          # [] or list of ints for peak months
        moonSpawner,          # True/False: lunar-linked spawning
        spawnPhase            # list of phase names, e.g. ["Waning gibbous"]
    )
    print(f"Spawning schedule built: starttime={starttime} sec, segments={len(events)}")
    # --- end spawning events setup ---



    #FOR TESTING
    #idx = random.sample(range(len(lon)), 1000)
    #lon = [lon[i] for i in idx]
    #lat = [lat[i] for i in idx]
    #dep = [dep[i] for i in idx]

    #set up particle set
    LarvalParticle = define_Larval_Class(typeSpawner)
    pset = ParticleSet.from_list(fieldset = fieldset,
                                 pclass = LarvalParticle,
                                 lon = lon, lat = lat,
                                 depth = dep, time = starttime,
                                 repeatdt = datetime.timedelta(minutes=dt).total_seconds(),
                                 siteID = range(1, (len(lon) + 1)),
                                 nowBuoyancy = [firstBuoy] * len(lon))

    print("Number of spawning sites:")
    print(len(lon))

    #functions to kernels
    #kernels = pset.Kernel(HandleBoundaries) + pset.Kernel(AdvectionRK4_3D_Lanaau) + pset.Kernel(RWDiffusion_Lanaau) + pset.Kernel(Buoyancy) + pset.Kernel(KillFish) + pset.Kernel(SettleFish) + pset.Kernel(GrowFish)
    #kernels = pset.Kernel(InitDepth) + pset.Kernel(GrowFishComplex) + pset.Kernel(KillFish) + pset.Kernel(SettleFish2) + pset.Kernel(AdvectionRK4_3D_Lanaau) + pset.Kernel(RWDiffusion_Lanaau) + pset.Kernel(Behavior) + pset.Kernel(Buoyancy) + pset.Kernel(HandleBoundaries) + pset.Kernel(KeepInOcean) + pset.Kernel(CheckOutOfBounds)
    #output_file = ParticleFile(name="/home/sshedd/working/emily_parcels/outputs/practice_run.zarr", outputdt=timedelta(hours = 24))
    #counter = 1
    #n_events = len(events)
    print("Reached point A: after pset creation")
    # Build kernel stack (add your physical/behavior kernels here)
    kernels = (
        pset.Kernel(InitDepth)
        + pset.Kernel(GrowFishComplex)
        + pset.Kernel(KillFish)
        + pset.Kernel(SettleFish2)
        + pset.Kernel(AdvectionRK4_3D_Lanaau)
        + pset.Kernel(KeepInOcean)        # v3 safety: handle surface errors
        + pset.Kernel(CheckOutOfBounds)   # v3 safety: delete OOB particles
    )
    print("Reached point B: kernels built")
    # --- Output setup (Zarr folder). Keep this simple and robust.
    from parcels import ParticleFile    # <- make sure this import exists at top or leave here
    from datetime import timedelta
    from pathlib import Path
    import shutil

    #out_path = Path("/home/sshedd/working/emily_parcels/outputs/practice_run2.zarr")
    #out_path.parent.mkdir(parents=True, exist_ok=True)
    #if out_path.exists():
    #    shutil.rmtree(out_path)

    #output_file = ParticleFile(name=str(out_path), outputdt=datetime.timedelta(hours=6))
    #output_file = ParticleFile(str(out_path), pset, outputdt=timedelta(hours=1))
    
    print("Skipping ParticleFile for debug run")
    output_file = None
    
    counter = 1
    n_events = len(events)
    print(f"Reached point C: n_events = {n_events}")
    
    for duration in events:
        
         print("On " + str(counter) + " out of " + str(n_events) + " steps")

         #if weʻre at an EVEN number, we donʻt spawn
         if (counter % 2) == 0:
             pset.repeatdt = None

         #if weʻre at an ODD number, we do spawn
         else:
             pset.repeatdt = datetime.timedelta(minutes=dt).total_seconds()
             #pset.repeat_starttime = runningTime

         #pset.execute(kernels, runtime=duration, dt=datetime.timedelta(minutes=dt), output_file=output_file)
         pset.execute(kernels, runtime=duration, dt=datetime.timedelta(minutes=dt), output_file=output_file)

         #pset.execute(kernels, runtime=duration, dt=timedelta(minutes=dt), output_file=output_file)
                      #recovery={StatusCode.ErrorThroughSurface: submergeParticle,
                                #StatusCode.ErrorOutOfBounds: deleteParticle})

         #pset.execute(kernels, runtime=duration, dt=timedelta(minutes=dt), output_file=output_file,
         #            recovery={StatusCode.ErrorOutOfBounds: deleteParticle})
         counter += 1
 # flush to disk
#    output_file.close()

main()
