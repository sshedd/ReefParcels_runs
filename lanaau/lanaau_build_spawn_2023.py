import netCDF4 as nc
import numpy as np
import random
import datetime
from lanaau_temporal import getModelDay, getDate, getMoon, getTimeOfDay

def getSpawningSites(spawn_file, depth_file, maxDepth, n, typeSpawner):
  #read in files
  spawn_1km = nc.Dataset(spawn_file)
  depth_1km = nc.Dataset(depth_file)
  
  #set how many times we release a particle from a cell per time step
  n_particles = n
              
  #generate random spawn sites within each habitat cell
  coordx = []; coordy = []
  spawn_array_1km = spawn_1km["spawn"]
  for x in range(0, (spawn_array_1km.shape[1])):
    for y in range(0, (spawn_array_1km.shape[0])):
      if spawn_array_1km[y,x] != "--":
        for i in range(n_particles):
          coordx.append(x)
          coordy.append(y)
    
  n = len(coordx)    
  dx = (max(spawn_1km["lon"]) - min(spawn_1km["lon"])) / (spawn_array_1km.shape[1] * 2)
  dy = (max(spawn_1km["lat"]) - min(spawn_1km["lat"])) / (spawn_array_1km.shape[0] * 2)
  lon = spawn_1km["lon"][coordx].tolist()
  lat = spawn_1km["lat"][coordy].tolist()
  newlat_1km = (lat + np.random.uniform(low = -dy, high = dy, size = n)).tolist()
  newlon_1km = (lon + np.random.uniform(low = -dx, high = dx, size = n)).tolist()

  #randomly assign starting depths
  rand_depths_1km = []
  
  if typeSpawner == 0: #broadcast - max depth to surface
    for i in range(n):
      try:
        rd = random.randrange(0, maxDepth+1, 1)
      except ValueError:
        rd = maxDepth #for opihi at surface
        
      depth = depth_1km["maxdepth"][coordy[i], coordx[i]].tolist()
      if (not np.isnan(depth)) and (rd > depth):
        rd = depth
      rand_depths_1km.append(rd)
    
  elif typeSpawner == 1: #brood - max depth or ocean floor, whichever shallower
    for i in range(n):
      depth = depth_1km["maxdepth"][coordy[i], coordx[i]].tolist()
      
      if depth > maxDepth:
        rand_depths_1km.append(maxDepth)
      else:
        rand_depths_1km.append(depth)

  return(newlat_1km, newlon_1km, rand_depths_1km)
  

def getSpawningEvents(monthSpawner, spawnMonths, timedSpawner, timeSpawnBegin, \
  timeSpawnEnd, peakSpawner, peakSpawnMonths, moonSpawner, spawnPhase):
    
  startTimestamp = getDate(0) #model starts
  endTimestamp = getDate(851) #model stops
  
  calendar = [startTimestamp + datetime.timedelta(days=x) for x in range((endTimestamp - startTimestamp).days + 1)]
  
  spawningdays = []
  peakspawningdays = []
  addSpawn = True
  addPeak = True
  
  for tmpDay in calendar:
    if monthSpawner and tmpDay.month not in spawnMonths: #if we donʻt spawn this month
      addSpawn = False
    if (not peakSpawner) or (peakSpawner and tmpDay.month not in peakSpawnMonths): #if we donʻt peak spawn this month
      addPeak = False
    if (addSpawn and moonSpawner and getMoon(tmpDay) not in spawnPhase): #if we donʻt spawn this moon phase
      addSpawn = False
      addPeak = False
      
    if addSpawn == True:
      spawningdays.append(tmpDay)
    if addPeak == True:
      peakspawningdays.append(tmpDay)
      
    addSpawn = True
    addPeak = True
    
  #moon, month, no peak, no time of day
  #C. exarata
  if moonSpawner and monthSpawner and not peakSpawner and not timedSpawner:
    #generate one random spawning event per site per spawning day
    spawn_starts = [x + datetime.timedelta(hours = random.randint(1, 23)) for x in spawningdays]
    spawn_stops = [x + datetime.timedelta(hours = 1) for x in spawn_starts]
    
  #moon, month, no peak, time of day
  #P. meandrina
  if moonSpawner and monthSpawner and not peakSpawner and timedSpawner:
    #add time to spawning days
    spawn_starts = [x + timeSpawnBegin for x in spawningdays]
    spawn_stops = [x + timeSpawnEnd for x in spawningdays]
    
  #no moon, no month, peak, no time of day
  #Panulirus spp.
  if not moonSpawner and not monthSpawner and peakSpawner and not timedSpawner:
    
    print("panulirus spp.")
    
    spawn_starts = []
    spawn_stops = []
    counter = 0
    
    peak_dur = datetime.timedelta(hours = 2)
    non_peak_dur = datetime.timedelta(minutes = 12)
    
    for x in spawningdays:
      if x in peakspawningdays: #one per site every five days
        if (counter % 5) == 0:
          tmp_start = datetime.timedelta(hours = random.randint(1, 22))
          spawn_starts.append(x + tmp_start)
          spawn_stops.append(x + tmp_start + peak_dur)
        
      else:
        if (counter % 10) == 0:
          tmp_start = datetime.timedelta(hours = random.randint(1, 22))
          spawn_starts.append(x + tmp_start)
          spawn_stops.append(x + tmp_start + non_peak_dur)   
        
      counter += 1
    
  #no moon, no month, peak, time of day
  #Octopus cyanea
  if not moonSpawner and not monthSpawner and peakSpawner and timedSpawner:
    
    print("octopus cyanea")
    
    spawn_starts = []
    spawn_stops = []
    counter = 0
    
    #non-peak duration
    nonPeakEnd = timeSpawnBegin + (((timeSpawnEnd + datetime.timedelta(hours = 24)) - timeSpawnBegin) / 10)
    
    for x in spawningdays:
      #during peak, generate one event per site every five days
      if x in peakspawningdays:
        if (counter % 5) == 0:
          spawn_starts.append(x + timeSpawnBegin)
          spawn_stops.append(x + datetime.timedelta(hours = 24) + timeSpawnEnd) #NOTE - increased by 1 day bc spanning midnight
          
      #during non-peak, 1/10 particles every ten days
      else:
        if (counter % 10) == 0:
          spawn_starts.append(x + timeSpawnBegin)
          spawn_stops.append(x + nonPeakEnd) 
          
      counter += 1
      
  #no moon, month, peak, time of day    
  #C. perspicillatus, S. rubroviolaceus, P. porphyreus, C. strigosus
  if not moonSpawner and monthSpawner and peakSpawner and timedSpawner:
    
    print("c perspicillatus")
    print("s rubroviolaceus")
    print("p porphyreus")
    print("c strigosus")
    print("c. melampygus")
    
    spawn_starts = []
    spawn_stops = []
    counter = 0
    
    #non-peak duration
    nonPeakEnd = timeSpawnBegin + (((timeSpawnEnd + datetime.timedelta(hours = 24)) - timeSpawnBegin) / 10)
    
    for x in spawningdays:
      #during peak, generate one spawning event per site every five days
      if x in peakspawningdays:
        if (counter % 5) == 0:
          #c. melampygus near dusk and dawn
          if timeSpawnBegin == datetime.timedelta(hours = 18) and timeSpawnEnd == datetime.timedelta(hours = 6):
            spawn_starts.append(x + timeSpawnBegin)
            spawn_stops.append(x + (timeSpawnBegin + datetime.timedelta(hours = 3)))
            spawn_starts.append(x + datetime.timedelta(hours = 24) + (timeSpawnEnd - datetime.timedelta(hours = 3)))
            spawn_stops.append(x + datetime.timedelta(hours = 24) + timeSpawnEnd)
            
          else:
            spawn_starts.append(x + timeSpawnBegin)
            if timeSpawnBegin > timeSpawnEnd:
              spawn_stops.append(x + datetime.timedelta(hours = 24) + timeSpawnEnd)
            else:
              spawn_stops.append(x + timeSpawnEnd)
      
      #during non-peak, 1/10 particles every ten days
      else:
        if (counter % 10) == 0:
          #c. melampygus near dusk and dawn
          if timeSpawnBegin == datetime.timedelta(hours = 18) and timeSpawnEnd == datetime.timedelta(hours = 6):
            spawn_starts.append(x + timeSpawnBegin)
            spawn_stops.append(x + (timeSpawnBegin + (datetime.timedelta(hours = 3) / 10)))
            spawn_starts.append(x + datetime.timedelta(hours = 24) + (timeSpawnEnd - (datetime.timedelta(hours = 3) / 10)))
            spawn_stops.append(x + datetime.timedelta(hours = 24) + timeSpawnEnd)
            
          else:
            spawn_starts.append(x + timeSpawnBegin)
            if timeSpawnBegin > nonPeakEnd:
              spawn_stops.append(x + datetime.timedelta(hours = 24) + nonPeakEnd)
            else:
              spawn_stops.append(x + nonPeakEnd)
          
      counter += 1

  #moon, month, peak, time of day
  #P. sexfilis, C. strigosus, C. melampygus
  if moonSpawner and monthSpawner and peakSpawner and timedSpawner:

    print("p. sexfilis")
    print("c. ignoblis")

    spawn_starts = []
    spawn_stops = []
    counter = 0
    
    #non-peak duration
    nonPeakEnd = timeSpawnBegin + (((timeSpawnEnd + datetime.timedelta(hours = 24)) - timeSpawnBegin) / 10)
    
    for x in spawningdays:
      #during peak, one spawning event per day
      if x in peakspawningdays:
        spawn_starts.append(x + timeSpawnBegin)
        if timeSpawnBegin > timeSpawnEnd:
          spawn_stops.append(x + datetime.timedelta(hours = 24) + timeSpawnEnd)
        else:
          spawn_stops.append(x + timeSpawnEnd)
      #during non-peak, 1/10 number of particles every other day
      else:
        if (counter % 2) == 0:
          spawn_starts.append(x + timeSpawnBegin)
          spawn_stops.append(x + nonPeakEnd)  
          
      counter += 1
  
  #get rid of starts/ends that go past the end of the model
  spawn_starts = [x for x in spawn_starts if x < endTimestamp]
  spawn_stops = [x for x in spawn_stops if x < endTimestamp]
  if len(spawn_starts) - len(spawn_stops) == 1:
    spawn_stops.append(endTimestamp)
  
  #build list of [spawning duration, no spawning duration, spawning duration, etc.]
  duration_spawning = []
  lastStop = ""
    
  for n in range(len(spawn_starts)):
    this_start = spawn_starts[n]
    this_stop = spawn_stops[n]
    
    if lastStop != "":
      #add no spawning duration
      duration_spawning.append(this_start - lastStop)
      
    #add spawning durations
    duration_spawning.append(this_stop - this_start)
    
    lastStop = this_stop
      
  #add extra chunk of time to the end
  #starttime = (spawningdays[0] - startTimestamp).total_seconds() #won't work if time-specific! try:
  starttime = (spawn_starts[0] - startTimestamp).total_seconds()
  end_chunk = (endTimestamp - startTimestamp) - sum(duration_spawning, datetime.timedelta()) - (spawn_starts[0] - startTimestamp)
  if end_chunk != 0:
    duration_spawning.append(end_chunk)
        
  return(starttime, duration_spawning)
