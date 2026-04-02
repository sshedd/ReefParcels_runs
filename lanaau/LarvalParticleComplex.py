from parcels import JITParticle, Variable
import numpy as np
from operator import attrgetter

def define_Larval_Class(spawnType, ptype=JITParticle, ):
  class LarvalParticle(ptype):

    #broadcast (0) or brood (1)
    spawnType = Variable('spawnType', dtype=np.float32, to_write=False)

    ##set up whether particle is pre settlement (0), settling (1), settled (2), or past (3)
    settleStage = Variable('settleStage', initial = 0, dtype=np.float32, to_write=True)

    ##set up whether particle is dead (1) or not (0)
    isDead = Variable('isDead', initial = 0, dtype=np.float32, to_write=True)
    ##set up how a particle died - not dead (0), OOB (1), beached (2), random (3), max PLD (4)
    howDead = Variable('howDead', initial = 0, dtype=np.float32, to_write=True)

    ##location variables
    siteID = Variable('siteID', initial = attrgetter('siteID'), dtype=np.float32, to_write=True)

    #life history variables
    birthday = Variable('birthday', initial=attrgetter('time'), dtype=np.float32, to_write=True) #date of birth (seconds)
    settleday = Variable('settleday', dtype=np.float32, to_write=True) #date of settlement 
    age = Variable('age', initial=0, dtype=np.float32, to_write=True) #age in days
    isBeached = Variable('isBeached', initial=0, dtype=np.float32, to_write=False) #whether or not we are beached
    timeBeached = Variable('timeBeached', initial=0, dtype=np.float32, to_write=False) #time beached in hours
    lastBeached = Variable('lastBeached', initial=0, dtype=np.float32, to_write=False) #set up last timestep particle was bleached
    whereBeachedLat = Variable('whereBeachedLat', initial=0, dtype=np.float32, to_write=False)
    whereBeachedLon = Variable('whereBeachedLon', initial=0, dtype=np.float32, to_write=False)
    OOB = Variable('OOB', initial=0, dtype=np.float32, to_write=False) #whether or not weʻre OOB (1) or not (0)
    inSettleZone = Variable('inSettleZone', initial=1, dtype=np.float32, to_write=False) #whether or not we're currently in a settle zone (1) or not (0)
    
    #development variables
    length = Variable('length', initial=0, dtype=np.float32, to_write=True)
    lastTheta = Variable('lastTheta', initial=-100, dtype=np.float32, to_write=True)
    thisTheta = Variable('thisTheta', initial=-100, dtype=np.float32, to_write=True)
  
    #advection variables
    distance = Variable('distance', initial=0, dtype=np.float32, to_write=True) #total distance travelled
    prev_lon = Variable('prev_lon', initial=attrgetter('lon'), dtype=np.float64, to_write=True) #previous longitude
    prev_lat = Variable('prev_lat', initial=attrgetter('lat'), dtype=np.float64, to_write=True) #previous latitude
    prev_depth = Variable('prev_depth', initial=attrgetter('depth'), dtype=np.float64, to_write=True) #previous depth
    init_lon = Variable('init_lon', initial=attrgetter('lon'), dtype=np.float32, to_write=True) #initial longitude
    init_lat = Variable('init_lat', initial=attrgetter('lat'), dtype=np.float32, to_write=True) #initial latitude
    init_depth = Variable('init_depth', initial=attrgetter('depth'), dtype=np.float32, to_write=True)
    final_lon = Variable('final_lon', initial=0, dtype=np.float32, to_write=True)
    final_lat = Variable('final_lat', initial=0, dtype=np.float32, to_write=True)
    final_depth = Variable('final_depth', initial=0, dtype=np.float32, to_write=True)
    
    #set up complex larval stage: egg (0), larva1 (1), larva2 (2), larva3 (3), larva4 (4), juvenile (10)
    complexStage = Variable('complexStage', dtype=np.float32, to_write=True)
    preferredDepth = Variable('preferredDepth', initial=0, dtype=np.float32, to_write=True)
    nowSpeedMin = Variable('nowSpeedMin', initial=0, dtype=np.float32, to_write=True)
    nowSpeedMax = Variable('nowSpeedMax', initial=0, dtype=np.float32, to_write=True)
    nowGrowthType = Variable('nowGrowthType', initial=0, dtype=np.float32, to_write=False)
    nowVertical = Variable('nowVertical', initial=0, dtype=np.float32, to_write=False)
    nowHorizontal = Variable('nowHorizontal', initial=0, dtype=np.float32, to_write=False)
    nowBuoyancy = Variable('nowBuoyancy', initial=attrgetter('nowBuoyancy'), dtype=np.float32, to_write=False)
    tempSpeedMin = Variable('tempSpeedMin', initial=0, dtype=np.float32, to_write=False)
    tempSpeedMax = Variable('tempSpeedMax', initial=0, dtype=np.float32, to_write=False)
    tempSpeed_b = Variable('tempSpeed_b', initial=0, dtype=np.float32, to_write=False)
    
    def __init__(self, *args, **kwargs):
      """Custom initialisation function which calls the base
      initialisation and adds the instance variable p"""
      super(LarvalParticle, self).__init__(*args, **kwargs)
      self.spawnType = spawnType

      if spawnType == 0:
        complexStage = Variable('complexStage', initial = 0, dtype=np.float32, to_write=True) #only used in complex model
        
      elif spawnType == 1:
        complexStage = Variable('complexStage', initial = 1, dtype=np.float32, to_write=True) #only used in complex model
            
  return LarvalParticle
