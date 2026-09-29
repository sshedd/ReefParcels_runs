from parcels import ParcelsRandom as ParcelsRandom
from random_von_mises import random_von_Mises
import math

def HandleBoundaries(particle, fieldset, time):
  
  if particle.isDead == 0 and not particle.settleStage == 2:
  
    plat = particle.lat
    plon = particle.lon

    #check if we're out of bounds and remove if we are
    if plat <= fieldset.latmin or plat >= fieldset.latmax or plon <= fieldset.lonmin or plon >= fieldset.lonmax:
      particle.OOB = 1
    
    #if not OOB, check for beaching
    #maybe just if last coords = this coords; could also check UVW
    #run this as very last step
    #maybe a while loop?
    else:
      if ((particle.lon == particle.prev_lon) and (particle.lat == particle.prev_lat)) or particle.isBeached == 1:
      
        pdp = particle.depth
        plb = particle.lastBeached
  
        Nc = fieldset.F_N[time, pdp, plat, plon]
        Wc = fieldset.F_W[time, pdp, plat, plon]
  
        #if beached, give it a gentle "push" away from land and towards the sea
        #same scale as diffusion
        D1 = 0.2; R1 = 6378137
        push_lat = (((ParcelsRandom.uniform(0., 1.) * math.sqrt(D1 * particle.dt)) / R1) * (180 / math.pi))
        push_lon = (((ParcelsRandom.uniform(0., 1.) * math.sqrt(D1 * particle.dt)) / (R1 * math.cos(math.pi * particle.lat / 180))) * (180 / math.pi))
        
        if Nc != 0: #sea to the north
          particle.lat += push_lat
        else:
          Sc = fieldset.F_S[time, pdp, plat, plon] #test for sea to the south
          if Sc != 0: 
            particle.lat -= push_lat
            
        if Wc != 0: #sea to the west
          particle.lon += push_lon
        else:
          Ec = fieldset.F_E[time, pdp, plat, plon] #test for sea to the east
          if Ec != 0: 
            particle.lon -= push_lon
  
        #update beaching status
        #keep track of when + where particle last beached
        particle.lastBeached = time
        particle.whereBeachedLat = plat
        particle.whereBeachedLon = plon
  
        #track time beached
        if (time - plb) == particle.dt:
          particle.timeBeached += particle.dt
        else:
          particle.timeBeached = particle.dt
  
        #test to see if we were successful
        (u, v, w) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        if u == 0 and v == 0:
          particle.isBeached = 1
        else:
          particle.isBeached = 0
  
      else:
        particle.isBeached = 0
    
def AdvectionRK4_3D_Lanaau(particle, fieldset, time):
    """Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.
    Function needs to be converted to Kernel object before execution"""
    
    #only run for living, not settled, non beached fish
    
    #print("advection")
    
    if particle.isDead == 0 and not particle.settleStage == 2 and particle.isBeached == 0 and particle.OOB == 0:
    #   
    #   advect = True
    #   
    #   #have to force interpolation on finer grid around coasts - hacky but works
    #   if particle.lon > 200 and particle.lon < 205.25 and particle.lat > 18.5 and particle.lat < 21.65:
    # 
    #     u1 = fieldset.U1[time, particle.depth, particle.lat, particle.lon]
    #     v1 = fieldset.V1[time, particle.depth, particle.lat, particle.lon]
    #     w1 = fieldset.W1[time, particle.depth, particle.lat, particle.lon]
    #     
    #     lon1 = particle.lon + u1*.5*particle.dt
    #     lat1 = particle.lat + v1*.5*particle.dt
    #     dep1 = particle.depth + w1*.5*particle.dt
    #     
    #     u2 = fieldset.U1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     v2 = fieldset.V1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     w2 = fieldset.W1[time + .5 * particle.dt, dep1, lat1, lon1]
    #     lon2 = particle.lon + u2*.5*particle.dt
    #     lat2 = particle.lat + v2*.5*particle.dt
    #     dep2 = particle.depth + w2*.5*particle.dt
    #     
    #     u3 = fieldset.U1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     v3 = fieldset.V1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     w3 = fieldset.W1[time + .5 * particle.dt, dep2, lat2, lon2]
    #     lon3 = particle.lon + u3*particle.dt
    #     lat3 = particle.lat + v3*particle.dt
    #     dep3 = particle.depth + w3*particle.dt
    #     
    #     u4 = fieldset.U1[time + particle.dt, dep3, lat3, lon3]
    #     v4 = fieldset.V1[time + particle.dt, dep3, lat3, lon3]
    #     w4 = fieldset.W1[time + particle.dt, dep3, lat3, lon3]
    #             
    #     #test for beaching
    #     if u4 == 0 and v4 == 0 and w4 == 0:
    #       advect = False
    #       particle.isBeached = 1
    #    
    #  #otherwise, use nested fields as usual
    #  else:
      advect = True
      n_valid_u = 0.
      n_valid_v = 0.
      n_valid_w = 0.

      (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
      lon1 = particle.lon + u1*.5*particle.dt
      lat1 = particle.lat + v1*.5*particle.dt
      dep1 = particle.depth + w1*.5*particle.dt
      if u1 != 0:
        n_valid_u += 1.
      if v1 != 0:
        n_valid_v += 1.
      if w1 != 0:
        n_valid_w += 1.
      
      (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
      lon2 = particle.lon + u2*.5*particle.dt
      lat2 = particle.lat + v2*.5*particle.dt
      dep2 = particle.depth + w2*.5*particle.dt
      if u2 != 0:
        n_valid_u += 1.
      if v2 != 0:
        n_valid_v += 1.
      if w2 != 0:
        n_valid_w += 1.
      
      (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
      lon3 = particle.lon + u3*particle.dt
      lat3 = particle.lat + v3*particle.dt
      dep3 = particle.depth + w3*particle.dt
      if u3 != 0:
        n_valid_u += 1.
      if v3 != 0:
        n_valid_v += 1.
      if w3 != 0:
        n_valid_w += 1.      
      
      (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
      if u4 != 0:
        n_valid_u += 1.
      if v4 != 0:
        n_valid_v += 1.
      if w4 != 0:
        n_valid_w += 1.
        
      #test for beaching
      if u4 == 0 and v4 == 0:
        advect = False
        particle.isBeached = 1
      
      #if we're not beached
      if advect:
        
        #exclude 0-velocity (land) from calculation
        #particle.lon += ((u1 + 2*u2 + 2*u3 + u4) / n_valid_u * particle.dt)
        #particle.lat += ((v1 + 2*v2 + 2*v3 + v4) / n_valid_v * particle.dt)
        #particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)
        

        # Preserve original zero-velocity handling, but avoid 0/0 -> NaN.
        # If all four samples for one component are zero, displacement
        # in that component is zero.
        if n_valid_u > 0:
          particle.lon += ((u1 + 2*u2 + 2*u3 + u4) / n_valid_u * particle.dt)

        if n_valid_v > 0:
          particle.lat += ((v1 + 2*v2 + 2*v3 + v4) / n_valid_v * particle.dt)

        if n_valid_w > 0:
          particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)
        #out1 = ((u1 + 2*u2 + 2*u3 + u4) / n_valid_u * particle.dt)
        #out2 = ((v1 + 2*v2 + 2*v3 + v4) / n_valid_v * particle.dt)
        
        #if move would place particle through surface or seafloor, don't make that move
        #test_dep = particle.depth + ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)
        #move = True
      
        #if test_dep <= 0:
        #  move = False
        #  particle.depth = 0.1
        #else:
        #  ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
        #  if test_dep > ocean_dep:
        #    move = False
        #    particle.depth = ocean_dep - 0.1
          
        ##set new depth
        #if move == True:
        #  particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / n_valid_w * particle.dt)

#TODO
def Behavior(particle, fieldset, time):
  
    #print("behavior")

    #only run for living, not egg, not settled, non beached fish
    if particle.isDead == 0 and not particle.complexStage == 0 and not particle.settleStage == 2 and particle.isBeached == 0 and particle.OOB == 0:
      #horizontal behavior
      if particle.nowHorizontal != 0 and particle.inSettleZone == 0: #avoids weird behavior with larvae swimming onto land
        #set fieldset.k = 3
        #set fieldset.B (detection distance) = 4
        
        #if we're swimming for the first time, set initial angle
        if particle.lastTheta == -100:
          #test_angle = fieldset.angle[time, particle.depth, particle.lat, particle.lon]
          #particle.lastTheta = ParcelsRandom.uniform(0, (2*math.pi))
          particle.lastTheta = ParcelsRandom.uniform(-math.pi, math.pi) #TESTING - or 0 to 2pi?

        else:
          particle.lastTheta = particle.thisTheta 
        
        C = fieldset.distance[time, particle.depth, particle.lat, particle.lon] #pull from distance-interpolated field
        
        #check to see if we're beyond detection distance, set no sensing ability
        if C > fieldset.B:
          sensingAbility = 0
        else:
          sensingAbility = 1 - (C / fieldset.B)
        
        #angle to closest reef, relative to field
        thetaGoal = fieldset.angle[time, particle.depth, particle.lat, particle.lon] #pull from angle-interpolated field
        
        if C <= 2 and thetaGoal == 0: 
          turningAngle = particle.lastTheta #avoid turning 180 when already within habitat
          phi = 0
        else:
          phi = math.atan2(math.sin(thetaGoal - particle.lastTheta), math.cos(thetaGoal - particle.lastTheta))
          turningAngle = -sensingAbility * phi #double check whether sensingAbility positive or negative
          #turningAngle = (-sensingAbility) * (particle.lastTheta - thetaGoal)
        
        mu = turningAngle
        kappa = fieldset.k
        
        a = 1.0 + math.sqrt(1.0 + 4.0 * kappa**2)
        b = (a - math.sqrt(2.0 * a)) / (2.0 * kappa)
        r = (1.0 + b**2) / (2.0 * b)
        
        # Pseudo-random numbers sampled from a uniform distribution [0,1]
        done = False
        
        while done == False:
        
          rand1 = ParcelsRandom.uniform(0., 1.)
          rand2 = ParcelsRandom.uniform(0., 1.)
          rand3 = ParcelsRandom.uniform(0., 1.)
          
          z = math.cos(math.pi * rand1)
          f = (1.0 + r * z) / (r + z)
          c = kappa * (r - f)
          
          if (rand3 - 0.5) >= 0:
            sign = 1
          else:
            sign = -1
          
          if (((c * (2.0 - c) - rand2) > 0.0) or ((math.log(c/rand2) + 1.0 - c) > 0.0)):
          
              done = True
              theta = math.fmod(sign * math.acos(f) + (mu - particle.lastTheta), (2*math.pi))   
              nowTheta = theta # *****triple***** CHECK
              #print(nowTheta)

        #get u/v and move fish accordingly
        swimmingSpeed = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax)
        
        #convert speed from m/s to decimal degrees/s (approximate)
        R = 6378137
        latSwimmingSpeed = (swimmingSpeed / R) * (180 / math.pi)
        lonSwimmingSpeed = (swimmingSpeed / (R * math.cos(math.pi * particle.lat / 180)))*(180 / math.pi) 
        
        #does the move put us on land? if so, try not swimming as far
        counter = 0
        tmpLonSwimmingSpeed = lonSwimmingSpeed
        tmpLatSwimmingSpeed = latSwimmingSpeed
        
        while counter < 3:
          uOrient = tmpLonSwimmingSpeed * math.cos(nowTheta)
          vOrient = tmpLatSwimmingSpeed * math.sin(nowTheta)
          test_lon = particle.lon + (uOrient * particle.dt)
          test_lat = particle.lat + (vOrient * particle.dt)  
          (u_test, v_test, w_test) = fieldset.UVW[time, particle.depth, test_lat, test_lon]
          
          #if non-zero, move and break loop
          if u_test != 0 and v_test != 0:
            if particle.nowHorizontal == 1:
              particle.lon += (uOrient * particle.dt)
              particle.lat += (vOrient * particle.dt)
            elif particle.nowHorizontal == 2:
              #get time of day
              secondsThisDay = math.fmod(time, 86400.)
              if secondsThisDay >= 21600 and secondsThisDay < 64800:
                particle.lon += (uOrient * particle.dt)
                particle.lat += (vOrient * particle.dt)  
            break
            
          #if zero, decrease distance
          else:
            tmpLonSwimmingSpeed = tmpLonSwimmingSpeed / 2
            tmpLatSwimmingSpeed = tmpLatSwimmingSpeed / 2
          
          counter += 1
        
        #OLD BELOW
        #if move puts particle on land, don't make it
        #uOrient = lonSwimmingSpeed * math.cos(nowTheta) #speed in decimal degrees, theta in radians
        #vOrient = latSwimmingSpeed * math.sin(nowTheta) #speed in decimal degrees, theta in radians
        #test_lon = particle.lon + (uOrient * particle.dt)
        #test_lat = particle.lat + (vOrient * particle.dt)
        #(u_test, v_test, w_test) = fieldset.UVW[time, particle.depth, test_lat, test_lon]
        
        #OLD BELOW
        #if u_test != 0 and v_test != 0:
        #  if particle.nowHorizontal == 1:
        #    particle.lon += (uOrient * particle.dt)
        #    particle.lat += (vOrient * particle.dt)
        #  elif particle.nowHorizontal == 2:
        #    #get time of day
        #    secondsThisDay = math.fmod(time, 86400.)
        #    if secondsThisDay >= 21600 and secondsThisDay < 64800:
        #      particle.lon += (uOrient * particle.dt)
        #      particle.lat += (vOrient * particle.dt)              
        
        #save value
        particle.thisTheta = nowTheta
        
      #vertical behavior DONE
      if particle.nowVertical != 0:
        
        #get current depth and change with swimming
        d_now = particle.depth
        d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
      
        #phototaxis DONE
        if particle.nowVertical == 1 or particle.nowVertical == 2 or particle.nowVertical == 3:
          
          remain_at_depth = False
        
          #get time of day
          daytime = False
          secondsThisDay = math.fmod(time, 86400.)
          
          if secondsThisDay >= 21600 and secondsThisDay < 64800:
            daytime = True
          
          #positive phototaxis
          if particle.nowVertical == 1:
            if daytime: #swim up
              particle.depth -= d_delta
              
            else: #remain at depth
              remain_at_depth = True
          
          #negative phototaxis
          elif particle.nowVertical == 2:
            if daytime: #remain at depth
              remain_at_depth = True
            else: #swim up
              particle.depth -= d_delta
          
          # #negative phototaxis during full moon
          # elif particle.nowVertical == 3:
          #   
          #   #I am so sad that I'm hardcoding this but I will try to fix it later
          #   #TODO: some tracker module that gives current moon phase given initial phase
          #   #can have this track time of day too
          #   
          #   if not daytime:
          #     
          #     moon_c = 0
          #     time_test = floor(time / (60*60*24))
          #     moon_swimming = False
          #     
          #     while moon_c < len(fullMoons):
          #       if time_test == fullMoons[moon_c]:
          #         moon_swimming = True
          #       moon_c += 1
          #       
          #     if moon_swimming == True:
          #       particle.depth -= d_delta
          #     else:
          #       remain_at_depth = True
          #       
          #   else:
          #     remain_at_depth = True
        
        #depth constraint
        if (particle.nowVertical == 4) or (remain_at_depth == True):
          
          #get difference between current depth and preferred depth
          d_pref = particle.preferredDepth
          d_dif = math.fabs(d_now - d_pref)

          #are we far enough away to move?
          if d_dif > 1:
            if d_now < d_pref: #above preferred depth, swim down
              if d_delta >= d_dif:
                particle.depth = d_pref
              else:
                particle.depth += d_delta
                
            elif d_now > d_pref: #below preferred depth, swim up
              if d_delta >= d_dif:
                particle.depth = d_pref
              else:
                particle.depth -= d_delta
            
        #swimming down
        if (particle.nowVertical == 5):
          if d_now > maxAdultDepth:
            particle.depth -= d_delta
          else:  
            particle.depth += d_delta
          
        #swimming up
        if (particle.nowVertical == 6):
          particle.depth -= d_delta
            
        #finally, make sure we're in bounds
        if particle.depth <= 0:
          particle.depth = 1
        else:
          ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
          if particle.depth >= ocean_dep and ocean_dep != 0:
            particle.depth = ocean_dep - 1
        
        
          #if within depth bounds, swim up, down, or stay
          #else:
          #  vert_chance = ParcelsRandom.uniform(0, 1)
          #  if vert_chance < (1./3.):
          #    d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
          #    particle.depth += d_delta
          #  elif vert_chance > (2./3.):
          #    d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
          #    particle.depth -= d_delta
 
def Behavior_vertonly(particle, fieldset, time):

    #only run for living, not egg, not settled, non beached fish
    if particle.isDead == 0 and not particle.complexStage == 0 and not particle.settleStage == 2 and particle.isBeached == 0 and particle.OOB == 0:
      # #horizontal behavior
      # if particle.nowHorizontal != 0:
      #   #set fieldset.k = 3
      #   #set fieldset.B (detection distance) = 5
      #   
      #   #biased correlated random walk
      #   
      #   #if we're swimming for the first time, set initial angle
      #   if particle.lastTheta == -100:
      #     #test_angle = fieldset.angle[time, particle.depth, particle.lat, particle.lon]
      #     particle.lastTheta = ParcelsRandom.uniform(0, (2*math.pi))
      # 
      #   else:
      #     particle.lastTheta = particle.thisTheta 
      #   
      #   C = fieldset.distance[time, particle.depth, particle.lat, particle.lon] #pull from distance-interpolated field
      #   sensingAbility = 1 - (C / fieldset.B)
      #   
      #   thetaGoal = fieldset.angle[time, particle.depth, particle.lat, particle.lon] #pull from angle-interpolated field
      #   turningAngle = (-sensingAbility) * (particle.lastTheta - thetaGoal)
      #   
      #   mu = turningAngle
      #   kappa = fieldset.k
      #   
      #   #nextTheta = random_von_Mises(fieldset.k, turningAngle, lastTheta)
      #   a = 1.0 + math.sqrt(1.0 + 4.0 * kappa**2)
      #   b = (a - math.sqrt(2.0 * a)) / (2.0 * kappa)
      #   r = (1.0 + b**2) / (2.0 * b)
      #   
      #   # Pseudo-random numbers sampled from a uniform distribution [0,1]
      #   done = False
      #   
      #   while done == False:
      #   
      #     rand1 = ParcelsRandom.uniform(0., 1.)
      #     rand2 = ParcelsRandom.uniform(0., 1.)
      #     rand3 = ParcelsRandom.uniform(0., 1.)
      #     
      #     z = math.cos(math.pi * rand1)
      #     f = (1.0 + r * z) / (r + z)
      #     c = kappa * (r - f)
      #     
      #     if (rand3 - 0.5) >= 0:
      #       sign = 1
      #     else:
      #       sign = -1
      #     
      #     if (((c * (2.0 - c) - rand2) > 0.0)  or ((math.log(c/rand2) + 1.0 - c) > 0.0)):
      #     
      #         done = True
      #         theta = math.fmod((sign * math.acos(f) + mu), (2*math.pi))   
      #         nowTheta = theta + particle.lastTheta
      # 
      #   #get u/v and move fish accordingly
      #   swimmingSpeed = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax)
      #   
      #   #convert speed from m/s to decimal degrees/s (approximate)
      #   R = 6378137
      #   latSwimmingSpeed = (swimmingSpeed / R) * (180 / math.pi)
      #   lonSwimmingSpeed = (swimmingSpeed / (R * math.cos(math.pi * particle.lat / 180)))*(180 / math.pi) 
      #   
      #   #if move puts particle on land, don't make it
      #   uOrient = lonSwimmingSpeed * math.cos(nowTheta) #speed in decimal degrees, theta in radians
      #   vOrient = latSwimmingSpeed * math.sin(nowTheta) #speed in decimal degrees, theta in radians
      #   test_lon = particle.lon + (uOrient * particle.dt)
      #   test_lat = particle.lat + (vOrient * particle.dt)
      #   (u_test, v_test, w_test) = fieldset.UVW[time, particle.depth, test_lat, test_lon]
      #   
      #   if u_test != 0 and v_test != 0:
      #     particle.lon += (uOrient * particle.dt)
      #     particle.lat += (vOrient * particle.dt)
      #   
      #   #save value
      #   particle.thisTheta = nowTheta
        
      #vertical behavior DONE
      if particle.nowVertical != 0:
      
        #phototaxis DONE
        if particle.nowVertical == 1 or particle.nowVertical == 2 or particle.nowVertical == 3:
          
          #get current depth and change with swimming
          d_now = particle.depth
          d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
          remain_at_depth = False
        
          #get time of day
          daytime = False
          secondsThisDay = math.fmod(time, 86400.)
          
          if secondsThisDay >= 21600 and secondsThisDay < 64800:
            daytime = True
          
          #positive phototaxis
          if particle.nowVertical == 1:
            if daytime: #swim up
              particle.depth -= d_delta
              
            else: #remain at depth
              remain_at_depth = True
          
          #negative phototaxis
          elif particle.nowVertical == 2:
            if daytime: #remain at depth
              remain_at_depth = True
            else: #swim up
              particle.depth -= d_delta
          
          #negative phototaxis during full moon
          elif particle.nowVertical == 3:
            if not daytime and floor(time / (60*60*24)): #in fieldset.fullMoons: #swim up TODO moons
              particle.depth -= d_delta
            else: #remain at depth
              remain_at_depth = True
        
        #depth constraint DONE
        if (particle.nowVertical == 4) or (remain_at_depth == True):
          
          #if above preferred depth, swim down
          if d_now < particle.preferredDepth:
            particle.depth += d_delta
          
          #if below preferred depth, swim up
          elif d_now > particle.preferredDepth:
            particle.depth -= d_delta
            
        #finally, make sure we're in bounds
        if particle.depth <= 0:
          particle.depth = 1
        else:
          ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
          if particle.depth >= ocean_dep and ocean_dep != 0:
            particle.depth = ocean_dep - 1
        
        
          #if within depth bounds, swim up, down, or stay
          #else:
          #  vert_chance = ParcelsRandom.uniform(0, 1)
          #  if vert_chance < (1./3.):
          #    d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
          #    particle.depth += d_delta
          #  elif vert_chance > (2./3.):
          #    d_delta = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax) * particle.dt
          #    particle.depth -= d_delta
 
def Buoyancy(particle, fieldset, time):
    
    #only run for living, not settled, non beached fish
    if particle.isDead == 0 and not particle.settleStage == 2 and particle.isBeached == 0 and particle.OOB == 0:
   
      #add buoyancy
      b1 = particle.nowBuoyancy
      particle.depth += (b1 * particle.dt)
      
      #test for breaking the surface or ocean floor
      if particle.depth <= 0:
        particle.depth = 1
      else:
        ocean_dep = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
        if particle.depth >= ocean_dep and ocean_dep != 0:
          particle.depth = ocean_dep - 1
         
def RWDiffusion_Lanaau(particle, fieldset, time):
  #only run for living, not settled, non beached fish
  if particle.isDead == 0 and not particle.settleStage == 2 and particle.isBeached == 0 and particle.OOB == 0:
  
    D = 0.2; R = 6378137 #D = diffusion constant, R = radius of earth in meters
    particle.lat += ((ParcelsRandom.uniform(-1., 1.) * math.sqrt(D * particle.dt)) / R) * (180 / math.pi)
    particle.lon += ((ParcelsRandom.uniform(-1., 1.) * math.sqrt(D * particle.dt)) / (R * math.cos(math.pi * particle.lat / 180))) * (180 / math.pi)

def deleteParticle(particle, fieldset, time):
    particle.OOB = 1
    particle.delete()
    
def submergeParticle(particle, fieldset, time):
    particle.depth = 1 
    
def KillFish(particle, fieldset, time):
  
  #only run for living larval fish
  if not particle.settleStage == 2 and particle.isDead == 0:
    
    #check if weʻre OOB
    if particle.OOB == 1:
      particle.isDead = 1
      particle.howDead = 1
  
    #check if we're over PLD
    if particle.settleStage == 3:
      particle.isDead = 1
      particle.howDead = 4
  
    #check if we've been beached for a long time (24 hrs)
    if particle.timeBeached >= (24*60*60) and particle.isDead == 0:
      particle.isDead = 1
      particle.howDead = 2
  
    #kill by random instantaneous mortality
    if particle.isDead == 0:
        if particle.complexStage == 0: #egg
            mort = fieldset.eggMort
        elif particle.complexStage > 0: #larvae
            mort = fieldset.larvaeMort
          
    if ParcelsRandom.random() <= mort:
        particle.isDead = 1
        particle.howDead = 3
  
    #if particle is dead, remove particle
    if particle.isDead:
      particle.delete() 

def SettleFish(particle, fieldset, time):
  if particle.isDead == 0 and particle.settleStage == 1: #if particle is living and currently settling
    
    settleTest = fieldset.settle[time, particle.depth, particle.lat, particle.lon]
  
    #test if we're somewhere we can settle
    #try out coin-flip settlement?
    #if (settleTest) > 0 and (ParcelsRandom.uniform(0., 1.) >= 0.5):
    particle.settleStage = 2 #if we do, mark as settled and note location
    particle.final_lon = particle.lon
    particle.final_lat = particle.lat
    particle.final_depth = particle.depth
    particle.settleday = time

def InitDepth(particle, fieldset, time):
  
    #print("init")
  
    #check - did we just hatch? if so randomize depth based on adult distribution
    #hacky but easiest way to do this
    if particle.age < 0.0083:
        biodepth = fieldset.maxAdultDepth
        phydepth = fieldset.seafloor[time, particle.depth, particle.lat, particle.lon]
        if phydepth < biodepth:
          maxdepth = phydepth
        else:
          maxdepth = biodepth
            
        #pick a random depth
        init_depth = ParcelsRandom.uniform(0., maxdepth)
        particle.depth = init_depth
     
def GrowFishComplex(particle, fieldset, time):
  
    #print("grow")
  
    #only run for living, not juvenile fish
    if particle.isDead == 0 and not particle.complexStage == 10:

        #AGING
        #advance age (in days)
        mydt = particle.dt
        age_check = particle.age
        particle.age += (mydt / (60*60*24))
        
        #GROWTH
        if particle.complexStage != 0: #if not an egg, grow
          particle.length = (fieldset.growth_m * particle.age) + fieldset.growth_c
    
        #DEVELOPMENT
        #check if we hatch (larval stage zero)
        if ((particle.complexStage == 0) and (particle.age >= fieldset.hatchTime)) or ((age_check < 0.0083) and (particle.spawnType == 1)):
          #egg hatches
          particle.complexStage = 1
          
          #what behaviors for this stage?
          particle.nowBuoyancy = fieldset.zeroBuoyancy
          particle.nowHorizontal = fieldset.zeroHorizontal
          particle.nowVertical = fieldset.zeroVertical
          particle.nowGrowthType = fieldset.zeroGrowthType
          particle.preferredDepth = ParcelsRandom.uniform(fieldset.zeroMinDepth, fieldset.zeroMaxDepth)
          particle.tempSpeedMin = fieldset.zeroSpeedMin
          particle.tempSpeedMax = fieldset.zeroSpeedMax
          particle.tempSpeed_b = fieldset.zeroSpeed_b
        
        #check if we develop (larval stage one)
        elif fieldset.devOneTime != -100 and particle.complexStage == 1 and particle.age >= fieldset.devOneTime:
          #larvae develops from stage 1 to stage 2
          particle.complexStage = 2
          
          #what behaviors for this stage?
          particle.nowBuoyancy = fieldset.oneBuoyancy
          particle.nowHorizontal = fieldset.oneHorizontal
          particle.nowVertical = fieldset.oneVertical
          particle.nowGrowthType = fieldset.oneGrowthType
          particle.preferredDepth = ParcelsRandom.uniform(fieldset.oneMinDepth, fieldset.oneMaxDepth)
          particle.tempSpeedMin = fieldset.oneSpeedMin
          particle.tempSpeedMax = fieldset.oneSpeedMax
          particle.tempSpeed_b = fieldset.oneSpeed_b
          
        #check if we develop (larval stage two)
        elif fieldset.devTwoTime != -100 and particle.complexStage == 2 and particle.age >= fieldset.devTwoTime:
          #larvae develops from stage 2 to stage 3
          particle.complexStage = 3
          
          #what behaviors for this stage?
          particle.nowBuoyancy = fieldset.twoBuoyancy
          particle.nowHorizontal = fieldset.twoHorizontal
          particle.nowVertical = fieldset.twoVertical
          particle.nowGrowthType = fieldset.twoGrowthType
          particle.preferredDepth = ParcelsRandom.uniform(fieldset.twoMinDepth, fieldset.twoMaxDepth)
          particle.tempSpeedMin = fieldset.twoSpeedMin
          particle.tempSpeedMax = fieldset.twoSpeedMax
          particle.tempSpeed_b = fieldset.twoSpeed_b
          
        #check if we develop (larval stage three)
        elif fieldset.devThreeTime != -100 and particle.complexStage == 3 and particle.age >= fieldset.devThreeTime:
          #larvae develops from stage 3 to stage 4
          particle.complexStage = 4

          #what behaviors for this stage?
          particle.nowBuoyancy = fieldset.threeBuoyancy
          particle.nowHorizontal = fieldset.threeHorizontal
          particle.nowVertical = fieldset.threeVertical
          particle.nowGrowthType = fieldset.threeGrowthType
          particle.preferredDepth = ParcelsRandom.uniform(fieldset.threeMinDepth, fieldset.threeMaxDepth)
          particle.tempSpeedMin = fieldset.threeSpeedMin
          particle.tempSpeedMax = fieldset.threeSpeedMax  
          particle.tempSpeed_b = fieldset.threeSpeed_b
        
        elif fieldset.devFourTime != -100 and particle.complexStage == 4 and particle.age >= fieldset.devFourTime:
          #larvae develops from stage 4 to stage 5
          particle.complexStage = 5
          
          #what behaviors for this stage?
          particle.nowBuoyancy = fieldset.fourBuoyancy
          particle.nowHorizontal = fieldset.fourHorizontal
          particle.nowVertical = fieldset.fourVertical
          particle.nowGrowthType = fieldset.fourGrowthType
          particle.preferredDepth = ParcelsRandom.uniform(fieldset.fourMinDepth, fieldset.fourMaxDepth)
          particle.tempSpeedMin = fieldset.fourSpeedMin
          particle.tempSpeedMax = fieldset.fourSpeedMax  
          particle.tempSpeed_b = fieldset.fourSpeed_b
        
        #check if we've settled
        if particle.settleStage == 2:
            particle.complexStage = 10
      
        #check if we need to update settleStage (past PLD)
        if particle.age >= fieldset.maxPLD and particle.settleStage != 2:
            particle.settleStage = 3
        elif particle.settleStage == 0 and fieldset.minPLD <= particle.age and particle.age <= fieldset.maxPLD:
            particle.settleStage = 1
      
        #advance approx distance travelled (Euclidian)
        lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
        lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
        particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
        
        if (particle.settleStage == 0) or (particle.settleStage == 1):
            particle.prev_lon = particle.lon  # Set the stored values for next iteration
            particle.prev_lat = particle.lat
            particle.prev_depth = particle.depth
            
        #set up stage-appropriate swimming speeds
        if (particle.nowHorizontal != 0) or (particle.nowVertical != 0): #if we swim
          if particle.nowGrowthType == 0: #age-dependent speed
            particle.nowSpeedMin = (particle.tempSpeedMin * particle.age) + particle.tempSpeed_b
            particle.nowSpeedMax = (particle.tempSpeedMax * particle.age) + particle.tempSpeed_b           
          elif particle.nowGrowthType == 1: #length-dependent speed
            particle.nowSpeedMin = (particle.tempSpeedMin * particle.length) + particle.tempSpeed_b
            particle.nowSpeedMax = (particle.tempSpeedMax * particle.length) + particle.tempSpeed_b
          elif particle.nowGrowthType == 2: #absolute speed
            particle.nowSpeedMin = particle.tempSpeedMin
            particle.nowSpeedMax = particle.tempSpeedMax           
            

