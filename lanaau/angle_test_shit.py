import math
import numpy as np

C=3.5 # current distance 
B = 4 #max sensing distance
thetaGoal = math.pi #angle to reef
lastTheta = 0
nowTheta = 0
kappa = 50

while C > 0.0:
  lastTheta = nowTheta #last angle
  thetaGoal = math.pi #angle to reef
  
  #check to see if we're beyond detection distance, set no sensing ability
  if C > B:
    sensingAbility = 0
  else:
    sensingAbility = 1 - (C / B)
  
  if C <= 2 and thetaGoal == 0: 
    turningAngle = lastTheta #avoid turning 180 when already within habitat
    phi = 0
  else:
    phi = math.atan2(math.sin(thetaGoal - lastTheta), math.cos(thetaGoal - lastTheta))
    turningAngle = -sensingAbility * phi #double check whether sensingAbility positive or negative
    #turningAngle = (-sensingAbility) * (particle.lastTheta - thetaGoal)
  
  mu = turningAngle
  
  a = 1.0 + math.sqrt(1.0 + 4.0 * kappa**2)
  b = (a - math.sqrt(2.0 * a)) / (2.0 * kappa)
  r = (1.0 + b**2) / (2.0 * b)
  
  # Pseudo-random numbers sampled from a uniform distribution [0,1]
  done = False
  
  while done == False:
  
    rand1 = np.random.rand()
    rand2 = np.random.rand()
    rand3 = np.random.rand()
    
    z = math.cos(math.pi * rand1)
    f = (1.0 + r * z) / (r + z)
    c = kappa * (r - f)
    
    if (rand3 - 0.5) >= 0:
      sign = 1
    else:
      sign = -1
    
    if (((c * (2.0 - c) - rand2) > 0.0) or ((math.log(c/rand2) + 1.0 - c) > 0.0)):
    
        done = True
        theta = math.fmod(sign * math.acos(f) + (mu - lastTheta), (2*math.pi))   
        nowTheta = theta # *****triple***** CHECK
  
  print(theta, phi, mu, thetaGoal,lastTheta,sensingAbility)
  C = C - 0.25






#STOP HERE

lt = particle.lastTheta
print("phi = ")
print(phi)
print("mu = ")
print(mu)
print("theta_goal = ")
print(thetaGoal)
print("theta_last = ")
print(lt)
print("nowTheta = ")
print(nowTheta)
print("d_t = ")
print(sensingAbility)
print("")

#get u/v and move fish accordingly
swimmingSpeed = ParcelsRandom.uniform(particle.nowSpeedMin, particle.nowSpeedMax)

#convert speed from m/s to decimal degrees/s (approximate)
R = 6378137
latSwimmingSpeed = (swimmingSpeed / R) * (180 / math.pi)
lonSwimmingSpeed = (swimmingSpeed / (R * math.cos(math.pi * particle.lat / 180)))*(180 / math.pi) 

#if move puts particle on land, don't make it
uOrient = lonSwimmingSpeed * math.cos(nowTheta) #speed in decimal degrees, theta in radians
vOrient = latSwimmingSpeed * math.sin(nowTheta) #speed in decimal degrees, theta in radians
test_lon = particle.lon + (uOrient * particle.dt)
test_lat = particle.lat + (vOrient * particle.dt)
(u_test, v_test, w_test) = fieldset.UVW[time, particle.depth, test_lat, test_lon]

if u_test != 0 and v_test != 0:
  particle.lon += (uOrient * particle.dt)
  particle.lat += (vOrient * particle.dt)

#save value
particle.thisTheta = nowTheta
