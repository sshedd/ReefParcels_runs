import math

def random_von_Mises(mu, kappa, thetaLast):
    
    """
        from https://dlwhittenbury.github.io/ds-1-sampling-and-visualising-the-von-mises-distribution.html
        modified to work in base Python
        rand_von_Mises(N,mu,kappa)
        ==========================
    
        Generates theta an Nx1 array of samples of a von Mises distribution
        with mean direction mu and concentration kappa.
    
        INPUT:
        
            * mu - mean direction (float)
            * kappa - concentration (float)
            * thetaLast - previous angle (float)
    
        OUTPUT:
    
            * theta - an Nx1 array of samples of a von Mises distribution with mean
            direction mu and concentration kappa.
    
         References:
         ===========
    
         Algorithm first given in
    
         [1] D. J. Best and N. I. Fisher, Efficient Simulation of the von Mises
         Distribution, Applied Statistics, 28, 2, 152--157, (1979).
    
         Also given in the following textbook/monograph
    
         [2] N. I. Fisher, Statistical analysis of circular data, Cambridge University Press, (1993).
    
    """
    
    # Checks
    # =======
    
    #  mu should be a real scalar. It can wrap around the circle, so it can be negative, positive and also
    #  outside the range [0,2*pi].
    if (type(mu) is not float) and (type(mu) is not int):
        raise TypeError("mu must be a real scalar number.")
    
    # kappa should be positive real scalar
    if (type(kappa) is not float) and (type(kappa) is not int):
        raise TypeError("kappa must be a positive float.")
    if kappa < 0:
        raise Exception("kappa must be a positive float.")
    
    #  SPECIAL CASE
    # ==============
    
    #  As kappa -> 0 one obtains the uniform distribution on the circle
    #float_epsilon = np.finfo(float).eps
    #if kappa <= float_epsilon:
    #    theta = 2.0 * math.pi * ParcelsRandom.random.rand(N,1) # [0,1] -> [0,2*pi]
    #    return theta
    
    # MAIN BODY OF ALGORITHM
    # =======================
    
    # Used same notation as Ref.~[2], p49
    
    a = 1.0 + math.sqrt(1.0 + 4.0 * kappa**2)
    b = (a - math.sqrt(2.0 * a)) / (2.0 * kappa)
    r = (1.0 + b**2) / (2.0 * b)
    
    # Pseudo-random numbers sampled from a uniform distribution [0,1]
    done = False
    
    while done == False:
    
      U1 = ParcelsRandom.uniform(0., 1.)
      U2 = ParcelsRandom.uniform(0., 1.)
      U3 = ParcelsRandom.uniform(0., 1.)
      
      z = math.cos(math.pi * U1)
      f = (1.0 + r * z) / (r + z)
      c = kappa * (r - f)
      
      #theta = 0
      
      if (U3 - 0.5) >= 0:
        sign = 1
      else:
        sign = -1
      
      if (((c * (2.0 - c) - U2) > 0.0)  or ((math.log(c/U2) + 1.0 - c) > 0.0)):
      
          done = True
          theta = (sign * math.acos(f) + mu) % (2*math.pi)     
          thetaNow = theta + thetaLast

    return thetaNow
