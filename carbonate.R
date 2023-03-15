carbonate = function(TEMP,alk,DIC,pH=8.) {
  
  # Carbonate equilibrium system based on 
  # P. Tans, "Why Carbon Dioxide from Fossil Fuel Burning Won't Go Away" 
  #  In: MacAladay, J. (ed),  "Environmental Chemistry," 1996. pp. 
  
  # Note that we are ignoring contributions of PO4 and SiO4 to alkalinity, which may be 
  # important regionally but probably lead to errors in pCO2 < 1 ppmv globally
  # See Follows et al (2006) for a discussion of this approximation
  
  # Default values if the user doesn't specify at runtime are for preindustrial global surface
  
  # Tans defines pCO2 via Henry's Law in atmospheres
  # I'm converting to ppmv for clarity 
  # 
  # Variable definitions
  
  # INPUT VARIABLES
  #       T      #  Temperature (Celsius)
  #       TA     #  titration alkalinity (equivalents/kg)
  #       DIC    #  total dissolved inorganic carbon (mol/kg)
  
  # OUTPUT VARIABLES
  #       pH     #  pH of the water (-log(H))
  #       pCO2    #  partial pressure of CO2 (Pascals)
  #       iter   iteration count
  
  # CHEMICAL & THERMODYNAMIC COEFFICIENTS  
  #       K0     #  Henry's Law constant
  #       K1     #  first dissociation coefficient for H2CO3
  #       K2     #  second dissociation coefficient for H2CO3
  #       Kb     #  dissociation constant for boric acid
  #       S      #  Salinity (g/kg)
  #       Boron  #  total boron (mol/kg)
  #       H0     #  "old" concentration of hydrogen ion (mol/kg)
  
  # INTERNAL VARIABLES USED TO CALCULATE pCO2 and pH
  #       H      #  current concentration of hydrogen ion (mol/kg)
  #       CA     #  carbonate alkalinity (eq / l)
  #       CO2aq  #  concentration of aqueous CO2
  #       diff.H  #  difference in successive estimates of H+ (mol/kg)
  #       tiny.diff.H   (1e-15) test for convergence 
  #       a      #  first term in quadratic for eq 12
  #       b      #  second term in quadratic for eq 12
  #       c      # third term in quadratic for eq 12
  
  # Convert input units to mks
  TEMP <- TEMP + 273.15 # temperature from Celsius to Kelvin
  alk <- alk * 1.e-6    # microequivalents to equivalents
  DIC <- DIC * 1.e-6    # micromoles to moles
  
  # Set values of prescribed constants
  S <- 34.78                 # Salinity in ppt
  Boron <- 1.179e-5 * S      # Total Boron mole/kg as a fraction of salinity
  
  # Carbonate and boric acid equilibrium constants as functions of temp and S
  K0 <- exp(-60.2409 + 9345.17/TEMP + 23.3585*log(TEMP/100) 
            + S * (0.023517 - 0.00023656*TEMP +0.0047036*(TEMP/100)^2) )
  
  K1 <- exp(2.18867 - 2275.036/TEMP - 1.468591 * log(TEMP) 
            + (-0.138681 - 9.33291/TEMP) * sqrt(S) + 0.0726483*S    
            - 0.00574938 * S ^1.5)
  
  K2 <- exp(-0.84226 - 3741.1288/TEMP -1.437139 * log(TEMP)
            + (-0.128417 - 24.41239/TEMP)*sqrt(S) + 0.1195308 * S   
            - 0.0091284 * S ^1.5 )
  
  Kb <- exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S 
              + 1.726 * S ^1.5 - 0.0993*S^2) / TEMP                
             + (148.0248 + 137.194 * sqrt(S) + 1.62247 * S)       
             + (-24.4344 - 25.085 * sqrt(S) - 0.2474 * S) * log(TEMP)
             + 0.053105 * sqrt(S) * TEMP)      
  
  # Iterate for H and CA by repeated solution of eqs 13 and 12
  H <- 10.^(-pH)                 # initial guess from arg list      
  diff.H <- H     
  tiny.diff.H <- 1.e-15 
  
  iter <- 0
  
  while (diff.H > tiny.diff.H) {     # iterate until H converges
    
    H.old <- H                      # remember old value of H
    
    # solve Tans' equation 13 for carbonate alkalinity from TA
    CA <- alk - (Kb/(Kb+H)) * Boron      
    
    # solve quadratic for H (Tans' equation 12)
    a <- CA
    b <- K1 * (CA - DIC)
    c <- K1 * K2 * (CA - 2 * DIC)
    H <- (-b + sqrt(b^2 - 4. * a * c) ) / (2. * a)  
    
    # How different is new estimate from previous one?
    diff.H <- abs(H - H.old)
    iter <- iter + 1
    
  }
  
  # Now solve for CO2 from equation 11 and pCO2 from eq 4
  CO2aq <- CA / (K1/H + 2*K1*K2/H^2)  # Eq 11
  pCO2 <- CO2aq / K0 * 1.e6           # Eq 4 (converted to ppmv)
  pH <- -log10(H)
  
  return(list(TEMP=TEMP-273.15, alk=alk*1.e6, DIC=DIC*1.e6, pH=pH, pCO2=pCO2, iter=iter))      
}
