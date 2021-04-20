#Ros Balbi
rosB<-function(Ta,W,sav,d,m,u,elev,slop){
	# Ta: air temperature (C)
	# W: total fuel load (Kg/m2)
	# sav: Surface area to volume ratio (1/m)
	# d: fuel bed depth (m)
	# m: fuel moisture (%Drymass)
	# u: wind speed (m/s)
	# v: wind direction (degrees)
	# slop: topo slope (degrees)
	# asp: topo asp (degrees)
	# Ta=30;W=0.7;sav=4000;d=0.7;m=30;u=1;v=30;slop=5;elev=0;asp=30
	###############################################################################################
	# 02.define the constants
	###############################################################################################
	kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
	kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
	kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
	kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
	kSecInDay <- 86400  # number of seconds in a day
	kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
	kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
	kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
	kB <- 5.670374e-8   # Stefan-Boltzmann constant https://physics.nist.gov/cgi-bin/cuu/Value?sigma
	kH<- 1.74e7 	     # fuel heat of combustion of the pyrolysis gases, J/kg(Awad et al., 2020)
	chi_o<- 0.3         #radiant heat fraction (Awad et al., 2020)
	st<- 17              # stoichiometric coefficient, 17 (Balbi et al., 2020); 9 (Awad et al., 2020)
	tau_o <- 75591      # flame residence time parameter, s/m (Balbi et al., 2020)
	n<- 4 # optical depth parameter
	dp <- 512 # fuel particle density 398(Balbi et al., 2020), kg/m3  ;512.7915(Scott and Burgan, 2005)
	Cpf <- 1172 # fuel specific heat capacity J·kg–1·K–1 Nelson, 2000
	Cpw <- 4180 # water specific Heat (Balbi et al., 2020) 
	Ti<- 600 #ignition temperature ,K  600(Balbi et al., 2020) 673 K (Awad et al., 2020) 
	Kd <-130 # drag coefficient s/m (Balbi et al., 2020) 
	r00 <- 2.5e-5 # model coefficient (Balbi et al., 2020) 
	#ef<- 0.9# flame emmisivity, Fons, 1946;Rossi et al., 2010
	#ef<- 0.227 # flame emmisivity, Fons, 1946;Rossi et al., 2010
	#theta_f<- 20*pi/180 # flame angle 20 deg
	Tvap<-373 # vaporization temperature(Balbi et al., 2020) 
	###############################################################################################
	# 01. Define some functions
	###############################################################################################
	# ************************************************************************
	# Name:     specific_heat
	# Inputs:   double (tc), air temperature, degrees C
	# Returns:  double, specific heat of moist air, J/kg/K
	# Features: This function calculates the spefic heat of moist air
	# Ref:      Tsilingris (2008), Thermophysical and transport properties of
	#           humid air at temperature range between 0 and 100 Â°C, Energy
	#           Conversion and Management, vol. 49, pp. 1098--1110.
	# ************************************************************************
	specific_heat <- function(tc) {
		
		tc<-ifelse(tc < 0,0,tc)
		cp <- 1.0045714270 +
		(2.050632750e-3)*tc -
		(1.631537093e-4)*tc*tc +
		(6.212300300e-6)*tc*tc*tc -
		(8.830478888e-8)*tc*tc*tc*tc +
		(5.071307038e-10)*tc*tc*tc*tc*tc
		cp <- (1e3)*cp
		
		return(cp)
	}
	# ************************************************************************
	# Name:     enthalpy_vap
	# Inputs:   double (tc), air temperature, degrees C
	# Returns:  double, J/kg
	# Features: This function calculates the temperature-dependent enthalpy
	#           of vaporization (latent heat of vaporization)
	# Ref:      Eq. 8, Henderson-Sellers (1984), A new formula for latent heat
	#             of vaporization of water as a function of temperature, Quarterly
	#             Journal of the Royal Meteorological Society, vol. 110, pp. 1186--
	#             1190.
	# ************************************************************************
	enthalpy_vap <- function(tc) {
		1.91846e6*((tc + 273.15)/(tc + 273.15 - 33.91))^2
	}
	elev2pres <- function(z) {
		kPo*(1 - kL*z/kTo)^(kG*kMa/(kR*kL))
	}
	p_dryair<-function(elev,t){
		#calc dry air density [kg/m3] using pv=nRT
		(elev2pres(elev)*kMa)/(kR*(t+273.15))
	}
	
	###############################################################################################
	# 1. calc energy required for ignition
	###############################################################################################
	Tak<-Ta+273.15
	# moisture to fraction
	m<-m/100
	#air specific heat (J/kg/K)
	Cpa<-specific_heat(Ta)
	#mean packing ratio for fuel complex
	beta<-(W/d)/dp
	#air density
	da<-p_dryair(elev,Ta)
	#LAI fuel
	lai<-beta*sav*d/2
	# enthalpy vap (J/kg)
	h<-enthalpy_vap(Ta)
	#energy req for ignition (J/kg)
	q<-Cpf*(Ti-Tak) + m*(h+Cpw*(Tvap-Tak))
	#Flame temperature (K)
	Tf<-Tak + ((kH*(1-chi_o))/((st+1)*Cpa))
	#u0	Upward gas velocity at flame body mid-height on flat terrain (m/s)
	u_o<-2*((st+1)/tau_o)*(Tf/Tak)*(dp/da)*pmin(lai,2*pi)
	# flame height
	h_fl<-u_o^2/(kG*((Tf/Tak)-1))
	###############################################################################################
	# 2. calc  ROS component for radiative heat (m/s)
	###############################################################################################
	# 
	Rb<-pmin(lai/pi,1)*(kB*Tf^4/(beta*dp*q))
	###############################################################################################
	# 3.  calc  ROS component for convective heat (m/s)
	###############################################################################################
	Rc<-(sav*kH/(q*tau_o))*pmin(d,(2*pi)/(sav*beta))*(((d*tan(slop*pi/180))/(2*d+h_fl))+(u*exp(-1*Kd*sqrt(beta)*Rb)/u_o))
	###############################################################################################
	# 3.  calc  ROS component for flame radiation heat (m/s)
	###############################################################################################
	#calc radiative factor
	A<-pmin(lai/2*pi,1) * (chi_o*kH)/(4*q)
	# calc flame tilt angle (radians)
	tilt_fl<-atan(tan(slop*pi/180)+(u/u_o))
	Rr<-A*Rb*((1+sin(tilt_fl)-cos(tilt_fl))/(1+(Rb*cos(tilt_fl)/(sav*r00))))
	###############################################################################################
	# 3.  calc  ROS first approximation (m/s)
	###############################################################################################
	R<-Rb+Rc+Rr
	Rc<-(sav*kH/(q*tau_o))*pmin(d,(2*pi)/(sav*beta))*(((d*tan(slop*pi/180))/(2*d+h_fl))+(u*exp(-1*Kd*sqrt(beta)*R)/u_o))
	Rr<-A*R*((1+sin(tilt_fl)-cos(tilt_fl))/(1+(R*cos(tilt_fl)/(sav*r00))))
	###############################################################################################
	# 4.  iterate
	###############################################################################################
	for(i in 1:100){
		R<-Rb+Rc+Rr
		Rc<-(sav*kH/(q*tau_o))*pmin(d,(2*pi)/(sav*beta))*(((d*tan(slop*pi/180))/(2*d+h_fl))+(u*exp(-1*Kd*sqrt(beta)*R)/u_o))
		Rr<-A*R*((1+sin(tilt_fl)-cos(tilt_fl))/(1+(R*cos(tilt_fl)/(sav*r00))))
	}
	R
	
	
}
