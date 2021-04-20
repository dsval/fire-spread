#physics ROS

rosPhys<-function(Ta,W,sav,d,m,u,v,slop,asp){
	# Ta: air temperature (C)
	# W: total fuel load (Kg/m2)
	# sav: Surface area to volume ratio (1/m)
	# d: fuel bed depth (m)
	# m: fuel moisture (%Drymass)
	# u: wind speed (m/s)
	# v: wind direction (degrees)
	# slop: topo slope (degrees)
	# asp: topo asp (degrees)
	# Ta=30;W=0.7;sav=4000;d=0.7;m=30;u=1;v=30;slop=5;asp=30
	###############################################################################################
	# 02.define the constants
	###############################################################################################
	#
	
	kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
	kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
	kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
	kSecInDay <- 86400  # number of seconds in a day
	kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
	kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
	kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
	kB <- 5.670374e-8   # Stefan-Boltzmann constant https://physics.nist.gov/cgi-bin/cuu/Value?sigma
	kH<- 17.4e3 	     # fuel heat of combustion of the pyrolysis gases, kJ/kg(Awad et al., 2020)
	chi_o<- 0.3         #radiant heat fraction (Awad et al., 2020)
	st<- 9              # stoichiometric coefficient
	tau_o <- 75591      # flame residence time parameter, s/m
	n<- 4 # optical depth parameter
	dp <- 16.01876 # fuel particle density, kg/m3 (Scott and Burgan, 2005)
	Cpf <- 1172 # fuel specific heat capacity J·kg–1·K–1 Nelson, 2000
	Ti<- 673 #ignition temperature ,  K (Awad et al., 2020)
	ef<- 0.9# flame emmisivity, Fons, 1946;Rossi et al., 2010
	#ef<- 0.227 # flame emmisivity, Fons, 1946;Rossi et al., 2010
	theta_f<- 20*pi/180 # flame angle 20 deg
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
	p_dryair<-function(elev,t){
		#calc dry air density [kg/m3] using pv=nRT
		kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
		kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
		(elev2pres(elev)*kMa)/(kR*(t+273.15))
	}
	###############################################################################################
	# 1. calc extinction moisture (% of dry mass)
	###############################################################################################
	# temp to kelvin
	Tak<-Ta+273.15
	# moisture to fraction
	m<-m/100
	#angles to rad
	#v=v*pi/180;slop=slop*pi/180;asp=asp*pi/180
	#air specific heat (J/kg/K)
	Cpa<-specific_heat(Ta)
	#mean packing ratio for fuel complex
	beta<-(W/d)/dp
	# enthalpy vap (J/kg)
	h<-enthalpy_vap(Ta)
	#Flame temperature (K)
	Tf<-Tak + ((kH*1000*(1-chi_o))/((st+1)*Cpa))
	# flame residence time
	#tau<-tau_o/sav
	#LAI fuel
	lai<-beta*sav*d/2
	#air cooling effect
	ab<-1-(0.47/(2*lai)^(1/2))
	#extinction moisture (% of dry fuel)(Awad et al., 2020)
	Mx<-((tau_o*kB*Tf^4*ab)/(n*dp*h))-((Cpf*(Ti-Tak))/h)
	###############################################################################################
	# 2. calc wind and slope factor
	###############################################################################################
	# McRae wind factor (Sharples, 2008)
	wf<-(1-tan(theta_f)*tan(slop)*cos(v-asp))^-1
	##wind coefficient rothermel
	# to ft
	sigma.tot<-sav/3.281
	C=7.47*exp(-0.133*sigma.tot^.55)
	B=0.02526*sigma.tot^.54
	E=0.715*exp(-3.59*10^(-4)*sigma.tot)
	#optimum packing ratio
	beta.op=3.348*sigma.tot^(-0.8189)
	rpr=beta/beta.op #relative packing ratio
	fw=C*(u*54.6806649)^B*rpr^(-E)
	#slope coefficient
	fs=5.275*beta^(-0.3)*(slop/100)^2
	###############################################################################################
	# 3. calc rate of spread (m/s)
	###############################################################################################
	heat_source<-ef*kB*Tf^4*(1+fw+fs)*ab
	ros<-(1/(beta*dp))*(ef*kB*Tf^4*ab/(Cpf*(Ti-Tak)+m*h))*(1+fw+fs)
	list(ros=ros,Mx=Mx,HS=heat_source)
	
}
