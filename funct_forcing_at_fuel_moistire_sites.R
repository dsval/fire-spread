#funct forcing at fuel moistire sites
library(xts)
getForcing_fuelsites<-function(site,lat,lon,moisture_samples){
	#site=sites$Site[59];lat=sites$lat[59];lon=sites$lon[59]
	df<-subset(moisture_samples,moisture_samples$Site==site)
	df$Date<-as.Date(df$Date)
	df<-df[order(df$Date),]
	yend<-as.numeric(format(df$Date[length(df$Date)],'%Y'))
	ystart<-as.numeric(format(df$Date[1],'%Y'))
	if(yend>2020){
		yend<-2020
	}
	
	#download daymet data
	daymet_data<- try(daymetr::download_daymet("Daymet",
		lat = lat,
		lon = lon,
		start = ystart,
		end = yend,
		internal = TRUE))
	if(class(daymet_data)=="try-error"){
		df<-data.frame( yearmon= format(df$Date,'%Y-%m'),   Date=df$Date , X= df$X,GACC=df$GACC, 
					State =df$State, Group=df$Group , Site=df$Site ,Fuel= df$Fuel , 
					Percent=df$Percent,   P.mm.d= rep(NA,length(df$Percent)),   
			swin.Wm2.d=rep(NA,length(df$Percent)), swe.mm.d= rep(NA,length(df$Percent)),  tmax.C.d=rep(NA,length(df$Percent)),  tmin.C.d=rep(NA,length(df$Percent)),   vp.Pa.d=rep(NA,length(df$Percent)) ,  
			photop.h.d=rep(NA,length(df$Percent)), tmean.C.d=rep(NA,length(df$Percent)),  P.mm.m= rep(NA,length(df$Percent)),    swin.Wm2.m=rep(NA,length(df$Percent)), swe.mm.m=rep(NA,length(df$Percent)), 
			tmax.C.m=rep(NA,length(df$Percent)), tmin.C.m=rep(NA,length(df$Percent)),   vp.Pa.m=rep(NA,length(df$Percent)),    photop.h.m=rep(NA,length(df$Percent)), tmean.C.m=rep(NA,length(df$Percent)))
	}else{
		daymet_data[[7]]$dates<-as.Date(paste0(daymet_data[[7]]$year,'-',daymet_data[[7]]$yday),format='%Y-%j')
		#get daily forcing
		forcing<-subset(daymet_data[[7]],daymet_data[[7]]$dates %in% df$Date)
		forcing$srad..W.m.2.<-(forcing$srad..W.m.2.*forcing$dayl..s.)/86400
		forcing$photop<-(forcing$dayl..s./86400)*24
		forcing$tmean<-(forcing$tmax..deg.c.+forcing$tmin..deg.c.)/2
		forcing<-forcing[,4:12]
		names(forcing)<-c('P.mm.d','swin.Wm2.d','swe.mm.d',
			'tmax.C.d','tmin.C.d','vp.Pa.d','dates','photop.h.d','tmean.C.d')
		df<-merge(df,forcing,by.x='Date',by.y='dates')
		#get monthly forcing
		forcing_m<-xts(daymet_data[[7]][,1:9],daymet_data[[7]]$dates)
		forcing_m$srad..W.m.2.<-(forcing_m$srad..W.m.2.*forcing_m$dayl..s.)/86400
		forcing_m$photop<-(forcing_m$dayl..s./86400)*24
		forcing_m$tmean<-(forcing_m$tmax..deg.c.+forcing_m$tmin..deg.c.)/2
		prec_m<-apply.monthly(forcing_m$prcp..mm.day.,sum,na.rm=T)
		forcing_m<-apply.monthly(forcing_m,mean,na.rm=T)
		forcing_m$prcp..mm.day.<-prec_m
		forcing_m<-as.data.frame(forcing_m[,4:11])
		names(forcing_m)<-c('P.mm.m','swin.Wm2.m','swe.mm.m',
			'tmax.C.m','tmin.C.m','vp.Pa.m','photop.h.m','tmean.C.m')
		forcing_m$yearmon<-format(time(prec_m),'%Y-%m')
		df$yearmon<-format(df$Date,'%Y-%m')		
		df<-merge(df,forcing_m,by='yearmon')
	}	
		
		
	
	df
}
