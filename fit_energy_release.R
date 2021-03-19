#fit energy release
library(firebehavioR)
library(Rothermel)
#### load fire models
data(SFM_metric, fuelMoisture)
fuelMoisture$fmLitter <- NULL
data (firexp)
fuelModels<-SFM_metric[14:53,]
############################################################################
## 1. fuel load
############################################################################
## get total load
fuelModels$total_load<-fuelModels$Load_1h+fuelModels$Load_10h+fuelModels$Load_100h+fuelModels$Load_Live_Herb+fuelModels$Load_Live_Woody
############################################################################
## 2. fuel structure
############################################################################
## calc weighted mean sav
weights_sav<-cbind(	fuelModels$Load_1h/fuelModels$total_load,
	fuelModels$Load_10h/fuelModels$total_load,
	fuelModels$Load_100h/fuelModels$total_load,
	fuelModels$Load_Live_Herb/fuelModels$total_load,
	fuelModels$Load_Live_Woody/fuelModels$total_load)

fuelModels$savmean<-apply(cbind(fuelModels[,7:11],weights_sav),1,FUN=function(x){weighted.mean(x=x[1:5],w=x[6:10])})
############################################################################
## 3. test
############################################################################

wind<-seq(0.1,10,0.1)
#test single model all moisture scenarios, single value wind
test1<-ros(modeltype=fuelModels[1,1],w=fuelModels[1,2:6],s=fuelModels[1,7:11],delta=fuelModels[1,12],
	mx.dead=fuelModels[1,13],h=fuelModels[1,14:18],m=fuelMoisture[,1:5],u=1,slope=0.0)
#test single model all moisture scenarios + 100 wind values	
test2<-mapply(ros,u=wind,MoreArgs = list(modeltype=fuelModels[1,1],w=fuelModels[1,2:6],s=fuelModels[1,7:11],delta=fuelModels[1,12],	mx.dead=fuelModels[1,13],h=fuelModels[1,14:18],m=fuelMoisture[,1:5],slope=0.0),SIMPLIFY = F)
############################################################################
## 4. generate rothermel predictions, all the models*moisture scenarios*100 wind values
############################################################################

rothermel_sims<-list()
#loop through the models
for(i in 1:40){
	rothermel_sims[[i]]<-mapply(ros,u=wind,MoreArgs = list(modeltype=fuelModels[i,1],w=fuelModels[i,2:6],s=fuelModels[i,7:11],delta=fuelModels[i,12],mx.dead=fuelModels[i,13],h=fuelModels[i,14:18],m=fuelMoisture[,1:5],slope=0.0),SIMPLIFY = F)
}
rothermel.df<-lapply(rothermel_sims,FUN=function(x){do.call(rbind,lapply(x,as.data.frame))})
rothermel.df<-do.call(rbind,rothermel.df)
############################################################################
## 4. generate simplified inputs
############################################################################
moisture<-list()
df<-list()
for (i in 1:40){
	moisture[[i]]<-apply(fuelMoisture[,1:5],1,weighted.mean,w=weights_sav[i,])
	df[[i]]<-expand.grid(m=moisture[[i]],u=wind)
	df[[i]]$total_load<-rep(fuelModels$total_load[i],length(df[[i]][,1]))
	df[[i]]$sav_mean<-rep(fuelModels$savmean[i],length(df[[i]][,1]))
}
df<-do.call(rbind,df)
df$H_src<-rothermel.df$Heat.source..kW.m2.
############################################################################
## 4. fit logistic eqn 1st approx, assume 600 kW/m2 as max heat source
############################################################################
df$rsp_dummy<-log(df$H_src/(600-df$H_src))
df_lm<-subset(df,!is.infinite(df$rsp_dummy) & !is.na(df$rsp_dummy))
mod1_glm<-lm(rsp_dummy~total_load+sav_mean+m+u,data=df_lm)
coeff_mod1_glm<-coefficients(mod1_glm)
#############################################################################
#### 5. fit non linear eqn
#############################################################################
mod2_nls<-nls(H_src~600/(1+exp(k1+k2*total_load+k3*sav_mean+k4*m+k5*u)),
			data=df,
start = list(k1=coeff_mod1_glm[1],k2=coeff_mod1_glm[2],k3=coeff_mod1_glm[3],
			k4=coeff_mod1_glm[4],k5=coeff_mod1_glm[5])
	)

coeff_mod2_nls<-coefficients(mod2_nls)
