library(tidyverse)
library(lubridate)
library(knitr)
opts_chunk$set(echo=F,warning=F,message=F,fig.width = 11,fig.height = 5,cache=F)
library(boot)
library(raster)
library(tidyr)
library(ggplot2)
# MASS This fucks with dplyr, load later
library(ggmap)
library(dplyr)
library(chron)
library(gridExtra)
library(stringr)
library(R2jags)
library(maptools)
library(reshape2)
#New model to be run, flag turned off if just updating.
newModel<-T

migrators<-read.csv("./All_animals_cleaned.csv")
migrators$Ptt<-as.factor(migrators$Ptt)
summary(migrators)
#we have to leave out the tags in 2016 that did 24 on, 24 off
migrators<-migrators%>%filter(Ptt%in%c("112699", "121207", "121208", "121210", "123224", "123232", "123236", "131130", "131132", "166122", "166123", "166125", "166126", "166128"))
#test that it worked
mig<-migrators%>%
  group_by(Ptt)%>%
  filter(row_number()==1)
mig


migrators$timestamp<-as.character(migrators$timestamp)
migrators$timestamp<-as.POSIXct(migrators$timestamp,format="%Y-%m-%d %H:%M:%S",tz="GMT")

#first we will check if there are duplicate timestamps for any individuals
time<-migrators %>% group_by(Ptt,timestamp) %>% summarize(n=n()) %>% filter(n>1)
#looks like no. 
# Now lets test if there are any duplicate locations for any of the animals
location<-migrators%>%group_by(Ptt, Latitude, Longitude)%>%summarize(n=n())%>%filter(n>1)
#but thats ok we keep it
ungroup(migrators)

#remove empty timestamps
migrators<-migrators[!is.na(migrators$timestamp),]#there aren't any

migrators<-migrators%>%filter(Quality!="Z")
mig<-droplevels(migrators)

mig2<-mig

mig2$Date<-ymd(mig2$Date)

f_l<-mig2%>%group_by(Ptt)%>%filter(row_number()==1 | row_number()==n())




#standardize column names to match the simulation
#Create an animal tag.

mig2<-mig2%>%dplyr::select(1,3,5,7:9,21:24)
names(mig2) <- c("X", "Animal", "Date", "argos.lc", "y", "x", "timestamp", "Month", "Day", "Year" )
mig3<-mig2
#rearrange
mig3<-mig3%>%dplyr::select("X", "Animal", "argos.lc", "x", "y", "timestamp", "Month", "Year", "Date")

#holding mig 4 for later
mig4<-mig3


mig3$Month<-factor(mig3$Month,levels=month.name)

##Time is the beginning of the first point.
step_length=12

dumb1<-split(mig3,mig3$Animal)

#time diff function
timed<-function(d,step_length){
  d$j[1]<-0
  for (x in 2:nrow(d)){
    d$j[x]<-as.numeric(difftime(as.POSIXct(d$timestamp[x]),as.POSIXct(d$timestamp[x-1]),units="mins"))/(step_length*60)
  }
  
  #Split out track endings
  ends<-c(1,which(d$j>1),nrow(d))
  
  for(w in 2:length(ends)){
    d[ends[w-1]:ends[w],"Track"]<-w-1
  }
  
  #remove tracks that are shorter than 24 hrs
  track_time<-d %>% group_by(Track) %>% summarize(mt=difftime(max(as.POSIXct(timestamp)),min(as.POSIXct(timestamp)),units="hours")) %>% filter(mt>=24) %>% .$Track
  
  d<-d[d$Track %in% track_time,]
  
  #renumber the tracks
  d$Track<-as.numeric(as.factor(d$Track))
  return(d)
}

dumb2<-lapply(dumb1,timed,step_length=step_length)

#Format matrices for jags
mig3<-bind_rows(dumb2)

######recode whales
#mig3$Animal<-as.numeric(as.factor(mig3$Animal))

dumb3<-split(mig3,list(mig3$Animal,mig3$Track),drop=TRUE)

dumb4<-lapply(dumb3,function(x){
  #How many observations in each step length segment
  x$step<-as.numeric(cut(as.POSIXct(x$timestamp),paste(step_length,"hours")))
  return(x)
})

mig3<-bind_rows(dumb4)

#refactor animal
mig3$Ptt<-mig3$Animal
mig3$Animal<-as.numeric(as.factor(mig3$Animal))


mig3<-mig3 %>% group_by(Animal,Track,step) %>% mutate(jStep=1:n())

mig4<-mig3%>%mutate(j=ifelse(step == 1 & jStep == 1,0,j))  
splitmig<-split(mig4,list(mig4$Animal,mig4$Track),drop=TRUE)

splitmigsum<-lapply(splitmig,function(x){
  x$sumj<-cumsum(x$j)
  return(x)
})

mig5<-bind_rows(splitmigsum)

revtrunc <- function(x) { sign(x) * (x - floor(x)) } 
mig5$basej<-revtrunc(mig5$sumj)                    

mig5<-mig5%>%mutate(step_new=floor(sumj)+1)

#gonna have to redo jstep too
mig5<-mig5 %>% group_by(Animal,Track,step_new) %>% mutate(jStep_new=1:n())
ungroup(mig5)
mig6<-mig5%>%dplyr::select(Ptt, Animal, argos.lc, x, y, timestamp, Month, Year, Date, Track, step_new, jStep_new, basej)

names(mig6)[names(mig6) == 'step_new'] <- 'step'
names(mig6)[names(mig6)== 'jStep_new']<-'jStep'
names(mig6)[names(mig6)== 'basej']<-'j'
head(mig6)


mig6$Month<-factor(mig6$Month,levels=month.name)

#total number of steps per track/animal
steps_all<-mig6 %>% group_by(Animal,Track) %>% summarize(n=length(unique(step)))


#Cast time array
j<-acast(mig6,Animal~Track~step~jStep,value.var="j")

#how many observations per individual in each step. Basically the j step. 
mig6$step<-factor(mig6$step,levels=1:max(steps_all$n))
idx<-melt(table(mig6$Animal,mig6$Track,mig6$step))
colnames(idx)<-c("Animal","Track","step","jStep")
idx<-acast(data=idx,Animal~Track~step) # what this is saying is that
#we get 225 matrices, 1 for each step, rownames are animal number, column names are track #, values are the number of observations in that animal and track at that step

#month array
mig6$MonthF<-as.numeric(factor(mig6$Month,levels=month.name))

MonthA<-acast(mig6,Animal~Track~step,value.var="MonthF",fun.aggregate = min)
MonthA[!is.finite(MonthA)]<-NA

#Individuals
ind=length(unique(mig6$Animal))

#tracks per indivudal
tracks<-mig6 %>% group_by(Animal) %>% summarize(tracks=length(unique(Track))) %>% .$tracks

#steps per track
#resulting matrix - rows are Animal ID, column names are Track number, each spot is number of steps
steps<-acast(steps_all,Animal~Track,value.var="n")

#index array for merging
mig6$Index<-1:nrow(mig6)

#obs array
obs<-melt(mig6,measure.vars=c("x","y"))
obs<-acast(obs,Animal~Track~step~jStep~variable)

MonthIndex<-acast(mig6,Animal~Track~step,value.var="Month",fun=function(x) names(which.max(table(x))) )

#argos error class array
mig6$argos.lc<-factor(mig6$argos.lc,levels=c(3,2,1,0,"A","B"))
mig6$numargos<-as.numeric(mig6$argos.lc)
obs_class<-acast(mig6,Animal~Track~step~jStep,value.var="numargos")
write.csv(mig6,"./mig6noarrangedModel12_416.csv")


library(MASS)
#source jags file



#prior cov shape
#the left is the jags and the right is the R
R <- diag(c(1,1))
data=list(argos=obs,steps=steps,R=R,ind=ind,j=j,idx=idx,tracks=tracks,argos_class=obs_class,Month=MonthA,Months=max(MonthA,na.rm=T))

#parameters to track
pt<-c("theta","gamma","phi","alpha","state","gamma_mu","gamma_tau")

if(newModel){
  system.time(jagM<-jags.parallel(model.file = "./Bens_model.R",data=data,n.chains=2,parameters.to.save=pt,n.iter=270000,n.burnin=250000,n.thin=8,DIC=FALSE))
}

summar<-head(jagM$BUGSoutput$summary,44)
write.csv(jagM$BUGSoutput$summary, "./rhat.csv")
write.csv(summar,"./rhat_head.csv")

#delete jags objects
#rm(data)
#rm(argos)
#rm(obs)
#rm(j)
gc()

#bind chains
pc<-melt(jagM$BUGSoutput$sims.array)
str(pc)
#rm(jagM)
gc()
write.csv(pc, "./pc12again.csv")
colnames(pc)<-c("Draw","chain","par","value")

#extract parameter name
#what this does is it creates a news column called parameter that just matches the parameter name from par column

pc$parameter<-data.frame(str_match(pc$par,"(\\w+)"))[,-1]

#Extract index
#Now we are splitting pc along parameter lines into groups
splitpc<-split(pc,pc$parameter)

#single index
#I think this means we are taking the theta list only and we are adding the behavior
splitpc[c("theta")]<-lapply(splitpc[c("theta")],function(x){
  sv<-data.frame(str_match(x$par,"(\\w+)\\[(\\d+)]"))[,3]
  pc<-data.frame(x,Behavior=sv)
  return(pc)
})

##Double index
#single index
splitpc[c("alpha","gamma")]<-lapply(splitpc[c("alpha","gamma")],function(x){
  sv<-data.frame(str_match(x$par,"(\\w+)\\[(\\d+),(\\d+)]"))[,3:4]
  colnames(sv)<-c("Behavior","Month")
  pc<-data.frame(x,sv)
  return(pc)
})

#Three index
#Adding Animal, Track, step and behavior to Phi
splitpc[c("phi")]<-lapply(splitpc[c("phi")],function(x){
  #As matrices
  sv<-data.frame(str_match(x$par,"(\\w+)\\[(\\d+),(\\d+),(\\d+),(\\d+)]"))[,3:6]
  colnames(sv)<-c("Animal","Track","step","Behavior")
  pc<-data.frame(x,sv)
})


#State index
#adding Animal, Track and step to state
splitpc[c("state")]<-lapply(splitpc[c("state")],function(x){
  #As matrices
  sv<-data.frame(str_match(x$par,"(\\w+)\\[(\\d+),(\\d+),(\\d+)]"))[,3:5]
  colnames(sv)<-c("Animal","Track","step")
  pc<-data.frame(x,sv)
})

#bind all matrices back together
pc<-bind_rows(splitpc) 
# some warning about unequal factor levels
rm(splitpc)

#associate month with state
mm<-melt(MonthIndex)
colnames(mm)<-c("Animal","Track","step","MonthA")
pc<-merge(pc,mm,by=c("Animal","Track","step"),all.x=T)

#complete level matching, ugly but successful

pc[pc$parameter %in% c("alpha","gamma"),]$MonthA<-month.name[as.numeric(pc$Month[pc$parameter %in% c("alpha","gamma")])]
pc$MonthA<-factor(pc$MonthA,levels=month.name)

#plot all but phi and state
#I think this is the chains but we need a better way of visualizing



write.csv(pc,"./outdata_parse12.csv")