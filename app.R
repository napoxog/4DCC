t= 0:1500
ricker <- function (t=0, f=25) {
  s=1/f
  return (2/sqrt(3*s)/pi^0.25*(1-(t/s)^2)*exp(-t*t/2/s/s))
}
#w=1
#t=seq(-w,w,0.1)
#a=ricker(t)

seed=runif(1)

genRCtrace <- function (nEvents=10, scale=1, freq=25) {
  .Random.seed <- seed
  tEvents=runif(n = nEvents)*0.6+0.2
  tScals=runif(n = nEvents,min = -1, max = 1)
  tFreqs=runif(n = nEvents,min = freq/2, max = freq*2)
  odf=data.frame(tEvents,tScals,tFreqs)
  print(summary(odf))
  return(odf)
}

traceRC <- genRCtrace (3)
print(traceRC)
  
getRandomRicker <- function(t=seq(1,1000), traceRC = traceRC, scale=1) {
  len=length(t)
  
  num=nrow(traceRC)
  tat=min(t) + traceRC$tEvent*diff(range(t))
  scals=traceRC$tScals*scale/sqrt(num)
  freqs=traceRC$tFreqs
  
  a=sin(t/2/pi/1000)/10
  #print(summary(data.frame(tat,scals,freqs)))
  for(i in seq(1,len)) {
    #a[i]=0
    for(w in seq(1,num)) {
      a[i]=a[i]+ricker(t = (t[i]-tat[w])/1000,f = freqs[w]) * scals[w]
    }
  }
  return(a)
}

getFun <- function (t,f=5) 
{
  
  return(getRandomRicker(t = t, traceRC,scale = 1))
  #browser()
  t=t/1000*(2*pi*f)
  f1=0.5*sin(t)+sin(t)*cos(t/3)
  w1=sapply(cos(t/5),FUN=function(x) { max(0.2,x)})
  return (f1*w1)
}

y=getFun(t)
y2=getFun(t*0.90+20)

par(mfrow = c(2,1))
plot(x=t,y=y,t='l')
lines(x=t,y=y2,col='red')
title(paste("transform: t=t*0.990+20"))

cc <- ccf(x = y,y=y2,lag.max = length(t)/2,type = 'correlation',plot = F)
plot(x=cc$lag,y=cc$acf,t='l')
title(paste("max CC=",max(cc$acf),"@ DT=",
                                         cc$lag[which(cc$acf==max(cc$acf))]))

print(c(cc$lag[which(cc$acf==max(cc$acf))],max(cc$acf)))



dataSleipRAW <- read.csv2(file = "Sleipner_IL1840_XL1130.csv",sep = ',',quote = '"',dec = '.')
dataSleipDSE <- read.csv2(file = "Sleipner_IL1840_XL1130_DSE.csv",sep = ',',quote = '"',dec = '.')
dataSleip<-dataSleipDSE


#plot(t,getRandomRicker(t,10,1),t='l')

require(shiny)
require(shinydashboard)
require(shinythemes)
require(shinyWidgets)
require(signal)
require(spectral)
require(ggplot2)
require(gridExtra)
require(plotly)
#### GUI ####
ui <- dashboardPage(#skin = "black", 
  
  # Application title
  dashboardHeader(title = "4D alignment QC"),
  
  dashboardSidebar(collapsed = F,width = 340,
      sidebarMenu( id = "tabs",#,collapsed=F,width = 12,background = 'black'), 
                   menuItem(tabName = 'dt',text = h3('Full trace DT Analysis', .noWS="outside"), selected=F),#icon = icon('chart-area')),
                   materialSwitch('useReal',"Use Sleipner 4D data",right = T,value = F),
                   materialSwitch('useAutoWin',"Automatic window size",right = T,value = T),
                   materialSwitch('prgrsvApply',"Progressive alignment",right = T,value = T),
                   materialSwitch('strtchApply',"Apply de-Stretch",right = T,value = T),
                   materialSwitch('phaseApply',"Apply de-Phase",right = T,value = T),
                   materialSwitch('multiFactor',"Use 3-factror analysis",right = T,value = T),
                   menuItem(tabName = 'cc',text = h3('Local Cross-corelation QC'), selected = T),#icon = icon('chart-line')),
                   menuItem(tabName = 'mfs',text = h3('Multi-factor Slice QC'), selected=F),#icon = icon('chart-line')),
                   menuItem(tabName = 'prm',text = h3('Parameters'), startExpanded = T,#,icon = icon('gear')
                     menuItem(tabName='prm-pr',text = h4("Display:"), startExpanded = T, 
                       sliderInput("tWinLoc",label = "Window Location, %",min = 0, max =100, step = 1,value = 50),
                       sliderInput("maxShift",label = "Max shift displayed, ms",min = 10, max =200, step = 10,value = 100),
                       sliderInput("sliceID",label = "Stretch-factor slice, ",min = 0.5, max =1.5, step = 0.1,value = 1)
                     ),
                     menuItem(tabName='prm-cc',text = h4("Correlation parameters:"),
                        sliderInput("tWin",label = "Window size, %",min = 5, max =50, step = 1,value = 30),
                        sliderInput("nWin",label = "Number of esimations",min = 10, max =100, step = 10,value = 10),
                        sliderInput("corrLim",label = "Correlation threshold, %",min = 0, max =1, step = 0.1,value = 0.3),
                        checkboxInput('weiDT',"Weight DT by Cross-correlation",value = F),
                        checkboxInput('shapeCC',"Shape Cross-correlation with Gaussian shaper",value = T),
                        sliderInput("shaperRange",label = "Gaussian sigma range",min = 1, max =3, step = 0.1,value = 3),
                        sliderInput("stretchLags",label = "Stretch Analysis resolution",min = 5, max =100, step = 5,value = 30)
                     ),
                     menuItem(tabName='prm-sg',text = h4("Synthetics generation:"),
                       sliderInput("nEvent",label = "Events number",min = 1, max =100, step = 1,value = 37),
                       sliderInput("tFreq",label = "Base frequency",min = 5, max =100, step = 1,value = 35),
                       sliderInput("tScaler",label = "Scale T axis by",min = 0.8, max =1.2, step = 0.01,value = 1.11),
                       sliderInput("tShift",label = "Shift T axis by",min = -50, max =50, step = 2,value = -50)
                     ),
                     menuItem(tabName='prm-tr',text = h4("Transformations:"),
                       checkboxInput('useAbs',"Use absolute values",value = F),
                       checkboxInput('useSquared',"Use squared values",value = F),
                       checkboxInput('useEnvelop',"Use signal Envelop",value = F),
                       checkboxInput('useNorm',"Normalize inputs (transformed)",value = T),
                       checkboxInput('useFilter',"Apply BP Filter",value = F),
                       sliderInput("freqRange",label = "Filter frequencies gate, Hz",min = 0, max =60, step = 1,value = c(10,20))
                     )
                   )
      )
  ),
  dashboardBody(
    tabItems(
    tabItem(tabName = 'cc',
      plotlyOutput('plot',height = "1200")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'mfs',
            plotlyOutput('plotMF',height = "600")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'dt',
      plotlyOutput('plotDT', height = "1200")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'prm',
            plotOutput('plotDTprm', height = "1200")#,width = "10", height = "800")
    )
  ))
)

applyTaper <- function (x, size=1/8) {
  X=x
  taper=round(1:round((length(x)*size)))
  taperW=(1-cos(pi/(length(x)*size)*taper))/2
  #browser();
  X[taper]= taperW[taper] * x[taper]
  X[(length(x)-taper+1)]=taperW[taper] * x[(length(x)-taper+1)]
  return(X)
}

applyMatrixRotate <- function (data, angle) {
  ## convert degrees to radians
  angle <- angle * pi / 180
  
  ## build rotation matrix
  rot <- rbind(c(0, -sin(angle), cos(angle)),
               c(0, cos(angle), sin(angle)),
               c(1, 0, 0))
  
  ## homogenise input data sets
  # if(nrow(data) == 2) {
  #   data <- rbind(data,
  #                 rep(0, ncol(data)))
  # }
  #data=applyTaper(data)
  data <- rbind(data, rep(1, length(data)),rep(0, length(data)))
  #data <- rbind(data, data,data)
  ## rotate signal traces
  data_out <- t(t(data) %*% rot)
  
  ## prepare output data
  #data_out <- data_out[3:1,]
  #browser()
  #phase=90;w=c(1:(250/length(Xfft))); X_ph = Xfft*(cos(phase/180*pi) + -1i*sin(phase/180*pi)); x_ph = ifft(X_ph);xxx=applyMatrixRotate (x, phase);xxx1=xxx[1,];xxx2=xxx[2,]
  #plot(Re(x_ph),t='l',col='red',ylim=c(-2,2)); lines(x) ; lines (X,col='blue'); lines(xxx1,col='green');lines(xxx2,col='pink')
  return(data_out[1,])
}

applyHilbertRotate <- function(x, phaseDeg=180.) {
  #X=applyTaper(x)
  
  #print("phase")
  #print(summary(wave))
  waveH=hilbert(f=1000, x)
  #summary(Re(ht)-wave$a)
  #print(phase)
  
  phase=Re(waveH)*0 + phaseDeg/180.*pi
  #print(phase)
  #print(class(cos(phase)*Re(waveH)-sin(phase)*Im(waveH)))
  #wave$a <- wave$a - (cos(phase)*Re(waveH)-sin(phase)*Im(waveH))
  X <-   Re(cos(phase)*Re(waveH)-sin(phase)*Im(waveH))[,1]
  #browser()
  #print(summary(wave))
  #t=c(1:length(X))
  #t=seq(-length(X)/2,length(X)/2,length.out=length(X))
  #X=approx(x=t+(phaseDeg/10),y=X,xout = t,method = "linear",rule = 2)$y
  return(x)
}

applyPhaseRotation <- function(x, phase=180.) {
  #X = rep(0,length(x)*2)
  #range=round(0.01+c((length(X)/4+1):round((length(X)*3/4))))
  #X [range] = x[1:length(x)]
  X=applyTaper(x)
  #taper=round(1:round((length(x)/4)))
  #taperW=(1-cos(pi*4/(length(x))*taper))/2
  #X[taper]= taperW[taper] * x[taper]
  #X[(length(x)-taper+1)]=taperW[taper] * x[(length(x)-taper+1)]
  #X = rep(0,length(x))
  #X = x
  Xfft = fft(X);
  #w = (2*pi/length(Xfft)) * (1:(length(Xfft)));
  w=c(1:(250/length(Xfft)))  
  k = pi*phase/180.;
  X_ph = Xfft*exp(-1i*(w+k));
  x_ph = ifft(X_ph)/length(X_ph);
  applyMatrixRotate (X, phase) 
  applyHilbertRotate (X, phase)
  #phase=120;plot(Re(H(X)*exp(1i*pi*phase/180)),t='l',ylim=c(-3,3)); lines(X,col='red')
  browser()
  return(Re(x_ph))
}

dePhaseRange = seq(-90,90,20)

dePhase <- function(x,y)  {
  lento = length(x)
  strch = dePhaseRange
  nstr=length(strch)
  ccs=strch
  sdxy=strch
  ccs_ys=list()
  y1=y
  ccfLag=0.5
  #x=applyTaper(x)
  #applyPhaseRotation(x)
  #cc1=rep(0,to)
  #  print(paste("before",strch))
  y_spl=y#[round(lento*0.25):round(lento*0.75)]
  
  sdorig=sd(x-y,na.rm = T)
  ccforig=ccf(x = x,y=y,lag.max = length(x)*ccfLag,type = 'correlation',plot = F,demean = F)
  ccorig=max(ccforig$acf)
  ccs=ccorig # fill ccs with maxCCorig
  #if(ccorig<0.6) 
  #  return(list(cc=1,stretch=0,y=y))
  
  maxCClag=ccforig$lag[which(ccorig==ccforig$acf)]
  assym=(maxCClag/diff(range(ccforig$lag))*2)*10
  #print(c(maxCClag,range(ccforig$lag),assym))
  
  t=seq(-10+assym,10+assym,length.out=length(y_spl))
  #t=seq(-10,10,length.out=length(y_spl))
  
  #print(range(t))

  for(ns in c(1:nstr)){
    y1=applyHilbertRotate(y,phase = strch[ns])
    #y1=spline(x=t*strch[ns],y=y_spl,xout = t,method = "fmm")$y
    #browser()
    sdxy[ns]=sd(x-y1,na.rm = T)
    cc=ccf(x = x,y=y1,lag.max = length(x)*ccfLag,type = 'correlation',plot = F,demean = F)
    ccs[ns]=max(cc$acf)
    ccs_ys[ns]=list(y1)
    #if(sd(x)>0.5)     browser()
    #print(paste("Stretch:",round(strch[ns],digits=3),"; CC=",round(max(cc$acf),digits = 3)))
    
  }
  #browser()
  maxCC=max(ccs)
  maxsidx=match(maxCC,ccs)
  
  minSD=min(sdxy)
  minSDidx=match(minSD,sdxy)
  
  print(paste("Phase:",round(strch[maxsidx],3),", CC=",round(ccorig,3)," >> ",round(maxCC,3),", SD=",round(sdorig,3)," >> ",round(minSD,3)))
  #if(sd(x)>0.5) browser() 
  #t=c(1:length(x)); plot(x=t,y=x,t='l',col='red'); lines(x=t,y=y); lines(x=t*strch[maxsidx]+max(t)*(strch[maxsidx]-1)/2,y=y1,col='blue'); plot(x=strch,y=ccs,t='l')
  #if(maxCC<0.0) return(list(cc=1,stretch=1,y=y))
  
  adjusted=list(cc=ccs,sd=sdxy,stretch=strch[maxsidx],y=unlist(ccs_ys[maxsidx]))
  return(adjusted)
}
  

stretchAdjustRange=c(0.5,1.5)

deStretch <- function(x,y,nstr=11) {
  lento = length(x)
  strch = seq(stretchAdjustRange[1],stretchAdjustRange[2],length.out=nstr)
  ccs=strch
  sdxy=strch
  ccs_ys=list()
  y1=y
  ccfLag=1
  #x=applyTaper(x)
  #applyPhaseRotation(x)
  #cc1=rep(0,to)
  #  print(paste("before",strch))
  y_spl=y#[round(lento*0.25):round(lento*0.75)]

  sdorig=sd(x-y,na.rm = T)
  ccforig=ccf(x = x,y=y,lag.max = length(x)*ccfLag,type = 'correlation',plot = F,demean = F)
  ccorig=max(ccforig$acf)
  ccs=ccorig # fill ccs with maxCCorig
  #if(ccorig<0.6) return(list(cc=ccs,stretch=1,y=y))

  maxCClag=ccforig$lag[which(ccorig==ccforig$acf)]
  assym=(maxCClag/diff(range(ccforig$lag))*2)*10
  #print(c(maxCClag,range(ccforig$lag),assym))
  
  t=seq(-10+assym,10+assym,length.out=length(y_spl))
  #t=seq(-10,10,length.out=length(y_spl))
  
  #print(range(t))
  
  for(ns in c(1:nstr)){
    y1=approx(x=t/strch[ns],y=y_spl,xout = t,method = "linear",rule = 2)$y
    #y1=spline(x=t*strch[ns],y=y_spl,xout = t,method = "fmm")$y
    y1=applyTaper(y1)
    #browser()
    sdxy[ns]=sd(x-y1,na.rm = T)
    cc=ccf(x = x,y=y1,lag.max = length(x)*ccfLag,type = 'correlation',plot = F,demean = F)
    ccs[ns]=max(cc$acf)
    ccs_ys[ns]=list(y1)
    #if(sd(x)>0.5)     browser()
    #print(paste("Stretch:",round(strch[ns],digits=3),"; CC=",round(max(cc$acf),digits = 3)))
    
  }
  maxCC=max(ccs)
  maxsidx=match(maxCC,ccs)

  minSD=min(sdxy)
  minSDidx=match(minSD,sdxy)

  print(paste("Stretch:",round(strch[maxsidx],3),", CC=",round(ccorig,3)," >> ",round(maxCC,3),", SD=",round(sdorig,3)," >> ",round(minSD,3)))
  #print(paste("Stretch:",round(strch[minSDidx],3),", SD=",round(sdorig,3)," >> ",round(minSD,3)))
  #if(sd(x)>0.5) browser() 
  #t=c(1:length(x)); plot(x=t,y=x,t='l',col='red'); lines(x=t,y=y); lines(x=t*strch[maxsidx]+max(t)*(strch[maxsidx]-1)/2,y=y1,col='blue'); plot(x=strch,y=ccs,t='l')
  #if(maxCC<0.0) return(list(cc=1,stretch=1,y=y))
  #browser()
  
  adjusted=list(cc=ccs,sd=sdxy,stretch=strch[maxsidx],y=unlist(ccs_ys[maxsidx]))
  return(adjusted)
}


getDT3dim <- function(x,y) {
  ccfLag=1
  sdorig=sd(x-y,na.rm = T)
  ccforig=ccf(x = x,y=y,lag.max = length(x)*ccfLag,type = 'correlation',plot = F,demean = F)
  ccorig=max(ccforig$acf)
  ccs=ccorig # fill ccs with maxCCorig
  #if(ccorig<0.6) return(list(cc=ccs,stretch=1,y=y))
  
  nstr=15
  lento = length(x)*2-1
  strch = seq(stretchAdjustRange[1],stretchAdjustRange[2],length.out=nstr)
  ccs=strch
  phs = dePhaseRange
  nphs= length(phs)
  #browser()
  CCnrows=nstr*nphs*lento
  CCmatrix=data.frame(
    st=rep(0,CCnrows),
    ph=rep(0,CCnrows),
    dt=rep(0,CCnrows),
    cc=rep(0,CCnrows))
  sdxy=strch
  ccs_ys=list()
  #y1=y
  #y_spl=y#[round(lento*0.25):round(lento*0.75)]
  
  
  maxCClag=ccforig$lag[which(ccorig==ccforig$acf)]
  assym=(maxCClag/diff(range(ccforig$lag))*2)*10
  #print(c(maxCClag,range(ccforig$lag),assym))
  
  #t=seq(-10+assym,10+assym,length.out=length(y))
  t=seq(-10,10,length.out=length(y))
  
  #print(range(t))
  for(ns in c(1:nstr)){
    y1=approx(x=t/strch[ns],y=y,xout = t,method = "linear",rule = 2)$y
    #y1=spline(x=t*strch[ns],y=y_spl,xout = t,method = "fmm")$y
    #y1=applyTaper(y1)
    for(np in c(1:nphs)){
      ph=dePhaseRange[np]
      st=strch[ns]
      y1=applyHilbertRotate(y1,ph)
      ccy <- ccf(x = x, y = y1, lag.max = length(x), type = 'correlation',plot = F,demean = F)
      CCrow=(ns-1)*nphs*lento+(np-1)*lento
      CCrange=c(CCrow:(CCrow+lento-1))+1
      #browser()
      CCmatrix$ph[CCrange] = rep(ph,lento)
      CCmatrix$st[CCrange] = rep(st,lento)
      CCmatrix$dt[CCrange] = ccy$lag
      CCmatrix$cc[CCrange] = ccy$acf
    }  
    #browser()
  }
  maxCC=max(CCmatrix$cc)
  maxsidx=match(maxCC,CCmatrix$cc)
  #y1=approx(x=t/CCmatrix[maxsidx]$st,y=y,xout = t,method = "linear",rule = 2)$y
  #y1=applyTaper(y1)
  #y1=applyHilbertRotate(y1,CCmatrix[maxsidx]$ph)
  #y1=approx(x=t/CCmatrix[maxsidx]$st,y=y,xout = t,method = "linear",rule = 2)$y
  #browser()
#  minSD=min(sdxy)
#  minSDidx=match(minSD,sdxy)
  resStr=paste(names(CCmatrix),round(CCmatrix[maxsidx,],3),collapse = ' ')
  print(paste("3dim DT:",round(CCmatrix$dt[maxsidx],3),", CC=",round(ccorig,3)," >> ",round(maxCC,3),' //',resStr))
  #print(paste("Stretch:",round(strch[minSDidx],3),", SD=",round(sdorig,3)," >> ",round(minSD,3)))
  #if(sd(x)>0.5) browser() 
  #t=c(1:length(x)); plot(x=t,y=x,t='l',col='red'); lines(x=t,y=y); lines(x=t*strch[maxsidx]+max(t)*(strch[maxsidx]-1)/2,y=y1,col='blue'); plot(x=strch,y=ccs,t='l')
  #if(maxCC<0.0) return(list(cc=1,stretch=1,y=y))
  #istr=10; CCfirst=istr*nphs*lento; CCrange=c(CCfirst:(CCfirst + nphs*lento-1))+1; im = CCmatrix[CCrange,];plot_ly(data=im) %>% add_heatmap(x=~ph,y=~dt,z=~cc) %>% layout(title = paste('Slice for stretch factor ',strch[istr+1],' max CC=',round(max(im$cc),3)))
  #browser()
  
#  adjusted=list(cc=ccs,sd=sdxy,stretch=strch[maxsidx],y=unlist(ccs_ys[maxsidx]))
  #print(paste(names(CCmatrix),CCmatrix[maxsidx,]))
  return(list(ccs=CCmatrix,maxsidx=maxsidx))
}

getDT <- function (x,y, corrLim, shapeCC=1, useStretch, usePhase, nstr,smpDT=2) {
  #print(summary(x))
  #print(summary(y))
  #print(paste("getDT: ",sd(x),length(x),sd(y),length(y)))
  #
  #x=applyTaper(x)
  #y=applyTaper(y)
  #maxCC=getDT3dim(x,y)
  #return(c(maxCC$cc,maxCC$dt))
  
  stretchDT=1
  if(useStretch){
    str_y=deStretch(x,y,nstr)
    y=str_y$y
    stretchDT=str_y$stretch
  }
  if(usePhase){
    str_ph=dePhase(x,y)
    y=str_ph$y
  }
  ccy <- ccf(x = x, y = y,lag.max = length(x),type = 'correlation',plot = F,demean = F)
  ccy$lag=ccy$lag

  if(!is.na(shapeCC)) ccy$acf=gausswin(length(ccy$acf),shapeCC)*ccy$acf
  ccsy=max(ccy$acf)
  dtsy=0
  #print(summary(x))
  #print(summary(y))
  #print(corrLim)
  #print(ccsy)
  if(ccsy>corrLim) dtsy=ccy$lag[which(ccy$acf==ccsy)]*smpDT*stretchDT
 
  return(c(ccsy,dtsy))
}

getDTSS <- function (x,y, corrLim) {
  #print(paste("getDT: ",length(x),length(y)))
  ccy <- ccf(x = x,y=y,lag.max = length(x)/2,type = 'correlation',plot = F)
  ccy$lag=ccy$lag
  
  
  ccsy=max(ccy$acf)
  dtsy=0
  
  if(ccsy>corrLim) dtsy=ccy$lag[which(ccy$acf==ccsy)]
  
  return(c(ccsy,dtsy))
}

getFiltered <- function (x, freq = c(5,10)) {
  fN=1/0.002/2
  freq=freq/fN
  #browser()
  bf= butter(1,freq) 
  return(filtfilt(bf,x))
}

applyTransform <- function(x,useFilter=F,freqRange=F,useAbs=F,useSquared=F,useEnvelop=F,useNorm=F) {
  x=applyTaper(x)
  if(useFilter)  x<-getFiltered(x,freqRange)
  if(useAbs)     x<-abs(x)
  if(useSquared) x<-x**2
  if(useEnvelop) x<-Re(spectral::envelope(x))
  if(useNorm)    x<-(x-mean(x))/sd(x)
  
  return (x)
}

getDTvector <- function(n=100, df,tWin,useFilter,freqRange,useAbs,useSquared,useEnvelop, useNorm, corrLim, freq=NULL, shapeCC=1, useStretch, nstr, usePhase,multiFactor) {
  dtsy=rep(0,n)
  ccsy=rep(0,n)
  #dtsz=rep(0,n)
  #ccsz=rep(0,n)
  smpDT=2
  ts=rep(0,n)
  tWin=tWin/100*n
  maxDT=diff(range(df$t))/n/2
  for(tWinLoc in c(1:n)) {
    winRange=c(round(nrow(df)*min(1,max(0,(tWinLoc-tWin/2))/n)),
               round(nrow(df)*min(1,max(0,(tWinLoc+tWin/2))/n)))
    #print(paste("winRange: ",winRange))
    t=(winRange[1]:winRange[2])*smpDT
    x=applyTransform(df[winRange[1]:winRange[2],2],useFilter,freqRange,useAbs,useSquared,useEnvelop,useNorm)
    y=applyTransform(df[winRange[1]:winRange[2],3],useFilter,freqRange,useAbs,useSquared,useEnvelop,useNorm)
    
    #print(summary(data.frame(x,y)))
    if(multiFactor) {
      maxCC=getDT3dim(x,y)
      maxsidx=maxCC$maxsidx
      dty = c(maxCC$ccs$cc[maxsidx],maxCC$ccs$dt[maxsidx])
    }
    else dty = getDT(x=x,y=y,corrLim = corrLim, shapeCC = shapeCC, useStretch = useStretch,usePhase = usePhase, nstr = nstr,smpDT = smpDT)
    dty[2]=max(min(dty[2],maxDT),-maxDT)
    #browser()

    ccsy[tWinLoc]=dty[1]
    dtsy[tWinLoc]=dty[2]

    ts[tWinLoc]=smpDT*mean(winRange)
    
    #dtsz[tWinLoc]=ccz$lag[which(ccz$acf==ccsz[tWinLoc])]
  }
  return(data.frame(cbind(ts,ccsy,dtsy)))
}

#getDTvector(dataSleip,5)

#return()
synLen=1000

rctvs <- reactiveValues(traceRC = traceRC)
  
server <- function(input, output,session) {
  #browser()
  observeEvent(input$nEvent, { rctvs$traceRC <- genRCtrace (nEvents = input$nEvent,freq = input$tFreq)})
  observeEvent(input$tFreq,  { rctvs$traceRC <- genRCtrace (nEvents = input$nEvent,freq = input$tFreq)})

  observe({
  print("strated new session...")
  print(paste0(
              "   protocol: ", session$clientData$url_protocol, 
              "   hostname: ", session$clientData$url_hostname, 
              "   pathname: ", session$clientData$url_pathname, 
              "   port: ",     session$clientData$url_port,     
              "   search: ",   session$clientData$url_search   
  ))
  })
  
  output$plot <- renderPlotly({
    #### Local Data ####
    if(input$useReal) {
      winRange=c(floor(nrow(dataSleip)*max(1,(input$tWinLoc-input$tWin/2))/100),
                 floor(nrow(dataSleip)*min(99,(input$tWinLoc+input$tWin/2))/100))
      #print(paste("winRange: ",winRange))
      t=(winRange[1]:winRange[2])*2
      x=dataSleip[winRange[1]:winRange[2],3]
      y=dataSleip[winRange[1]:winRange[2],1]
      z=dataSleip[winRange[1]:winRange[2],2]
    } else {
      winRange=c(floor(synLen*max(1,(input$tWinLoc-input$tWin/2))/100),
                     floor(synLen*min(99,(input$tWinLoc+input$tWin/2))/100))
      #print(paste("winRange: ",winRange))
      trc<-rctvs$traceRC
      subRange=(winRange[1]:winRange[2])
      t=subRange*2
      x=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      y=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      z=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
    }
    #print(winRange)
    df=data.frame(cbind(t=t,x=x,y=y,z=z))
    #print(summary(df))

    df$xx=applyTransform(df$x,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)
    df$yy=applyTransform(df$y,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)
    df$zz=applyTransform(df$z,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)
    #### Local Calc  ####
    #print(summary(df))
    #getDTvector(t,y,y2,winRange)
    
    a4syy=list()
    a4szz=list()
    a4pyy=list()
    a4pzz=list()
    
    if(input$strtchApply){
      a4syy=deStretch(df$xx,y=df$yy,input$stretchLags)
      a4szz=deStretch(df$xx,y=df$zz,input$stretchLags)
      df$yy=a4syy$y
      df$zz=a4szz$y
    }
    
    if(input$phaseApply){
      a4pyy=dePhase(df$xx,y=df$yy)
      a4pzz=dePhase(df$xx,y=df$zz)
      df$yy=a4pyy$y
      df$zz=a4pzz$y
    }
    
    cc <- ccf(x = df$xx,y=df$yy,lag.max = nrow(df),type = 'correlation',plot = F,demean = F)
    cc$lag=cc$lag*2

    ccz <- ccf(x = df$xx,y=df$zz,lag.max = nrow(df),type = 'correlation',plot = F,demean = F)
    ccz$lag=ccz$lag*2
    
    if(input$shapeCC) cc$acf=gausswin(length(cc$acf),input$shaperRange)*cc$acf
    if(input$shapeCC) ccz$acf=gausswin(length(ccz$acf),input$shaperRange)*ccz$acf
    
    if(input$useReal) title="Sleipner_IL1840_XL1130"
    else title=paste("transform: t=t*",input$tScaler,"+",input$tShift)
    
    dt=cc$lag[which(cc$acf==max(cc$acf))]
    dtz=ccz$lag[which(ccz$acf==max(ccz$acf))]
    # par(mfrow = c(2,1))
    # plot (x=df$t,y=df$x,t='l',ylim=range(df[,-1]))
    # lines(df$t,y=df$y,col='red')
    # lines(df$t,y=df$z,col='blue')
    # lines(df$t,y=df$xx,col='black',lwd=2)
    # lines(df$t,y=df$yy,col='red',lwd=2)
    # lines(df$t,y=df$zz,col='blue',lwd=2)
    # lines(x=rep(mean(range(df$t)),2),y=c(-10,10),col='pink')
    # 
    # if(input$useReal) title("Sleipner_IL1840_XL1130")
    # else title(paste("transform: t=t*",input$tScaler,"+",input$tShift))
    # 
    # plot (x=cc$lag,y=cc$acf,t='l',col='red')
    # lines(x=ccz$lag,y=ccz$acf,col='blue')
    # lines(x=c(0,0),y=c(-2,2),col='pink')
    # 
    # title(main=paste("max CC=",c(round(digits=3,max(cc$acf)),round(digits=3,max(ccz$acf))),"@ DT=",c(dt,dtz)/2),#," ~",signif (digits = 3,100*dt/nrow(df)*2),"%"),
    #       sub=paste("DT estimation error = ", round(digits=4,100*((input$tShift-dt)/(input$tShift))),"%"))
    #### Local Plot ####
    p4=p3=plot_ly() %>% add_text(x=0,y=0,text = 'NA')
    p1 = plot_ly(x=df$t,legendgroup = "Trace") %>% 
      layout() %>%
      add_lines(legendgroup = "traces",legendgrouptitle = list(text = "<b>Traces in window</b>"),y=df$x, color=I("gray70"), name = 'base') %>%
      add_lines(legendgroup = "traces",y=df$y, color=I("pink"),name = 'mon 1') %>%
      add_lines(legendgroup = "traces",y=df$z, color=I("lightblue"), name = 'mon 2') %>%
      add_lines(legendgroup = "traces",y=df$xx, color=I("black"), name = 'base transformed') %>%
      add_lines(legendgroup = "traces",y=df$yy, color=I("red"),name = 'mon 1 transformed') %>%
      add_lines(legendgroup = "traces",y=df$zz, color=I("blue"), name = 'mon 2 transformed') 
    p2 = plot_ly(x=cc$lag[,1,]*2,legendgroup = "CC") %>% 
      layout(yaxis=list(range=c(-1,1))) %>%
      add_lines(legendgroup = "cc",legendgrouptitle = list(text = "<b>Correlation in window</b>"),y=cc$acf, color=I("red"),name = 'R mon 1') %>%
      add_lines(legendgroup = "cc",y=ccz$acf, color=I("blue"),name = 'R mon 2') %>%
      add_lines(legendgroup = "cc",y=gausswin(length(cc$acf),input$shaperRange), color=I("lightgreen"),name = 'R shaper')
    if(input$strtchApply) {
    p3 = plot_ly(x=seq(stretchAdjustRange[1],stretchAdjustRange[2],length.out=length(a4syy$cc)),legendgroup = "CCstr") %>% 
      layout(yaxis=list(range=c(-1,1))) %>%
      add_lines(legendgroup = "CCstr",legendgrouptitle = list(text = "<b>Correlation vs stretch</b>"),y=a4syy$cc, color=I("red"),name = 'R mon 1') %>%
      add_lines(legendgroup = "CCstr",y=a4szz$cc, color=I("blue"),name = 'R mon 2') 
    }
    if(input$phaseApply) {
    p4 = plot_ly(x=dePhaseRange,legendgroup = "CCph") %>% 
      layout(yaxis=list(range=c(-1,1))) %>%
      add_lines(legendgroup = "CCph",legendgrouptitle = list(text = "<b>Correlation vs phase</b>"),y=a4pyy$cc, color=I("red"),name = 'R mon 1') %>%
      add_lines(legendgroup = "CCph",y=a4pzz$cc, color=I("blue"),name = 'R mon 2') 
    }
    
    subplot(p1,p2,p3,p4,shareY=T,nrows = 4) %>% layout(title=title,showlegend = T,legend=list(tracegroupgap=200))
  })
  #### MFS Data ####
  output$plotMF <- renderPlotly({
    if(input$useReal) {
      winRange=c(floor(nrow(dataSleip)*max(1,(input$tWinLoc-input$tWin/2))/100),
                 floor(nrow(dataSleip)*min(99,(input$tWinLoc+input$tWin/2))/100))
      #print(paste("winRange: ",winRange))
      t=(winRange[1]:winRange[2])*2
      x=dataSleip[winRange[1]:winRange[2],3]
      y=dataSleip[winRange[1]:winRange[2],1]
      z=dataSleip[winRange[1]:winRange[2],2]
    } else {
      winRange=c(floor(synLen*max(1,(input$tWinLoc-input$tWin/2))/100),
                 floor(synLen*min(99,(input$tWinLoc+input$tWin/2))/100))
      #print(paste("winRange: ",winRange))
      trc<-rctvs$traceRC
      subRange=(winRange[1]:winRange[2])
      t=subRange*2
      x=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      y=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      z=getRandomRicker(traceRC = trc,t = seq(1,synLen)*2)[subRange]
    }
    #print(winRange)
    df=data.frame(cbind(t=t,x=x,y=y,z=z))
    #print(summary(df))
    
    df$xx=applyTransform(df$x,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)
    df$yy=applyTransform(df$y,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)
    df$zz=applyTransform(df$z,input$useFilter,input$freqRange,input$useAbs,input$useSquared,input$useEnvelop,input$useNorm)

    CCmatrix1 = getDT3dim(df$xx,df$yy)
    CCmatrix2 = getDT3dim(df$xx,df$zz)
    
    nphs=length(dePhaseRange)
    lento=length(df$xx)*2-1
    strch=seq(stretchAdjustRange[1],stretchAdjustRange[2],length.out=15)
    istr=match(min(abs(strch-input$sliceID)),abs(strch-input$sliceID))
    CCfirst=istr*nphs*lento; 
    print(paste("istr=",istr))
    CCrange=c(CCfirst:(CCfirst + nphs*lento-1))+1; 
    im1 = CCmatrix1$ccs[CCrange,];
    im2 = CCmatrix2$ccs[CCrange,];
    p1=plot_ly(data=im1) %>% add_heatmap(x=~ph,y=~dt,z=~cc)  
      #layout(title = paste('Slice for stretch factor ',strch[istr+1],' max CC=',round(max(im1$cc),3)))
    p2=plot_ly(data=im2) %>% add_heatmap(x=~ph,y=~dt,z=~cc) 
      #layout(title = paste('Slice for stretch factor ',strch[istr+1],' max CC=',round(max(im2$cc),3)))
    
    subplot(p1,p2,shareY=T,nrows = 1) %>% layout(title=paste('Slice for stretch factor ',round(strch[istr+1],3),' max CC=',round(max(im2$cc),3)),showlegend = T,legend=list(tracegroupgap=200))
  })
  
  output$plotDT <- renderPlotly({
    #### Glob Data ####
    if(input$useReal) {
      t=seq(1,nrow(dataSleip))*2
      x=dataSleip[,3]
      y=dataSleip[,1]
      z=dataSleip[,2]
    } else {
      t=seq(1,synLen)*2
      trc<-rctvs$traceRC
      x=getRandomRicker(traceRC = trc,t = t)
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      y=getRandomRicker(traceRC = trc,t = t)
      trc$tEvents=trc$tEvents*input$tScaler+input$tShift/1000
      z=getRandomRicker(traceRC = trc,t = t)
    }
    #print(winRange)
    df=data.frame(cbind(t=t,x=x,y=y,z=z))
    
    spec=spec.fft(df$x,df$t/1000/pi)
    medfreq=median(spec$fx[spec$fx>0][amax(spec$PSD[spec$fx>0 & spec$PSD>mean(spec$PSD)])])-
      sd(spec$fx[spec$fx>0][amax(spec$PSD[spec$fx>0 & spec$PSD>mean(spec$PSD)])])
    print(paste("med freq, Hz = ",medfreq))
    if(input$useAutoWin) {
      t=seq(-10,10, length.out=1000)
      s=ricker(t=t,f = medfreq)
      sdt=max(20,min(500,diff(range(t[abs(s)>1e-3]))*1000))
      sdt=max(20,min(500,4000/medfreq))
      sdt_prc=round(100*sdt/diff(range(df$t)))
      updateSliderInput(session = session,"tWin",value = sdt_prc)
    } else {
      sdt_prc=input$tWin
      sdt=input$tWin*diff(range(df$t))/100
    }
    print(paste("win estimate, ms = ",sdt," (", sdt_prc, "%)"))
    nDT=input$nWin
    #browser()
    #print(summary(df))
    #### Glob Calc ####
    dtsy <- getDTvector(df = df[,-4], n = nDT,
                        tWin = sdt_prc, #input$tWin,
                        useFilter = input$useFilter,
                        freqRange = input$freqRange,
                        useAbs = input$useAbs,
                        useSquared = input$useSquared,
                        useEnvelop = input$useEnvelop,
                        useNorm = input$useNorm,
                        corrLim = input$corrLim,
                        useStretch = input$strtchApply,
                        usePhase = input$phaseApply,
                        nstr = input$stretchLags,
                        multiFactor = input$multiFactor,
                        shapeCC = if(input$shapeCC) input$shaperRange else NA,
                        freq = if(input$useFilter) input$freqRange else NULL)
    if(input$prgrsvApply) {
      q=df$z
      shft_y = round((df$t - approx(dtsy$ts,dtsy$dtsy,xout = df$t,method = "linear")$y)/2)
      dt_range=range(shft_y,na.rm = T)
      shft_y[is.na(shft_y)]=round(seq(dt_range[1],dt_range[2],length.out=nrow(df)))[is.na(shft_y)]
      df$z=df$z[shft_y]
      #browser()
      dtsz <- getDTvector(df = df[,-3], n = nDT,
                          tWin = sdt_prc, #input$tWin,
                          useFilter = input$useFilter,
                          freqRange = input$freqRange,
                          useAbs = input$useAbs,
                          useSquared = input$useSquared,
                          useEnvelop = input$useEnvelop,
                          useNorm = input$useNorm,
                          corrLim = input$corrLim,
                          useStretch = input$strtchApply,
                          usePhase = input$phaseApply,
                          nstr = input$stretchLags,
                          multiFactor = input$multiFactor,
                          shapeCC = if(input$shapeCC) input$shaperRange else NA,
                          freq = if(input$useFilter) input$freqRange else NULL)
      dtsz$dtsy=dtsz$dtsy+dtsy$dtsy
      df$z=q
    } else {
      dtsz <- getDTvector(df = df[,-3], n = nDT,
                          tWin = sdt_prc, #input$tWin,
                          useFilter = input$useFilter,
                          freqRange = input$freqRange,
                          useAbs = input$useAbs,
                          useSquared = input$useSquared,
                          useEnvelop = input$useEnvelop,
                          useNorm = input$useNorm,
                          corrLim = input$corrLim,
                          useStretch = input$strtchApply,
                          usePhase = input$phaseApply,
                          nstr = input$stretchLags,
                          multiFactor = input$multiFactor,
                          shapeCC = if(input$shapeCC) input$shaperRange else NA,
                          freq = if(input$useFilter) input$freqRange else NULL)
    }
    colnames(dtsz) <- c('ts', 'ccsz','dtsz')
    dts=cbind(dtsy,dtsz[,-1])
    #browser()
    if(input$useReal) title="Sleipner_IL1840_XL1130"
    else title=paste("transform: t=t*",input$tScaler,"+",input$tShift)
    spline_method="linear"
    # if(input$weiDT){
    #   shft_y = round((df$t - spline(dts$ts,dts$dtsy*dts$ccsy**2,xout = df$t,method = spline_method)$y)/2)
    #   shft_z = round((df$t - spline(dts$ts,dts$dtsz*dts$ccsz**2,xout = df$t,method = spline_method)$y)/2)
    # } else {
    #   shft_y = round((df$t - spline(dts$ts,dts$dtsy,xout = df$t,method = spline_method)$y)/2)
    #   shft_z = round((df$t - spline(dts$ts,dts$dtsz,xout = df$t,method = spline_method)$y)/2)
    # }
    # ddRange_y=range(shft_y,na.rm = T)
    # ddRange_z=range(shft_z,na.rm = T)
    # shft_y[is.na(shft_y)]=0#seq(ddRange_y[1],ddRange_y[2],nrow(df))[is.na(shft_y)]
    # shft_z[is.na(shft_z)]=0#seq(ddRange_z[1],ddRange_z[2],nrow(df))[is.na(shft_z)]
    # 
    # shft_y = sapply(shft_y,FUN = function(x) {max(1,min(nrow(df),x))})
    # shft_z = sapply(shft_z,FUN = function(x) {max(1,min(nrow(df),x))})
    
    if(input$weiDT){
      dts$dtsy=dts$dtsy*dts$ccsy**2
      dts$dtsz=dts$dtsz*dts$ccsz**2
    }
    # appliedDT_yt = spline(x=dts$ts-dts$dtsy,y=dts$ts,xout = df$t,method = spline_method)$y
    # appliedDT_zt = spline(x=dts$ts-dts$dtsz,y=dts$ts,xout = df$t,method = spline_method)$y
    
    appliedDT_yt = approx(x=dts$ts-dts$dtsy,y=dts$ts,xout = df$t,method = 'linear',rule=2)$y
    appliedDT_y = approx(x=appliedDT_yt,df$y,xout = df$t,method = spline_method)$y
    appliedDT_zt = approx(x=dts$ts-dts$dtsz,y=dts$ts,xout = df$t,method = 'linear',rule=2)$y
    appliedDT_z = approx(x=appliedDT_zt,df$z,xout = df$t,method = spline_method)$y
    
    print(paste("sd Y:",sd(appliedDT_y-df$y, na.rm = T)))
    print(paste("sd Z:",sd(appliedDT_z-df$z, na.rm = T)))
    # ytdiff = c(diff(appliedDT_yt),0)
    # ytdiff[ytdiff<-20]=10
    # ytdiff[ytdiff>20]=10
    # appliedDT_yt = diffinv(ytdiff)[-1]+appliedDT_yt[1]
    # 
    # ztdiff = c(diff(appliedDT_zt),0)
    # ztdiff[ztdiff<-20]=10
    # ztdiff[ztdiff>20]=10
    # appliedDT_zt = diffinv(ztdiff)[-1]+appliedDT_zt[1]
    
    # plot (x=dts$ts,y=dts$ts-dts$dtsz,t='l',col='blue'); lines(x=dts$ts,y=dts$ts-dts$dtsy,col='red');plot((appliedDT_zt),t='l',col='blue'); lines((appliedDT_yt),col='red')
    # plot(x=appliedDT_zt,y=df$z, t='l',col='blue');lines(x=appliedDT_yt,y=df$y,col='red')
    #browser()
    #print(summary(cbind(shft_y,shft_z)))
    
    #twin=c(1:100)/100*max(df$t)
    
    # par(mfrow = c(3,1))
    # plot (x=df$t,y=df$x-1,t='l',ylim=c(-4,4))#range(df[,-1]))
    # lines(x=df$t,y=df$y,col='red')
    # lines(x=df$t,y=df$z+1,col='blue')
    # lines(x=rep(input$tWinLoc/100*nrow(df)*2,2),y=c(-10,10),col='pink')
    # if(input$useReal) title("Sleipner_IL1840_XL1130")
    # else title(paste("transform: t=t*",input$tScaler,"+",input$tShift))
    # 
    # plot (dts$dtsy,t='l', col='red',ylim=c(-input$maxShift,input$maxShift))#range(dts$dtsy,dts$dtsz))
    # lines(dts$dtsy*dts$ccsy**2,col='red',lwd=2)
    # lines(dts$dtsz,col='blue')
    # lines(dts$dtsz*dts$ccsz**2,col='blue',lwd=2)
    # lines(x=c(0,nrow(dts)),y=c(0,0),col='green')
    # lines(x=rep(input$tWinLoc,2),y=c(-1000,1000),col='pink')
    # 
    # plot (dts$ccsy,t='l',ylim=range(0,1),col='red')
    # lines(dts$ccsz,col='blue')
    # lines(dts$ccsy**2,col='red',lwd=2)
    # lines(dts$ccsz**2,col='blue',lwd=2)
    # lines(x=c(0,nrow(dts)),y=c(input$corrLim,input$corrLim),col='green')
    # lines(x=rep(input$tWinLoc,2),y=c(-10,10),col='pink')
    #browser()
    #### Glob Plot ####
    half_win = sdt_prc/100*nrow(df)
    win_pos = input$tWinLoc/100*nrow(df)*2
    #browser()
    winBox=data.frame(x=c(win_pos-half_win,win_pos-half_win,win_pos+half_win,win_pos+half_win),
                      y=c(0,1,1,0))
    #print(winBox)
    p1 = plot_ly(x=df$t) %>% 
      layout(yaxis=list(range=c(-4,4))) %>%
      add_polygons(legendgroup = "decor",legendgrouptitle = list(text = "<b>Decorations</b>"),
                   x=winBox$x,y=(winBox$y-0.5)*8,name="Processing window", color=I("lightgreen")) %>%
      add_lines(legendgroup = "traces",legendgrouptitle = list(text = "<b>Traces</b>"),y=df$x-1, color=I("black"), name = 'base') %>%
      add_lines(legendgroup = "mon1",legendgrouptitle = list(text = "<b>Montor 1</b>"),y=df$y, color=I("red"),name = 'trace') %>%
      add_lines(legendgroup = "mon2",legendgrouptitle = list(text = "<b>Montor 2</b>"),y=df$z+1, color=I("blue"), name = 'trace') %>%
      #add_lines(legendgroup = "mon1",y=-1.1+df$y[shft_y], color=I("pink"), name = 'aligned') %>%
      #add_lines(legendgroup = "mon2",y=-1.2+df$z[shft_z], color=I("lightblue"), name = 'aligned') %>%
      add_lines(legendgroup = "mon1",y=-1.1+appliedDT_y, color=I("pink"), name = 'aligned') %>%
      add_lines(legendgroup = "mon2",y=-1.2+appliedDT_z, color=I("lightblue"), name = 'aligned') %>%
      add_lines(legendgroup = "decor",x=rep(input$tWinLoc/100*nrow(df)*2,2),y=c(-10,10),color=I('green'),name = 'Window loc.')
    p2 = plot_ly(x=dts$ts) %>% 
      layout(yaxis=list(range=c(-input$maxShift,input$maxShift)),xaxis=list(range=range(df$t))) %>%
      add_polygons(legendgroup = "decor",legendgrouptitle = list(text = "<b>Decorations</b>"),
                   x=winBox$x,y=(winBox$y-0.5)*input$maxShift*2,name="Processing window", color=I("lightgreen")) %>%
      add_lines(legendgroup = "mon1",y=dts$dtsy, color=I("red"),name = 'Delta T') %>%
      add_lines(legendgroup = "mon2",y=dts$dtsz, color=I("blue"),name = 'Delta T') %>%
      add_lines(legendgroup = "mon1",y=dts$dtsy*dts$ccsy**2, color=I("pink"),name = 'Delta T * R<sup>2') %>%
      add_lines(legendgroup = "mon2",y=dts$dtsz*dts$ccsz**2, color=I("lightblue"),name = 'Delta T * R<sup>2') %>%
      add_lines(legendgroup = "decor",x=rep(input$tWinLoc/100*nrow(df)*2,2),y=c(-input$maxShift,input$maxShift),color=I('green'),name = 'Window loc.',showlegend = F)
    p3 = plot_ly(x=dts$ts) %>% 
      layout(yaxis=list(range=c(0,1)),xaxis=list(range=range(df$t))) %>%
      add_polygons(legendgroup = "decor",legendgrouptitle = list(text = "<b>Decorations</b>"),
                   x=winBox$x,y=winBox$y,name="Processing window", color=I("lightgreen")) %>%
      add_lines(legendgroup = "mon1",y=dts$ccsy, color=I("red"),name = 'R') %>%
      add_lines(legendgroup = "mon2",y=dts$ccsz, color=I("blue"),name = 'R') %>%
      add_lines(legendgroup = "mon1",y=dts$ccsy**2, color=I("pink"),name = 'R<sup>2</sup>') %>%
      add_lines(legendgroup = "mon2",y=dts$ccsz**2, color=I("lightblue"),name = 'R<sup>2</sup>') %>%
      add_lines(legendgroup = "decor",x=range(df$t),y=c(input$corrLim,input$corrLim), color=I("black"),name = 'Correlation limit') %>%
      add_lines(legendgroup = "decor",x=rep(win_pos,2),y=c(0,1),color=I('green'),name = 'Window loc.',showlegend = F)
    
    # wf=waterfall(y=df$x,x=df$t/1000,nf = 1)
    
    
    
    # p4 = plot_ly(x=df$t) %>%
    #  layout(xaxis=list(range=c(0,250))) %>%
    #  add_lines(legendgroup = "decor",legendgrouptitle = list(text = "<b>Decorations</b>"),
    #            x=spec$fx[spec$fx>0],y=spec$PSD[spec$fx>0])
    #browser()
    subplot(p1,p2,p3,shareX=T,shareY=T,nrows = 3) %>% 
      layout(title=title,showlegend = T,legend=list(tracegroupgap=0))
    
  })
}

shinyApp(ui = ui, server = server,)