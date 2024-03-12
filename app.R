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
  tEvents=runif(n = nEvents)*0.8
  tScals=runif(n = nEvents,min = -1, max = 1)
  tFreqs=runif(n = nEvents,min = freq/2, max = freq*2)
  odf=data.frame(tEvents,tScals,tFreqs)
  #print(summary(odf))
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



dataSleip <- read.csv2(file = "Sleipner_IL1840_XL1130.csv",sep = ',',quote = '"',dec = '.')



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

ui <- dashboardPage(#skin = "black", 
  
  # Application title
  dashboardHeader(title = "Flattening QC"),
  
  dashboardSidebar(collapsed = F,width = 340,
      sidebarMenu( id = "tabs",#,collapsed=F,width = 12,background = 'black'), 
                   menuItem(tabName = 'dt',text = 'Full trace DT Analysis',icon = icon('chart-area')), selected=T,
                   menuItem(tabName = 'cc',text = 'Cross-corelation in window',icon = icon('chart-line')),
                   materialSwitch('useReal',"Use Sleipner data",value = F),
                   materialSwitch('prgrsvApply',"Use Progressive alignment",value = T),
                   menuItem(tabName = 'prm',text = 'Parameters',icon = icon('gear'), startExpanded = T,
                     menuItem(tabName='prm-pr',text = "Display:", startExpanded = T, 
                       sliderInput("tWinLoc",label = "Window Location, %",min = 0, max =100, step = 1,value = 50),
                       sliderInput("maxShift",label = "Max shift displayed, ms",min = 10, max =200, step = 10,value = 30)
                     ),
                     menuItem(tabName='prm-cc',text = "Correlation parameters:",
                        sliderInput("tWin",label = "Window size, %",min = 5, max =50, step = 1,value = 10),
                        sliderInput("corrLim",label = "Correlation threshold, %",min = 0, max =1, step = 0.1,value = 0.3),
                        checkboxInput('weiDT',"Weight DT by Cross-correlation",value = F),
                        checkboxInput('shapeCC',"Shape Cross-correlation with Gaussian shaper",value = T),
                        sliderInput("shaperRange",label = "Gaussian sigma range",min = 1, max =3, step = 0.1,value = 1)
                     ),
                     menuItem(tabName='prm-sg',text = "Synthetics generation:",
                       sliderInput("nEvent",label = "Events number",min = 1, max =100, step = 1,value = 8),
                       sliderInput("tFreq",label = "Base frequency",min = 5, max =100, step = 1,value = 15),
                       sliderInput("tScaler",label = "Scale T axis by",min = 0.8, max =1.2, step = 0.01,value = 1.05),
                       sliderInput("tShift",label = "Shift T axis by",min = -50, max =50, step = 2,value = -10)
                     ),
                     menuItem(tabName='prm-tr',text = "Transformations:",
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
      plotlyOutput('plot',height = "800")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'dt',
      plotlyOutput('plotDT', height = "1200")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'prm',
            plotOutput('plotDTprm', height = "1200")#,width = "10", height = "800")
    )
  ))
)

getDT <- function (x,y, corrLim, shapeCC=1) {
  #print(summary(x))
  #print(summary(y))
  #print(paste("getDT: ",sd(x),length(x),sd(y),length(y)))
  ccy <- ccf(x = x,y=y,lag.max = length(x)/2,type = 'correlation',plot = F)
  ccy$lag=ccy$lag

  if(!is.na(shapeCC)) ccy$acf=gausswin(length(ccy$acf),shapeCC)*ccy$acf
  ccsy=max(ccy$acf)
  dtsy=0
  #print(summary(x))
  #print(summary(y))
  #print(corrLim)
  #print(ccsy)
  if(ccsy>corrLim) dtsy=ccy$lag[which(ccy$acf==ccsy)]*2

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

getDTvector <- function(df,tWin,useAbs,useSquared,useEnvelop, useNorm, corrLim, freq=NULL, shapeCC=1) {

  n=100
  dtsy=rep(0,n)
  ccsy=rep(0,n)
  dtsz=rep(0,n)
  ccsz=rep(0,n)
  ts=rep(0,n)
  for(tWinLoc in c(0:n)) {
    winRange=c(floor(nrow(df)*max(0,(tWinLoc-tWin/2))/100),
               floor(nrow(df)*min(100,(tWinLoc+tWin/2))/100))
    #print(paste("winRange: ",winRange))
    t=(winRange[1]:winRange[2])*2
    x=df[winRange[1]:winRange[2],2]
    y=df[winRange[1]:winRange[2],3]
    #z=df[winRange[1]:winRange[2],4]
    
    if(!is.null(freq)) {
      x=getFiltered(x,freq)
      y=getFiltered(y,freq)
      #z=getFiltered(z,freq)
    }
    if(useAbs) {
      x=abs(x)
      y=abs(y)
      #z=abs(z)
    } 
    if(useSquared) {
      x=x**2
      y=y**2
      #z=z**2
    } 
    if(useEnvelop) {
      x=Re(envelope(x))
      y=Re(envelope(y))
      #z=Re(envelope(z))
    } 
    if(useNorm) {
      x=(x-mean(x))/sd(x)
      y=(y-mean(y))/sd(y)
      #z=(z-mean(z))/sd(z)
    }
    #print(summary(data.frame(x,y)))

    dty = getDT(x,y,corrLim,shapeCC)
    #dtz = getDT(x,z,corrLim,shapeCC)
    
    ccsy[tWinLoc]=dty[1]
    dtsy[tWinLoc]=dty[2]

    #ccsz[tWinLoc]=dtz[1]
    #dtsz[tWinLoc]=dtz[2]
    
    ts[tWinLoc]=mean(t)
    
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
    df$xx=(df$x)
    df$yy=(df$y)
    df$zz=(df$z)
    
    if(input$useFilter){
      df$xx=getFiltered(df$xx,input$freqRange)
      df$yy=getFiltered(df$yy,input$freqRange)
      df$zz=getFiltered(df$zz,input$freqRange)
    }
    if(input$useAbs) {
      df$xx=abs(df$xx)
      df$yy=abs(df$yy)
      df$zz=abs(df$zz)
    } 
    if(input$useSquared) {
      df$xx=df$xx**2
      df$yy=df$yy**2
      df$zz=df$zz**2
    } 
    if(input$useEnvelop) {
      df$xx=Re(spectral::envelope(df$xx))
      df$yy=Re(spectral::envelope(df$yy))
      df$zz=Re(spectral::envelope(df$zz))
    }
    if(input$useNorm) {
      df$xx=(df$xx-mean(df$xx))/sd(df$xx)
      df$yy=(df$yy-mean(df$yy))/sd(df$yy)
      df$zz=(df$zz-mean(df$zz))/sd(df$zz)
    }
    #print(summary(df))
    #getDTvector(t,y,y2,winRange)
    cc <- ccf(x = df$xx,y=df$yy,lag.max = nrow(df)/2,type = 'correlation',plot = F)
    cc$lag=cc$lag*2

    ccz <- ccf(x = df$xx,y=df$zz,lag.max = nrow(df)/2,type = 'correlation',plot = F)
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
    
    p1 = plot_ly(x=df$t,legendgroup = "Trace") %>% 
      layout() %>%
      add_lines(legendgroup = "traces",legendgrouptitle = list(text = "<b>Traces in window</b>"),y=df$x, color=I("gray70"), name = 'base') %>%
      add_lines(legendgroup = "traces",y=df$y, color=I("pink"),name = 'mon 1') %>%
      add_lines(legendgroup = "traces",y=df$z, color=I("lightblue"), name = 'mon 2') %>%
      add_lines(legendgroup = "traces",y=df$xx, color=I("black"), name = 'base transformed') %>%
      add_lines(legendgroup = "traces",y=df$yy, color=I("red"),name = 'mon 1 transformed') %>%
      add_lines(legendgroup = "traces",y=df$zz, color=I("blue"), name = 'mon 2 transformed') 
    p2 = plot_ly(x=cc$lag[,1,],legendgroup = "CC") %>% 
      layout(yaxis=list(range=c(-1,1))) %>%
      add_lines(legendgroup = "cc",legendgrouptitle = list(text = "<b>Correlation in window</b>"),y=cc$acf, color=I("red"),name = 'R mon 1') %>%
      add_lines(legendgroup = "cc",y=ccz$acf, color=I("blue"),name = 'R mon 2') %>%
      add_lines(legendgroup = "cc",y=gausswin(length(cc$acf),input$shaperRange), color=I("lightgreen"),name = 'R shaper')
      
    subplot(p1,p2,shareY=T,nrows = 2) %>% layout(title=title,showlegend = T,legend=list(tracegroupgap=300))
  })
  
  output$plotDT <- renderPlotly({
    
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
    #print(summary(df))
    dtsy <- getDTvector(df = df[,-4],
                        tWin = input$tWin,
                        useAbs = input$useAbs,
                        useSquared = input$useSquared,
                        useEnvelop = input$useEnvelop,
                        useNorm = input$useNorm,
                        corrLim = input$corrLim,
                        shapeCC = if(input$shapeCC) input$shaperRange else NA,
                        freq = if(input$useFilter) input$freqRange else NULL)
    if(input$prgrsvApply) {
      q=df$z
      shft_y = round((df$t - approx(dtsy$ts,dtsy$dtsy,xout = df$t,method = "linear")$y)/2)
      dt_range=range(shft_y,na.rm = T)
      shft_y[is.na(shft_y)]=round(seq(dt_range[1],dt_range[2],length.out=nrow(df)))[is.na(shft_y)]
      df$z=df$z[shft_y]
      #browser()
      dtsz <- getDTvector(df = df[,-3],
                          tWin = input$tWin,
                          useAbs = input$useAbs,
                          useSquared = input$useSquared,
                          useEnvelop = input$useEnvelop,
                          useNorm = input$useNorm,
                          corrLim = input$corrLim,
                          shapeCC = if(input$shapeCC) input$shaperRange else NA,
                          freq = if(input$useFilter) input$freqRange else NULL)
      dtsz$dtsy=dtsz$dtsy+dtsy$dtsy
      df$z=q
    } else {
      dtsz <- getDTvector(df = df[,-3],
                          tWin = input$tWin,
                          useAbs = input$useAbs,
                          useSquared = input$useSquared,
                          useEnvelop = input$useEnvelop,
                          useNorm = input$useNorm,
                          corrLim = input$corrLim,
                          shapeCC = if(input$shapeCC) input$shaperRange else NA,
                          freq = if(input$useFilter) input$freqRange else NULL)
    }
    colnames(dtsz) <- c('ts', 'ccsz','dtsz')
    dts=cbind(dtsy,dtsz[,-1])
    #browser()
    if(input$useReal) title="Sleipner_IL1840_XL1130"
    else title=paste("transform: t=t*",input$tScaler,"+",input$tShift)
    
    if(input$weiDT){
      shft_y = round((df$t - approx(dts$ts,dts$dtsy*dts$ccsy**2,xout = df$t,method = "linear")$y)/2)
      shft_z = round((df$t - approx(dts$ts,dts$dtsz*dts$ccsz**2,xout = df$t,method = "linear")$y)/2)
    } else {
      shft_y = round((df$t - approx(dts$ts,dts$dtsy,xout = df$t,method = "linear")$y)/2)
      shft_z = round((df$t - approx(dts$ts,dts$dtsz,xout = df$t,method = "linear")$y)/2)
    }
    
    shft_y[is.na(shft_y)]=c(1:nrow(df))[is.na(shft_y)]
    shft_z[is.na(shft_z)]=c(1:nrow(df))[is.na(shft_z)]
    
    shft_y = sapply(shft_y,FUN = function(x) {max(1,min(nrow(df),x))})
    shft_z = sapply(shft_z,FUN = function(x) {max(1,min(nrow(df),x))})
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
    half_win = input$tWin/100*nrow(df)
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
      add_lines(legendgroup = "mon1",y=-1.1+df$y[shft_y], color=I("pink"), name = 'aligned') %>%
      add_lines(legendgroup = "mon2",y=-1.2+df$z[shft_z], color=I("lightblue"), name = 'aligned') %>%
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
    
    
    subplot(p1,p2,p3,shareX=T,shareY=T,nrows = 3) %>% 
      layout(title=title,showlegend = T,legend=list(tracegroupgap=0))
    
  })
}

shinyApp(ui = ui, server = server,)