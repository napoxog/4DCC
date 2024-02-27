t= 0:1500
getFun <- function (t,f=5) 
{
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



dataSleip <- read.csv2(file = "F:\\R\\Projects\\CC\\Sleipner_IL1840_XL1130.csv",sep = ',',quote = '"',dec = '.')


require(shiny)
require(shinydashboard)
require(shinythemes)
require(signal)
require(spectral)

ui <- dashboardPage(#skin = "black", 
  
  # Application title
  dashboardHeader(title = "Flattening QC"),
  
  dashboardSidebar(collapsed = F,width = 340,
      sidebarMenu( id = "tabs",#,collapsed=F,width = 12,background = 'black'), 
                   menuItem(tabName = 'cc',text = 'Cross-Correlation',icon = icon('gear')),
                   menuItem(tabName = 'dt',text = 'Delta T',icon = icon('folder')), selected=T,
                   menuItem(tabName = 'prm',text = 'Parameters',icon = icon('gear'), startExpanded = T,
                     menuItem(tabName='prm-sg',text = "Synthetics generation:",
                       sliderInput("tFreq",label = "Base frequency",min = 5, max =100, step = 5,value = 5),
                       sliderInput("tScaler",label = "Scale T axis by",min = 0.8, max =1.2, step = 0.01,value = 1),
                       sliderInput("tShift",label = "Shift T axis by",min = -50, max =50, step = 2,value = 0)
                     ),
                     menuItem(tabName='prm-tr',text = "Transformations:",
                       checkboxInput('useAbs',"Use absolute values",value = F),
                       checkboxInput('useSquared',"Use squared values",value = F),
                       checkboxInput('useEnvelop',"Use signal Envelop",value = F),
                       checkboxInput('useNorm',"Normalize",value = T)
                     ),
                     menuItem(tabName='prm-pr',text = "Processing:", startExpanded = T, 
                       sliderInput("tWinLoc",label = "Window Location, %",min = 0, max =100, step = 1,value = 50),
                       sliderInput("tWin",label = "Window size, %",min = 5, max =50, step = 1,value = 10),
                       sliderInput("corrLim",label = "Correlation threshold, %",min = 0, max =1, step = 0.1,value = 0.3),
                       checkboxInput('useReal',"Use Sleipner data",value = T),
                       checkboxInput('useFilter',"ApplyFilter",value = T),
                       sliderInput("freqRange",label = "Filter frequency, Hz",min = 0, max =60, step = 1,value = c(10,20))
                     )
                   )
      )
  ),
  dashboardBody(
    tabItems(
    tabItem(tabName = 'cc',
      plotOutput('plot',height = "800")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'dt',
      plotOutput('plotDT', height = "1200")#,width = "10", height = "800")
    ),
    tabItem(tabName = 'prm',
            plotOutput('plotDTprm', height = "1200")#,width = "10", height = "800")
    )
  ))
)

getDT <- function (x,y, corrLim) {
  ccy <- ccf(x = x,y=y,lag.max = length(x)/2,type = 'correlation',plot = F)
  ccy$lag=ccy$lag*2

  
  ccsy=max(ccy$acf)
  dtsy=0
  
  if(ccsy>corrLim) dtsy=ccy$lag[which(ccy$acf==ccsy)]

  return(c(ccsy,dtsy))
}

getDTSS <- function (x,y, corrLim) {
  ccy <- ccf(x = x,y=y,lag.max = length(x)/2,type = 'correlation',plot = F)
  ccy$lag=ccy$lag*2
  
  
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

getDTvector <- function(df,tWin,useAbs,useSquared,useEnvelop, useNorm, corrLim, freq=NULL) {

  dtsy=rep(0,100)
  ccsy=rep(0,100)
  dtsz=rep(0,100)
  ccsz=rep(0,100)
  for(tWinLoc in c(0:100)) {
    winRange=c(floor(nrow(df)*max(1,(tWinLoc-tWin/2))/100),
               floor(nrow(df)*min(100,(tWinLoc+tWin/2))/100))
    t=(winRange[1]:winRange[2])*2
    x=df[winRange[1]:winRange[2],3]
    y=df[winRange[1]:winRange[2],1]
    z=df[winRange[1]:winRange[2],2]
    
    if(!is.null(freq)) {
      x=getFiltered(x,freq)
      y=getFiltered(y,freq)
      z=getFiltered(z,freq)
    }
    if(useAbs) {
      x=abs(x)
      y=abs(y)
      z=abs(z)
    } 
    if(useSquared) {
      x=x**2
      y=y**2
      z=z**2
    } 
    if(useEnvelop) {
      x=Re(envelope(x))
      y=Re(envelope(z))
      y=Re(envelope(y))
    } 
    if(useNorm) {
      x=(x-mean(x))/sd(x)
      y=(y-mean(y))/sd(y)
      z=(z-mean(z))/sd(z)
    }
    
    dty = getDT(x,y,corrLim)
    dtz = getDT(x,z,corrLim)
    
    ccsy[tWinLoc]=dty[1]
    dtsy[tWinLoc]=dty[2]

    ccsz[tWinLoc]=dtz[1]
    dtsz[tWinLoc]=dtz[2]
    
    #dtsz[tWinLoc]=ccz$lag[which(ccz$acf==ccsz[tWinLoc])]
  }
  return(data.frame(cbind(ccsy,ccsz,dtsy,dtsz)))
}

#getDTvector(dataSleip,5)

#return()
  
server <- function(input, output,session) {
  output$plot <- renderPlot({
    if(input$useReal) {
      winRange=c(floor(nrow(dataSleip)*max(1,(input$tWinLoc-input$tWin/2))/100),
                 floor(nrow(dataSleip)*min(100,(input$tWinLoc+input$tWin/2))/100))
      t=(winRange[1]:winRange[2])*2
      x=dataSleip[winRange[1]:winRange[2],3]
      y=dataSleip[winRange[1]:winRange[2],1]
      z=dataSleip[winRange[1]:winRange[2],2]
    } else {
      synLen=1500
      winRange=c(floor(synLen*max(1,(input$tWinLoc-input$tWin/2))/100),
                 floor(synLen*min(100,(input$tWinLoc+input$tWin/2))/100))
      t=(winRange[1]:winRange[2])*2
      x=getFun(t,f = input$tFreq)
      y=getFun(t*input$tScaler+input$tShift,f = input$tFreq)
      z=getFun(t*input$tScaler-input$tShift,f = input$tFreq)
    }
    #print(winRange)
    df=data.frame(cbind(t=t,x=x,y=y,z=z))
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
    print(summary(df))
    #getDTvector(t,y,y2,winRange)
    cc <- ccf(x = df$xx,y=df$yy,lag.max = nrow(df)/2,type = 'correlation',plot = F)
    cc$lag=cc$lag*2

    ccz <- ccf(x = df$xx,y=df$zz,lag.max = nrow(df)/2,type = 'correlation',plot = F)
    ccz$lag=ccz$lag*2
    
    par(mfrow = c(2,1))
    plot(x=df$t,y=df$x,t='l',ylim=range(df[,-1]))
    lines(df$t,y=df$y,col='red')
    lines(df$t,y=df$z,col='blue')
    lines(df$t,y=df$xx,col='black',lwd=2)
    lines(df$t,y=df$yy,col='red',lwd=2)
    lines(df$t,y=df$zz,col='blue',lwd=2)
    lines(x=rep(mean(range(df$t)),2),y=c(-10,10),col='pink')
    title(paste("transform: t=t*",input$tScaler,"+",input$tShift))

    plot(x=cc$lag,y=cc$acf,t='l',col='red')
    lines(x=ccz$lag,y=ccz$acf,col='blue')
    lines(x=c(0,0),y=c(-2,2),col='pink')
    
    dt=cc$lag[which(cc$acf==max(cc$acf))]
    title(main=paste("max CC=",round(digits=3,max(cc$acf)),"@ DT=",dt," ~",signif (digits = 3,100*dt/nrow(df)*2),"%"),
          sub=paste("DT estimation error = ", round(digits=4,100*((input$tShift-dt)/(input$tShift))),"%"))
    
  })
  
  output$plotDT <- renderPlot({
    
    if(input$useReal) {
      t=seq(0,nrow(dataSleip))*2
      x=dataSleip[,3]
      y=dataSleip[,1]
      z=dataSleip[,2]
    } else {
      synLen=1000
      t=seq(0,synLen)*2
      x=getFun(t,f = input$tFreq)
      y=getFun(t*input$tScaler+input$tShift,f = input$tFreq)
      z=getFun(t*input$tScaler-input$tShift,f = input$tFreq)
    }
    #print(winRange)
    df=data.frame(cbind(x=x,y=y,z=z))
    dts <- getDTvector(df = df,
                       tWin = input$tWin,
                       useAbs = input$useAbs,
                       useSquared = input$useSquared,
                       useEnvelop = input$useEnvelop,
                       useNorm = input$useNorm,
                       corrLim = input$corrLim,
                       freq = if(input$useFilter) input$freqRange else NULL)
    
    
    par(mfrow = c(3,1))
    plot(df$x-1,t='l',ylim=range(df))
    lines(df$y,col='red')
    lines(df$z+1,col='blue')
    lines(x=rep(input$tWinLoc/100*nrow(df),2),y=c(-10,10),col='pink')
    
    plot(dts$dtsy,t='l', col='red',ylim=c(-50,50))#range(dts$dtsy,dts$dtsz))
    lines(dts$dtsz,col='blue')
    lines(x=c(0,nrow(dts)),y=c(0,0),col='green')
    lines(x=rep(input$tWinLoc,2),y=c(-1000,1000),col='pink')
    
    plot(dts$ccsy,t='l',ylim=range(dts$ccsy,dts$ccsz),col='red')
    lines(dts$ccsz,col='blue')
    lines(x=c(0,nrow(dts)),y=c(input$corrLim,input$corrLim),col='green')
    lines(x=rep(input$tWinLoc,2),y=c(-10,10),col='pink')
    #ggplot()+geom_raster(data = matrix(dataSleip))
  })
}

shinyApp(ui = ui, server = server,)