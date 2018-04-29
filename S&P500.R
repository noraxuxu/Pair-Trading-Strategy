setwd("C://r4")

packs<-c("quantmod","usdm","portes","dplyr","zoo","tseries","IBrokers","beepr","httpuv","gmailr","plotrix","Quandl","data.table")
install.packages(packs)

library(dplyr)
library(quantmod)
library(usdm)
library(portes)
library(zoo)
library(tseries)
library(beepr)
library(gmailr)
library(httpuv)
library(IBrokers)
library(plotrix)
library(Quandl)
library(data.table)

context.stocks = c('SPY','XLF','ABT', 'ABBV', 'ACN', 'ACE', 'ADBE', 'ADT', 'AAP', 'AES', 'AET', 'AFL', 'AMG', 'A', 'GAS', 'APD', 'ARG', 'AKAM', 'AGN', 'ALXN', 'ALLE', 'ADS', 'ALL', 'ALTR', 'MO', 'AMZN', 'AEE', 'AAL', 'AEP', 'AXP', 'AIG', 'AMT', 'AMP', 'ABC', 'AME', 'AMGN', 'APH', 'APC', 'ADI', 'AON', 'APA', 'AIV', 'AMAT', 'ADM', 'AIZ', 'T', 'ADSK', 'ADP', 'AN', 'AZO', 'AVGO', 'AVB', 'AVY', 'BHI', 'BLL', 'BAC', 'BK', 'BCR', 'BXLT', 'BAX', 'BBT', 'BDX', 'BBBY', 'BRK.B', 'BBY', 'BLX', 'HRB', 'BA', 'BWA', 'BXP', 'BSK', 'BMY', 'BRCM', 'BF.B', 'CHRW', 'CA', 'CVC', 'COG', 'CAM', 'CPB', 'COF', 'CAH', 'HSIC', 'KMX', 'CCL', 'CAT', 'CBG', 'CBS', 'CELG', 'CNP', 'CTL', 'CERN', 'CF', 'SCHW', 'CHK', 'CVX', 'CMG', 'CB', 'CI', 'XEC', 'CINF', 'CTAS', 'CSCO', 'C', 'CTXS', 'CLX', 'CME', 'CMS', 'COH', 'KO', 'CCE', 'CTSH', 'CL', 'CMCSA', 'CMA', 'CSC', 'CAG', 'COP', 'CNX', 'ED', 'STZ', 'GLW', 'COST', 'CCI', 'CSX', 'CMI', 'CVS', 'DHI', 'DHR', 'DRI', 'DVA', 'DE', 'DLPH', 'DAL', 'XRAY', 'DVN', 'DO', 'DTV', 'DFS', 'DISCA', 'DISCK', 'DG', 'DLTR', 'D', 'DOV', 'DOW', 'DPS', 'DTE', 'DD', 'DUK', 'DNB', 'ETFC', 'EMN', 'ETN', 'EBAY', 'ECL', 'EIX', 'EW', 'EA', 'EMC', 'EMR', 'ENDP', 'ESV', 'ETR', 'EOG', 'EQT', 'EFX', 'EQIX', 'EQR', 'ESS', 'EL', 'ES', 'EXC', 'EXPE', 'EXPD', 'ESRX', 'XOM', 'FFIV', 'FB', 'FAST', 'FDX', 'FIS', 'FITB', 'FSLR', 'FE', 'FSIV', 'FLIR', 'FLS', 'FLR', 'FMC', 'FTI', 'F', 'FOSL', 'BEN', 'FCX', 'FTR', 'GME', 'GPS', 'GRMN', 'GD', 'GE', 'GGP', 'GIS', 'GM', 'GPC', 'GNW', 'GILD', 'GS', 'GT', 'GOOGL', 'GOOG', 'GWW', 'HAL', 'HBI', 'HOG', 'HAR', 'HRS', 'HIG', 'HAS', 'HCA', 'HCP', 'HCN', 'HP', 'HES', 'HPQ', 'HD', 'HON', 'HRL', 'HSP', 'HST', 'HCBK', 'HUM', 'HBAN', 'ITW', 'IR', 'INTC', 'ICE', 'IBM', 'IP', 'IPG', 'IFF', 'INTU', 'ISRG', 'IVZ', 'IRM', 'JEC', 'JBHT', 'JNJ', 'JCI', 'JOY', 'JPM', 'JNPR', 'KSU', 'K', 'KEY', 'GMCR', 'KMB', 'KIM', 'KMI', 'KLAC', 'KSS', 'KRFT', 'KR', 'LB', 'LLL', 'LH', 'LRCX', 'LM', 'LEG', 'LEN', 'LVLT', 'LUK', 'LLY', 'LNC', 'LLTC', 'LMT', 'L', 'LOW', 'LYB', 'MTB', 'MAC', 'M', 'MNK', 'MRO', 'MPC', 'MAR', 'MMC', 'MLM', 'MAS', 'MA', 'MAT', 'MKC', 'MCD', 'MHFI', 'MCK', 'MJN', 'MMV', 'MDT', 'MRK', 'MET', 'KORS', 'MCHP', 'MU', 'MSFT', 'MHK', 'TAP', 'MDLZ', 'MON', 'MNST', 'MCO', 'MS', 'MOS', 'MSI', 'MUR', 'MYL', 'NDAQ', 'NOV', 'NAVI', 'NTAP', 'NFLX', 'NWL', 'NFX', 'NEM', 'NWSA', 'NEE', 'NLSN', 'NKE', 'NI', 'NE', 'NBL', 'JWN', 'NSC', 'NTRS', 'NOC', 'NRG', 'NUE', 'NVDA', 'ORLY', 'OXY', 'OMC', 'OKE', 'ORCL', 'OI', 'PCAR', 'PLL', 'PH', 'PDCO', 'PAYX', 'PNR', 'PBCT', 'POM', 'PEP', 'PKI', 'PRGO', 'PFE', 'PCG', 'PM', 'PSX', 'PNW', 'PXD', 'PBI', 'PCL', 'PNC', 'RL', 'PPG', 'PPL', 'PX', 'PCP', 'PCLN', 'PFG', 'PG', 'PGR', 'PLD', 'PRU', 'PEG', 'PSA', 'PHM', 'PVH', 'QRVO', 'PWR', 'QCOM', 'DGX', 'RRC', 'RTN', 'O', 'RHT', 'REGN', 'RF', 'RSG', 'RAI', 'RHI', 'ROK', 'COL', 'ROP', 'ROST', 'RLC', 'R', 'CRM', 'SNDK', 'SCG', 'SLB', 'SNI', 'STX', 'SEE', 'SRE', 'SHW', 'SIAL', 'SPG', 'SWKS', 'SLG', 'SJM', 'SNA', 'SO', 'LUV', 'SWN', 'SE', 'STJ', 'SWK', 'SPLS', 'SBUX', 'HOT', 'STT', 'SRCL', 'SYK', 'STI', 'SYMC', 'SYY', 'TROW', 'TGT', 'TEL', 'TE', 'TGNA', 'THC', 'TDC', 'TSO', 'TXN', 'TXT', 'HSY', 'TRV', 'TMO', 'TIF', 'TWX', 'TWC', 'TJK', 'TMK', 'TSS', 'TSCO', 'RIG', 'TRIP', 'FOXA', 'TSN', 'TYC', 'UA', 'UNP', 'UNH', 'UPS', 'URI', 'UTX', 'UHS', 'UNM', 'URBN', 'VFC', 'VLO', 'VAR', 'VTR', 'VRSN', 'VZ', 'VRTX', 'VIAB', 'V', 'VNO', 'VMC', 'WMT', 'WBA', 'DIS', 'WM', 'WAT', 'ANTM', 'WFC', 'WDC', 'WU', 'WY', 'WHR', 'WFM', 'WMB', 'WEC', 'WYN', 'WYNN', 'XEL', 'XRX', 'XLNX', 'XL', 'XYL', 'YHOO', 'YUM', 'ZBH', 'ZION', 'ZTS')  

SP500_EOD = paste("EOD/",context.stocks,".11",sep="")
write.csv(SP500_EOD,"SP500_EOD.csv",row.names = F)
write.csv(context.stocks,"SP500.csv",row.names = F)

SP500_EOD = read.csv(file="SP500_EOD.csv",head=TRUE,stringsAsFactors=FALSE,sep=";")
stocks = unique(SP500_EOD$x)

SP500 = read.csv(file="SP500.csv",head=TRUE,stringsAsFactors=FALSE,sep=";")
stockNames = unique(SP500$x)

#entryPoint=2.2
#exit_signal_1=0.5
#exit_signal_2=2.7
#lookbackPeriod=15  see line 144 to set lookback Period

#Retrieve historic prices from Quandl package
Quandl.api_key("xxxxxxxxxxxxxxx")  # quandl licence key
closing = Quandl(stocks, start_date="2017-04-01", end_date="2017-11-17") 
beep(4)
row.names(closing) = closing[,1]
closing = closing[,-1]
colnames(closing) = stockNames
closing2 = closing[, apply(closing,2,function(x) !any(is.na(x)))]
closingprice = closing2
stockNames = colnames(closingprice)

#Add real-time data
#real_time_price = getQuote(stockNames)$Last
#time = Sys.time()
#closingprice_new = rbind(closing,real_time_price)
#row.names(closingprice_new)[dim(closingprice_new)[1]] = as.character(time)
#closingprice = closingprice_new

#Calculate daily return
return = c()
for (i in 1:dim(closingprice)[2]){
  #Compute the Daily Returns
  prices=closingprice[,i]
  n <- length(prices)
  dailyReturn <- log(prices[-1]/prices[-n])
  return=cbind(return,dailyReturn)
}
return = as.data.frame(return)
colnames(return) = stockNames


x1= return[,1]
x2= return[,2]
df=as.data.frame(cbind(x1, x2))
df=cbind(rownames(df),df)
colnames(df)=c("Date","x1","x2")
options(warn=-1) 
x2= log( return[,2])
x2[is.na(x2)] <- 0
x2[which(x2==Inf)] <- 0
x2[which(x2==-Inf)] <- 0 
#
# Now, after filtering the contents of x2 to replace na, Inf, and - Inf with 0, 
# we map x1 and the log of x2 to the dataframe df
#
df=as.data.frame(cbind(x1, x2))
#
# Next, declare four vectors to be used in the next part of the script
#
betaCoeff=c()
t_val=c()
p_val=c()
p_val1=c()
stockSpecificReturn=c()
for (i in (1+dim(df)[2]):(length(stockNames)))
{
  reg=lm(return[,i] ~ x1+ x2)
  plot(main=stockNames[i], reg$residuals)
  acf(main=stockNames[i], reg$residuals)
  betaCoeff=rbind(betaCoeff,reg$coefficients)
  
  stationaryTest= adf.test(reg$residuals,alternative = "stationary")
  p_val=rbind(p_val,sapply(stationaryTest, "[[", 1) [4])
  
  #Stock Specific Return
  stockSpecificReturn= cbind(stockSpecificReturn, return[,i] - (x1+ x2))
  
  stationaryTest_Stock= adf.test((return[,i] - (x1+ x2))[complete.cases(return[,i] - (x1+ x2))],alternative = "stationary")
  p_val1=rbind(p_val1,sapply(stationaryTest_Stock, "[[", 1) [4])
  
  #Get t-value
  t_stat= summary(reg)$coef[,"t value"]
  t_val=rbind(t_val,t_stat)
  #fit=return[,i]-fitted(reg)
}

#Beta Coefficient Dataframe 
rownames(betaCoeff)=stockNames[(1+dim(df)[2]):(length(stockNames))]
betaCoeff=as.data.frame(betaCoeff)

#T_Val Dataframe 
rownames(t_val)=stockNames[(1+dim(df)[2]):(length(stockNames))]
t_val=as.data.frame(t_val)

#Remove the Stocks Having T stat Absolute Value less than 1.95
#print("4. Remove stocks having T stat less than 1.95")
#Sys.time()
#
stockUsed=c()
for (i in 1:dim(t_val)[1])
{
  if (abs(t_val[i,2])>=1.95 & abs(t_val[i,3])>=1.95)
  {
    stockUsed= rbind(stockUsed,rownames(t_val)[i])
  }
}

#Generate pairs from stocks above
pairs=combn(stockUsed,2)



#Analysis on certain days before the most recent day (check potential trading)
lookingforward_summary = c()
A<-0
B<-0
divergence = 2.25
lookbackPeriod = 15
day_prior = 7
positions=c()

for(i in 1:dim(pairs)[2]){
  A<-round(i/200)
  if(A>B){B<-A;
  print(B*200)}
  
  #Name of Stock Pairs
  stock1_Name= pairs[1,i]
  stock2_Name= pairs[2,i]
  
  #Index of Stocks
  getIndex1=which(colnames(closingprice)==stock1_Name)
  getIndex2=which(colnames(closingprice)==stock2_Name)
  
  #Prices of Stocks
  stockPrice1=closingprice[1:(dim(closingprice)[1]),getIndex1]
  stockPrice2=closingprice[1:(dim(closingprice)[1]),getIndex2]
  
  #Returns of Stocks
  returnStock1=return[,getIndex1]
  returnStock2=return[,getIndex2]
  
  #Index of the Stock in Beta Coefficient Matrix
  
  stock1Index=which(row.names(betaCoeff)==stock1_Name)
  stock2Index=which(row.names(betaCoeff)==stock2_Name)
  
  #Calculate the Common Factors for Both Stocks
  #commonFactorStock1= betaCoeff$x1[stock1Index] *x1 + betaCoeff$x2[stock1Index]*x2  
  #commonFactorStock2= betaCoeff$x1[stock2Index] *x1 + betaCoeff$x2[stock2Index]*x2
  
  #Specific Returns
  #stock1_Specific= stockSpecificReturn[,getIndex1-2]
  #stock2_Specific= stockSpecificReturn[,getIndex2-2]
  
  #lookbackPeriod we can change these numbers
  
  #Look at the historical data to get the optimal hedge ratio
  
  #Get stock prices within the look back period
  stock1_look = stockPrice1[(dim(closingprice)[1]-lookbackPeriod-day_prior):(dim(closingprice)[1]-day_prior)]
  stock2_look = stockPrice2[(dim(closingprice)[1]-lookbackPeriod-day_prior):(dim(closingprice)[1]-day_prior)]
  
  #regr= lm(commonFactorStock1[t:(t+lookbackPeriod)] ~ commonFactorStock2[t:(t+lookbackPeriod)])
  regr= lm( stock1_look ~ stock2_look)
  #regr= lm(stock1_Specific[t:(t+lookbackPeriod)]~ stock2_Specific[t:(t+lookbackPeriod)])
  
  #Regression Details
  #summary(regr)
  #acf(reg$residuals)
  
  #Hedge Ratio
  hedgeRatio=regr$coefficients[2]
  #pairs_hedgeRatio=rbind(pairs_hedgeRatio, c(stock1,stock2,hedgeRatio))
  
  #Spread With Common Factors Difference
  #spread= commonFactorStock1[(t+lookbackPeriod):(t+2*lookbackPeriod)]-hedgeRatio*commonFactorStock2[(t+lookbackPeriod):(t+2*lookbackPeriod)]
  
  #Spread with Specific Factors
  #spread= stock1_Specific[(t+lookbackPeriod):(t+2*lookbackPeriod)]-hedgeRatio*stock1_Specific[(t+lookbackPeriod):(t+2*lookbackPeriod)]
  
  #Spread with Stock Price Difference
  spread= stock1_look-hedgeRatio*stock2_look
  
  #Scale the spread
  spread_scale = scale(spread)
  
  #Mean of Spread
  spreadMean=mean(spread)
  
  #Std of Spread
  spreadStd= sd(spread)
  
  #Z-score
  z_score=(spread-spreadMean)/spreadStd
  
  #If z_score is not beyond the desirable standard deviation, then quit and check next pair
  if(abs(tail(z_score,n=1)) < divergence){
    next
  }else{
    #Two days before the most recent day
    dates=row.names(closingprice)[(dim(closingprice)[1])-day_prior]
    
    #Compute the Z Score
    z_score_current = tail(z_score,n=1)
    positions = data.frame(z_score=z_score_current, PositionLong= ifelse(z_score_current<=-divergence,1,0),
                           PositionShort= ifelse(z_score_current>=divergence,-1,0))
    positions$OverallPosition=positions$PositionLong+positions$PositionShort
    
    #Get the Date and Overall Position
    positions=positions[c(1,4)]
    
    positions=cbind(dates,stockPrice1[(dim(closingprice)[1])-day_prior],stockPrice2[(dim(closingprice)[1])-day_prior],hedgeRatio,positions)
    
    colnames(positions)=c("Dates", "stock1Price","stock2Price","HedgeRatio","z_score","Strategy")
    
    positions = mutate(positions, Stock1 = stock1_Name, Stock2 = stock2_Name)
    positions = positions[c("Dates","Stock1","stock1Price","Stock2","stock2Price",
                            "z_score","HedgeRatio","Strategy")]
    
    #Add columns Long/Short, #Shares, Position
    positions$STK1_longshort = ifelse(positions$Strategy==1,"Long","Short")
    positions$STK2_longshort = ifelse(positions$Strategy==1,"Short","Long")
    
    positions$STK1_Position = 10000
    positions$STK1_Share = positions$STK1_Position / positions$stock1Price
    
    positions$STK2_Share = positions$STK1_Share * positions$HedgeRatio
    positions$STK2_Position = positions$stock2Price * positions$STK2_Share
    
    positions$Position_diff = ifelse(positions$STK1_longshort=="Long", positions$STK1_Position - positions$STK2_Position,
                                     positions$STK2_Position - positions$STK1_Position)
    
    positions = positions[c("Dates","Stock1","stock1Price","STK1_longshort","STK1_Share","STK1_Position",
                            "Stock2","stock2Price","STK2_longshort","STK2_Share","STK2_Position",
                            "z_score","HedgeRatio","Position_diff","Strategy")]
    
    positions$p_diff_sd = ifelse(positions$STK1_longshort=="Long", 
                                 sd(positions$STK1_Position - stock2_look * positions$HedgeRatio * positions$STK1_Position /stock1_look),
                                 sd(stock2_look * positions$HedgeRatio * positions$STK1_Position /stock1_look -  positions$STK1_Position))
    
    #Add current positions to summary
    lookingforward_summary = unique(rbind(lookingforward_summary, positions))
    
    #Calculate the Convergence and Divergence Spread
    #convergencePoint= convergence_signal*spreadStd+spreadMean
    #divergencePoint= divergence_signal*spreadStd+spreadMean
    
    lookingforward_summary = lookingforward_summary %>% arrange(abs(Position_diff))
    
    
    #calculate the standard deviation for the position difference 
    
  }
}
rownames(lookingforward_summary) = 1:dim(lookingforward_summary)[1]
write.csv(lookingforward_summary,paste("C:\\rOutput\\SP500 Good Pairs Selected",dates,divergence,"SD",lookbackPeriod,"days", day_prior, "prior days.csv"),row.names = FALSE)


#Plotting
dates=row.names(closingprice)[(dim(closingprice)[1])]
pdf(file=paste("C:\\rOutput\\SP500.TradingPlot",dates,divergence,"SD",lookbackPeriod,"days", day_prior, "prior days.pdf"), onefile=T)

for(i in 1:dim(lookingforward_summary)[1]){
  
  stock1=lookingforward_summary$Stock1[i]
  stock2=lookingforward_summary$Stock2[i]
  
  #Index of Stocks
  getIndex1=which(colnames(closingprice)==stock1)
  getIndex2=which(colnames(closingprice)==stock2)
  
  #Prices of Stocks
  stockPrice1=closingprice[1:(dim(closingprice)[1]),getIndex1]
  stockPrice2=closingprice[1:(dim(closingprice)[1]),getIndex2]
  
  stock1_look = stockPrice1[(dim(closingprice)[1]-lookbackPeriod-day_prior):(dim(closingprice)[1]-day_prior)]
  stock2_look = stockPrice2[(dim(closingprice)[1]-lookbackPeriod-day_prior):(dim(closingprice)[1]-day_prior)]
  
  #Regression between Stock1 and Stock2
  regr= lm(stock1_look ~ stock2_look)
  #Hedge Ratio
  hedgeRatio=regr$coefficients[2]
  
  #Spread with Stock Price Difference
  spread= stock1_look-hedgeRatio*stock2_look
  #Scale the spread
  spread_scale = scale(spread)
  #Mean of Spread
  spreadMean=mean(spread)
  #Std of Spread
  spreadStd= sd(spread)
  
  #Calculate new spread 
  stock1_look_new = stockPrice1[(dim(closingprice)[1]-day_prior+1):(dim(closingprice)[1])]
  stock2_look_new = stockPrice2[(dim(closingprice)[1]-day_prior+1):(dim(closingprice)[1])]
  
  spread_new = stock1_look_new - hedgeRatio*stock2_look_new
  spreadNew = (spread_new - spreadMean)/spreadStd
  
  #Combine old and new
  spread_scale = c(spread_scale, spreadNew)
  
  #Plot the trading
  lookback_date = as.Date(rownames(closingprice)[(dim(closingprice)[1]-lookbackPeriod-day_prior):(dim(closingprice)[1])])
  tableinplot1 = data.frame(Strategy=lookingforward_summary$Strategy[i], HedgeRatio=round(lookingforward_summary$HedgeRatio[i], 4), LookbackPeriod=paste(lookbackPeriod,"days"),
                            position_diff=round(lookingforward_summary$Position_diff[i],4), diff_SD=round(lookingforward_summary$p_diff_sd[i],4))
  tableinplot2 = data.frame(STK_Name=c(lookingforward_summary$Stock1[i],lookingforward_summary$Stock2[i]), Prices=c(lookingforward_summary$stock1Price[i],lookingforward_summary$stock2Price[i]),
                            Short_Long=c(lookingforward_summary$STK1_longshort[i],lookingforward_summary$STK2_longshort[i]))
  rownames(tableinplot2) = c("STK1","STK2")
  
  plot(lookback_date, spread_scale, type="b",main=paste(stock1, "_vs_", stock2,"(",dates,")"),
       col=ifelse(lookback_date==as.Date(row.names(closingprice)[length(rownames(closingprice))-day_prior]),"red","black"),pch=ifelse(lookback_date==dates, 16, 1), cex=ifelse(lookback_date==dates, 1.5, 1))
  abline(h=0,col="blue",lty=2,lwd=2)
  addtable2plot("topleft", table=tableinplot1, display.rownames=F, bty="o", cex=0.6, hlines=T, vlines=T)
  addtable2plot("bottomleft", table=tableinplot2, display.rownames=T, bty="o", cex=0.6, hlines=T, vlines=T)
}
dev.off()

