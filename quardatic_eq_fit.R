# Fitting 2nd order quardatic equation through time and parameter set
# coeff.final = data.frame(variable=as.factor(rep(c("k","Y","af","as","sf"),length(vol))), 
#                    c1=numeric(no.var*length(vol)), c2=numeric(no.var*length(vol)), c3=numeric(no.var*length(vol)))

for (j in 1:length(vol)) {
  coeff = data.frame(volume = as.factor(rep(vol[j],no.var)), variable=as.factor(c("k","Y","af","as","sf")), 
                     c1=numeric(no.var), c2=numeric(no.var), c3=numeric(no.var))
  parameter = subset(summary.param,volume==vol[j])
  keeps <- c('Date', 'variable','Parameter')
  parameter = parameter[ , keeps, drop = FALSE]
  parameter = dcast(parameter, Date ~ variable, value.var="Parameter")
  keeps <- c("Date","k","Y","af","as","sf")
  parameter = parameter[ , keeps, drop = FALSE]
  
  for (i in 1:no.var) {
    model.fit = data.frame(parameter$Date,parameter[,i+1])
    names(model.fit) = c("Date","parameter")
    model.fit$Date = as.Date(model.fit$Date)
    model.fit$Date = model.fit$Date[]-model.fit$Date[1]+1
    # model.fit$Date = model.fit$Date[]-model.fit$Date[1]+ceiling(0.5*(model.fit$Date[2]-model.fit$Date[1]))
    
    x = model.fit$Date
    y = model.fit$parameter
    fit <- lm( y~poly(x,2) )
    coefficients(fit) # model coefficients
    coeff[i,(ncol(coeff)-2):ncol(coeff)] = coefficients(fit)
    plot(x,y,xlim=c(0,121),ylim=c(0,1),main=paste(coeff$variable[i]))
    xx <- seq(1,121, length.out=121)
    yy = predict(fit, data.frame(x=xx))
    lines(xx, yy, col='blue')
  }
  if (j == 1) {
    coeff.final = coeff
  }
  if (j > 1) {
    coeff.final = rbind(coeff.final,coeff)
  }
}

c3 = dcast(coeff.final, volume ~ variable, value.var="c3")
c3[nrow(c3)+1, 2:ncol(c3)] = apply(c3[,2:ncol(c3)],2,min)
c3[nrow(c3)+1, 2:ncol(c3)] = apply(c3[,2:ncol(c3)],2,max)
c3 = c3[(nrow(c3)-1):nrow(c3),2:ncol(c3)]
dimnames(c3)[[1]] <- c("Min", "Max")

write.csv(c3, file = "rawdata/c3.csv", row.names = FALSE)


