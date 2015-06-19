
modelMat = designCombn(matrix(1,1,len(grep('X',colnames(dieDesign))))) 
modelMat = rbind(rep(0,ncol(modelMat)),modelMat)
modelMat = cbind(matrix(1,nrow(modelMat),ncol= len(grep('X',colnames(dieDesign),invert=T))),modelMat)

IC = matrix(NA,nrow(modelMat),2)
colnames(IC) = c('BIC','AIC')
rownames(IC) = seq(1,nrow(modelMat))
for(i in 1:nrow(modelMat)){
  
  whr1 = which(modelMat[i,] == 1)
  Mdesign = dieDesign[,whr1]
  fit = glmer(dieResp ~ Mdesign - 1 + (1|dieIndiv), family = binomial)
  IC[i,1] = BIC(fit) #beyondS and beyondN = 1 out of range 
  IC[i,2] = AIC(fit) 
  
}

num1 = a.n(names(sort(apply(apply(IC,2,rank),1,mean)))[1])

whr1 = which(modelMat[num1,] == 1)
Mdesign = dieDesign[,whr1]
fit = glmer(dieResp ~ Mdesign - 1 + (1|dieIndiv), family = binomial)
summary(fit)

coes = coefficients(summary(fit))[,1]
coeSD = coefficients(summary(fit))[,2]

prob = function(y) 1/(1+exp(-y))

dummyMat = designCombn(matrix(1,1,len(coes)-1)) 
dummyMat = rbind(rep(0,ncol(dummyMat)),dummyMat)
dummyMat = cbind(rep(1,nrow(dummyMat)),dummyMat)
colnames(dummyMat) = colnames(Mdesign)

ys = dummyMat %*% coes
probDb = prob(ys)
names(probDb) = apply(dummyMat,1,function(x) paste(names(x)[which(x==1)],collapse='|'))
probDb = probDb[grep('Struct',names(probDb),invert=T)]
probDb = probDb[grep('BeyondS\\|BeyondN',names(probDb),invert=T)]
round(rev(sort(probDb)) *100,2)

treatByInd = list()
for(i in 1:nrow(dummyMat)){
  treatByInd[[i]] = matchAllCol(Mdesign,dummyMat[i,]) 
  if(len(matchAllCol(Mdesign,dummyMat[i,])) == 0) treatByInd[[i]] = NA
}
#rans = unlist(ranef(fit))[fit@pp$Zt@i] #fit@pp$Zt@i #individuals
rans = unlist(ranef(fit))[match(dieIndiv,rownames(fit@pp$Zt))]
samples = 10000#samples = unlist(lapply(treatByInd,length))

nsam = 10000
probsMat = matrix(NA,nsam,len(ys))
colnames(probsMat) = apply(dummyMat,1,function(x) paste(names(x)[which(x==1)],collapse='|'))
whr = grepl('Struct',colnames(probsMat)) == FALSE & grepl('BeyondS\\|BeyondN',colnames(probsMat)) == FALSE
probsMat = probsMat[,whr]
dummyMatT = dummyMat[whr,]
for(i in 1:nsam){ 
  tmpB = mvrnorm(1,coes,sqrtMat(vcov(fit)))
  tmp = dummyMatT %*% tmpB + rans[unlist(lapply(treatByInd,sample,1)[whr])]
  probsMat[i,] = prob(rnorm(len(tmp),tmp,sd=fit@theta))
  progress(i,nsam,100)
}
results = apply(probsMat,2,quantile,c(.025,.5,.975))*100

resp = rbind(results[,c(1,2,3,7,8,9)],results[,c(4,5,6,10,11,12)])
rownames(resp) = paste(rep('db|',3),rownames(tmp),sep="")

mods = c('Int|BeyondS','Int','Int|BeyondN',
         'Int|Duke|BeyondS','Int|Duke','Int|Duke|BeyondN')
names(mods) = c('HF|S','HF','HF|N',
                'DF|S','DF','DF|N')
probsMat = probsMat[,match(mods,colnames(probsMat))]
names(probsMat) = names(mods)
plot(density(probsMat[,1]*100))
t.test(x = probsMat[,1], y = probsMat[,2],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
t.test(x = probsMat[,1], y = probsMat[,3],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
t.test(x = probsMat[,2], y = probsMat[,3],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)

boxplot(c(probsMat[,1],probsMat[,2],probsMat[,3])~rep(c('1.S','2.IN','3.N'),each = nsam),
        notch=T,horizontal = F)
title('HF')

t.test(x = probsMat[,4], y = probsMat[,5],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
t.test(x = probsMat[,4], y = probsMat[,6],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)
t.test(x = probsMat[,5], y = probsMat[,6],alternative = c("two.sided"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95)

boxplot(c(probsMat[,4],probsMat[,5],probsMat[,6])~rep(c('1.S','2.IN','3.N'),each = nsam)
        ,notch=T,horizontal = F)
title('DF')

overlap(probsMat[,1],probsMat[,2])
overlap(probsMat[,2],probsMat[,3])
overlap(probsMat[,4],probsMat[,5])
overlap(probsMat[,5],probsMat[,6])

overlap(probsMat[,2],probsMat[,1])
overlap(probsMat[,3],probsMat[,2])
overlap(probsMat[,5],probsMat[,4])
overlap(probsMat[,6],probsMat[,5])




####growth

modelMat = designCombn(matrix(1,1,len(grep('X',colnames(design))))) 
modelMat = rbind(rep(0,ncol(modelMat)),modelMat)
modelMat = cbind(matrix(1,nrow(modelMat),ncol= len(grep('X',colnames(design),invert=T))),modelMat)

IC = matrix(NA,nrow(modelMat),2)
colnames(IC) = c('BIC','AIC')
rownames(IC) = seq(1,nrow(modelMat))
for(i in 1:nrow(modelMat)){
  
  whr = which(modelMat[i,] == 1)
  Mdesign = design[,whr]
  fit = lmer(Grow ~ Mdesign - 1 + (1 |Indiv),REML=T)
  IC[i,1] = BIC(fit) #beyondS and beyondN = 1 out of range 
  IC[i,2] = AIC(fit) 
  #IC[i,3] = fit@devcomp$cmp['REML']
  
}

num1 = a.n(names(sort(apply(apply(IC,2,rank),1,mean)))[1])

whr1 = which(modelMat[num1,] == 1)
Mdesign = design[,whr1]
fit = lmer(Grow ~ Mdesign - 1 + (1 |Indiv),REML=T)
summary(fit)

coes = coefficients(summary(fit))[,1]
coeSD = coefficients(summary(fit))[,2]

dummyMat = designCombn(matrix(1,1,len(coes)-1)) 
dummyMat = rbind(rep(0,ncol(dummyMat)),dummyMat)
dummyMat = cbind(rep(1,nrow(dummyMat)),dummyMat)
colnames(dummyMat) = colnames(Mdesign)


ys = dummyMat %*% coes
probDb = ys
names(probDb) = apply(dummyMat,1,function(x) paste(names(x)[which(x==1)],collapse='|'))
probDb = probDb[grep('struct',names(probDb),invert=T)]
probDb = probDb[grep('beyond.S\\|beyond.N',names(probDb),invert=T)]
probDb = probDb[grep('diam',names(probDb),invert=T)]
round(rev(sort(probDb)),2)

tmp = Mdesign
tmp[,'diam'] = dummyMat[,'diam'] =  0
tmp[,'struct'] = dummyMat[,'struct']  = 0
treatByInd = list()
for(i in 1:nrow(dummyMat)){
  treatByInd[[i]] = matchAllCol(tmp,dummyMat[i,]) 
  if(len(matchAllCol(tmp,dummyMat[i,])) == 0) treatByInd[[i]] = NA
}
#rans = unlist(ranef(fit))[fit@pp$Zt@i] #fit@pp$Zt@i #individuals

rans = unlist(ranef(fit))[match(Indiv,rownames(fit@pp$Zt))]
samples = 10000#samples = unlist(lapply(treatByInd,length))

mods = c('int|beyond.S','int','int|beyond.N',
         'int|gap|beyond.S|gapXbeyond.S','int|gap','int|gap|beyond.N',
         'int|duke|beyond.S|dukeXbeyond.S','int|duke','int|duke|beyond.N',
         'int|gap|duke|beyond.S|dukeXbeyond.S|gapXbeyond.S','int|gap|duke','int|gap|duke|beyond.N')
names(mods) = c('HF|S','HF','HF|N','HF|S|Gap','HF|Gap','HF|N|Gap',
                'DF|S','DF','DF|N','DF|S|Gap','DF|Gap','DF|N|Gap')

#tapply(rans,apply(tmp,1,function(x) paste(x,collapse="")),mean)
nsam = 10000
probsMat = matrix(NA,nsam,len(ys))
colnames(probsMat) = apply(dummyMat,1,function(x) paste(names(x)[which(x==1)],collapse='|'))
whr = match(mods,colnames(probsMat))
probsMat=probsMat[,whr]
dummyMatT = dummyMat[whr,]
betaSam = matrix(NA,nsam,5)
for(i in 1:nsam){
  tmpB = mvrnorm(1,coes,sqrtMat(vcov(fit)))
  betaSam[i,] = c(tmpB[c(6,7)],tmpB[6]+tmpB[8],tmpB[6]+tmpB[8]+tmpB[9],tmpB[6]+tmpB[9])
  tmp = dummyMatT %*% tmpB + rans[unlist(lapply(treatByInd,sample,1)[whr])]
  tmp[tmp < 0] = 0
  probsMat[i,] = tnorm(len(tmp),0,Inf,tmp,fit@theta)
  progress(i,nsam,100)
}
results = apply(probsMat,2,quantile,c(.025,.5,.975),na.rm=T)
colnames(results) = names(mods)
apply(betaSam,2,quantile,c(.025,.5,.975),na.rm=T)
tmp = rbind(results[,c(1,2,3,7,8,9)],results[,c(4,5,6,10,11,12)])
rownames(tmp) = paste(c(rep('gr|',3),rep('gr|G|',3)),rownames(tmp),sep="")
resp = rbind(resp,tmp)

par(mfrow=c(2,2))
boxplot(c(probsMat[,match(mods,colnames(probsMat))][,1],probsMat[,match(mods,colnames(probsMat))][,2],probsMat[,match(mods,colnames(probsMat))][,3]) ~ 
          rep(c('1.S','2.IN','3.N'),each = 10000),notch=T)
title('HF')
boxplot(c(probsMat[,match(mods,colnames(probsMat))][,4],probsMat[,match(mods,colnames(probsMat))][,5],probsMat[,match(mods,colnames(probsMat))][,6])~
          rep(c('1.S','2.IN','3.N'),each = 10000),notch=T)
title('HF Gap')
boxplot(c(probsMat[,match(mods,colnames(probsMat))][,7],probsMat[,match(mods,colnames(probsMat))][,8],probsMat[,match(mods,colnames(probsMat))][,9])~
          rep(c('1.S','2.IN','3.N'),each = 10000),notch=T)
title('DF')
boxplot(c(probsMat[,match(mods,colnames(probsMat))][,10],probsMat[,match(mods,colnames(probsMat))][,11],probsMat[,match(mods,colnames(probsMat))][,12])~
          rep(c('1.S','2.IN','3.N'),each = 10000),notch=T)
title('DF Gap')



resp = read.csv(paste(path,'/SimpleLimits/response.csv',sep=''),header=T)
rownames(resp) = resp[,1]
resp = resp[,-1]


allA.lo = t(cbind(t(resp[4,1:3]),t(resp[4,4:6]),t(resp[7,1:3]),t(resp[7,4:6])))
allA = t(cbind(t(resp[5,1:3]),t(resp[5,4:6]),t(resp[8,1:3]),t(resp[8,4:6])))
allA.hi = t(cbind(t(resp[6,1:3]),t(resp[6,4:6]),t(resp[9,1:3]),t(resp[9,4:6])))

allR.lo = t(cbind(t(resp[1,1:3]),t(resp[1,4:6]),t(resp[1,1:3]),t(resp[1,4:6])))
allR = t(cbind(t(resp[2,1:3]),t(resp[2,4:6]),t(resp[2,1:3]),t(resp[2,4:6])))
allR.hi = t(cbind(t(resp[3,1:3]),t(resp[3,4:6]),t(resp[3,1:3]),t(resp[3,4:6])))


par(mfrow=c(2,1),mar=c(3, 4, 2, 2), oma=c(0,0,0,0),xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
kk=1
for(k in c('allA','allR')){
  
  co = get(k); co.lo = get(paste(k,'.lo',sep='')); co.hi = get(paste(k,'.hi',sep=''))
  numR = nrow(co)
  numC = ncol(co)
  num = (numR+1)*ncol(co) + 1
  step = 1/num
  ys = step
  minX = min(c(co.lo,co.hi))
  maxX = max(c(co.lo,co.hi))
  minY = 0; maxY = 1
  xTit = ''; yTit=''
  
  Title = c('Growth','Dieback','Germination','Mortaility')[kk]
  kk=1+kk
  plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, 
       main='', xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), 
       family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
  for(i in seq(0,ncol(co),by=2)){
    rect(minX-10,i*numR*step,
         maxX+10,i*numR*step +  numR*step,
         col = rgb(.95,.95,.95), border=NA
    )
    
  }
  xs = rep(seq(step,numC*step,by = step),nrow(co)) + rep(seq(0,1-step*numC,step*numR),each=ncol(co))
  title(bquote(bold(.(Title))),line=.9)
  axis(4,at = xs,labels=rep(c('South','Within','North'),nrow(co)),tick=F,cex.axis=.55,hadj = .6)
  axis(1,tick=F,cex = .9);#axis(3,tick=F,cex = .9)
  axis(2,at = seq(numR*step,1,by = numR*step),
       labels=c('HF','DF','HF|Gap','DF|Gap'),
       tick=F,cex.axis=.9,hadj = 0,padj=1.1,adj = 0,line=-.5)
  
  abline(h = seq(0,1,numR*step),lty=3)
  
  abline(v=0,lty=2,col='grey')
  pty1 = c(0,1,2,5)
  pty2 = c(22,21,24,23)
  for(j in 1:nrow(co)){
    if(j %% 2 == 0){b.col = 'white'}else{b.col =  rgb(.95,.95,.95)}
    for(i in 1:numC){
      
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col= b.col,lwd = 8)
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col=c(3,2,1)[i],lwd = 1.5)
      points(co[j,i],ys,col=c(3,2,1)[i],cex=.75,pch=pty2[i],bg=c(3,2,1)[i])
      
      ys = ys + step  
    }
    ys = ys + step  
  }
}

Survival = list.files(paste(path,'/SimpleLimits',sep=''),full.names=T)
load(Survival[grep("fittedSurv",Survival)])

