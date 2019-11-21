# Created by Adam Parker, August 2019
# Generates the gifs showing growing beeswarms.
# Modified by Dorothy Bishop, Sept 2019
# Can generate datasets that are not plotted, but which are used to test strategies


# load packages ================================================================================================
library("ggplot2")
library("gganimate")
library("beeswarm")
library("dplyr")
library("animation")
library("FSA")
library("data.table")
library("ggbeeswarm")
library("effsize")
library("BayesFactor")

source("split_vio.R")

# simulate data ================================================================================================
set.seed(5)#make data reproducible
options(scipen = 999) #turn off sci notation
gifmake<-0 #set to 1 for creating gifs for study; 0 if making lots of datasets for testing strategies
Elist<-c(3,5,8,0) #effect size x 10

nsets<-c(6,6,6,18) #samples when doing gifs
if(gifmake==0){
nsets<-c(60,60,60,180) #N sets at each effect size - here doing lots so can check accuracy rates of differnt strategies
}

N_size <- c(20, 40, 80, 160, 320, 640, 20000) #sample sizes to use
allN<-length(N_size)-1
thistrial<-0 #initialise counter

#make dataframe to save observed eff size and bayes factor for each N for each run
mybf<-data.frame(matrix(nrow=(sum(nsets)*allN),ncol=7))
names(mybf)<-c('run','lettercode','Nsize','trueES','obsES','t','BF')

#make dataframe to save simulated data for each run
N= 10000
mydf<-data.frame(matrix(nrow=(N*2),ncol=2))
colnames(mydf)<-c('Group','Score')
mydf$Group <- 'Experimental'
mydf$Group[1:N] <- 'Control'
mydf$col<-2 #red and blue dots
mydf$col[1:N]<-4

#start loop here
for (effsize in 1:length(Elist)){
  E<-Elist[effsize]/10
  ntrial<-nsets[effsize]
  for (i in 1:ntrial){
    thistrial<-thistrial+1 #global counter across effect sizes, to use in mybf dataframe
  #for each effect size, create giant population dataset in mydf, with half Exp and half Control

mydf$Score <- rnorm((N*2),mean=E,sd=1) #take random normal sample with mean of E
mydf$Score[1:N] <- rnorm(N,mean=0,sd=1) #take random normal sample with mean of zero
mydf$count <- seq(1, 20000, by = 2) # assign all odd numbers
mydf$count[1:N] <- seq(2, 20000, by = 2) # assign control even numbers




#identify range of rows to write to for mybf
bfstartrow<-(thistrial-1)*allN+1
bfendrow<-thistrial*allN
mybf[bfstartrow:bfendrow,1]<-thistrial
mybf[bfstartrow:bfendrow,2]<-LETTERS[i]
mybf[bfstartrow:bfendrow,3]<-N_size[1:(length(N_size)-1)]
mybf[bfstartrow:bfendrow,4]<-E


# create data set which includes accumulating data for each level
datalist = list()
jcounter<-0
  for (j in N_size){
    jcounter<-jcounter+1
      mydfL <- subset(mydf, count <= (j))
      mydfL$n <- j
      datalist[[j]] <- mydfL
      if(j<max(N_size)){
        mybf[(bfstartrow-1+jcounter),5]<- (-1)*cohen.d(mydfL$Score ~ as.factor(mydfL$Group))$estimate
        #multiply by -1 bcs factor takes control first, so all eff sizes are neg
        
        bf = ttestBF(x = mydfL$Score[1:(j/2)],y=mydfL$Score[((j/2)+1):j], paired=FALSE)
        myt<-t.test(x=mydfL$Score[((j/2)+1):j],y = mydfL$Score[1:(j/2)], paired=FALSE)
        mybf[(bfstartrow-1+jcounter),6]<-myt$statistic
        mybf[(bfstartrow-1+jcounter),7]<-as.data.frame(bf)$bf
      }
  }
# bind the data    
big_data = do.call(rbind, datalist)  
big_data$iter <- as.factor(big_data$n)

big_data$cnd_sample <- big_data$n/2
big_data$cnd_sample <- as.factor(big_data$cnd_sample)

# Create gif ====================================================================================================
if(gifmake==1){
# create plot for animation
p <- ggplot(subset(big_data, iter != 20000), aes(x= cnd_sample, y= Score, color= Group)) +
  geom_beeswarm(cex = 1.5, dodge.width = .2) + 
  xlab("Sample size per condition") + ylim(-3,3)
p <- p + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top") +
  theme(text = element_text(size=15))
#p


# set animation values
anim <- p + transition_states(n,
              transition_length = 10,
              state_length = 60) +
            enter_grow() + exit_shrink()
anim <- animate(anim, fps=1, duration=12) #creates gif with 12 frames,2 of each N
# print animation
#anim

# save

filename<-paste0('gifs_pngs/E',Elist[effsize],'_',LETTERS[i],'.gif')
anim_save(filename, anim)
# answer plot ==================================================================================================

truedat <- big_data

# Maintain all visible data for n= 320
truedat <-
  truedat %>%
  mutate(Score2 = ifelse(iter == 20000, Score, -1000000000))

truedat <-
  truedat %>%
  mutate(iter2 = ifelse(iter == 20000, 320, as.numeric(cnd_sample)*10))
truedat$iter2 <- as.factor(truedat$iter2)

# create means for horiztonal lines
d_dat <- subset(truedat, c(iter2=="320"))
meanC<- subset(truedat, c(Group=="Control", iter=="320"))
meanE<- subset(truedat, c(Group=="Experimental", iter=="320"))
meanC<-mean(meanC$Score, na.rm= TRUE)
meanE<-mean(meanE$Score, na.rm= TRUE)

# create answer plot
q <- ggplot(truedat, aes(x= iter2, y= Score2, color= Group, fill= Group)) +
  geom_hline(yintercept = 0, linetype="dotdash", color = "#F8766D", size=1) + 
  geom_hline(yintercept = 0, linetype="longdash", color = "#00BFC4", size=1) + 
  geom_split_violin(trim= TRUE, alpha= 0.5) +
  geom_boxplot(width = 0.1, position = position_dodge(width=0.25), outlier.alpha = 0) + 
  xlab("Sample size per condition") + ylab("Score") + ylim(-3,3)
q <- q + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top")  +
  theme(text = element_text(size=15))
q

filename<-paste0('gifs_pngs/E',Elist[effsize],'_',LETTERS[i],'.png')
ggsave(filename , dpi= 300, width= 7, height = 7)  
}
# # ensure that test matches full data
# testDF <- big_data[big_data$iter == "20000",]
# test_res <- t.test(testDF$Score~testDF$Group)$p.value
# if (test_res < .05) {
#   "Yes, there is a statistical difference (p< .05)"
# } else {
#   "Oh no, this is only a numerical, statistically non-significant difference (p> .05)"
# }



# Clear df ======================================================================================================

#remove(list = ls())
  }
}
#logBF is more useful, so

mybf$BF<-log(mybf$BF)

colnames(mybf)[7]<-c('logBFH1')
csvname<-'bigBFdata.csv'
if(gifmake==1){csvname<-'gifs_pngs/BFdata.csv'}
#NB when log scale used, the evidence for null is just the same as for H1, but negative
write.csv(mybf,csvname,row.names=FALSE)


#Checking computation of Bayes factor for diff between means
xall = mydfL$Score[1:(j/2)] #just using last block of data 
yall = mydfL$Score[((j/2)+1):j]
myN<-5
#take first N of each
x<-xall[1:myN]
y<-yall[1:myN]
myt<-t.test(x,y,alternative='greater')$statistic
bf=ttestBF(x=x,y=y,paired=FALSE)
bf<-as.data.frame(bf)$bf #computed bf for these data
bft<-ttest.tstat(myt, myN, myN, nullInterval = NULL, rscale = "medium",
                 complement = FALSE, simple = FALSE)
heff<-2
diff<-mean(x)-mean(y)
sdiff<-1/sqrt(myN-1)
score0<-(0-diff)/sdiff
score1<-(heff-diff)/sdiff
p0<-pt(score0,lower.tail=TRUE)
p1<-pt(score1,lower.tail=TRUE)
p0<-pcauchy(score0,lower.tail=TRUE)
p1<-pcauchy(score1,lower.tail=TRUE)
#d0<-dt(score0,df=(myN*2-1))
#d1<-dt(score1,df=(myN*2-1))

compBF<-p1/p0
#compBFd<-d1/d0
bf
bft$bf
compBF
log(compBF)

#########2nd attempt 
mydfa<-mybf[1:20,  ]
mydfa$bf<-exp(mydfa$logBFH1)
mydfa$bft<-NA
mydfa$pt0<-mydfa$pt1<-mydfa$OR<-mydfa$logOR<-NA
estes<-1 #this agrees with the other bayes estimate
for (i in 1:20){
  mydfa$bft[i]<-ttest.tstat(mydfa$t[i], mydfa$Nsize[i]/2, mydfa$Nsize[i]/2, nullInterval = NULL, rscale = "medium",
                     complement = FALSE, simple = FALSE)$bf
  mydfa$pt0[i]<-pt(mydfa$t[i],df=(mydfa$Nsize[i]-2),lower.tail=TRUE)
  mydfa$pt1[i]<-pt((estes-mydfa$t[i]),df=(mydfa$Nsize[i]-2),lower.tail=TRUE)
  mydfa$OR[i]<-mydfa$pt0[i]/mydfa$pt1[i]
  mydfa$logOR[i]<-log(mydfa$OR[i])
  
  
}

