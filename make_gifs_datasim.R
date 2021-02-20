# Created by Adam Parker, August 2019
# Generates the gifs showing growing beeswarms.
# Modified by Dorothy Bishop, Sept 2019 to add spreadsheet creation
# Option for larger spreadsheet with datasets that are not plotted, but which are used to test strategies

# Gorilla Game is here:  https://gorilla.sc/admin/project/7839#

set.seed(10) #we need gifs to be reproducible! Only change seed if you want to make new gifs - 
# and keep a note of this seed if you do change, in case we need to make identical gifs

options(scipen = 999) #turn off scientific notation
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
library(gtools) # for chr


################################################################################
#function for split violin plots
################################################################################
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}
################################################################################


################################################################################
# function to convert number to letter pair for values > 26 (like Excel columns)
################################################################################
num2alpha<-function(n){
  i<-as.integer((n-1)/26)
  pre<-''
  if(i>0){
    pre<-chr(64+i)
  }
  post<-chr(64+n%%26)
if(post=='@'){post<-'Z'}
  myletters<-paste0(pre,post)
  return(myletters)
}
################################################################################

################################################################################
# function for computing log likelihood of True vs Null using pnorm
################################################################################
makeLL <- function(ES,ESobs,Nsample){

  grain=.025 # grain determines step size when computing probabilities : we take values of ESobs+/- grain/2
#ES is effect size of true effect
#ESobs is observed effect size

# Nsample is sample size per group

  se<-1/sqrt(Nsample)
  
  Nnull <- Nsample*(pnorm((ESobs+(grain/2)),0,se,lower.tail=T)-pnorm((ESobs-(grain/2)),0,se,lower.tail=T))
  #use range of observed effect for each value of myint +/- .05
  Nexp <-  Nsample*(pnorm((ESobs+(grain/2)),ES,se,lower.tail=T)-pnorm((ESobs-(grain/2)),ES,se,lower.tail=T))
  
  LL <- log(Nexp/Nnull)
  return(LL)
}

################################################################################


# simulate data ================================================================================================

gifmake<-1 #set to 1 for creating gifs for study; 0 if making lots of datasets for testing effectiveness of different resonsestrategies
#Elist<-c(3,5,8,0) #effect size x 10 - use this if you want ES of .3, .5 and .8
#Here we'll just have one effect size, but script can handle multiple ESs. Always include zero so we have null cases included.

Elist<-c(3,0) #comment out if you want all 3 effect sizes

#nsets<-c(10,10,10,30) # N samples at each effect size for when creating gifs for Gorilla
nsets<-c(40,40)  #comment out if you want all 3 effect sizes - here we have the N for true and null simulations
blocklength<-20 #sum of nsets should be divisible exactly by blocklength - this is N trials per block
#When we create the gorilla spreadsheet, a divider will be put in at each block

# for testing ================================================================================================
dotest<-0 #set to 1 to just create a small number of gifs for testing; reset to zero to make full set
if(dotest==1){
nsets<-c(2,2) #for testing
blocklength<-2}
#set dotest to zero once testing complete.
# ================================================================================================

#Set parameters for making many sets so can check accuracy rates of different strategies : we won't get gifs for these
if(gifmake==0){
nsets<-c(180,180) #N sets at each effect size - here doing lots
}

N_size <- c(20, 40, 80, 160, 320,640, 20000) #sample sizes to use - this is total for both groups, so halve for sample size per condition
allN<-length(N_size)-1 #ignore final N_size which is used for feedback graph
thistrial<-0 #initialise counter 

#make dataframe to save observed eff size, LL and t for each N for each run
mybf<-data.frame(matrix(nrow=(sum(nsets)*allN),ncol=14))
 names(mybf)<-c('run','lettercode','Nsize','level','trueES','obsES','t','LL','BF','BF2','meanE','sdE','meanC','sdC')
 #make spreadsheet for Gorilla
myspread<-data.frame(matrix(nrow=(3+length(nsets)+sum(nsets)),ncol=58)) #rows need to allow for intro and breaks
names(myspread)<-c('randomise_blocks','randomise_trials','display','ANSWER','image1','image1f','ES','EarningCorrect','EarningWrong',
                   'skip_ref',
                   'meanC1','sdC1','meanE1','sdE1',
                   'meanC2','sdC2','meanE2','sdE2',
                   'meanC3','sdC3','meanE3','sdE3',
                   'meanC4','sdC4','meanE4','sdE4',
                   'meanC5','sdC5','meanE5','sdE5',
                   'meanC6','sdC6','meanE6','sdE6',
                   'ObsE1','ObsE2','ObsE3','ObsE4','ObsE5','ObsE6',
                   'LL1','LL2','LL3','LL4','LL5','LL6',
                   't1','t2','t3','t4','t5','t6',
                   'BF1','BF2','BF3','BF4','BF5','BF6')
spreadrowcount <- 2 #keep count of row number for writing to spreadsheet - will start at row 3, after demo and instructions
myspread$display<-'beeswarms' #default
#add basic info to the first two rows of spreadsheet
myspread$display[1]<-'instructions'
myspread$display[2]<-'demo'
myspread$image1[1:2]<-'E0_A.gif'
myspread$image1f[1:2]<-'E0_A.png'
myspread$EarningCorrect<-10
myspread$EarningWrong<- (-10)
myspread$ANSWER <-'Blue=Pink' #default
myblock<-1
blockcount<-0


#make dataframe to save simulated data for each run
N= N_size[length(N_size)]/2 #biggest sample size to simulate, half will be True and half Null
mydf<-data.frame(matrix(nrow=(N*2),ncol=2))
colnames(mydf)<-c('Group','Score')
mydf$Group <- 'Experimental'
mydf$Group[1:N] <- 'Control'
mydf$col<-2 #red and blue dots
mydf$col[(N+1):(N*2)]<-4

#start loop here
for (effsize in 1:length(Elist)){
  E<-Elist[effsize]/10 #true effect size in large population
  ntrial<-nsets[effsize] 
  for (i in 1:ntrial){
    thistrial<-thistrial+1 #global counter across effect sizes, to use in mybf dataframe
    
    spreadrowcount<-spreadrowcount+1 #row number for gorilla spreadsheet
    blockcount<-blockcount+1 #counter that just resets when length of block is exceeded
    if(blockcount>blocklength)
       {blockcount<-0 #reset counter and add a row that specifies a break
         myblock<-myblock+1
         myspread[spreadrowcount,]<-NA
         myspread$display[spreadrowcount]<-'break'
         spreadrowcount<-spreadrowcount+1 }
    myspread$skip_ref[spreadrowcount]<-myblock #this indicates which block
    myspread$ES[spreadrowcount]<-E
    if(E>0)
    {myspread$ANSWER[spreadrowcount] <-'Blue>Pink'}
    
    
  #for each effect size, create giant population dataset in mydf, with half Exp and half Control

mydf$Score <- rnorm((N*2),mean=E,sd=1) #take random normal sample with mean of E
mydf$Score[1:N] <- rnorm(N,mean=0,sd=1) #take random normal sample with mean of zero
mydf$count <- seq(1, 2*N, by = 2) # assign all odd numbers
mydf$count[1:N] <- seq(2, 2*N, by = 2) # assign control even numbers


#identify range of rows to write to for mybf
bfstartrow<-(thistrial-1)*allN+1
bfendrow<-thistrial*allN
mybf[bfstartrow:bfendrow,1]<-thistrial
mybf[bfstartrow:bfendrow,2]<-num2alpha(i) #use previously defined formula (A-Z and then AA AB etc; allows us to go beyond Z)
mybf[bfstartrow:bfendrow,3]<-N_size[1:(length(N_size)-1)]/2 #this is N per condition (ie half N_size)
mybf[bfstartrow:bfendrow,4]<-1:(length(N_size)-1) #level - ie not actual number, but code for sample size from 1;6
mybf[bfstartrow:bfendrow,5]<-E



# create data set which includes accumulating data for each level
datalist = list()
jcounter<-0
  for (j in (N_size)){ #j/2 is n per group
    jcounter<-jcounter+1
      mydfL <- subset(mydf, count <= (j))
      mydfL$n <- j/2
      datalist[[j]] <- mydfL
      if(j<(max(N_size))){
        ESobs<-(-1)*effsize::cohen.d(mydfL$Score~as.factor(mydfL$Group))$estimate
        mybf[(bfstartrow-1+jcounter),6]<- ESobs
        #multiply by -1 bcs factor takes control first, so all eff sizes are neg
        # NB need to specify cohen.d from effsize package - otherwise conflict with psych package
        
        bf = ttestBF(x = mydfL$Score[1:(j/2)],y=mydfL$Score[((j/2)+1):j], paired=FALSE,nullInterval=c(-Inf,0)) # Bayes Factor from t-test
        #I have difficulty understanding this. It is on raw rather than log scale - directional test gives 2 BFs?
        myt<-t.test(mydfL$Score ~ mydfL$Group)
        mybf[(bfstartrow-1+jcounter),7]<-(-myt$statistic)
        mybf[(bfstartrow-1+jcounter),8]<-makeLL(Elist[1]/10,ESobs,(j/2))
        mybf[(bfstartrow-1+jcounter),9]<-as.data.frame(bf)$bf[1] #this is the only way to get into the bf output!
        mybf[(bfstartrow-1+jcounter),10]<-as.data.frame(bf)$bf[2] #2 BFs generated for directional test
        mybf[(bfstartrow-1+jcounter),11]<-myt$estimate[2]
        mybf[(bfstartrow-1+jcounter),13]<-myt$estimate[1]
        mybf[(bfstartrow-1+jcounter),12]<-sd(mydfL$Score[((j/2)+1):j])
        mybf[(bfstartrow-1+jcounter),14]<-sd(mydfL$Score[1:(j/2)])
      }
  }
#for completeness can add p-values to mybf
mybf$p<-1-pt(mybf$t,(2*mybf$Nsize-2))

# add means and sds to spreadsheet
myspread[spreadrowcount,seq(11,34,4)]<-mybf$meanC[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(12,34,4)]<-mybf$sdC[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(13,34,4)]<-mybf$meanE[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(14,34,4)]<-mybf$sdE[bfstartrow:bfendrow]
myspread[spreadrowcount,35:40]<-mybf$obsES[bfstartrow:bfendrow]
myspread[spreadrowcount,41:46]<-mybf$LL[bfstartrow:bfendrow]
myspread[spreadrowcount,47:52]<-mybf$t[bfstartrow:bfendrow]
myspread[spreadrowcount,53:58]<-mybf$BF[bfstartrow:bfendrow]
# bind the data    
big_data = do.call(rbind, datalist)  
big_data$iter <- as.factor(big_data$n)

big_data$cnd_sample <- big_data$n
big_data$cnd_sample <- as.factor(big_data$cnd_sample)

# Create gif ====================================================================================================
if(gifmake==1){
# create plot for animation
  plotdat<-subset(big_data, iter != (max(N_size)/2))
p <- ggplot(plotdat, aes(x= cnd_sample, y= Score, color= Group)) +
  geom_beeswarm(cex = 1.5, dodge.width = .2) + 
  xlab("Sample size per condition") + ylim(-3,3)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.75)  #add line corresponding to mean

p <- p + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top") +
  theme(text = element_text(size=15))
#p
ggsave(paste0('gifs_pngs/static jpg versions/myplot',E,'_',num2alpha(i),'.jpg'),width= 7, height = 7,p)

# set animation values
anim <- p + transition_states(n,
              transition_length = 10,
              state_length = 60) +
            enter_grow() + exit_shrink()
anim <- animate(anim, fps=1, duration=12) #creates gif with 12 frames,2 of each N
# print animation
#anim

# save

filename<-paste0('gifs_pngs/E',Elist[effsize],'_',num2alpha(i),'.gif')
anim_save(filename, anim)
# answer plot ==================================================================================================

truedat <- big_data

# Maintain all visible data for n= 320
truedat <-
  truedat %>%
  mutate(Score2 = ifelse(iter == (max(N_size)/2), Score, -1000000000))

truedat <-
  truedat %>%
  mutate(iter2 = ifelse(iter == (max(N_size)/2), 320, as.numeric(cnd_sample)*10))
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
  geom_hline(yintercept = E, linetype="longdash", color = "#00BFC4", size=1) + 
  geom_split_violin(trim= TRUE, alpha= 0.5) +
  geom_boxplot(width = 0.1, position = position_dodge(width=0.25), outlier.alpha = 0) + 
  xlab("Sample size per condition") + ylab("Score") + ylim(-3,3)
q <- q + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top")  +
  theme(text = element_text(size=15))
q

filename<-paste0('gifs_pngs/E',Elist[effsize],'_',LETTERS[i],'.png')
ggsave(filename , dpi= 300, width= 7, height = 7)  
}



csvname<-'bigBFdata.csv'
if(gifmake==1){csvname<-'gifs_pngs/BFdata.csv'}
#NB when log scale used, the evidence for null is just the same as for H1, but negative
write.csv(mybf,csvname,row.names=FALSE)

myspread$image1f[spreadrowcount]<-paste0('E',Elist[effsize],'_',num2alpha(i),'.png')
myspread$image1[spreadrowcount]<-paste0('E',Elist[effsize],'_',num2alpha(i),'.gif')
myspread$randomise_trials[spreadrowcount]<-1
  }
}
myspread$display[spreadrowcount]<-'end'
write.csv(myspread,'gorilla_spreadsheets/spreadsheet1.csv',row.names=FALSE)

#Some plots to check distributions

#bf has some crazy big values that mess with graphs, so we'll get rid of that
w<-which(mybf$BF > 50000)
mybf$BF[w]<-NA
plot(mybf$t,mybf$obsES,col=as.factor(mybf$Nsize))
abline(h=0)
abline(v=0)
plot(log(mybf$BF),mybf$obsES,col=as.factor(mybf$Nsize))
abline(h=0)
abline(v=0)

plot(mybf$LL,mybf$obsES,col=as.factor(mybf$Nsize),pch=(13+10*mybf$trueES),xlab='Log Likelihood',ylab='Observed Effect size')
abline(h=0)
abline(v=0)
legend(-20, 1.4, legend=(N_size[1:6]/2),
       col=1:6,  pch=15,cex=0.8)
text(-15,1.5,'N per group',cex=.8)

legend(-30, 1.4, legend=c(0,.3),
        pch=c(13,16),cex=0.8)
text(-28,1.5,'True effect size',cex=.8)

plot(mybf$level,mybf$t,col=as.factor(mybf$Nsize))
plot(mybf$level,mybf$LL,col=as.factor(mybf$Nsize))
plot(mybf$level,mybf$BF,col=as.factor(mybf$Nsize))
plot(mybf$level,mybf$obsES,col=as.factor(mybf$Nsize))
plot(log(mybf$BF),mybf$LL,col=as.factor(mybf$Nsize))
abline(h=0)
abline(v=0)   #The first BF from the direction BF t-test is similar to my computed LL, but scale is compacted

# mybf$BFratio <- mybf$BF/mybf$BF2
# plot(log(mybf$BFratio),mybf$LL,col=as.factor(mybf$Nsize))
# abline(h=0)
# abline(v=0)  #ratio of the 2 BF doesn't appear meaningful here