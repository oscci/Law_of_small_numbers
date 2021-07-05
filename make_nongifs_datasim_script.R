# Based on make_gifs_datasim_script.
# Dorothy Bishop, 5th July 2021, based on Adam Parker script from 2019
# For new version of task, we need pngs or jpgs of individual beeswarms at each sample size.
# So script is same as before, except that instead of creating gifs from the data, we make many individual plots. 
# The final plot with all display sizes is retained.

# The spreadsheet also needs to be modified for the new version of game

# Gorilla Game is here:  https://gorilla.sc/admin/project/7839# - NEED TO UPDATE THIS

## Explanation

# For the Gorilla game, we need jpgs that show how a dataset changes as additional datapoints are sampled. 
# We also create the spreadsheet that will be used to control display of trials 
# in the Gorilla game, and we will save the data behind the gifs so we can check how participants' responses 
# relate to the observed effect size, etc.

# Stimuli will be presented in blocks, each of which contains a mixture of two effect sizes: 
# either zero (null) or E (true effect). Here we used E = .3, but this can be modified easily.

set.seed(10) #we need stimuli to be reproducible! Only change seed if you want to make new set of stimuli - 
# and keep a note of this seed if you do change, in case we need to revert to originals

options(scipen = 999) #turn off scientific notation
# load packages ================================================================================================
library("ggplot2")
library("gganimate")
library("beeswarm")
library("dplyr")

library("FSA")
library("data.table")
library("ggbeeswarm")
library("effsize")
library("BayesFactor")
library(gtools) # for chr
library(permute) #for shuffling effect sizes


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
# Function to convert number to letter pair for values > 26 (like Excel columns)
# These are used for labelling the stimuli presented on a given trial
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
# This is not used in the game, but may be useful for exploring response strategies
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

jpgmake<- 1 #set to 1 for creating jpgs for study; set to 2 to just make gorilla worksheet, 0 to make many spreadsheets
#0 if making lots of datasets for testing effectiveness of different response strategies

#Elist<-c(3,5,8,0) #effect size x 10 - use this if you want ES of .3, .5 and .8
#Here we'll just have one effect size, but script can handle multiple ESs. Always include zero so we have null cases included.

Elist<-c(3,0) #comment out if you want all 3 effect sizes

#can here specify earnings for correct or error responses
corrEarn <- 8
wrongEarn <- (-4)

nsets<-c(40,40)  # here we have the N stimulus sets for true and null simulations
blocklength<-20 #sum of nsets should be divisible exactly by blocklength - this is N trials per block
#When we create the gorilla spreadsheet, a divider will be put in at each block
Nblock<-sum(nsets)/blocklength
# for testing ================================================================================================
dotest<-1 #set to 1 to just create a small number of jpgs for testing; reset to zero to make full set
if(dotest==1){
  nsets<-c(4,4) #for testing
  blocklength<-4}
#set dotest to zero once testing complete.
# ================================================================================================

#Set parameters for making many sets so can check accuracy rates of different strategies : we won't get gifs for these
if(jpgmake==0){
  nsets<-c(180,180) #N sets at each effect size - here doing lots
}
# ================================================================================================

N_size <- c(20, 40, 80, 160, 320,640, 20000) #sample sizes to use - this is total for both groups, so halve for sample size per condition
allN<-length(N_size)-1 #ignore final N_size which is used for feedback graph only
thistrial<-0 #initialise counter 

#make dataframe to save observed eff size, LL and t for each N for each run
dfsummary<-data.frame(matrix(nrow=(sum(nsets)*allN),ncol=14))
names(dfsummary)<-c('run','lettercode','Nsize','level','trueES','obsES','t','LL','BF','BF2','meanE','sdE','meanC','sdC')
#make spreadsheet for Gorilla
myspread<-data.frame(matrix(nrow=(2+Nblock+sum(nsets)),ncol=64)) #rows need to allow for intro and breaks
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
                   'BF1','BF2','BF3','BF4','BF5','BF6',
                   'olap1','olap2','olap3','olap4','olap5','olap6')
myspread$display<-'beeswarms' #default
#add basic info to the first two rows of spreadsheet
myspread$display[1]<-'instructions'
myspread$display[2]<-'demo'
myspread$image1[1:2]<-'E5_demo.gif' #predetermined gifs for the demo
myspread$image1f[1:2]<-'E5_demo.png'
myspread$EarningCorrect<-NA
myspread$EarningWrong<- NA
myspread$ANSWER <-'Blue=Pink' #default



#make dataframe to save simulated data for each run
N= N_size[length(N_size)]/2 #biggest sample size to simulate, half will be True and half Null
mydf<-data.frame(matrix(nrow=(N*2),ncol=2))
colnames(mydf)<-c('Group','Score')
mydf$Group <- 'Experimental'
mydf$Group[1:N] <- 'Control'
mydf$col<-2 #red and blue dots
mydf$col[(N+1):(N*2)]<-4

myblock<-1
blockcount<--1
spreadrowcount <- 2 #keep count of row number for writing to spreadsheet - will start at row 3, after demo and instructions

#start loop here
for (effsize in 1:length(Elist)){
  E<-Elist[effsize]/10 #true effect size in large population
  ntrial<-nsets[effsize] 
  for (i in 1:ntrial){
    thistrial<-thistrial+1 #global counter across effect sizes, to use in dfsummary dataframe
    
    spreadrowcount<-spreadrowcount+1 #row number for gorilla spreadsheet
    blockcount<-blockcount+1 #counter that just resets when length of block is exceeded
    if(blockcount==blocklength)
    {blockcount<-0 #reset counter and add a row that specifies a break
          
          myspread[spreadrowcount,]<-NA #NA in all columns
          if(myblock < Nblock){
          myspread$display[spreadrowcount]<-'break'
          }
         spreadrowcount<-spreadrowcount+1 #advance to next row
          myblock<-myblock+1
    }
           
    myspread$skip_ref[spreadrowcount]<-myblock #this indicates which block
    myspread$ES[spreadrowcount]<-E
    if(E>0)
    {myspread$ANSWER[spreadrowcount] <-'Blue>Pink'}
    
    
    #for each effect size, create giant population dataset in mydf, with half Exp and half Control
    
    mydf$Score <- rnorm((N*2),mean=E,sd=1) #take random normal sample with mean of E
    mydf$Score[1:N] <- rnorm(N,mean=0,sd=1) #take random normal sample with mean of zero
    mydf$count <- seq(1, 2*N, by = 2) # assign all odd numbers
    mydf$count[1:N] <- seq(2, 2*N, by = 2) # assign control even numbers
    
    
    #identify range of rows to write to for dfsummary
    bfstartrow<-(thistrial-1)*allN+1
    bfendrow<-thistrial*allN
    dfsummary[bfstartrow:bfendrow,1]<-thistrial
    dfsummary[bfstartrow:bfendrow,2]<-num2alpha(i) #use previously defined formula (A-Z and then AA AB etc; allows us to go beyond Z)
    dfsummary[bfstartrow:bfendrow,3]<-N_size[1:(length(N_size)-1)]/2 #this is N per condition (ie half N_size)
    dfsummary[bfstartrow:bfendrow,4]<-1:(length(N_size)-1) #level - ie not actual number, but code for sample size from 1;6
    dfsummary[bfstartrow:bfendrow,5]<-E
    
    
    
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
        dfsummary[(bfstartrow-1+jcounter),6]<- ESobs
        #multiply by -1 bcs factor takes control first, so all eff sizes are neg
        # NB need to specify cohen.d from effsize package - otherwise conflict with psych package
        
        bf = ttestBF(x = mydfL$Score[1:(j/2)],y=mydfL$Score[((j/2)+1):j], paired=FALSE,nullInterval=c(-Inf,0)) # Bayes Factor from t-test
        #I have difficulty understanding this. It is on raw rather than log scale - directional test gives 2 BFs?
        myt<-t.test(mydfL$Score ~ mydfL$Group)
        dfsummary[(bfstartrow-1+jcounter),7]<-(-myt$statistic)
        dfsummary[(bfstartrow-1+jcounter),8]<-makeLL(Elist[1]/10,ESobs,(j/2))
        dfsummary[(bfstartrow-1+jcounter),9]<-as.data.frame(bf)$bf[1] #this is the only way to get into the bf output!
        dfsummary[(bfstartrow-1+jcounter),10]<-as.data.frame(bf)$bf[2] #2 BFs generated for directional test
        dfsummary[(bfstartrow-1+jcounter),11]<-myt$estimate[2]
        dfsummary[(bfstartrow-1+jcounter),13]<-myt$estimate[1]
        dfsummary[(bfstartrow-1+jcounter),12]<-sd(mydfL$Score[((j/2)+1):j])
        dfsummary[(bfstartrow-1+jcounter),14]<-sd(mydfL$Score[1:(j/2)])

      }
    }
    
    
    # add means and sds to spreadsheet
    myspread[spreadrowcount,seq(11,34,4)]<-dfsummary$meanC[bfstartrow:bfendrow]
    myspread[spreadrowcount,seq(12,34,4)]<-dfsummary$sdC[bfstartrow:bfendrow]
    myspread[spreadrowcount,seq(13,34,4)]<-dfsummary$meanE[bfstartrow:bfendrow]
    myspread[spreadrowcount,seq(14,34,4)]<-dfsummary$sdE[bfstartrow:bfendrow]
    myspread[spreadrowcount,35:40]<-dfsummary$obsES[bfstartrow:bfendrow]
    myspread[spreadrowcount,41:46]<-dfsummary$LL[bfstartrow:bfendrow]
    myspread[spreadrowcount,47:52]<-dfsummary$t[bfstartrow:bfendrow]
    myspread[spreadrowcount,53:58]<-dfsummary$BF[bfstartrow:bfendrow]
    
    myspread$EarningCorrect[spreadrowcount]<-corrEarn
    myspread$EarningWrong[spreadrowcount]<-wrongEarn
    
    # bind the data    
    big_data = do.call(rbind, datalist)  
    big_data$iter <- as.factor(big_data$n)
    
    big_data$cnd_sample <- big_data$n
    big_data$cnd_sample <- as.factor(big_data$cnd_sample)
    
    # Create jpgs ====================================================================================================
    if(jpgmake==1){
      # create plots
      plotdat<-subset(big_data, iter != (max(N_size)/2))
      #plot modified to make points semi-transparent with alpha=.4, so lines are visible
      p <- ggplot(plotdat, aes(x= cnd_sample, y= Score, color= Group)) +
        geom_beeswarm(cex = 1.5, dodge.width = .15,alpha=.4) + 
        xlab("Sample size per condition") + ylim(-3,3)+
        stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                     geom = "crossbar", width = .75)  #add line corresponding to mean
      
      p <- p + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top") +
        theme(text = element_text(size=15))
      #p
      ggsave(paste0('gorilla_jpgs/bees',E,'_',num2alpha(i),'.jpg'),width= 7, height = 7,p)
      
      #now save beeswarms for individual sample sizes
      for (n in N_size[1:6]){ #omit last one
        databit<-plotdat[plotdat$cnd_sample==(n/2),]
        px <-  ggplot(databit, aes(x= cnd_sample, y= Score, color= Group)) +
          geom_beeswarm(cex = 1.5, dodge.width = .15,alpha=.4) + 
          xlab("Sample size per condition") + ylim(-3,3)+
          stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                       geom = "crossbar", width = .75)  #add line corresponding to mean
        
        px <- px + theme_bw() + theme(legend.title = element_blank()) + theme(legend.position="top") +
          theme(text = element_text(size=15))
        px
     
      ggsave(paste0('gorilla_jpgs/individual_plots/bees',E,'_',num2alpha(i),'_N',n,'.jpg'),width= 4, height = 7,px)
      }
      
    }
    
 
 #Save original data that figs were made from 
  csvname<-'gorilla_jpgs/bigdata.csv'

  #for completeness can add p-values to dfsummary
  dfsummary$p<-1-pt(dfsummary$t,(2*dfsummary$Nsize-2))
  #NB when log scale used, the evidence for null is just the same as for H1, but negative
  write.csv(dfsummary,csvname,row.names=FALSE)
  
  myspread$image1f[spreadrowcount]<-paste0('E',Elist[effsize],'_',num2alpha(i),'.jpg')
  myspread$randomise_trials[spreadrowcount]<-1
}

}
myspread$display[(spreadrowcount+1)]<-'end'

# Shuffle rows of spreadsheet aiming for equal N 0 and E effect sizes in each block
trialrows <- which(myspread$randomise_trials==1)
triallist <- trialrows[shuffle(trialrows)]
myspread1 <- myspread
myspread1[triallist,]<-myspread1[trialrows,]
myspread1$skip_ref[triallist] <- myspread$skip_ref[triallist] #skip_ref is name for block (not sure why!); restore this from original so that blocks correctly numbered

#add column that specifies array size at which absolute LL exceeds 3.68 (equal to 40:1 chance of null:true or vv)
# This can be used for bonus.
myspread1$bonus <- NA
c <-which(colnames(myspread1)=='LL1')
for (ii in 1:nrow(myspread1)){
  if(!is.na(myspread1$LL1[ii])){ #ignore rows with blanks (block dividers)
     w<-first(which(abs(myspread1[i,c:(c+5)])>3.68)) #find first instance of LL>3.68
      if(length(w)>0)
        {myspread1$bonus[ii] <- w}
       if(length(w)==0)
       {myspread1$bonus[ii] <- 6} #if no instance, then set bonus for 6
  }
}

write.csv(myspread1,'gorilla_spreadsheets/spreadsheet_postRR.csv',row.names=FALSE)

