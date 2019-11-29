# Created by Adam Parker, August 2019
# Generates the gifs showing growing beeswarms.
# Modified by Dorothy Bishop, Sept 2019 to add spreadsheet creation
# Option for larger spreadsheet with datasets that are not plotted, but which are used to test strategies

# Gorilla Game is here:  https://gorilla.sc/admin/project/7839#

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
source("split_vio.R")


# function to convert number to letter pair for values > 26
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

# simulate data ================================================================================================

gifmake<-1 #set to 1 for creating gifs for study; 0 if making lots of datasets for testing strategies
Elist<-c(3,5,8,0) #effect size x 10
Elist<-c(3,0) #comment out if you want all 3 effect sizes

nsets<-c(10,10,10,30) # N samples at each effect size for when creating gifs for Gorilla
nsets<-c(30,30)  #comment out if you want all 3 effect sizes
blocklength<-15 #sum of nsets should be divisible exactly by blocklength

# for testing ================================================================================================
dotest<-0
if(dotest==1){
nsets<-c(4,4,4,12) #for testing
blocklength<-6}
#set dotest to zero once testing complete.
# ================================================================================================

if(gifmake==0){
nsets<-c(60,60,60,180) #N sets at each effect size - here doing lots so can check accuracy rates of differnt strategies
}

N_size <- c(20, 40, 80, 160, 320, 640, 20000) #sample sizes to use
allN<-length(N_size)-1
thistrial<-0 #initialise counter 

#make dataframe to save observed eff size and bayes factor for each N for each run
mybf<-data.frame(matrix(nrow=(sum(nsets)*allN),ncol=11))
 names(mybf)<-c('run','lettercode','Nsize','trueES','obsES','t','BF','meanC','sdC','meanE','sdE')
 #make spreadsheet for Gorilla
myspread<-data.frame(matrix(nrow=(2+length(nsets)+sum(nsets)),ncol=34)) #rows need to allow for intro and breaks
names(myspread)<-c('randomise_blocks','randomise_trials','display','ANSWER','image1','image1f','ES','EarningCorrect','EarningWrong',
                   'skip_ref',
                   'meanC1','sdC1','meanE1','sdE1',
                   'meanC2','sdC2','meanE2','sdE2',
                   'meanC3','sdC3','meanE3','sdE3',
                   'meanC4','sdC4','meanE4','sdE4',
                   'meanC5','sdC5','meanE5','sdE5',
                   'meanC6','sdC6','meanE6','sdE6')
spreadrowcount <- 2 #keep count of row number for writing to spreadsheet - will start at row 3, after demo and instructions
myspread$display<-'beeswarms' #default
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
    
    spreadrowcount<-spreadrowcount+1
    blockcount<-blockcount+1
    if(blockcount>blocklength)
       {blockcount<-0
         myblock<-myblock+1
         myspread[spreadrowcount,]<-NA
         myspread$display[spreadrowcount]<-'break'
         spreadrowcount<-spreadrowcount+1 }
    myspread$skip_ref[spreadrowcount]<-myblock
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
mybf[bfstartrow:bfendrow,2]<-num2alpha(i) #use previously defined formula (allows us to go beyond Z)
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
        mybf[(bfstartrow-1+jcounter),5]<- (-1)*effsize::cohen.d(mydfL$Score~as.factor(mydfL$Group))$estimate
        #multiply by -1 bcs factor takes control first, so all eff sizes are neg
        # NB need to specify cohen.d from effsize package - otherwise conflict with psych package
        
        bf = ttestBF(x = mydfL$Score[1:(j/2)],y=mydfL$Score[((j/2)+1):j], paired=FALSE)
        myt<-t.test(x=mydfL$Score[((j/2)+1):j],y = mydfL$Score[1:(j/2)], paired=FALSE)
        mybf[(bfstartrow-1+jcounter),6]<-myt$statistic
        mybf[(bfstartrow-1+jcounter),7]<-as.data.frame(bf)$bf
        mybf[(bfstartrow-1+jcounter),8]<-mean(mydfL$Score[((j/2)+1):j])
        mybf[(bfstartrow-1+jcounter),10]<-mean(mydfL$Score[1:(j/2)])
        mybf[(bfstartrow-1+jcounter),9]<-sd(mydfL$Score[((j/2)+1):j])
        mybf[(bfstartrow-1+jcounter),11]<-sd(mydfL$Score[1:(j/2)])
      }
  }

# add means and sds to spreadsheet
myspread[spreadrowcount,seq(11,34,4)]<-mybf$meanC[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(12,34,4)]<-mybf$sdC[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(13,34,4)]<-mybf$meanE[bfstartrow:bfendrow]
myspread[spreadrowcount,seq(14,34,4)]<-mybf$sdE[bfstartrow:bfendrow]

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

myspread$image1f[spreadrowcount]<-paste0('E',Elist[effsize],'_',LETTERS[i],'.png')
myspread$image1[spreadrowcount]<-paste0('E',Elist[effsize],'_',LETTERS[i],'.gif')
myspread$randomise_trials[spreadrowcount]<-1
  }
}
myspread$display[spreadrowcount]<-'end'
write.csv(myspread,'gorilla_bits/spreadsheet1.csv',row.names=FALSE)