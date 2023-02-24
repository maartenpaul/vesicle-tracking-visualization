#required packages
library(plyr)
library(tidyverse)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)

#run these functions before running analysis
msd_analyze_data_BJ <- function(directory,condition_list,track_file_name="tracks.simple.filtered.txt",extension=""){
  segments_all <- list()
  msd_fit_all <- list()
  tracks_all <- list()
  
  for (i in 1:length(condition_list)){ #loop over the different conditions of the acquired data which is organized in sub-folders
    segments <- list()
    tracks <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    for(j in 1:length(filelist)){
      tracks_simple <- read_tsv(file.path(filelist[j],track_file_name),col_names = F) #import tab separate file
      names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error") #adjust column names
      segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(dirname(file.path(filelist[j],track_file_name)))) #analyze the segments angle and displacements
      write.table(segments[[j]],file = file.path(dir,paste0(basename(filelist[j]),"_segment_stats.txt")))
      tracks[[j]] <- data.frame(TRACK_STAT(tracks_simple),"cellID"=basename(dirname(file.path(filelist[j],track_file_name)))) #analyze the segments angle and displacements
      write.table(tracks[[j]],file = file.path(dir,paste0(basename(filelist[j]),"_tracks_stats.txt")))
     # make plot per cell
      angles_data <- segments[[j]]$angle[segments[[j]]$angle>=0]
      angles <- hist(c(angles_data,abs(360-angles_data)),breaks = seq(0,360,15),main = paste(filelist[j]),xlab = "angle")
      
      angles <- data.frame('mids'=angles$mids,'density'=angles$density)
      
      histogram <- ggplot(angles,aes(x = mids,y=density))+geom_bar(stat='identity',width=15,fill='white',color="black") +ylim(c(0,.01))+ggtitle(paste(basename(dirname(dirname(filelist[j]))),basename(dirname(filelist[j])),basename(filelist[j]),sep = "_"))
      ggsave(histogram,filename = file.path(filelist[j],paste0("histogram",".pdf")))
      ggsave(histogram,filename = file.path(filelist[j],paste0("histogram",".png")))
      
      rose_plot <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=15) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
        theme(legend.position = "none",text = element_text(size=15))+xlab("")+ylab("")+ggtitle(paste(basename(dirname(dirname(filelist[j]))),basename(dirname(filelist[j])),basename(filelist[j]),sep = "_"))+theme(axis.text.x=element_blank(), #remove x axis labels
                                                                                                                           axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                           axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                           axis.ticks.y=element_blank() ) #remove y axis ticks
      ggsave(rose_plot,filename = file.path(filelist[j],paste0("radialplot",".pdf")))
      ggsave(rose_plot,filename = file.path(filelist[j],paste0("radialplot",".png")))
      
    }
    segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))
    tracks_all[[basename(dir)]] <- data.frame(ldply(tracks),"condition"=basename(dir))
    
 
  }
  #save data to the folder
  save(segments_all,file=file.path(directory,paste0("segments_all",extension,".Rdata")))
  save(tracks_all,file=file.path(directory,paste0("tracks_all",extension,".Rdata")))
  

}
track_stat <- function(x,framerate=30,pxsize=100){
  
  x$X <- (x$X*pxsize)/1000
  x$Y <- (x$Y*pxsize)/1000
  out <- ddply(x,.variables = "track",.fun= function(x) {
    speed <- 0
    for (i in 2:nrow(x)){
      speed <- speed + (x$X[i]-x$X[i-1])^2+(x$Y[i]-x$Y[i-1])^2/((x$frame[i]-x$frame[i-1])*framerate)
    }
    speed <- speed/nrow(x)
    
    coord <- cbind(x$X,x$Y)
    # use pricipal component analysis on X and Y coordinates to get eigenvectors: major and minor axis
    D <- princomp(coord)
    angle <- atan2(D$loadings[2,1],D$loadings[1,1])
    #calculate convex hull and futher statistics
    y <- chull(coord)
    area <- 0 #set to zero because of error, parameter not used
    perimeter <- 0 #set to zero because of error, parameter not used
    D_chull <- 0 #set to zero because of error, parameter not used
    
    #return(data.frame("sd"=((sd(x$X)+sd(x$Y))/2)*2.35,"N"=nrow(x),"channel"=1))
    # return(data.frame(,"N"=nrow(x),"channel"=1))
    return(data.frame("N"=nrow(x),"meanX"=mean(x$X),"meanY"=mean(x$Y),"meanspeed"=speed ,
                      "sd"=((sd(x$X)+sd(x$Y))/2),"sdpri"=((D$sdev[1]+D$sdev[2])/2),"major"=D$sdev[1],"minor"=D$sdev[2],
                      "width"=(max(D$scores[,1])-min(D$scores[,1])),"ratio"=(D$sdev[1]/D$sdev[2]),"angle"=angle,
                      "chull_area"=area,"chull_perimeter"=perimeter,"chull_major"=D_chull,"chull_minor"=D_chull,meanInt=mean(x$intensity),medianInt=median(x$intensity)))
  })
  
  return(out)
}

TRACK_STAT <- function(x,framerate=30,pxsize=100){
  UseMethod("TRACK_STAT")
}

TRACK_STAT.default <- function(x,framerate=30,pxsize=100){
  stop("TRACK_STAT requires data frame")
}

TRACK_STAT.data.frame <-  function(x,framerate=30,pxsize=100){
  track_stat(x,framerate,pxsize)
  
}

TRACK_STAT.list <-  function(x,framerate=30,pxsize=100){
  llply(x,function(x){
    TRACK_STAT(x,framerate,pxsize)
  })
}


#input variables
root_dir <- ""

datasets <- c("Hypoxia/MSD24","Hypoxia/MSD25","Normoxia/MSD22","Normoxia/MSD23")

conditions <- c("Hypoxia","Normoxia")

channels <- data.frame(channel=c("a","b"),names=c("GFP-Clathrin","EGF-TMR"))

#process data
for (dir in datasets) {
  directory <- file.path(root_dir,dir)
  condition_list <- list.dirs(directory,full.names = F,recursive = F)
  msd_analyze_data_BJ(directory,condition_list,track_file_name="tracks.simple.filtered.txt")
}



#make plots per replicate per condition per timepoint
for (dir in datasets) {
  segments_file <- file.path(root_dir,dir,"segments_all.Rdata")
  load(segments_file)
  segments_all <- ldply(segments_all)
  
  assy_coef_results <- segments_all %>%
    group_by(cellID, condition) %>%
    filter(angle>=0)%>%
    dplyr::summarize(assy_coef=length(angle[angle>=150])/length(angle[angle<=30]))
    
  write_tsv(assy_coef_results,file = file.path(root_dir,paste0(dir,"_assy_coef.txt")))
  
  for (i in 1:9){
    for(channel in c("a","b")){
      inputdata <- segments_all %>%
        dplyr::filter(cellID==i&str_detect(condition,channel)&angle>=0)
      
      angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),main = paste(dir,i,channels$name[which(channels$channel==channel)]),xlab = "angle")
      
      angles <- data.frame('mids'=angles$mids,'density'=angles$density)
      
      histogram <- ggplot(angles,aes(x = mids,y=density))+geom_bar(stat='identity',width=15,fill='white',color="black") +ylim(c(0,.01))+ggtitle(paste(dir,i,channels$name[which(channels$channel==channel)]))
      ggsave(histogram,filename = file.path(root_dir,dir,paste0("histogram_",i,"_",channels$name[which(channels$channel==channel)],".pdf")))
      ggsave(histogram,filename = file.path(root_dir,dir,paste0("histogram_",i,"_",channels$name[which(channels$channel==channel)],".png")))
      
      rose_plot <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=15) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
        theme(legend.position = "none",text = element_text(size=15))+xlab("")+ylab("")+ggtitle(paste(dir,i,channels$name[which(channels$channel==channel)]))+theme(axis.text.x=element_blank(), #remove x axis labels
                                                                                             axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                             axis.text.y=element_blank(),  #remove y axis labels
                                                                                             axis.ticks.y=element_blank() ) #remove y axis ticks
      ggsave(rose_plot,filename = file.path(root_dir,dir,paste0("radialplot_",i,"_",channels$name[which(channels$channel==channel)],".pdf")))
      ggsave(rose_plot,filename = file.path(root_dir,dir,paste0("radialplot_",i,"_",channels$name[which(channels$channel==channel)],".png")))
      
    }

  }

}

###make plots per condition per timepoint

for (condition in conditions) {
  dirs <- list.dirs(file.path(root_dir,condition),full.names = F,recursive = F)
  dirs <- dirs[-1]
  if (exists("segments_all2")){
    rm(segments_all2)
  }
  for (dir in dirs){
    segments_file <- file.path(root_dir,condition,dir,"segments_all.Rdata")
    load(segments_file)
    segments_all <- ldply(segments_all)
    if (exists("segments_all2")){
      rbind(segments_all2,segments_all)
    } else {
      segments_all2 <- segments_all
    }
} 

  for (i in 1:9){
    for(channel in c("a","b")){
      inputdata <- segments_all2 %>%
        dplyr::filter(cellID==i&str_detect(condition,channel)&angle>0)
      
      angles <- hist(c(inputdata$angle,abs(360-inputdata$angle)),breaks = seq(0,360,15),main = paste(dir,i,channels$name[which(channels$channel==channel)]),xlab = "angle")
      
      angles <- data.frame('mids'=angles$mids,'density'=angles$density)
      
      
      histogram <- ggplot(angles,aes(x = mids,y=density))+geom_bar(stat='identity',width=15,fill='white',color="black") +ylim(c(0,.01))+ggtitle(paste(dir,i,channels$name[which(channels$channel==channel)]))
      ggsave(histogram,filename = file.path(root_dir,paste0("histogram_",condition,"_",i,"_",channels$name[which(channels$channel==channel)],".pdf")))
      ggsave(histogram,filename = file.path(root_dir,paste0("histogram_",condition,"_",i,"_",channels$name[which(channels$channel==channel)],".png")))
      
      rose_plot <- ggplot(angles,aes(x = mids,y=density,fill=density))+geom_bar(stat='identity',width=15) +ylim(c(0,.01))+scale_x_continuous(breaks=seq(0, 350, 45))+coord_polar(start = 0.5*pi,direction = 1)+
        theme(legend.position = "none",text = element_text(size=15))+xlab("")+ylab("")+ggtitle(paste(dir,i,channel))+theme(axis.text.x=element_blank(), #remove x axis labels
                                                                                                                           axis.ticks.x=element_blank(), #remove x axis ticks
                                                                                                                           axis.text.y=element_blank(),  #remove y axis labels
                                                                                                                           axis.ticks.y=element_blank() ) #remove y axis ticks
      ggsave(rose_plot,filename = file.path(root_dir,paste0("radialplot_",condition,"_",i,"_",channels$name[which(channels$channel==channel)],".pdf")))
      ggsave(rose_plot,filename = file.path(root_dir,paste0("radialplot_",condition,"_",i,"_",channels$name[which(channels$channel==channel)],".png")))
      
    }
  }
}


###plot assymetry coefficients
rm(assycoef_all)
for (dir in datasets) {
  assycoef  <- read_tsv(file.path(root_dir,paste0(dir,"_assy_coef.txt")))
  dir_split <- unlist(strsplit(dir,"/"))
  assycoef$protein <- do.call(rbind,strsplit(assycoef$condition, split = "(?<=[0-9])\\s*(?=[a-zA-Z])", perl = TRUE))[,2] #substr(assycoef$condition,start = 2,stop = 2)[[1]][2]
  assycoef$condition <- dir_split[1]
  assycoef$replicate <- dir_split[2]
  
  
  if(exists("assycoef_all")){
    assycoef_all<-rbind(assycoef_all,assycoef)
  } else {
    assycoef_all <- assycoef
  }
  
 
}
p <- assycoef_all %>%
  ggplot(aes(x=as.character(cellID),y=(assy_coef),color=c(protein)))+geom_boxplot()+facet_grid(.~condition)+xlab("time")+ggtitle("Assymetry coefficient estimated per cell")
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_all.pdf"))
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_all.png"))

p <- assycoef_all %>%
  filter(condition=="Hypoxia")%>%
  ggplot(aes(x=as.character(cellID),y=(assy_coef),color=protein))+geom_boxplot()+facet_grid(condition~replicate)+xlab("time")
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_hypoxia.pdf"))
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_hypoxia.png"))

p <- assycoef_all %>%
  filter(condition=="Normoxia")%>%
  ggplot(aes(x=as.character(cellID),y=(assy_coef),color=protein))+geom_boxplot()+facet_grid(condition~replicate)+xlab("time")
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_normoxia.pdf"))
ggsave(p,filename = file.path(root_dir,"assyymetry_coef_normoxia.png"))

#make plots of mean intensity of tracks
if (exists("tracks_all2")){
  rm(tracks_all2)
}
for (condition in conditions) {
  dirs <- list.dirs(file.path(root_dir,condition),full.names = F,recursive = F)
  #dirs <- dirs[-1]
  for (dir in dirs){
    tracks_file <- file.path(root_dir,condition,dir,"tracks_all.Rdata")
    load(tracks_file)
    tracks_all <- ldply(tracks_all)
    dir_split <- unlist(strsplit(dir,"/"))
    tracks_all$protein <- do.call(rbind,strsplit(as.character(tracks_all$condition), split = "(?<=[0-9])\\s*(?=[a-zA-Z])", perl = TRUE))[,2] #substr(assycoef$condition,start = 2,stop = 2)[[1]][2]
    tracks_all$condition <- condition
    tracks_all$replicate <- dir
    
    

    if (exists("tracks_all2")){
      tracks_all2<- rbind(tracks_all2,tracks_all)
    } else {
      tracks_all2 <- tracks_all
    }
  } 
}


p <- tracks_all2 %>%
  ggplot(aes(x=as.character(cellID),y=meanInt,color=c(protein)))+geom_boxplot()+facet_grid(.~condition)+xlab("time")+ggtitle("Mean track intensity")
p
ggsave(p,filename = file.path(root_dir,"mean_int_all.pdf"))
ggsave(p,filename = file.path(root_dir,"mean_int_all.png"))

p <- tracks_all2 %>%
  ggplot(aes(x=as.character(cellID),y=medianInt,color=c(protein)))+geom_boxplot()+facet_grid(.~condition)+xlab("time")+ggtitle("Median track intensity")
p
ggsave(p,filename = file.path(root_dir,"median_int_all.pdf"))
ggsave(p,filename = file.path(root_dir,"median_int_all.png"))

p <- tracks_all2 %>%
  filter(condition=="Hypoxia")%>%
  ggplot(aes(x=as.character(cellID),y=meanInt,color=protein))+geom_boxplot()+facet_grid(condition~replicate)+xlab("time")
p
ggsave(p,filename = file.path(root_dir,"mean_int_hypoxia.pdf"))
ggsave(p,filename = file.path(root_dir,"mean_int_hypoxia.png"))

p <- tracks_all2 %>%
  filter(condition=="Normoxia")%>%
  ggplot(aes(x=as.character(cellID),y=meanInt,color=protein))+geom_boxplot()+facet_grid(condition~replicate)+xlab("time")
p
ggsave(p,filename = file.path(root_dir,"mean_int_normoxia.pdf"))
ggsave(p,filename = file.path(root_dir,"mean_int_normoxia.png"))

p <- tracks_all2 %>%
  group_by(cellID,protein,condition,replicate)%>%
  summarize(meanN=mean(N))%>%
  ggplot(aes(x=as.character(cellID),y=meanN,color=c(protein)))+geom_boxplot()+facet_grid(.~condition)+xlab("time")+ggtitle("Mean track intensity")
p
ggsave(p,filename = file.path(root_dir,"length_all.pdf"))
ggsave(p,filename = file.path(root_dir,"length_int_all.png"))
###

p <- tracks_all2 %>%
  group_by(cellID,protein,condition,replicate)%>%
  summarize(meanInt=mean(meanInt))%>%
  ggplot(aes(x=as.character(cellID),y=meanInt,color=c(protein)))+geom_boxplot()+facet_grid(.~condition)+xlab("time")+ggtitle("Mean track intensity")
p


for(x in datasets) {
  dataset <- unlist(strsplit(x,"/"))
  
  tracks_all2%>%
    dplyr::filter(condition==dataset[1],replicate==dataset[2])%>%
    dplyr::select(.id,track,N,meanX,meanY,meanspeed,angle,meanInt)%>%
    write_tsv(file = paste0(file.path(root_dir,x),"_track_stats.txt"))
}

