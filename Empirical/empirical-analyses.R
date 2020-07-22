library(tsvr)
library(tidyverse)
library(cowplot)
library(grid)

# Function that calculates SE
calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

###########################
# Plantago and Microseris #
##########################

# read in the Plantago-Microseris data
# data subset such that both species were at cover > 3 before gopher (year 0)
# were fully disturbed the next year (year 1) and then were not fully disturbed for at least 8 years
plmi <- read_csv("Empirical/plantago-microseris.csv")

#### Graph the aggregate underlying data #### 

# aggregate by year and species
plmi2 <- plmi %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 

# make the over time panel
a <- ggplot(subset(plmi2, rowdif < 18), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) 


#### Analyze each subplots with the tsvr #### 

## calculate the VRs within each quadID
outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plmi$quadID)

for (i in 1:length(quads)){
  
  # subset by rep
  subber <- subset(plmi, quadID == quads[i]) %>%
    tbl_df()
  
  # select species and fill 0s
  subber2 <- subber %>%
    select(year, species, cover) %>%
    spread(species, cover, fill = 0)
  
  # transpose the data
  d <- t(as.matrix(subber2[,2:dim(subber2)[2]]))
  
  # create a dataframe with replicate info to append to
  subdat <-subber %>%
    select(quadID, treatment) %>%
    unique()
  
  # calculate tsvr
  res0 <- vreq_classic(d)
  
  # append classic
  subdat$classicVR <- res0[[3]]
  
  # aggregate short vs long
  res <-tsvreq_classic(d)
  
  aggresLong<-aggts(res,res$ts[res$ts>=4])
  aggresShort<-aggts(res,res$ts[res$ts<4])
  
  # append short and long
  subdat$longVR <- aggresLong[[3]]
  subdat$shortVR <- aggresShort[[3]]
  
  # append to external dataframe
  siteout<-rbind(siteout, subdat)
  
}

# Calculate the mean and std of the tsvr values
siteout2 <- siteout %>%
  gather(metric, value, classicVR:shortVR) %>%
  group_by(metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))
siteout2$metric2 <- c("Classic", "Long 
                      Timescale", "Short 
                      Timescale")
siteout2$facorder <- c(3,2,1)


# plot the average and se of the tsvr
b <- ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2)  +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "lightseagreen", "darkgreen",  "blue"))  + 
  labs(x = "", y="Variance Ratio")



#######################
# Plantago and Elymus #
#######################
# read in the Plantago-Elymus data
# data subset such that both species were at cover > 3 before gopher (year 0)
# were fully disturbed the next year (year 1) and then were not fully disturbed for at least 8 years
plel <- read_csv("Empirical/plantago-elymus.csv")

#### Graph the aggregate underlying data #### 

# aggregate by year and species
plel2 <- plel %>%
  group_by(rowdif, quadID, species, treatment, minyear) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 

# make the over time panel
a2 <- ggplot(subset(plel2, rowdif < 10), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) 


#### Analyze each subplots with the tsvr #### 

## calculate the VRs within each quadID
outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plel$quadID)

for (i in 1:length(quads)){
  
  # subset by rep
  subber <- subset(plel, quadID == quads[i]) %>%
    tbl_df() 
  
  # select species and fill 0s
  subber2 <- subber %>%
    select(year, species, cover) %>%
    spread(species, cover, fill = 0)
  
  # transpose the data
  d <- t(as.matrix(subber2[,2:dim(subber2)[2]]))
  
  # create a dataframe with replicate info to append to
  subdat <-subber %>%
    select(quadID, treatment) %>%
    unique()
  
  # calculate tsvr
  res0 <- vreq_classic(d)
  subdat$classicVR <- res0[[3]]
  
  # append classic
  res <-tsvreq_classic(d)
  
  # aggregate short vs long
  aggresLong<-aggts(res,res$ts[res$ts>=4])
  aggresShort<-aggts(res,res$ts[res$ts<4])
  
  # append short and long
  subdat$longVR <- aggresLong[[3]]
  subdat$shortVR <- aggresShort[[3]]
  
  # append to external dataframe
  siteout<-rbind(siteout, subdat)
  
}

# Calculate the mean and std of the tsvr values
siteout2 <- siteout %>%
  gather(metric, value, classicVR:shortVR) %>%
  group_by(metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))
siteout2$metric2 <- c("Classic", "Long 
                      Timescale", "Short 
                      Timescale")
siteout2$facorder <- c(3,2,1)


# plot the average and se of the
b2 <- ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + ylim(0.75,1.25) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "lightseagreen", "darkgreen",  "blue"))  + 
  labs(x = "", y="Variance Ratio")


# put it all together
c <- plot_grid(a + theme(legend.position = "none") + annotate("text", x=0, y =32, label="a)", size = 5) +
                 scale_color_manual(values = c("black","darkgrey")) + labs(x="Time", y = "Percent Cover") + 
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75))+
                 annotate("text", x = 12, y = 24, label = "Plantago erecta",fontface = 'italic', color = "darkgrey") +
                 annotate("text", x = 6.75, y = 14, label = "Microseris douglasii",fontface = 'italic', color = "black") +
                 scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18)) + 
                 labs(x=""),
               b  + annotate("text", x=.6, y =1.25, label="b)", size = 5) +
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)) + labs (y=""),
               align = c("hv"))

c2 <- plot_grid(a2 + theme(legend.position = "none") + annotate("text", x=0, y =32, label="c)", size = 5) +
                  scale_color_manual(values = c("darkorchid1", "orange3")) + labs(x="Time", y = "Percent Cover") + 
                  theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)) + 
                  scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
                  annotate("text", x = 6, y = 31, label = "Plantago erecta",fontface = 'italic', color = "darkorchid1") +
                  annotate("text", x = 6, y = 10, label = "Elymus glaucus",fontface = 'italic', color = "orange3"),
                b2  + annotate("text", x=.6, y =1.25, label="d)", size = 5) +
                  theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)),
                align = c("hv"))

pdf("jr-empirical2.pdf", width = 10, height = 8)
plot_grid(c, c2, nrow = 2)
dev.off()
