
library(tidyverse)
library(factoextra)
library(cluster)
library(gridExtra)

# Export dataframe for selection only differential expressed genes
## Upload selected file

# Up regulated
setwd("~/RNA_122022/Kmeans")

## UPC3_vs_C4
dat <- read.csv("norm_UPC3_vs_C4.csv")
dat0 <- read.csv("norm_UPC3_vs_C4.csv")
class(dat)

dat

# set the first column as an index

dat <- dat %>% column_to_rownames(., var = "GeneID")

dat


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat <- as_tibble(scale(dat)) 
head (dat)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

fig1 <- fviz_nbclust(dat, kmeans, method = "silhouette") +
        labs(title = "")
fig1

##Elbow Method
set.seed(123)

fig2 <- fviz_nbclust(dat , kmeans, method = "wss") +
            labs(title = "")
fig2

###Gap Statistic
set.seed(123)

gap_stat <- clusGap(dat, FUN = kmeans, nstart = 2,
                    K.max = 10, B = 100, iter.max = 20)
?kmeans

fig3 <- fviz_gap_stat(gap_stat) +
                labs(title = "")
fig3

## Combine all results in one figure
grid.arrange(fig1, fig2, fig3, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

k2 <- kmeans(dat, centers = 2, nstart = 30)
k3 <- kmeans(dat, centers = 3, nstart = 30)
k4 <- kmeans(dat, centers = 4, nstart = 30)
k5 <- kmeans(dat, centers = 5, nstart = 30)


# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = dat) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = dat) + ggtitle("k = 3")
p3 <- fviz_cluster(k5, geom = "point",  data = dat) + ggtitle("k = 5")
p4 <- fviz_cluster(k4, geom = "point",  data = dat) + ggtitle("k = 4")

print(k2)
print(k3)
print(k5)

grid.arrange(p1, p2, p4, p3, nrow = 2)

# Compute k-means clustering with k = 3
set.seed(123)
final <- kmeans(dat, 3, nstart = 30, iter.max = 30)
print(final)
## I found k =2 to be the best
set.seed(123)
final_2 <- kmeans(dat, 2, nstart = 30, iter.max = 30)
print(final_2)
#Visualize final clusters

Up_C4vsC3 <- fviz_cluster(final_2, geom = "point",  data = dat) +
                      labs(title = "Up_C4/C4_vs_C3", tag="A") +
               theme(plot.title = element_text(hjust= 0.5),
                     panel.background = element_blank(),
                     axis.line = element_line())

Up_C4vsC3

Up_LEAFvsCot_day <- fviz_cluster(finalA, geom = "point",  data = dat2) +
  labs(title = "Up_L/L_vs_CD", tag="A") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())



#####

# Calulate cluster level means
dat
ClustMeans <- dat[1:2] %>%
  mutate(Cluster = final_2$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeans
##################################################
# Visualize   

ClustMeans %>% 
  gather(Samples, Clust, SP:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
Clu_1 <- dat0[1] %>%
  mutate(Cluster = final_2$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

Clu_1
write.csv(Clu_1,"Clust_1_UpC3vsC4.csv")

#
# Extract transcripts cluster 2
Clu_2 <- dat0[1] %>%
  mutate(Cluster = final_2$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

Clu_2

write.csv(Clu_2,"Clust_2_C3LvsC4.csv")

###############################################################################
###############################################################################

## Dow C3_vs_C4

dat1 <- read.csv("norm_DowC3_vs_C4.csv")
dat1_0 <- read.csv("norm_DowC3_vs_C4.csv")
class(dat1)

dat1

# set the first column as an index

dat1 <- dat1 %>% column_to_rownames(., var = "GeneID")

dat1


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat1 <- as_tibble(scale(dat1)) 
head (dat1)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

fig1_1 <- fviz_nbclust(dat1, kmeans, method = "silhouette") +
  labs(title = "")
fig1_1

##Elbow Method
set.seed(123)

fig2_2 <- fviz_nbclust(dat1 , kmeans, method = "wss") +
  labs(title = "")
fig2_2

###Gap Statistic
set.seed(123)

gap_stat1 <- clusGap(dat1, FUN = kmeans, nstart = 2,
                    K.max = 10, B = 100)
?kmeans

fig3_3 <- fviz_gap_stat(gap_stat1) +
  labs(title = "")
fig3_3

## Combine all results in one figure
grid.arrange(fig1_1, fig2_2, fig3_3, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

k2_2 <- kmeans(dat1, centers = 2, nstart = 30)
k3_3 <- kmeans(dat1, centers = 3, nstart = 30)
k4_4 <- kmeans(dat1, centers = 4, nstart = 30)
k5_5 <- kmeans(dat1, centers = 5, nstart = 30)


# plots to compare
p1_1 <- fviz_cluster(k2_2, geom = "point", data = dat1) + ggtitle("k = 2")
p2_2 <- fviz_cluster(k3_3, geom = "point",  data = dat1) + ggtitle("k = 3")
p3_3 <- fviz_cluster(k4_4, geom = "point",  data = dat1) + ggtitle("k = 4")
p4_4 <- fviz_cluster(k5_5, geom = "point",  data = dat1) + ggtitle("k = 5")

print(k2_2)
print(k3_3)
print(k5_5)

grid.arrange(p1_1, p2_2, p3_3, p4_4, nrow = 2)

# Compute k-means clustering with k = 3
## I found k =2 to be the best
set.seed(123)
final_3 <- kmeans(dat1, 2, nstart = 30, iter.max = 30)
print(final_3)
#Visualize final clusters

Dow_C4vsC3 <- fviz_cluster(final_3, geom = "point",  data = dat1) +
  labs(title = "Up_C3/C4_vs_C3", tag="B") +
  theme(plot.title = element_text(hjust= 0.5))

Dow_C4vsC3

############################
# Combining C3 and C4 results
# Save the legend of one of the sample
legend1 <- get_legend(Up_C4vsC3 + theme(legend.position="bottom"))

## Remove legends of all samples
Up_C4vsC31 <- Up_C4vsC3 + theme(legend.position="none")
Dow_C4vsC31 <- Dow_C4vsC3 + theme(legend.position="none")


# Combine them with a common legend and # Change the legend position

grid.arrange(Up_C4vsC31, Dow_C4vsC31, bottom=legend1$grobs[[1]], ncol=2)

##############################

# Calulate cluster level means
dat1
ClustMeans1 <- dat1[1:2] %>%
  mutate(Cluster = final_3$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeans1
##################################################
# Visualize   

ClustMeans1 %>% 
  gather(Samples, Clust, SP:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
DowClu_1 <- dat1_0[1] %>%
  mutate(Cluster = final_3$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

DowClu_1
write.csv(DowClu_1,"Clust_1_DowC3vsC4.csv")

#
# Extract transcripts cluster 2
DowClu_2 <- dat1_0[1] %>%
  mutate(Cluster = final_3$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

DowClu_2

write.csv(DowClu_2,"Clust_2_DowC3LvsC4.csv")
##########################################################################
##########################################################################

## norm_UPLeaf_vs_Coty_day

dat2 <- read.csv("norm_UPcountsLeaf_vs_Coty_day.csv")
dat2_0 <- read.csv("norm_UPcountsLeaf_vs_Coty_day.csv")
class(dat2)

dat2

# set the first column as an index

dat2 <- dat2 %>% column_to_rownames(., var = "GeneID")

dat2


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat2 <- as_tibble(scale(dat2)) 
head (dat2)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figA <- fviz_nbclust(dat2, kmeans, method = "silhouette") +
  labs(title = "")
figA

##Elbow Method
set.seed(123)

figA1 <- fviz_nbclust(dat2 , kmeans, method = "wss") +
  labs(title = "")
figA1

###Gap Statistic
set.seed(123)

gap_statA <- clusGap(dat2, FUN = kmeans, nstart = 2,
                    K.max = 10, B = 50, iter.max = 20)
?kmeans

figA2 <- fviz_gap_stat(gap_statA) +
  labs(title = "")
figA2

## Combine all results in one figure
grid.arrange(figA, figA1, figA2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kA1 <- kmeans(dat2, centers = 2, nstart = 30)
kA2 <- kmeans(dat2, centers = 3, nstart = 30)
kA3 <- kmeans(dat2, centers = 4, nstart = 30)
kA4 <- kmeans(dat2, centers = 5, nstart = 30)


# plots to compare
A1 <- fviz_cluster(kA1, geom = "point", data = dat2) + ggtitle("k = 2")
A2 <- fviz_cluster(kA2, geom = "point",  data = dat2) + ggtitle("k = 3")
A3 <- fviz_cluster(kA3, geom = "point",  data = dat2) + ggtitle("k = 4")
A4 <- fviz_cluster(kA4, geom = "point",  data = dat2) + ggtitle("k = 5")

print(kA1)
print(kA2)
print(kA3)

grid.arrange(A1, A2, A3, A4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalA <- kmeans(dat2, 2, nstart = 30, iter.max = 30)
print(finalA)
#Visualize final clusters

Up_LEAFvsCot_day <- fviz_cluster(finalA, geom = "point",  data = dat2) +
  labs(title = "Up_L/L_vs_CD", tag="A") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

Up_LEAFvsCot_day


#####

# Calulate cluster level means

ClustMeansA <- dat2[1:3] %>%
  mutate(Cluster = finalA$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansA
##################################################
# Visualize   

ClustMeansA %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
UpLcsCClu_1 <- dat2_0[1] %>%
  mutate(Cluster = finalA$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

UpLcsCClu_1

write.csv(UpLcsCClu_1,"UPLeaf_vs_Coty_dayClu1.csv")

#
# Extract transcripts cluster 2
UpLcsCClu_2 <- dat2_0[1] %>%
  mutate(Cluster = finalA$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

UpLcsCClu_2

write.csv(UpLcsCClu_2,"UpLvsCDClu_2.csv")

##########################################################################
##########################################################################

## norm_DowLeaf_vs_Coty_day

dat3 <- read.csv("norm_DowcountsLeaf_vs_Coty_day.csv")
dat3_0 <- read.csv("norm_DowcountsLeaf_vs_Coty_day.csv")
class(dat3)

dat3

# set the first column as an index

dat3 <- dat3 %>% column_to_rownames(., var = "GeneID")

dat3


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat3 <- as_tibble(scale(dat3)) 
head (dat3)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figB <- fviz_nbclust(dat3, kmeans, method = "silhouette") +
  labs(title = "")
figB

##Elbow Method
set.seed(123)

figB1 <- fviz_nbclust(dat3 , kmeans, method = "wss") +
  labs(title = "")
figB1

###Gap Statistic
set.seed(123)

gap_statB2 <- clusGap(dat3, FUN = kmeans, nstart = 2,
                     K.max = 10, B = 50)
?kmeans

figB2 <- fviz_gap_stat(gap_statB2) +
  labs(title = "")
figB2

## Combine all results in one figure
grid.arrange(figB, figB1, figB2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kB1 <- kmeans(dat3, centers = 2, nstart = 30)
kB2 <- kmeans(dat3, centers = 3, nstart = 30)
kB3 <- kmeans(dat3, centers = 4, nstart = 30)
kB4 <- kmeans(dat3, centers = 5, nstart = 30)


# plots to compare
B1 <- fviz_cluster(kB1, geom = "point", data = dat3) + ggtitle("k = 2")
B2 <- fviz_cluster(kB2, geom = "point",  data = dat3) + ggtitle("k = 3")
B3 <- fviz_cluster(kB3, geom = "point",  data = dat3) + ggtitle("k = 4")
B4 <- fviz_cluster(kB4, geom = "point",  data = dat3) + ggtitle("k = 5")

print(kB1)
print(kB2)
print(kB3)

grid.arrange(B1, B2, B3, B4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalB <- kmeans(dat3, 2, nstart = 30, iter.max = 30)
print(finalB)
#Visualize final clusters

Dow_LEAFvsCot_day <- fviz_cluster(finalB, geom = "point",  data = dat3) +
  labs(title = "Up_CD/L_vs_CD", tag="B") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

Dow_LEAFvsCot_day

#####

# Calulate cluster level means

ClustMeansB <- dat3[1:3] %>%
  mutate(Cluster = finalB$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansB
##################################################
# Visualize   

ClustMeansB %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
DowLcsCClu_1 <- dat3_0[1] %>%
  mutate(Cluster = finalB$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

DowLcsCClu_1

write.csv(DowLcsCClu_1,"DowClust_1_Leaf_vs_coty_day.csv")

#
# Extract transcripts cluster 2
DowLcsCClu_2 <- dat3_0[1] %>%
  mutate(Cluster = finalB$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

DowLcsCClu_2

write.csv(DowLcsCClu_2,"DowClust_2_Leaf_vs_coty_day.csv")

###############################################################################
###############################################################################


## norm_Up_CotyNight_vs_Coty_day

dat4 <- read.csv("norm_UPcountsConightVSday.csv")
dat4_0 <- read.csv("norm_UPcountsConightVSday.csv")
class(dat4)

dat4

# set the first column as an index

dat4 <- dat4 %>% column_to_rownames(., var = "GeneID")

dat4


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat4 <- as_tibble(scale(dat4)) 
head (dat4)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figC <- fviz_nbclust(dat4, kmeans, method = "silhouette") +
  labs(title = "")
figC

##Elbow Method
set.seed(123)

figC1 <- fviz_nbclust(dat4 , kmeans, method = "wss") +
  labs(title = "")
figC1

###Gap Statistic
set.seed(123)

gap_statC2 <- clusGap(dat4, FUN = kmeans, nstart = 2,
                      K.max = 10, B = 50)
?kmeans

figC2 <- fviz_gap_stat(gap_statC2) +
  labs(title = "")
figC2

## Combine all results in one figure
grid.arrange(figC, figC1, figC2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kC1 <- kmeans(dat4, centers = 2, nstart = 30)
kC2 <- kmeans(dat4, centers = 3, nstart = 30)
kC3 <- kmeans(dat4, centers = 4, nstart = 30)
kC4 <- kmeans(dat4, centers = 5, nstart = 30)


# plots to compare
C1 <- fviz_cluster(kC1, geom = "point", data = dat4) + ggtitle("k = 2")
C2 <- fviz_cluster(kC2, geom = "point",  data = dat4) + ggtitle("k = 3")
C3 <- fviz_cluster(kC3, geom = "point",  data = dat4) + ggtitle("k = 4")
C4 <- fviz_cluster(kC4, geom = "point",  data = dat4) + ggtitle("k = 5")

print(kC1)
print(kC2)
print(kC3)

grid.arrange(C1, C2, C3, C4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalC <- kmeans(dat4, 2, nstart = 30, iter.max = 30)
print(finalC)
#Visualize final clusters

Up_CNvsCD <- fviz_cluster(finalC, geom = "point",  data = dat4) +
  labs(title = "Up_CN/CN_vs_CD", tag="C") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

Up_CNvsCD

#####

# Calulate cluster level means

ClustMeansC <- dat4[1:3] %>%
  mutate(Cluster = finalC$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansC
##################################################
# Visualize   

ClustMeansC %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
UP_CDvsCNClu_1 <- dat4_0[1] %>%
  mutate(Cluster = finalC$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

UP_CDvsCNClu_1

write.csv(UP_CDvsCNClu_1,"Clust_1_UP_CDvsCNClu_1.csv")

#
# Extract transcripts cluster 2
UP_CDvsCNClu_2 <- dat4_0[1] %>%
  mutate(Cluster = finalC$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

UP_CDvsCNClu_2

write.csv(UP_CDvsCNClu_2,"Clust_2_UP_CDvsCNClu_2.csv")
############################################################################
############################################################################

## norm_Dow_CotyNight_vs_Coty_day

dat5 <- read.csv("norm_DowcountsConightVSday.csv")
dat5_0 <- read.csv("norm_DowcountsConightVSday.csv")
class(dat5)

dat5

# set the first column as an index

dat5 <- dat5 %>% column_to_rownames(., var = "GeneID")

dat5


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat5 <- as_tibble(scale(dat5)) 
head (dat5)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figD <- fviz_nbclust(dat5, kmeans, method = "silhouette") +
  labs(title = "")
figD

##Elbow Method
set.seed(123)

figD1 <- fviz_nbclust(dat5 , kmeans, method = "wss") +
  labs(title = "")
figD1

###Gap Statistic
set.seed(123)

gap_statD2 <- clusGap(dat5, FUN = kmeans, nstart = 2,
                      K.max = 10, B = 50, iter.max= 20)
?kmeans

figD2 <- fviz_gap_stat(gap_statD2) +
  labs(title = "")
figD2

## Combine all results in one figure
grid.arrange(figD, figD1, figD2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kD1 <- kmeans(dat5, centers = 2, nstart = 30)
kD2 <- kmeans(dat5, centers = 3, nstart = 30)
kD3 <- kmeans(dat5, centers = 4, nstart = 30)
kD4 <- kmeans(dat5, centers = 5, nstart = 30)


# plots to compare
D1 <- fviz_cluster(kD1, geom = "point", data = dat5) + ggtitle("k = 2")
D2 <- fviz_cluster(kD2, geom = "point",  data = dat5) + ggtitle("k = 3")
D3 <- fviz_cluster(kD3, geom = "point",  data = dat5) + ggtitle("k = 4")
D4 <- fviz_cluster(kD4, geom = "point",  data = dat5) + ggtitle("k = 5")

print(kD1)
print(kD2)
print(kD3)
print(kD4)
grid.arrange(D1, D2, D3, D4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalD <- kmeans(dat5, 2, nstart = 30, iter.max = 30)
print(finalD)
#Visualize final clusters

Dow_CNvsCD <- fviz_cluster(finalD, geom = "point",  data = dat5) +
  labs(title = "Up_CD/CN_vs_CD", tag="D") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

Dow_CNvsCD

#####

# Calulate cluster level means

ClustMeansD <- dat5[1:3] %>%
  mutate(Cluster = finalD$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansD
##################################################
# Visualize   

ClustMeansD %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
Dow_CDvsCNClu_1 <- dat5_0[1] %>%
  mutate(Cluster = finalD$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

Dow_CDvsCNClu_1

write.csv(Dow_CDvsCNClu_1,"Clust_1_Dow_CDvsCNClu_1.csv")

#
# Extract transcripts cluster 2
Dow_CDvsCNClu_2 <- dat5_0[1] %>%
  mutate(Cluster = finalD$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

Dow_CDvsCNClu_2

write.csv(Dow_CDvsCNClu_2,"Clust_2_Dow_CDvsCNClu_2.csv")

############################################################################
############################################################################

#############################################################################
#############################################################################

## norm_UpLeaf_vs_Coty_night

dat6 <- read.csv("norm_UPcountsLeaf_vs_Coty_night.csv")

dat6_0 <- read.csv("norm_UPcountsLeaf_vs_Coty_night.csv")
dat6

#######
# set the first column as an index

dat6 <- dat6 %>% column_to_rownames(., var = "GeneID")

dat6


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat6 <- as_tibble(scale(dat6)) 
head (dat6)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figE <- fviz_nbclust(dat6, kmeans, method = "silhouette") +
  labs(title = "")
figE

##Elbow Method
set.seed(123)

figE1 <- fviz_nbclust(dat6 , kmeans, method = "wss") +
  labs(title = "")
figE1

###Gap Statistic
set.seed(123)

gap_statE2 <- clusGap(dat6, FUN = kmeans, nstart = 2,
                      K.max = 10, B = 50)
?kmeans

figE2 <- fviz_gap_stat(gap_statE2) +
  labs(title = "")
figE2

## Combine all results in one figure
grid.arrange(figE, figE1, figE2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kE1 <- kmeans(dat6, centers = 2, nstart = 30)
kE2 <- kmeans(dat6, centers = 3, nstart = 30)
kE3 <- kmeans(dat6, centers = 4, nstart = 30)
kE4 <- kmeans(dat6, centers = 5, nstart = 30)


# plots to compare
E1 <- fviz_cluster(kE1, geom = "point", data = dat6) + ggtitle("k = 2")
E2 <- fviz_cluster(kE2, geom = "point",  data = dat6) + ggtitle("k = 3")
E3 <- fviz_cluster(kE3, geom = "point",  data = dat6) + ggtitle("k = 4")
E4 <- fviz_cluster(kE4, geom = "point",  data = dat6) + ggtitle("k = 5")

print(kE1)
print(kE2)
print(kE3)
print(kE4)
grid.arrange(E1, E2, E3, E4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalE <- kmeans(dat6, 2, nstart = 30, iter.max = 30)
print(finalE)
#Visualize final clusters

UP_Leaf_vs_CN <- fviz_cluster(finalE, geom = "point",  data = dat6) +
  labs(title = "Up_L/L_vs_CN", tag="E") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

UP_Leaf_vs_CN

#####

# Calulate cluster level means

ClustMeansE <- dat6[1:3] %>%
  mutate(Cluster = finalE$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansE
##################################################
# Visualize   

ClustMeansE %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
Up_LEAFvsCNClu_1 <- dat6_0[1] %>%
  mutate(Cluster = finalE$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

Up_LEAFvsCNClu_1

write.csv(Up_LEAFvsCNClu_1,"Clust_1_Up_LEAFvsCN.csv")

#
# Extract transcripts cluster 2
Up_LEAFvsCNClu_2 <- dat6_0[1] %>%
  mutate(Cluster = finalE$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

Up_LEAFvsCNClu_2

write.csv(Up_LEAFvsCNClu_2,"Clust_2_Up_LEAFvsCN.csv")


################################################################
################################################################

## norm_DowLeaf_vs_Coty_night

dat7 <- read.csv("norm_DowcountsLeaf_vs_Coty_night.csv")

dat7_0 <- read.csv("norm_DowcountsLeaf_vs_Coty_night.csv")
dat7

#######
# set the first column as an index

dat7 <- dat7 %>% column_to_rownames(., var = "GeneID")

dat7


##As we don’t want the clustering algorithm to depend to an arbitrary variable unit, we start by scaling/standardizing 
dat7 <- as_tibble(scale(dat7)) 
head (dat7)

##### Determine the number of cluster
##Silhouette method
set.seed(123)

figF <- fviz_nbclust(dat7, kmeans, method = "silhouette", iter.max = 20) +
  labs(title = "")
figF

##Elbow Method
set.seed(123)

figF1 <- fviz_nbclust(dat7 , kmeans, method = "wss") +
  labs(title = "")
figF1

###Gap Statistic
set.seed(123)

gap_statF2 <- clusGap(dat7, FUN = kmeans, nstart = 2,
                      K.max = 10, B = 50, iter.max = 20)
?kmeans

figF2 <- fviz_gap_stat(gap_statF2) +
  labs(title = "")
figF2

## Combine all results in one figure
grid.arrange(figF, figF1, figF2, nrow = 2)

# Compute Kmeans for K = 2 ,3 and  5
set.seed(123)

kF1 <- kmeans(dat7, centers = 2, nstart = 30)
kF2 <- kmeans(dat7, centers = 3, nstart = 30)
kF3 <- kmeans(dat7, centers = 4, nstart = 30)
kF4 <- kmeans(dat7, centers = 5, nstart = 30)


# plots to compare
F1 <- fviz_cluster(kF1, geom = "point", data = dat7) + ggtitle("k = 2")
F2 <- fviz_cluster(kF2, geom = "point",  data = dat7) + ggtitle("k = 3")
F3 <- fviz_cluster(kF3, geom = "point",  data = dat7) + ggtitle("k = 4")
F4 <- fviz_cluster(kF4, geom = "point",  data = dat7) + ggtitle("k = 5")

print(kF1)
print(kF2)
print(kF3)
print(kF4)
grid.arrange(F1, F2, F3, F4, nrow = 2)


# Compute k-means clustering with k = 2
set.seed(123)
finalF <- kmeans(dat7, 2, nstart = 30, iter.max = 30)
print(finalF)
#Visualize final clusters

Dow_Leaf_vs_CN <- fviz_cluster(finalF, geom = "point",  data = dat7) +
  labs(title = "Up_CN/L_vs_CN", tag="F") +
  theme(plot.title = element_text(hjust= 0.5),
        panel.background = element_blank(),
        axis.line = element_line())

Dow_Leaf_vs_CN

###############################################################
## Combining all results of leaf, CN and CD
# Save the legend of one of the sample
library(cowplot)
legend <- get_legend(Up_LEAFvsCot_day + 
                       theme(legend.position="bottom"))
## Remove legends of all samples
Up_LEAFvsCot_day1 <- Up_LEAFvsCot_day + theme(legend.position="none")

Up_LEAFvsCot_day1 <- Up_LEAFvsCot_day + theme(legend.position="none")
Dow_LEAFvsCot_day1 <- Dow_LEAFvsCot_day + theme(legend.position="none")
Up_CNvsCD1 <- Up_CNvsCD + theme(legend.position="none")
Dow_CNvsCD1 <- Dow_CNvsCD + theme(legend.position="none")
UP_Leaf_vs_CN1 <- UP_Leaf_vs_CN + theme(legend.position="none")
Dow_Leaf_vs_CN1 <- Dow_Leaf_vs_CN + theme(legend.position="none")

# Combine them with a common legend and # Change the legend position

grid.arrange(Up_LEAFvsCot_day1, Dow_LEAFvsCot_day1, Up_CNvsCD1, Dow_CNvsCD1,
             UP_Leaf_vs_CN1, Dow_Leaf_vs_CN1, bottom=legend$grobs[[1]], ncol=3)

###############################################################

# Calulate cluster level means

ClustMeansF <- dat7[1:3] %>%
  mutate(Cluster = finalF$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

ClustMeansF
##################################################
# Visualize   

ClustMeansF %>% 
  gather(Samples, Clust, CD:SS)%>%
  ggplot(aes(x=Samples , y = Clust, fill = Samples)) + geom_col(width = 0.5) + 
  facet_grid(.~ Cluster)+ 
  scale_fill_brewer(palette = "Accent")+
  ylab("Relative expression") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Extracting transcripts
# Extract transcripts cluster 1
Dow_LEAFvsCNClu_1 <- dat7_0[1] %>%
  mutate(Cluster = finalF$cluster) %>%
  filter(Cluster == 1) %>% ungroup() %>% count(GeneID)

Dow_LEAFvsCNClu_1

write.csv(Dow_LEAFvsCNClu_1,"Clust_1_Dow_LEAFvsCN.csv")

#
# Extract transcripts cluster 2
Dow_LEAFvsCNClu_2 <- dat7_0[1] %>%
  mutate(Cluster = finalF$cluster) %>%
  filter(Cluster == 2) %>% ungroup() %>% count(GeneID)

Dow_LEAFvsCNClu_2

write.csv(Dow_LEAFvsCNClu_2,"Clust_2_Dow_LEAFvsCN.csv")