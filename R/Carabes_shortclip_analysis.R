library(rhdf5); library(tidyverse); library(viridis)
library(gganimate); library(MoveR)

df = H5Fopen('data/data-raw/labels.v000.analysis.h5') 

# 1. Data import ----
df = H5Fopen('data/labels.v000.analysis.h5') 

list_df <- h5dump(df)
list_df$tracks
# 5824 frames
# 3 nodes (1: Head, 2: Thorax, Abdomen)
# X,Y
# 4 IDs

# Store all tracks
x <- list_df$tracks

# Extract head X,Y by ID
ID1.h <- x[1:5824,1,1:2,1]
ID2.h <- x[1:5824,1,1:2,2]
ID3.h <- x[1:5824,1,1:2,3]
ID4.h <- x[1:5824,1,1:2,4]
par(mfrow=c(2,2))
plot(ID1.h)
plot(ID2.h)
plot(ID3.h)
plot(ID4.h)

# Extract thorax X,Y by ID
ID1.t <- x[1:5824,2,1:2,1]
ID2.t <- x[1:5824,2,1:2,2]
ID3.t <- x[1:5824,2,1:2,3]
ID4.t <- x[1:5824,2,1:2,4]
par(mfrow=c(2,2))
plot(ID1.t)
plot(ID2.t)
plot(ID3.t)
plot(ID4.t)

# Extract abdomen X,Y by ID
ID1.a <- x[1:5824,3,1:2,1]
ID2.a <- x[1:5824,3,1:2,2]
ID3.a <- x[1:5824,3,1:2,3]
ID4.a <- x[1:5824,3,1:2,4]
par(mfrow=c(2,2))
plot(ID1.a)
plot(ID2.a)
plot(ID3.a)
plot(ID4.a)


# Combine in 1 dataframe
ID.h = data.frame(
  ID = factor(c(rep("ID1", 5824),
                rep("ID2", 5824),
                rep("ID3", 5824),
                rep("ID4", 5824))),
  part = 1,
  part.f = "head",
  t = rep(1:length(ID1.a[,1]), 4),
  X = c(ID1.h[,1],ID2.h[,1],ID3.h[,1],ID4.h[,1]),
  Y = c(ID1.h[,2],ID2.h[,2],ID3.h[,2],ID4.h[,2])
)

ID.t = data.frame(
  ID = factor(c(rep("ID1", 5824),
                rep("ID2", 5824),
                rep("ID3", 5824),
                rep("ID4", 5824))),
  part = 2,
  part.f = "thorax",
  t = rep(1:length(ID1.a[,1]), 4),
  X = c(ID1.t[,1],ID2.t[,1],ID3.t[,1],ID4.t[,1]),
  Y = c(ID1.t[,2],ID2.t[,2],ID3.t[,2],ID4.t[,2])
)


ID.a = data.frame(
  ID = factor(c(rep("ID1", 5824),
                rep("ID2", 5824),
                rep("ID3", 5824),
                rep("ID4", 5824))),
  part = 3,
  part.f = "abdomen",
  t = rep(1:length(ID1.a[,1]), 4),
  X = c(ID1.a[,1],ID2.a[,1],ID3.a[,1],ID4.a[,1]),
  Y = c(ID1.a[,2],ID2.a[,2],ID3.a[,2],ID4.a[,2])
)

df = rbind(ID.h, ID.t, ID.a)

# Export 
write.csv(file = "data/data-clean/carab_sleap_clean.csv", x = df, row.names = F)

# 2. Simple Figures ----

# Plot XY coordinates by ID
Fig1 = df %>% 
  ggplot(aes(x = X, y = Y, color = ID, group = "part")) +
  geom_point(alpha = .2) +
  scale_color_viridis(discrete = T) +
  facet_wrap(~ part) +
  coord_equal(ratio = 1) +
  theme_bw(14) +
  xlab("frame")

# Plot X-axes coordinates over time
Fig2 = df %>% 
  filter(ID == "ID1") %>% 
  ggplot(aes(x = t, y = X, group = ID, color = part.f)) +
  geom_point(alpha = .2) +
  scale_color_viridis(discrete = T) +
  # facet_wrap(~ID) +
  theme_bw(14) +
  xlab("frame")


# Animate!
Fig3 = df %>% 
  filter(part == 2) %>% # Show thorax tracks only
  ggplot(aes(x = X, y = Y, color = ID)) +
  geom_point(size = 4) +
  scale_color_viridis(discrete = T) +
  coord_equal(ratio = 1) +
  transition_time(t) +
  theme_bw(14) +
  xlab("frame") 

# Save figures
ggsave("Fig/Fig1.jpeg", plot = Fig1)
ggsave("Fig/Fig2.jpeg", plot = Fig2)
anim_save("Fig/Fig3.gif", plot = Fig3)

# 4.MoveR analysis ----
# Reimport cleaned data
tracksL <- readPlain("data/data-clean/carab_sleap_clean.csv", 
                     id = "ID", 
                     timeCol = "t", 
                     sep = ",")

# create a new identity variable considering indiv and body part
tracksL[["IDBody"]] <- paste(tracksL$identity, tracksL$part.f, sep = "_")

# convert it to a list of tracklet based on IDBody
tracks <- convert2Tracklets(tracksL, by = "IDBody")

# create a custom filter to remove NA
filterNA.x <- filterFunc(tracks,
                         toFilter = "x.pos",
                         customFunc = is.na)
# filter Na entries
tracksNoNa <- filterTracklets(tracks,
                              filter = filterNA.x)

# check the summary of the filtering process 
str(tracksNoNa[["SummaryFiltering"]])

# as filter is very conservative and split tracklets when data, here NA, are removed, we need to merge them again based on IDBody
# ! CAUTION ! HERE we are merging the tracklets again after filtering, it can result in biased computation of metrics as missing parts of the trajectories are
# replaced by a straight line.
tracksNoNaL <- convert2List(tracksNoNa[["CleanedTracklets"]])
tracksNoNa <- convert2Tracklets(tracksNoNaL, 
                                by = "IDBody")

# reorder the data according to the frame for each tracklet
tracksNoNa <- lapply(tracksNoNa, function(x) x[order(x[["frame"]]),])

# draw the tracklets
par(mfrow=c(1,1))
drawTracklets(tracksNoNa,
              timeCol = "frame",
              colId = "IDBody"#,
              # selTrack = c(1,2,3) # in case you want to display only some tracklets (here all body parts of ID1)
)

# compute speed over the different tracklets (IDBody)
tracksNoNa2 <- analyseTracklets(tracksNoNa,
                                customFunc = list(
                                  speed = function(x)
                                    speed(x, timeCol = "frame")
                                ))

# compute covariance between body part
## split the dataset according to body part
BP <- list()
for(i in c("head", "thorax", "abdomen")){
  temp <- tracksNoNa2[grep(i, names(tracksNoNa2))]
  BP[[i]] <- temp
}

## compute mean speed over time for each body parts 
SpeedRes <- list()
for (i in names(BP)) {
  SpeedRes[[i]] <- temporalTrend(
    BP[[i]],
    customFunc = list(
      speed_all =
        function(x)
          mean(x$speed, na.rm = T)
    ),
    timeCol = "frame",
    Tstep = 100, # sliding window over 100 frames
    sampling = 20,  # sampling step of 20 frames
    wtd = F
  )[["speed_all"]]
}

### plot the resulting speed trends
# plot the result
par(mfrow = c(1, 3))
for (p in seq_along(SpeedRes)) {
  plot(
    NULL,
    ylim = c(round(
      min(SpeedRes[[p]][, 1], na.rm = T),
      digits = 1
    ),
    round(
      max(SpeedRes[[p]][, 1], na.rm = T),
      digits = 1
    )),
    xlim = c(
      min(SpeedRes[[p]]$frame, na.rm = T),
      max(SpeedRes[[p]]$frame, na.rm = T)
    ),
    main = names(SpeedRes)[p],
    ylab = "speed (pixels/frame)",
    xlab = "Time (frames)"
  )
  lines(SpeedRes[[p]][, 1] ~ SpeedRes[[p]]$frame , col = "darkred")
}

## Compute the covariance between body parts (should add a function to MoveR but in the meantime the following should do the trick)

### remove the NA inserted at the beginning of each sliding window
SpeedResNoNa <- lapply(SpeedRes, function(x) na.omit(x))

### specify some starting parameter (the sliding window size in frame and the combination to do)
combination <- combn(names(SpeedResNoNa), 2)
Window = 10

### let's compute covariance between combinations of variables selected
CovRes = list()
for (j in seq(ncol(combination))) {
  ##### initialize an empty vector to retrieve the results
  VecL <- length(SpeedResNoNa[[combination[1, j]]]$speed_all)
  covResVec <-
    rep(NA, VecL)
  
  ##### loop over the dataset to perform covariance computation over the whole dataset
  for (i in 1:VecL) {
    start <- max(1, i - Window %/% 2)
    end <- min(length(x), i + Window %/% 2)
    if(i < Window %/% 2 | i > VecL - Window %/% 2){
      next
    }
    windowX <- SpeedResNoNa[[combination[1, j]]]$speed_all[start:end]
    windowY <- SpeedResNoNa[[combination[2, j]]]$speed_all[start:end]
    covResVec[i] <- cov(windowX, windowY, use = "complete.obs")
  }
  CovRes[[paste(combination[1, j], combination[2, j], sep = "_")]] <-
    covResVec
}

### plot the resulting covariances 
par(mfrow = c(1, 3))
for (p in seq_along(CovRes)) {
  plot(
    NULL,
    ylim = c(round(
      min(CovRes[[p]], na.rm = T),
      digits = 1
    ),
    round(
      max(CovRes[[p]], na.rm = T),
      digits = 1
    )),
    xlim = c(
      min(SpeedResNoNa[[1]]$frame, na.rm = T),
      max(SpeedResNoNa[[1]]$frame, na.rm = T)
    ),
    main = names(CovRes)[p],
    ylab = "covariance",
    xlab = "Time (frames)"
  )
  lines(CovRes[[p]] ~ SpeedResNoNa[[1]]$frame , col = "darkred")
}