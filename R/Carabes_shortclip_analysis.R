library(rhdf5); library(tidyverse); library(viridis)
library(gganimate)

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



df %>% 
  ggplot(aes(x = X, y = Y, color = ID, group = "part")) +
  geom_point(alpha = .2) +
  scale_color_viridis(discrete = T) +
  facet_wrap(~ part) +
  coord_equal(ratio = 1) +
  theme_bw(14) 


df %>% 
  filter(part == 2) %>% # Show thorax tracks onlys
  ggplot(aes(x = X, y = Y, color = ID)) +
  geom_point(size = 4) +
  scale_color_viridis(discrete = T) +
  coord_equal(ratio = 1) +
  transition_time(t) +
  theme_bw(14) 
