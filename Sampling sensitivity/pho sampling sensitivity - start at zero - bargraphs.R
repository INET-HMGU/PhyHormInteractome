# Stefan Altmann
# estimate sampling sensitivity of Y2H screen
# this script was executed using R3.6.1 and Windows 10

rm(list = ls())

setwd("~/../Desktop/Github/Sampling sensitivity")

library(openxlsx)

library(drc) # for fitting Michaelis Menten model
library(ggplot2) # for drawing

r.dir <- "Results/"
if(!dir.exists(r.dir)){
  dir.create(r.dir)
}

numbers.combined <- read.xlsx("Verified Interactions.xlsx", sheet = 1)

colnames(numbers.combined) <- c("numReps", "numInt")

mm.comb <- structure(list(S = c(3,2,2,2,1,1,1,0), v = c(rev(numbers.combined$numInt),0)), 
                     .Names = c("S", "v"), class = "data.frame", row.names = c(NA, -8L))

model.drm.comb <- drm (v ~ S, data = mm.comb, fct = MM.2())

print("Model")
print(summary(model.drm.comb))
print("Coefficients")
print(coef(model.drm.comb))

mml2.comb <- data.frame(S = seq(0, 10, length.out = 11))

mml2.comb$v <- predict(model.drm.comb, newdata = mml2.comb)

# save data of the model to xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "prediction")
openxlsx::writeData(wb, sheet = "prediction", mml2.comb, colNames = T, rowNames = T)

openxlsx::saveWorkbook(wb, file = paste(r.dir,"/Michaelis Menten Prediction.xlsx", sep = ""), overwrite = T)


# plot estimation for 10 repeats as line with measured values as points

ggplot(mm.comb, aes(x = S, y = v)) +
  theme_bw() +
  xlab("Repeats") +
  ylab("Interactions") +
  ggtitle("Sampling sensitivity PHO screen") +
  geom_point(alpha = 0.5) +
  geom_line(data = mml2.comb, aes(x = S, y = v), colour = "red")


ggsave(paste(r.dir,"PHI Sampling Sensitivity - Combined Screens.pdf", sep = ""), width = 6, height = 4)

### plot estimation for 10 repeats as barplot
### first 4 bars are measured values, following 7 grey bars are estimation

repeats <- seq(0,10)

meanValues <- c(0)
sdValues <- c(0)
for(i in 1:3){
  mean1 <- mean(mm.comb[which(mm.comb$S == i),2])
  sd1 <- sd(mm.comb[which(mm.comb$S == i),2])
  meanValues <- c(meanValues, mean1)
  sdValues <- c(sdValues, sd1)
}


my.means <- c(meanValues,mml2.comb$v[(2+3):11])

my.sds <- sdValues
my.sds[which(is.na(my.sds))] <- 0
for(k in (2+3):11){
  sd.temp <- (37.64418 * mml2.comb$v[k]/coef(model.drm.comb)[1])
  my.sds <- c(my.sds, sd.temp)
}
names(my.sds) <- NULL

my.type <- c(rep("measure",3+1),rep("estimation",10-3))

my.data <- data.frame(repeats, my.means, my.sds, my.type)

p <- ggplot(my.data, aes(x=repeats, y=my.means, fill=my.type)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=my.means-my.sds, ymax=my.means+my.sds), width=.0,
                    position=position_dodge(.9))+
  # scale_fill_manual(values=c('lightgray','darkorange'))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme_classic()
plot(p)
ggsave(paste(r.dir,"PHI Sampling Sensitivity - Combined Screens BARPLOT 2.pdf", sep = ""), width = 6, height = 4)


p <- p + geom_segment(x = 5.1, y = 529, xend = 5.1, yend = 529,
             arrow = arrow(length = unit(0.2, "cm")), color = "red", size = 0.5)
p <- p + geom_segment(x = 4.9, y = 529, xend = 4.9, yend = 529,
                 arrow = arrow(length = unit(-0.2, "cm")), color = "red", size = 0.5)
ggsave(paste(r.dir,"PHI Sampling Sensitivity - Combined Screens BARPLOT 2 incl value.pdf", sep = ""), width = 6, height = 4)
print(p)

q <- p + geom_point(data = mm.comb, aes(x = S, y = v),alpha = 0.5, fill = "red")
print(q)

ggsave(paste(r.dir,"PHI Sampling Sensitivity - Combined Screens BARPLOT 2 with points.pdf", sep = ""), width = 6, height = 4)



#### estimation as gray area with measured values as points
error.range <- c()
error.range <- c(error.range, 0)
for(i in 2:length(mml2.comb$S)){
  my.e <- 37.64418 * (mml2.comb$v[i]/615.95672)
  error.range <- c(error.range, my.e)
}

lower.values <- mml2.comb$v - error.range
upper.values <- mml2.comb$v + error.range

ids <- factor(c("1"))

values <- data.frame(
  id = ids,
  value = c(1)
)

positions <- data.frame(
  id = rep(ids,22),
  x = c(seq(0,10),rev(seq(0,10))),
  y = c(lower.values, rev(upper.values))
)

datapoly <- merge(values, positions, by = c("id"))

p <- ggplot(datapoly, aes(x = x, y = y, group=ids)) +
  geom_polygon(fill = 'gray') +
  geom_line(data=mml2.comb, aes(x=S, y=v))+
  geom_point(data = mm.comb, aes(x = S, y = v), color = 'darkorange')+
  theme_classic()

p <- p + geom_point(x=5, y=529, color = "blue")

print(p)

ggsave(paste(r.dir,"PHI Sampling Sensitivity - Range Polygon.pdf", sep = ""), width = 6, height = 4)

###############################
# estimation as boxplot
# measured values as points

p <- ggplot(mml2.comb, aes(x=S, y=v), fill = "#999999") + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=v-error.range, ymax=v+error.range), width=.0,
                position=position_dodge(.9))+
  # scale_fill_manual(values=c('lightgray','darkorange'))+
  # scale_fill_manual(values=c("#999999", "#E69F00"))+
  theme_classic()
p <- p + geom_point(mm.comb, mapping=aes(x = S, y = v),fill='blue', size = 2)
p <- p + geom_point(x=5, y=529, color = "red", size = 2)

print(p)


ggsave(paste(r.dir,"PHI Sampling Sensitivity - bars are estimation, points measured.pdf", sep = ""), width = 6, height = 4)






