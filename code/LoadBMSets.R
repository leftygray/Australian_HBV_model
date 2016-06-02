









##Load files
bgmort0to4 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/bgmort 0to4.csv", sep=",", header=FALSE)
bgmort5to14 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/bgmort 5to14.csv", sep=",", header=FALSE)
bgmort15to44 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/bgmort 15to44.csv", sep=",", header=FALSE)
bgmort45 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/bgmort 45plus.csv", sep=",", header=FALSE)
births <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/births.csv", sep=",", header=FALSE)
mig0to4chronic <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 0to4 chronic.csv", sep=",", header=FALSE)
mig0to4cleared <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 0to4 cleared.csv", sep=",", header=FALSE)
mig0to4sus <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 0to4 sus.csv", sep=",", header=FALSE)
mig5to14chronic <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 5to14 chronic.csv", sep=",", header=FALSE)
mig5to14cleared <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 5to14 cleared.csv", sep=",", header=FALSE)
mig5to14sus <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 5to14 sus.csv", sep=",", header=FALSE)
mig15to44chronic <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 15to44 chronic.csv", sep=",", header=FALSE)
mig15to44cleared <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 15to44 cleared.csv", sep=",", header=FALSE)
mig15to44sus <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 15to44 sus.csv", sep=",", header=FALSE)
mig45chronic <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 45plus chronic.csv", sep=",", header=FALSE)
mig45cleared <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 45plus cleared.csv", sep=",", header=FALSE)
mig45sus <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/mig 45plus sus.csv", sep=",", header=FALSE)
notifications_cum2014 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/notifications-cum2014.csv", sep=",", header=FALSE)
notifications2014 <- read.csv("/Volumes/public/SERP/Projects/BRISE/HBV/data/BM datasets/notifications2014.csv", sep=",", header=FALSE)

#merge all datasets with complete years 1951-2050
BMsets <- data.frame(bgmort0to4$V1, bgmort0to4$V2, bgmort5to14$V2, bgmort15to44$V2, bgmort45$V2, births$V2[1:100], mig0to4chronic$V2, mig0to4cleared$V2, mig0to4sus$V2, mig5to14chronic$V2, mig5to14cleared$V2, mig5to14sus$V2, mig15to44chronic$V2, mig15to44cleared$V2, mig15to44sus$V2, mig45chronic$V2, mig45cleared$V2, mig45sus$V2)

#rename column names
colnames(BMsets) <- c("Year", "bgmort0to4", "bgmort5to14", "bgmort15to44", "bgmort45", "births", "mig0to4chronic", "mig0to4cleared", "mig0to4sus", "mig5to14chronic", "mig5to14cleared", "mig5to14sus", "mig15to44chronic", "mig15to44cleared", "mig15to44sus", "mig45chronic", "mig45cleared", "mig45sus")

write.csv(BMsets, "/Volumes/public/SERP/Projects/BRISE/HBV/data/BMdataset.csv", row.names=FALSE)
