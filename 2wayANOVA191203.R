http://www.sthda.com/english/wiki/two-way-anova-test-in-r

[1] "AQP1.csv"  "DTX1.csv"  "DTX1.xlsx" "Hey1.csv"  "Hey1.xlsx"
setwd("C:/Users/user/Desktop/ABLstats")

AQP1=read.csv("AQP1.csv",header=T)
DTX1=read.csv("DTX1.csv",header=T)
Hey1=read.csv("Hey1.csv",header=T)

AQP1_26=AQP1[grep("01S-026",AQP1$patient),]
AQP1_26
AQP1_27=AQP1[grep("01S-027",AQP1$patient),]
AQP1_27

DTX1_26=DTX1[grep("01S-026",DTX1$patient),]
DTX1_26
DTX1_27=DTX1[grep("01S-027",DTX1$patient),]
DTX1_27

Hey1_26=Hey1[grep("01S-026",Hey1$patient),]
Hey1_26
Hey1_27=Hey1[grep("01S-027",Hey1$patient),]
Hey1_27

png(filename = "boxplotAQP1_26.png")
ggboxplot(AQP1_26, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotAQP1_27.png")
ggboxplot(AQP1_27, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotDTX1_26.png")
ggboxplot(DTX1_26, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotDTX1_27.png")
ggboxplot(DTX1_27, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotAQP1_26.png")
ggboxplot(AQP1_26, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotAQP1_27.png")
ggboxplot(AQP1_27, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()


png(filename = "boxplotHey1_26.png")
ggboxplot(Hey1_26, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()

png(filename = "boxplotHey1_27.png")
ggboxplot(Hey1_27, x="experiment", y="foldchange", color = "time", palette = c("#00AFBB", "#E7B800"))
dev.off()



AQP1_26.aov=aov(foldchange ~ time + experiment, data=AQP1_26)
AQP1_27.aov=aov(foldchange ~ time + experiment, data=AQP1_27)

DTX1_26.aov=aov(foldchange ~ time + experiment, data=DTX1_26)
DTX1_27.aov=aov(foldchange ~ time + experiment, data=DTX1_27)

Hey1_26.aov=aov(foldchange ~ time + experiment, data=Hey1_26)
Hey1_27.aov=aov(foldchange ~ time + experiment, data=Hey1_27)

summary(AQP1_26.aov)
summary(AQP1_27.aov)

summary(DTX1_26.aov)
summary(DTX1_27.aov)

summary(Hey1_26.aov)
summary(Hey1_27.aov)

