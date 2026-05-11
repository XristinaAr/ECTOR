## load data
temp.data <- read.table(file="temperature_growth_data.tsv",header=T, sep="\t")

head(temp.data)

##temp 10
t10.0 <- temp.data$temp10[temp.data$Time == 0]
t10.24 <- temp.data$temp10[temp.data$Time == 24]
t10.48 <- temp.data$temp10[temp.data$Time == 48]
t10.72 <- temp.data$temp10[temp.data$Time == 72]
t10.96 <- temp.data$temp10[temp.data$Time == 96]
t10.168 <- temp.data$temp10[temp.data$Time == 168]

t10.mean <- c(mean(t10.0),
              mean(t10.24),
              mean(t10.48),
              mean(t10.72),
              mean(t10.96),
              mean(t10.168))

t10.sd <- c(sd(t10.0),
            sd(t10.24),
            sd(t10.48),
            sd(t10.72),
            sd(t10.96),
            sd(t10.168))

##temp 14
t14.0 <- temp.data$temp14[temp.data$Time == 0]
t14.24 <- temp.data$temp14[temp.data$Time == 24]
t14.48 <- temp.data$temp14[temp.data$Time == 48]
t14.72 <- temp.data$temp14[temp.data$Time == 72]
t14.96 <- temp.data$temp14[temp.data$Time == 96]
t14.168 <- temp.data$temp14[temp.data$Time == 168]

t14.mean <- c(mean(t14.0),
              mean(t14.24),
              mean(t14.48),
              mean(t14.72),
              mean(t14.96),
              mean(t14.168))

t14.sd <- c(sd(t14.0),
            sd(t14.24),
            sd(t14.48),
            sd(t14.72),
            sd(t14.96),
            sd(t14.168))

##temp 20
t20.0 <- temp.data$temp20[temp.data$Time == 0]
t20.24 <- temp.data$temp20[temp.data$Time == 24]
t20.48 <- temp.data$temp20[temp.data$Time == 48]
t20.72 <- temp.data$temp20[temp.data$Time == 72]
t20.96 <- temp.data$temp20[temp.data$Time == 96]
t20.168 <- temp.data$temp20[temp.data$Time == 168]

t20.mean <- c(mean(t20.0),
              mean(t20.24),
              mean(t20.48),
              mean(t20.72),
              mean(t20.96),
              mean(t20.168))

t20.sd <- c(sd(t20.0),
            sd(t20.24),
            sd(t20.48),
            sd(t20.72),
            sd(t20.96),
            sd(t20.168))

##temp 24
t24.0 <- temp.data$temp24[temp.data$Time == 0]
t24.24 <- temp.data$temp24[temp.data$Time == 24]
t24.48 <- temp.data$temp24[temp.data$Time == 48]
t24.72 <- temp.data$temp24[temp.data$Time == 72]
t24.96 <- temp.data$temp24[temp.data$Time == 96]
t24.168 <- temp.data$temp24[temp.data$Time == 168]

t24.mean <- c(mean(t24.0),
              mean(t24.24),
              mean(t24.48),
              mean(t24.72),
              mean(t24.96),
              mean(t24.168))

t24.sd <- c(sd(t24.0),
            sd(t24.24),
            sd(t24.48),
            sd(t24.72),
            sd(t24.96),
            sd(t24.168))


##temp 26
t26.0 <- temp.data$temp26[temp.data$Time == 0]
t26.24 <- temp.data$temp26[temp.data$Time == 24]
t26.48 <- temp.data$temp26[temp.data$Time == 48]
t26.72 <- temp.data$temp26[temp.data$Time == 72]
t26.96 <- temp.data$temp26[temp.data$Time == 96]
t26.168 <- temp.data$temp26[temp.data$Time == 168]

t26.mean <- c(mean(t26.0),
              mean(t26.24),
              mean(t26.48),
              mean(t26.72),
              mean(t26.96),
              mean(t26.168))

t26.sd <- c(sd(t26.0),
            sd(t26.24),
            sd(t26.48),
            sd(t26.72),
            sd(t26.96),
            sd(t26.168))

time <- c(0,24,48,72,96,168)
jpeg(filename = "o_tauri_temp_growth_curves.jpeg",width = 470,height = 450)
plot(x = time,y=t20.mean,type = "l",lwd=2,ylim=c(0,12),col="white",xlab="Hours",ylab="Chl ug/ml",cex.lab=1.3)

polygon(x=c(time,rev(time)),
        y=c(t26.mean + t26.sd,rev(t26.mean - t26.sd)), lty="dashed",
        col = adjustcolor(colorRampPalette(c("blue", "red"))(4)[4], 
                          alpha.f=0.5))
lines(x=time, y=t26.mean, type="l",lwd=2,col=colorRampPalette(c("blue", "red"))(4)[4])


polygon(x=c(time,rev(time)),
        y=c(t20.mean + t20.sd,rev(t20.mean - t20.sd)), lty="dashed",
        col = adjustcolor(colorRampPalette(c("blue", "red"))(4)[3], 
                          alpha.f=0.5))
lines(x=time, y=t20.mean, type="l",lwd=2,col=colorRampPalette(c("blue", "red"))(4)[3])


polygon(x=c(time,rev(time)),
        y=c(t14.mean + t14.sd,rev(t14.mean - t14.sd)), lty="dashed",
        col = adjustcolor(colorRampPalette(c("blue", "red"))(4)[2], 
                          alpha.f=0.5))
lines(x=time, y=t14.mean, type="l",lwd=2,col=colorRampPalette(c("blue", "red"))(4)[2])

polygon(x=c(time,rev(time)),
        y=c(t10.mean + t10.sd,rev(t10.mean - t10.sd)), lty="dashed",
        col = adjustcolor(colorRampPalette(c("blue", "red"))(4)[1], 
                          alpha.f=0.5))
lines(x=time, y=t10.mean, type="l",lwd=2,col=colorRampPalette(c("blue", "red"))(4)[1])

legend("topleft",border = NA,legend = c("26ºC","20ºC","14ºC","10ºC"),lwd=3, col = rev(colorRampPalette(c("blue", "red"))(4)))
dev.off()
