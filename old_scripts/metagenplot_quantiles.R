genes = quantile1.genes
genes = quantile2.genes
genes = quantile3.genes
genes = quantile4.genes
genes = quantile5.genes




```{r all_k27_metageneplots}
#jpeg(filename = "images/h3k27me3_metageneplots.jpg",res = 250,width = 900,height = 900)
plot(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body.q1,spar = 0.5),     #quantile1
     type="l",col="blue",lwd=3,ylim=c(0,45),axes=F,xlab="",ylab="")
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body.q1,spar = 0.5),
      type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss.q1,spar = 0.5),
      type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes.q1,spar = 0.5),
      type="l",col="red",lwd=3),

lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body.q2,spar = 0.5),   #quantile2
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body.q2,spar = 0.5),
      type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss.q2,spar = 0.5),
      type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes.q2,spar = 0.5),
      type="l",col="red",lwd=3),

lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body.q3,spar = 0.5),   #quantile3
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body.q3,spar = 0.5),
            type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss.q3,spar = 0.5),
            type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes.q3,spar = 0.5),
            type="l",col="red",lwd=3),
      
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body.q4,spar = 0.5),   #quantile4
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body.q4,spar = 0.5),
                  type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss.q4,spar = 0.5),
                  type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes.q4,spar = 0.5),
                  type="l",col="red",lwd=3),    

lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.gene.body.q5,spar = 0.5),   #quantile5
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.internal.gene.body.q4,spar = 0.5),
            type="l",col="cyan",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tss.q5,spar = 0.5),
            type="l",col="green",lwd=3)
lines(smooth.spline(x=1:40,y=median.metageneplot.data.k27.tes.q5,spar = 0.5),
            type="l",col="red",lwd=3), 

axis(side = 1, at=c(0,10,30,40),labels = c("-2Kb","TSS","TES","2Kb"),las=2,lwd=2)
axis(side = 2,lwd=2)
#dev.off()
```