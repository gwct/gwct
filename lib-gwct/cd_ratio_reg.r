args=(commandArgs(TRUE))
#c = c(1,2,3,4,5,6,7,8,9)
#d = c(2,4,6,8,10,12,14,16,18)

#c = scan(file="C:/Users/Gregg/Desktop/Projects/dev/gwct-dev/c_test.txt")
#d = scan(file="C:/Users/Gregg/Desktop/Projects/dev/gwct-dev/d_test.txt")
script_outdir=args[1]

c = scan(file=paste(script_outdir, "/", args[2], sep=""))
d = scan(file=paste(script_outdir, "/", args[3], sep=""))

cd_reg = lm(c ~ d)
png(paste(script_outdir, "/", "cd_plot.png", sep=""),height=1024, width=1024)
par(mar=c(6,7,4,2)+0.1, mgp=c(4.5,1,0))
plot(d,c,xlab="Divergent substitutions",ylab="Convergent substitutions",pch=19,cex.lab=3,cex.axis=2,cex=1.5)
abline(cd_reg)
legend("topleft",bty="n",legend=paste("R^2 = ", format(summary(cd_reg)$adj.r.squared, digits=4)),cex=2)
dev.off()