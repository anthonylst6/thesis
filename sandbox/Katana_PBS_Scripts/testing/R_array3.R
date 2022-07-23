x <- rnorm(100)
plot(x)# Opening the graphical device
pdf("my_plot.pdf")

# Creating a plot
plot(rnorm(50))

# Closing the graphical device
dev.off() 

