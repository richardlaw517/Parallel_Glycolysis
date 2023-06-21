#CoCulture Linear Regression

data3<- read.csv("CoCulture_LinearRegression_N.csv", stringsAsFactors = TRUE)

#N11
intercept <- 1.014

Upshifts.N11R1.lm <- lm(I(N11R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N11R1.lm)

Upshifts.N11R2.lm <- lm(I(N11R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N11R2.lm)

Upshifts.N11R3.lm <- lm(I(N11R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N11R3.lm)

#N21
intercept <- 1.666

Upshifts.N21R1.lm <- lm(I(N21R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N21R1.lm)

Upshifts.N21R2.lm <- lm(I(N21R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N21R2.lm)

Upshifts.N21R3.lm <- lm(I(N21R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N21R3.lm)

#N12
intercept <- 0.381

Upshifts.N12R1.lm <- lm(I(N12R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N12R1.lm)

Upshifts.N12R2.lm <- lm(I(N12R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N12R2.lm)

Upshifts.N12R3.lm <- lm(I(N12R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.N12R3.lm)
