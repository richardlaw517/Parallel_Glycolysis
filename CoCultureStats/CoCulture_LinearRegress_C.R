#CoCulture Linear Regression

data3<- read.csv("CoCulture_LinearRegression_C.csv", stringsAsFactors = TRUE)

#C11
intercept <- 1                                 

Upshifts.C11R1.lm <- lm(I(C11R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C11R1.lm)

Upshifts.C11R2.lm <- lm(I(C11R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C11R2.lm)

Upshifts.C11R3.lm <- lm(I(C11R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C11R3.lm)

#C21
intercept <- 2

Upshifts.C21R1.lm <- lm(I(C21R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C21R1.lm)

Upshifts.C21R2.lm <- lm(I(C21R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C21R2.lm)

Upshifts.C21R3.lm <- lm(I(C21R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C21R3.lm)

#C12
intercept <- 0.5

Upshifts.C12R1.lm <- lm(I(C12R1 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C12R1.lm)

Upshifts.C12R2.lm <- lm(I(C12R2 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C12R2.lm)

Upshifts.C12R3.lm <- lm(I(C12R3 - intercept) ~ 0 + Upshifts, data = data3)
summary(Upshifts.C12R3.lm)
