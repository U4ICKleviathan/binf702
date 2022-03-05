sensitivity <- c(0,
                 0.567,
                 0.791,
                 0.866,
                 0.940,
                 1.000)
specificity <- c(1,
                 0.970,
                 0.818,
                 0.697,
                 0.182,
                 0.000)


plot(1-specificity, sensitivity, type="b",
     main="ROC for PACS film")



text((1-specificity)+c(.05, 
), sensitivity, paste("(",round(1-specificity, 2),", ", round(sensitivity, 2), ")", sep=""), cex=0.8) 