matchedControl <- read.csv("Final_Set_of_Matched_Controls.csv")

cat("Starting")
count <- 1
totalP <- 0


while (1) {
  count <- count + 1
  matchCopy <- matchedControl
  if ((count/100) == round(count/100,0)) {
    cat("Attempt number ", count, "\n")
  }
  
  oddAdjust <- sample(c(1,-1), 1)
  cancers <- filter(matchedControl, FinalDiagnosis == "Cancer")
  im <- filter(matchedControl, FinalDiagnosis == "IM")
  hg <- filter(matchedControl, FinalDiagnosis == "HG")
  ndbe <- filter(matchedControl, FinalDiagnosis == "NDBE")
  normal <- filter(matchedControl, FinalDiagnosis == "Normal")
  hv <- filter(matchedControl, FinalDiagnosis == "HV")
  
  samplesize <- c(nrow(cancers)/2, nrow(im) / 2, nrow(hg) / 2, nrow(ndbe) / 2, nrow(normal) / 2, nrow(hv) / 2)
  cancerChange <- sample(c(0.5, -0.5), 2, replace=F)
  samplesize[1:2] <- samplesize[1:2] - cancerChange
  
  for (i in 1:10) {
    jitterChange <- sample(1:4, 2, prob = c(1, 2, 4, 4), replace = F)
    samplesize[2 + jitterChange[1]] <- samplesize[2 + jitterChange[1]] + 0.5
    samplesize[2 + jitterChange[2]] <- samplesize[2 + jitterChange[2]] - 0.5
  }
  
  adjust <- which(samplesize != round(samplesize,0))
  adjustChange <- sample(c(0.5, -0.5), 2, replace = F)
  samplesize[adjust] <- samplesize[adjust] + adjustChange
  
  if (sum(samplesize) != 96 ) { stop( "Sample Size error\n")}
  
  Csampl <- sample_n(cancers, samplesize[1]) %>% select(Subject.number) %>% unlist() %>% as.character()
  IMsampl <- sample_n(im, samplesize[2]) %>% select(Subject.number) %>% unlist() %>% as.character()
  HGsampl <- sample_n(hg, samplesize[3]) %>% select(Subject.number) %>% unlist() %>% as.character()
  NDBEsampl <- sample_n(ndbe, samplesize[4]) %>% select(Subject.number) %>% unlist() %>% as.character()
  normsampl <- sample_n(normal, samplesize[5]) %>% select(Subject.number) %>% unlist() %>% as.character()
  hvsampl <- sample_n(hv, samplesize[6]) %>% select(Subject.number) %>% unlist() %>% as.character()
  
  matchCopy[,'Group'] <- 2
  
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% Csampl, 1))
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% IMsampl, 1))
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% HGsampl, 1))
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% NDBEsampl, 1))
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% normsampl, 1))
  matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% hvsampl, 1))
  
  
  #matchCopy <- mutate(matchCopy,  Group = replace(Group, Subject.number %in% NoCsampl, 1))
  
  
  CancerGenderP <- filter(matchCopy, FinalDiagnosis == "Cancer") %>% 
    summarise(pval = chisq.test(Gender, Group)$p.value)
  if (CancerGenderP < 0.05) { 
    #    cat("Gender fail\n")
    next()
  }
  
  genderP <- matchCopy %>% 
    summarise(pval = chisq.test(Gender, Group)$p.value)
  if (genderP < 0.05) { 
    #    cat("Gender fail\n")
    next()
  }
  
  DiagnosisP <- matchCopy %>% 
    summarise(pval = chisq.test(FinalDiagnosis, Group)$p.value)
  
  if (DiagnosisP < 0.05) { 
    #    cat("Diagnosis fail\n")
    next()
  }
  
  ageT <- t.test(age ~ Group, data = matchCopy)$p.value
  if (ageT < 0.05) {
    #cat ("Age T fail\n")
    next()
  }
  
  bmiT <- t.test(bmi ~ Group, data = matchCopy)$p.value
  if (bmiT < 0.05) {
    #cat("BMI t fail\n")
    next()
  }
  
  cigT <- t.test(CigsSmoked ~ Group, data = matchCopy)$p.value
  if (cigT < 0.05) {
    #cat("CIG t fail\n")
    next()
  }
  
  PPIT <- t.test(PPI ~ Group, data = matchCopy)$p.value
  if (PPIT < 0.05) {
    #cat("PPI t fail\n")
    next()
  }
  
  hbT <- t.test(hb ~ Group, data = matchCopy)$p.value
  if (hbT < 0.05) {
    #cat("Heartburn t fail\n")
    next()
  }
  
  newtotalP <- genderP * 4 + CancerGenderP * 4 + DiagnosisP * 5 + ageT * 3 + bmiT * 3 + cigT + PPIT + hbT
  if (newtotalP > totalP) {
    finalGroup <- matchCopy
    cat("We have a better match!!!\n")
    write.csv(finalGroup, "BestMatch.csv", row.names=F)
    totalP <- newtotalP
  }
}