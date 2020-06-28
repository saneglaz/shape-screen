  ##################################
  # 05_motif_detector.py.R
  #################################
  
  # Run motif_detector.py module
  list.train.chip <- list.files("gibbs_input/", pattern="positives_training_shape.dat")
  list.train.ctrl <- list.files("gibbs_input/", pattern="negatives_training_shape.dat")
  list.test.chip <- list.files("gibbs_input/", pattern="positives_validation_shape.dat")
  list.test.ctrl <- list.files("gibbs_input/", pattern="negatives_validation_shape.dat")
  path <- "/home/lbuo/chipistrello/myproject/"
  text <- NULL
  for(i in 1:length(list.train.chip)) {
    name <- stringr::str_split_fixed(list.train.chip[i], "[_]",2)[1,1]
    motif_detector.py <- paste("python2 ", getwd(), "/shape-motif/motif_detector.py", sep="")
    train.chip <- paste(getwd(), "/gibbs_input/", list.train.chip[i], sep="")
    train.ctrl <- paste(getwd(), "/gibbs_input/", list.train.ctrl[i], sep="")
    test.chip <- paste(getwd(), "/gibbs_input/", list.test.chip[i], sep="")
    test.ctrl <- paste(getwd(), "/gibbs_input/", list.test.ctrl[i], sep="")
    output <- paste(getwd(), "/motif_detector/", name, sep="")
    feature <- "ProT"
    train.chip.instance <- paste(getwd(), "/gibbs_out/", str_sub(list.train.chip[i],1,(str_length(list.train.chip[i])-4)), ".instance", sep="")
    train.ctrl.instance <- paste(getwd(), "/gibbs_out/", str_sub(list.train.ctrl[i],1,(str_length(list.train.ctrl[i])-4)), ".instance", sep="")
    test.chip.instance <- paste(getwd(), "/gibbs_out/", str_sub(list.test.chip[i],1,(str_length(list.test.chip[i])-4)), ".instance", sep="")
    test.ctrl.instance <- paste(getwd(), "/gibbs_out/", str_sub(list.test.ctrl[i],1,(str_length(list.test.ctrl[i])-4)), ".instance", sep="")
    
    x <- paste(motif_detector.py, 
               train.chip, 
               train.ctrl, 
               test.chip, 
               test.ctrl,
               output,
               feature,
               train.chip.instance,
               train.ctrl.instance,
               test.chip.instance,
               test.ctrl.instance,
               sep=" ")
    #system(x)
    text <- c(text, x)
  }
  writeLines(c("#!/bin/bash", text[1:30]), "scripts/motif_detector_1.sh")  
  writeLines(c("#!/bin/bash", text[31:60]), "scripts/motif_detector_2.sh")  
  writeLines(c("#!/bin/bash", text[61:90]), "scripts/motif_detector_3.sh")  
  
  # Run motif_detector scripts