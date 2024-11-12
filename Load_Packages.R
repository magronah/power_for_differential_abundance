pkgs <- scan("pkgs.txt", what = character(1))
lapply(pkgs, library, character.only = TRUE)

theme_set(theme_bw())
 

