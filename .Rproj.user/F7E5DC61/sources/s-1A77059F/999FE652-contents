convert_IDs <- function(dat, convert.IDs.numeric = FALSE) {

  x0 <- which(!duplicated(dat$family))
  x1 <- c(x0[-1], length(dat$family) + 1)
  # t gives the size of all families in the correct order
  t <- x1 - x0
  famnum <- rep.int(1:length(t), times = t)
  offset <- rep(0, nrow(dat))
  offset[1 + cumsum(t[-length(t)])] <- t[-length(t)]
  indiv <- 1:nrow(dat) - cumsum(offset)
  max.fam.size <- max(t)
  mother <- numeric(nrow(dat))
  father <- numeric(nrow(dat))
  mother[dat$mother == ""] <- -1
  father[dat$father == ""] <- -1
  for (i in 1:max.fam.size) {
    rootsi <- dat$indiv[indiv == i]
    mother[dat$mother %in% rootsi] <- i
    father[dat$father %in% rootsi] <- i
  }
  indiv <- as.character(indiv)
  mother <- as.character(mother)
  father <- as.character(father)
  mother[mother == "-1"] <- ""
  father[father == "-1"] <- ""
  dat$indiv <- indiv
  dat$mother <- mother
  dat$father <- father
  if (convert.IDs.numeric) {
    dat$indiv <- as.numeric(dat$indiv)
    dat$mother[dat$mother == ""] <- "0"
    dat$father[dat$father == ""] <- "0"
    dat$mother <- as.numeric(dat$mother)
    dat$father <- as.numeric(dat$father)
  }
  dat
}

