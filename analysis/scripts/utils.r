

# takes matrix same dimensions as cmat as n c.calc() above


stacker <- function(z, row = TRUE) {
  if (row == TRUE) {
    zz <- t(z[1, ])
    for (i in 2:length(z[, 1])) {
      zz <- rbind(zz, t(z[i, ]))
    }
  } else {
    zz <- paste(z[, 1])
    for (i in 2:length(z[1, ])) {
      zz <- c(zz, paste(z[, i]))
    }
    t(t(zz))
  }
  return(zz)
}
