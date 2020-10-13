Zeta.aux <- function(shape, qq, shift = 1) {
  
  
  
  LLL <- max(length(shape), length(qq))
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(qq   ) != LLL) qq    <- rep_len(qq,    LLL)
  
  if (any(qq < 12-1))
    warning("all values of argument 'q' should be 12 or more")
  aa <- qq
  
  
  B2 <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  kk <- length(B2)  # 8
  ans <- 1 / ((shape-1) * (shift + aa)^(shape-1)) +
    0.5 / (shift + aa)^shape
  res = NULL; res[1] = ans
  term <- (shape/2) / (shift + aa)^(shape+1)
  ans <- ans + term * B2[1]
  res[2] = ans
  for (mm in 2:kk) {
    term <- term * (shape+2*mm-3) *
      (shape+2*mm-2) / ((2*mm-1) * 2 * mm * (shift + aa)^2)
    ans <- ans + term * B2[mm]
    res[mm+1] = ans
  }
  cat(res, '\n')
  ifelse(aa - 1 <= qq, ans, rep(0, length(ans)))  # Handled above
}