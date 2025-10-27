test_that("aclhs1D sampling test", {
  df <- read.csv("csvs/data_1D.csv", header=T, na.strings="-999.25")
  ns <- 48
  w <- c(0.3,100,1)
  it <- 100
  vp <- aclhs.vario_params(8, 0, 90, 1)
  sd <- 1234

  function_answer <- aclhs(df=df, num_samples=ns, weights=w, iter=it,
                             vario_params=vp, seed=sd)

  temp_arr <- c(44,82,42,66,69,27,68,32,37,84,324,87,73,70,307,46,25,94,
                266,261,265,259,258,336,341,118,121,145,321,151,275,349,
                354,242,155,253,240,152,149,206,211,204,203,188,212,183,
                193,201)

  correct_answer <- NULL
  for (i in 1:length(temp_arr)) {
    correct_answer <- rbind(correct_answer, temp_arr[i])
  }

  expect_equal(function_answer, correct_answer)
})
