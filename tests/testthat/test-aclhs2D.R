test_that("aclhs2D sampling test", {
  df <- read.csv("csvs/data_2D.csv", header=T, na.strings="-999.25")
  ns <- 50
  w <- c(10, 1000, 0.001)
  it <- 100
  vp <- aclhs.vario_params(10, 0, 90, 1)
  sd <- 1234

  function_answer <- aclhs(df=df, num_samples=ns, weights=w, iter=it,
                             vario_params=vp, seed=sd)

  temp_arr <- c(544,31,258,805,744,810,344,133,343,590,469,558,514,108,71,435,
                132,99,719,211,383,421,200,255,520,318,326,269,53,209,408,18,
                359,443,355,148,297,527,820,646,737,870,872,49,618,577,892,723,
                487,457)

  correct_answer <- NULL
  for (i in 1:length(temp_arr)) {
    correct_answer <- rbind(correct_answer, temp_arr[i])
  }

  expect_equal(function_answer, correct_answer)
})
