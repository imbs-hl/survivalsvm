

svm <- survivalsvm(Surv(time, status) ~ ., veteran, gamma.mu = 0.1)

test_that("survival svm returns class survivalsvm", {
  expect_is(svm, "survivalsvm")
})
