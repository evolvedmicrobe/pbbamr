library(pbbamr)

test_that("combineDataFrames", {
  # Simple combine with a name vector
  a = data.frame(var1 = 1:3, var2 = 2:4)
  cnames = c("Go1", "Go2")
  l = list(a, a)
  all = combineConditions(l, cnames)
  expect_equal(nrow(all), 6)

  # Now let's try it with some new names
  newNames = c("Stop1", "Stop2")
  names(l) = newNames
  all = combineConditions(l)
  expect_equal(levels(all$Condition), newNames)

  # Now verify we didn't reset the default
  difNames = c("Slow1", "Slow2")
  names(l) <- difNames
  na = combineConditions(l)
  expect_equal(levels(na$Condition), difNames)

  all = combineConditions(l, cnames)
  expect_equal(levels(all$Condition), cnames)
})
