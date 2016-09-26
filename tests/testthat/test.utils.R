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

test_that("phredConversions", {
  # Simple combine with a name vector
  res  = toPhred(10 ^ (-(1:4)), 1e-3)
  expect_equal(res, c(10, 20, 30, 30))


  res = fromPhred(seq(10, 70, 10), 1e-5)
  expect_equal(res, c(1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-05, 1e-05))
})

test_that("holeNumberConversions", {
  hole = 50 * (2^16) + 95
  x = getHoleX(hole)
  expect_equal(x, 50)
  y = getHoleY(hole)
  expect_equal(y, 95)
})

test_that("loadReferencePath", {
  expected = normalizePath("./ReferenceSet/lambdaNEB/sequence/lambdaNEB.fasta")
  name = getReferencePath("./ReferenceSet/lambdaNEB/referenceset.xml")
  expect_equal(name, expected)
  expect_true(file.exists(name))
})

