context("treedata works")
test_that("treedata can handle matrix/dataframe input", {
  require(testthat)
  data(anolis)
  norn <- anolis$dat[,-1]
  rownames(norn) <- anolis$dat[,1]
  td0 <- make.treedata(anolis$phy, norn)
  td1 <- make.treedata(anolis$phy, as.matrix(anolis$dat))
  td2 <- make.treedata(anolis$phy, anolis$dat, name_column=1)
  jacknife <- sample(1:nrow(anolis$dat), 50, replace=FALSE)
  td3 <- make.treedata(anolis$phy, anolis$dat[jacknife,], name_column=1)
  td4 <- make.treedata(anolis$phy, as.matrix(anolis$dat[jacknife,]), name_column=1)
  
  ##Tests
  expect_equal(td0$dat$SVL, td1$dat$SVL,  tolerance=0.000001)
  expect_equal(td1$dat$SVL, td2$dat$SVL,  tolerance=0.000001)
  expect_equal(td1$dat$island, td2$dat$island)
  expect_equal(td3$dat$SVL, td4$dat$SVL,  tolerance=0.000001)
  expect_true(all(td3$phy$tip.label %in%  anolis$dat[jacknife,'X']))
  o3.1 <- order(td3$phy$tip.label)
  o3.2 <- order(anolis$dat[jacknife,'X'])
  expect_true(all(td3$dat$SVL[o3.1]==anolis$dat[jacknife,][o3.2, "SVL"]))
  
  expect_true(all(td4$phy$tip.label %in%  as.matrix(anolis$dat[jacknife,'X'])))
  o3.1 <- order(td3$phy$tip.label)
  o3.2 <- order(anolis$dat[jacknife,'X'])
  expect_true(all(td3$dat$SVL[o3.1]==as.matrix(anolis$dat[jacknife,][o3.2, "SVL"])))
  
})
