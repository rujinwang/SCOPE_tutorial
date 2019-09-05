context("foreach")
library(SCOPE)

test_that("Parallel computing works", {
  Gini = getGini(Y_sim)
  normObj.sim <- normalize_codex2_ns_noK(Y_qc =Y_sim,
                                         gc_qc = ref_sim$gc,
                                         norm_index = which(Gini<=0.12))

  Yhat.noK.sim=normObj.sim$Yhat
  beta.hat.noK.sim=normObj.sim$beta.hat
  ploidy.sim =  PreEst_ploidy(Y = Y_sim, Yhat = Yhat.noK.sim, ref = ref_sim)
  TimeParallel = system.time({normObj.scope.sim = normalize_scope_foreach(Y_qc = Y_sim, gc_qc = ref_sim$gc,
                                                           K = 1, ploidyInt = ploidy.sim,
                                                           norm_index = which(Gini<=0.12), T = 1:7,
                                                           beta0 = beta.hat.noK.sim, nCores = 2)})

  TimeSeq = system.time({normObj.scope.sim2 = normalize_scope(Y_qc = Y_sim, gc_qc = ref_sim$gc,
                                                    K = 1, ploidyInt = ploidy.sim,
                                                    norm_index = which(Gini<=0.12), T = 1:7,
                                                    beta0 = beta.hat.noK.sim)})

  expect_lte(TimeParallel[1], TimeSeq[1])
  expect_lte(TimeParallel[3], TimeSeq[3])
})
