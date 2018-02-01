model_xgb <- h2o.xgboost(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_xgb", nfolds = nfolds, keep_cross_validation_predictions = TRUE, fold_assignment = "Stratified")

h2o.performance(model_xgb, xval = T)
plot(h2o.performance(model_xgb, xval = T))

model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, distribution = "bernoulli", ntrees = 10, max_depth = 5, min_rows = 2, learn_rate = 0.2, nfolds = nfolds, fold_assignment = "Stratified", keep_cross_validation_predictions = TRUE, seed = 16052017)
h2o.performance(model_gbm, xval = T)
plot(h2o.performance(model_gbm, xval = T)
)
model_ensemble <- h2o.stackedEnsemble(x = features, y = target,
									  training_frame = train.set.h2o,
									  validation_frame = test.set.h2o,
									  model_id = "h2o_ensemble2",
									  base_models = list(model_gbm@model_id, model_drf@model_id, model_xgb@model_id))

summary(model_ensemble)
h2o.performance(model_ensemble, newdata = test.set.h2o)
h2o.performance(model_drf, newdata = test.set.h2o)
h2o.performance(model_gbm, newdata = test.set.h2o)
h2o.performance(model_xgb, newdata = test.set.h2o)

test.set.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.negative.set.rds")
test.set.h2o <- as.h2o(test.set.matched)
