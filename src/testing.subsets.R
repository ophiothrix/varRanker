# rm(list = ls())
# gc()
source("./src/model.check.test.unit.R")

models.summary <- read.csv("./reports/annotated.variant.summary.csv", row.names = 1, header = T)
models.summary$AUC <- NA
models.summary

training.set <- readRDS("./cache/all.variants.full.annotation.rds")
models.summary["all.variants.full.annotation", "AUC"] <- model.check.from.df(training.set)

training.set <- readRDS("./cache/imperfect.variants.full.annotation.rds")
models.summary["imperfect.variants.full.annotation", "AUC"] <- model.check.from.df(training.set)

training.set <- readRDS("./cache/perfect.variants.full.annotation.rds")
models.summary["perfect.variants.full.annotation", "AUC"] <- model.check.from.df(training.set)

training.set <- readRDS("./cache/all.variants.partial.annotation.rds")
models.summary["all.variants.partial.annotation", "AUC"] <- model.check.from.df(training.set)

training.set <- readRDS("./cache/imperfect.variants.partial.annotation.rds")
models.summary["imperfect.variants.partial.annotation", "AUC"] <- model.check.from.df(training.set)

training.set <- readRDS("./cache/perfect.variants.partial.annotation.rds")
models.summary["perfect.variants.partial.annotation", "AUC"] <- model.check.from.df(training.set)

write.csv(models.summary, "./reports/annotated.variant.summary.csv")
