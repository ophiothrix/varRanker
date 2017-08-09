library(ProjectTemplate)
create.project("../varRank", merge.strategy = "allow.non.conflict")

## Clean auxilary R files
system("rm ./diagnostics/1.R ./lib/helpers.R ./munge/01-A.R ./profiling/1.R ./src/eda.R ./tests/1.R")
