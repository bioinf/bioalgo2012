print("loading features...")
d <- read.csv("/labnas/students/adievsky/dievsky/bioalgo/1/final_release/train.csv")

d[["class"]] <- factor(d[["class"]])

library(nnet)
print("buildng model...")
model <- nnet( class ~ . , d, size = 4 )

save(model, file = "/labnas/students/adievsky/dievsky/bioalgo/1/final_release/model")
save(model, file = "/labnas/students/adievsky/dievsky/bioalgo/model")
save(model, file = "/labnas/students/ksenia/bioalgo/model")
save(model, file = "/labnas/students/ksenia/bioalgo/r/model")
save(model, file = "./model")
