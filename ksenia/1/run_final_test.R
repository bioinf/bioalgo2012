print("loading model")
load("./model" )

print("loading data to predict...")
data_test <- read.csv("/labnas/students/adievsky/dievsky/bioalgo/1/final_release/test.csv")

print("model:")
print(model)
p = predict( model, data_test, type = "class" )
write(p, file = "./predicted")

