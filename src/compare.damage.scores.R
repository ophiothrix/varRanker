positive.set <- readRDS("./cache/E126.annotated.fSkin_fibro.variants.rds")
negative.set <- readRDS("./cache/E126.annotated.fSkin_fibro.negative.set.rds")

pdf("./graphs/compare.damage.scores.pdf", 10, 6)
par(mfrow = c(1,2))
hist(positive.set$hocomoco.damage.score, freq = F, breaks = 100, col = "#00FF0088", ylim = c(0, 1), main = "Positive vs Negative HOCOMOCO")
hist(negative.set$hocomoco.damage.score, freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Positive set", "Negative set"), fill = c("#00FF0088", "#FF000088"), bty = "n")

hist(positive.set$hocomoco.damage.score, freq = F, breaks = 100, col = "#00FF0088", ylim = c(0, 1), main = "Positive vs Negative DNase only HOCOMOCO")
hist(negative.set$hocomoco.damage.score[negative.set$source == "dnase"], freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Positive set", "Negative set DNase only"), fill = c("#00FF0088", "#FF000088"), bty = "n")


hist(positive.set$jaspar.damage.score, freq = F, breaks = 100, col = "#00FF0088", ylim = c(0, 1), main = "Positive vs Negative jaspar")
hist(negative.set$jaspar.damage.score, freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Positive set", "Negative set"), fill = c("#00FF0088", "#FF000088"), bty = "n")

hist(positive.set$jaspar.damage.score, freq = F, breaks = 100, col = "#00FF0088", ylim = c(0, 1), main = "Positive vs Negative DNase only jaspar")
hist(negative.set$jaspar.damage.score[negative.set$source == "dnase"], freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Positive set", "Negative set DNase only"), fill = c("#00FF0088", "#FF000088"), bty = "n")
dev.off()
