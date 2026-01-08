
# R Data Analysis for Biologists

**Author:** Diya Sen  
**Audience:** PhD/Masters biology students (quantitative biology, ecology, plant sciences, biostatistics)  
**Dataset:** Built-in `iris` dataset (Fisher, 1936)  
**Format:** Each module contains a *Question*, *R Solution*, and *Interpretation*.

---

## 🚀 Setup: Load Data & Packages

```r
# Load dataset
data(iris)

# Optional: install and load required packages
packages <- c("ggplot2","ggfortify","MASS","cluster","randomForest","mgcv","effectsize","pwr","CCA")
installed <- packages %in% rownames(installed.packages())
if(any(!installed)) install.packages(packages[!installed], dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

# Ensure Species is a factor
iris$Species <- as.factor(iris$Species)
```

> **Tip:** If you cannot install some packages (e.g., `CCA`), skip the corresponding module or use alternatives.

---

## Module 1 — Getting Started: Loading & Inspecting Data (Beginner)

**Question:** What are the basic properties (dimension, variables, species count) of the iris dataset?

```r
dim(iris)         # rows, columns
str(iris)         # structure of variables
table(iris$Species)  # species counts
summary(iris)     # summary statistics
```

**Interpretation:** 150 samples, 4 numeric traits (*Sepal.Length*, *Sepal.Width*, *Petal.Length*, *Petal.Width*), and 3 species (*setosa*, *versicolor*, *virginica*).

---

## Module 2 — Basic Univariate Biology-Oriented Analysis (Beginner)

**Question:** What is the distribution of each morphological trait? Are any traits skewed?

```r
par(mfrow = c(2, 2))
hist(iris$Sepal.Length, main = "Sepal Length", col = "skyblue")
hist(iris$Sepal.Width,  main = "Sepal Width",  col = "lightgreen")
hist(iris$Petal.Length, main = "Petal Length", col = "salmon")
hist(iris$Petal.Width,  main = "Petal Width",  col = "orchid")
par(mfrow = c(1, 1))
```

**Interpretation:** Petal traits often show clearer separation across species; sepal width tends to be more variable.

---

## Module 3 — Normality Testing & Biological Implications (Beginner → Intermediate)

**Question:** Are morphological traits normally distributed? (Useful for selecting appropriate statistical models.)

```r
norm_tests <- lapply(iris[,1:4], shapiro.test)
norm_tests
```

**Interpretation:** If p < 0.05 → trait deviates from normality → consider nonparametric tests or transformations.

---

## Module 4 — Bivariate Relationships: Correlation & Interpretation (Intermediate)

**Question:** Which pairs of traits are most correlated? What does this imply about growth coordination?

```r
cor_matrix <- cor(iris[,1:4])
cor_matrix

pairs(iris[,1:4], col = iris$Species, pch = 19)
```

**Interpretation:** Petal length & width typically show the strongest correlation—suggesting coordinated floral development/allometry.

---

## Module 5 — Visualizing Species Differences (Intermediate)

**Question:** How do traits differ between species visually? Create species-stratified boxplots.

```r
ggplot(iris, aes(Species, Petal.Length, fill = Species)) +
  geom_boxplot(alpha = 0.8) +
  labs(title = "Petal Length by Species", x = "Species", y = "Petal Length (cm)") +
  theme_minimal()
```

**Interpretation:** Petal length separates species strongly, consistent with botanical taxonomy.

---

## Module 6 — Hypothesis Testing (Intermediate)

**Question:** Do species differ significantly in petal length? Test using ANOVA and post-hoc comparisons.

```r
aov_model <- aov(Petal.Length ~ Species, data = iris)
summary(aov_model)

# Post-hoc Tukey test
TukeyHSD(aov_model)
```

**Interpretation:** A significant ANOVA indicates at least one species differs; Tukey’s test identifies pairwise differences.

---

## Module 7 — MANOVA (PhD Entry Level)

**Question:** Do species differ in multivariate trait space (all four measurements together)?

```r
manova_model <- manova(as.matrix(iris[,1:4]) ~ Species)
summary(manova_model, test = "Pillai")
```

**Interpretation:** A significant Pillai’s trace shows species differ jointly across all traits (multivariate phenotype space).

---

## Module 8 — PCA (Dimensionality Reduction) (Graduate Level)

**Question:** How many major axes of morphological variation exist? Which traits dominate each principal component?

```r
pca <- prcomp(iris[,1:4], scale. = TRUE)
summary(pca)       # variance explained
pca$rotation       # loadings

# PCA biplot with loadings
autoplot(pca, data = iris, colour = "Species", loadings = TRUE, loadings.label = TRUE) +
  ggtitle("PCA of Iris Morphology")
```

**Interpretation:** PC1 typically reflects petal dimensions; PC2 often reflects sepal width. Species separation mainly along PC1.

---

## Module 9 — Linear Discriminant Analysis (LDA) (Graduate Level)

**Question:** Can we classify species accurately using their morphological traits? Which traits have strongest discriminant power?

```r
lda_model <- MASS::lda(Species ~ ., data = iris)
lda_model

pred <- predict(lda_model)$class
acc <- mean(pred == iris$Species)
acc

# Confusion matrix
table(Predicted = pred, True = iris$Species)
```

**Interpretation:** Accuracy is usually >95%. LD loadings indicate traits contributing most to separation.

---

## Module 10 — K-means Clustering & Cluster Validation (Advanced)

**Question:** If species labels were unknown, could unsupervised clustering recover the species structure?

```r
set.seed(123)
km <- kmeans(iris[,1:4], centers = 3, nstart = 50)

# Compare clusters to true species
table(Cluster = km$cluster, Species = iris$Species)

# Silhouette width (cluster validation)
d <- dist(iris[,1:4])
sil <- cluster::silhouette(km$cluster, d)
mean_sil <- mean(sil[,3])
mean_sil
```

**Interpretation:** High mean silhouette width indicates coherent clusters aligning with species; misclassification often occurs between *versicolor* and *virginica* due to overlap.

---

## Module 11 — Random Forest Feature Importance (PhD Level)

**Question:** Which traits are most important for predicting species identity under an ensemble learning method?

```r
set.seed(101)
rf <- randomForest::randomForest(Species ~ ., data = iris, importance = TRUE)
rf
randomForest::varImpPlot(rf, main = "Random Forest Feature Importance")
randomForest::importance(rf)
```

**Interpretation:** Petal length and petal width typically dominate importance, confirming their biological relevance for species differentiation.

---

## Module 12 — Nonlinear Modelling with GAM (Very Advanced)

**Question:** Does petal length scale nonlinearly with sepal length? (Developmental allometry.)

```r
gam_model <- mgcv::gam(Petal.Length ~ s(Sepal.Length), data = iris)
summary(gam_model)
plot(gam_model, shade = TRUE, residuals = TRUE)
```

**Interpretation:** A significant smooth term (edf > 1) suggests nonlinear scaling; relevant for modeling growth patterns.

---

## Optional — Effect Size & Power Analysis (PhD Design Consideration)

**Question:** What is the effect size of species differences in petal length, and is the sample size adequate?

```r
aov_model <- aov(Petal.Length ~ Species, data = iris)
eta <- effectsize::eta_squared(aov_model)
eta

# Convert eta^2 to Cohen's f for power analysis: f = sqrt(eta^2 / (1 - eta^2))
eta2 <- eta$Eta2[1]
cohen_f <- sqrt(eta2 / (1 - eta2))
cohen_f

pwr::pwr.anova.test(k = 3, f = cohen_f, sig.level = 0.05, power = 0.8)
```

**Interpretation:** Large effect sizes imply fewer samples needed to detect species differences; useful for planning phenotyping experiments.

---

## Capstone — Integrated Research Questions & Workflow

**Final Questions:**
1. What is the minimum set of traits required to distinguish species reliably?  
2. Are species separations linear or nonlinear?  
3. Do traits form developmental modules (e.g., sepal vs petal)?  
4. How much variance is explained by species identity vs intrinsic variation?

**Suggested R Workflow:**

```r
# 1) Multivariate separation
manova(as.matrix(iris[,1:4]) ~ Species)

# 2) Dimensionality & loadings
pca <- prcomp(iris[,1:4], scale. = TRUE)
summary(pca); pca$rotation

# 3) Supervised classification
lda_model <- MASS::lda(Species ~ ., data = iris)
mean(predict(lda_model)$class == iris$Species)

# 4) Unsupervised validation
set.seed(123); km <- kmeans(iris[,1:4], 3, nstart=50)
table(km$cluster, iris$Species)

# 5) Nonlinear relationships
gam_model <- mgcv::gam(Petal.Length ~ s(Sepal.Length), data = iris)
summary(gam_model)

# 6) Feature importance
rf <- randomForest::randomForest(Species ~ ., data = iris, importance = TRUE)
randomForest::importance(rf)
```

**Reporting Guidance:**
- Summarize MANOVA and PCA results to quantify multivariate separation.  
- Use LDA accuracy and confusion matrix to demonstrate predictive validity.  
- Discuss GAM smooths to justify nonlinear biological scaling.  
- Present Random Forest importance to prioritize morphological traits for future data collection.

---

## Reproducibility

```r
sessionInfo()
```

**Notes:** Set a random seed where appropriate; report versions of packages and R.

---

## References & Attribution
- Fisher, R.A. (1936). *The use of multiple measurements in taxonomic problems.* Annals of Eugenics.  
- R Core Team (2024). *R: A Language and Environment for Statistical Computing.* R Foundation for Statistical Computing.  
- Venables, W.N. & Ripley, B.D. (2002). *Modern Applied Statistics with S.* (for LDA via `MASS`).  
- Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R.* (for `mgcv`).

