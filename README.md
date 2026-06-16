# berm: Bootstrap-Enhanced Regularization Method

The `berm` package implements the Bootstrap-Enhanced Regularization Method (BERM), a statistical approach aiming to enhance the robustness and accuracy of variable selection and coefficient estimation in immunophenotyping datasets. These datasets are typically characterized by high multicollinearity and dependence, alongside substantially skewed distributions. By integrating bootstrapped confidence intervals with penalized regression techniques, BERM shows robust variable selection and precise coefficient estimation, effectively addressing the challenges posed by these complex data structures.

## Installation

You can install the released version of `berm` from Github with:

```{r}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("xiaorudong/berm")
```
## Example
Here is a basic example of how to use the `berm` package:

```{r}
library(berm)

set.seed(123)
mydata <- matrix(rnorm(100 * 10), ncol = 10)
colnames(mydata) <- paste("X", 1:10, sep = "_")
y <- rnorm(100)

# Fit a BERM model
res_berm <- berm(x = mydata, y = y)
summary(res_berm)

# Predict in a new dataset
newdata <- matrix(rnorm(100 * 10), ncol = 10)
colnames(newdata) <- paste("X", 1:10, sep = "_")
predict_res <- predict(object = res_berm, newdata = newdata)
head(predict_res)
```



## Tuning Parameters

### Default setting

By default, `unrestricted = FALSE`. Under this setting, `berm` is fixed at `alpha = 0.5` during the bootstrap variable selection step.

This default help to reduce computation time while maintaining stable performance. 


### Optional alpha tuning

Users may set `unrestricted = TRUE` to allow both `alpha` and `lambda` to be selected through cross-validation during each bootstrap iteration.

```r
res_berm <- berm(x = mydata, y = y, unrestricted = TRUE)
```

This option is useful when strong prior knowledge suggests an extreme sparsity pattern or coefficient structure. Examples include situations where only a very small number of predictors are expected to be associated with the outcome, or where a large proportion of predictors are expected to have non-zero effects.

Because tuning both `alpha` and `lambda` within every bootstrap iteration can substantially increase computation time, we recommend using the default setting (`unrestricted = FALSE`) unless additional flexibility is needed.


When `unrestricted = TRUE`, BERM performs cross-validation over:

```r
alpha_grid = seq(1, 0, length.out = 20)
lambda_grid = 10^seq(2, -3, by = -0.1)
```
Users may also provide custom `alpha_grid` and `lambda_grid` values based on their data characteristics or computational requirements.



## Reference

For methodological details and simulation studies, please see:

Dong X, Goyal A, Liang M, Brusko MA, Brusko TM, Bacher R. *Penalized Linear Models for Highly Correlated High-Dimensional Immunophenotyping Data*. arXiv. 2025.

https://arxiv.org/abs/2504.07771
