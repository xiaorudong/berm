# berm: Bootstrap-Enhanced Regularization Method

The `berm` package implements the Bootstrap-Enhanced Regularization Method (BERM), a statistical approach aiming to enhance the robustness and accuracy of variable selection and coefficient estimation in immunophenotyping datasets. These datasets are typically characterized by high multicollinearity and dependence, alongside substantially skewed distributions. By integrating bootstrapped confidence intervals with penalized regression techniques, BERM shows robust variable selection and precise coefficient estimation, effectively addressing the challenges posed by these complex data structures.

## Installation

You can install the released version of `berm` from Github with:

```{r}
devtools::install_github("xiaorudong/berm")
```
## Usage
Here is a basic example of how to use the berm package:

```{r}
library(berm)

set.seed(123)
mydata <- matrix(rnorm(100 * 10), ncol = 10)
colnames(mydata) <- paste("X", 1:10, sep = "_")
y <- rnorm(100)

# Fit BERM model
result <- BERM.model(x = mydata, y = y)
head(result$coef)
```