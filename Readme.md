# Introduction to R Programming

## Overview of R
- R is a programming language and environment for statistical computing and graphics.
- Applications include:
  - Data analysis
  - Statistical modeling
  - Machine learning
  - Data visualization
  - Bioinformatics and social sciences

---

## RStudio Environment
### Introduction
- RStudio is an integrated development environment (IDE) for R.
- Provides tools to write, debug, and visualize R code.

### Navigating the Interface
- **Source Pane**: Write and edit scripts.
- **Console**: Execute commands interactively.
- **Environment/History Pane**: View variables and command history.
- **Files/Plots/Packages/Help Pane**: Manage files, visualize plots, install packages, and access documentation.

### Writing and Running Scripts
- Scripts are saved as `.R` files.
- Run code using:
  - `Ctrl + Enter` (Windows/Linux)
  - `Cmd + Enter` (Mac)
- Output appears in the Console.

---

## Basics of R Syntax and Data Types
### Syntax and Conventions
- Case-sensitive language.
- Use `<-` or `=` for assignment.
- Functions are called with `function_name(arguments)`.

### Data Types
- **Numeric**: `x <- 10`
- **Character**: `y <- "Hello"`
- **Logical**: `z <- TRUE`
- **Factor**: `factor(c("Male", "Female", "Male"))`

### Variables and Operations
- Creating variables: `a <- 5`
- Arithmetic operations: `+`, `-`, `*`, `/`, `^`
- Logical operations: `&`, `|`, `!`
- Comparisons: `==`, `!=`, `<`, `>`, `<=`, `>=`

---

## Vectors, Lists, and Data Frames
### Vectors
- Create: `v <- c(1, 2, 3, 4)`
- Vectorized operations: `v * 2`
- Common functions: `sum(v)`, `mean(v)`, `length(v)`

### Lists
- Structure: Can hold different data types.
- Create: `list1 <- list(name="John", age=25, scores=c(90, 85, 88))`
- Access: `list1$name`

### Data Frames
- Fundamental structure for tabular data.
- Create:  
  ```R
  df <- data.frame(
    Name = c("Alice", "Bob"),
    Age = c(25, 30),
    Score = c(90, 85)
  )

