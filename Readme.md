# Introduction to R Programming

## Introduction to R
Overview of R and its applications.

Introduction to RStudio environment.

Navigating the RStudio interface.

Writing and running basic R scripts.

### Basics of R Syntax and Data Types & Variables and Basic Operations
Understanding R syntax and conventions. 
Introduction to R's data types: numeric, character, logical, and factor.
Working with different data types.
Creating and assigning variables.
Performing arithmetic operations in R.
Working with logical operations and comparison


### Vectors, Lists and Data Frames

Introduction to vectors in R.
Creating and manipulating vectors.
Vectorized operations.
Common functions for working with vectors (e.g., sum (), mean (), length ()).
Understanding lists and their structure.
Creating and accessing elements in a list.



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

