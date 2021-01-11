## -----------------------------------------------------------------------------
# +----------------------------------------------------------------------------+
# |                               the R code for the                           |
# |                                   Big R-book                               |
# |                           by Philippe J.S. De Brouwer                      |
# +----------------------------------------------------------------------------+
# | Copyright: (C) Philippe De Brouwer, 2020                                   |
# |            Owners of the book can use this code for private and commercial |
# |            use, except re-publishing.                                      |
# +----------------------------------------------------------------------------+
# | This document is a true copy of all the R code that appears in the book.   |
# | It is provided for you so that when you work with the book, you can copy   | 
# | and paste any code segment that you want to test without typing everything |
# | again.                                                                     |
# | No warranties or promises are given.                                       |
# | Have an intersting time and learn fast!                                    |
# | Philippe J.S. De Brouwer                                                   |
# +----------------------------------------------------------------------------+
# | For your convenience, and to help searching, each new code block is pre-   |
# | ceded by ## ---- and eventually followed by the caption of the plot that   |
# | is genderated by the code block.                                           |
# +----------------------------------------------------------------------------+
#

## ----eval=FALSE,echo=FALSE----------------------------------------------------
# Code to install all necessary packages:
p <- c(
"DiagrammeR",
"tidyverse",
"ggvis",
"pryr",
"rex",
"xlsx",
"RMySQL",
"MASS",
"carData",
"titanic",
"rpart",
"rpart.plot",
"ROCR",
"pROC",
"ggplot2",
"viridisLite",   # color schemes
"vioplot",
"grid",
"gridExtra",
"forecast",
"randomForest",
"RCurl",
"XML",
"stringr",
"plyr",
"reshape",
"sandwich",
"msm",
"quantmod",
"tm",
"SnowballC",
"RColorBrewer",
"wordcloud",
"quantmod",
"neuralnet",
"sodium",
"binr",
"diagram",
"latex2exp",
"ggfortify",
"cluster",
"flexdashboard",
"shiny",
"shinydashboard",
"mice",
"VIM",
"missForest",
"Hmisc",
"mi",
"ggrepel",
"plot3D",
"plotly",
"class",
"RColorBrewer",
"InformationValue",
"e1071",
"psych",
"nFactors",
"FactoMineR",
"knitr",
"devtools",
"roxygen2",
"zoo",
"xts",
"compiler",
"Rcpp",
"profr",
"proftools",
"SparkR",
"sparklyr",
"DBI",
"gpuR",
"parallel",
"Rcrawler"
)

system("cat /dev/null > bibTexPackagesNEW.bib")
for (i in 1:length(p)){
 print(paste("===",i, ": ",p[i],"==="))
 if (!require(p[i], character.only = TRUE)) {install.packages(p[i])}
 ref <- toBibtex(citation(p[i]))
 cat(ref, file="bibTexPackagesNEW.bib", append=TRUE, sep = "\n")
 }



## -----------------------------------------------------------------------------
# This is code
1+pi
Sys.getenv(c("EDITOR","USER","SHELL", "LC_NUMERIC"))


## ----fig.cap="An example showing the histogram of data generated from the normal distribution.\\label{fig:iq}"----
# generate 1000 random numbers between 0 and 100
x <- rnorm(1000, mean = 100, sd = 2)

# to illustrate previous, we show the histogram.
hist(x, col = "khaki3")

# This comment of the code is after the histogram.
# In rare cases the plot will be on the this page
# alone and this comment is the previous page.


## -----------------------------------------------------------------------------
# First, generate some data:
x <- c(1,2,3)

# Then calculate the mean:
mean(x)


## -----------------------------------------------------------------------------
mean(1:100)


## -----------------------------------------------------------------------------
# Code and especially the comments in it are part of 
# the normal flow of the text!


## ----eval=FALSE---------------------------------------------------------------
##  #addition
##  2 + 3
##  #product
##  2 * 3
##  #power
##  2**3
##  2^3
##  #logic
##  2 < 3
##  x <- c(1,3,4,3)
##  x.mean <- mean(x)
##  x.mean
##  y <- c(2,3,5,1)
##  x+y


## ----eval=FALSE---------------------------------------------------------------
## x <- scan()


## ----eval=FALSE---------------------------------------------------------------
## edit(x)


## ----eval=FALSE---------------------------------------------------------------
## my_function <- function(a,b)
## {
##   a +  b
## }


## -----------------------------------------------------------------------------
x.1 <- 5
x.1 + 3 -> .x
print(.x)


## -----------------------------------------------------------------------------
x.3 = 3.14
x.3


## -----------------------------------------------------------------------------
v1 <- c(1,2,3)
mean(v1, na.rm = TRUE)


## ----eval=FALSE---------------------------------------------------------------
## # List all variables
## ls()                   # hidden variable starts with dot
## ls(all.names = TRUE)   # shows all
## 
## # Remove a variable
## rm(x.1)                # removes the variable x.1
## ls()                   # x.1 is not there any more
## rm(list = ls())        # removes all variables
## ls()


## ----echo---------------------------------------------------------------------
# Booleans can be TRUE or FALSE:
x <- TRUE
class(x)

# Integers use the letter L (for Long integer):
x <- 5L
class(x)

# Decimal numbers, are referred to as 'numeric':
x <- 5.135
class(x)

# Complex numbers use the letter i (without multiplication sign):
x <- 2.2 + 3.2i
class(x)

# Strings are called 'character':
x <- "test"
class(x)


## -----------------------------------------------------------------------------
# Avoid this:
x <- 3L      # x defined as integer
x
x <- "test"  # R changes data type
x


## -----------------------------------------------------------------------------
# the function as.Data coerces its argument to a date
d <- as.Date(c("1852-05-12", "1914-11-5", "2015-05-01"))

# dates will work as expected
d_recent <- subset(d, d > as.Date("2005-01-01"))
print(d_recent)


## -----------------------------------------------------------------------------
x <- c(2, 2.5, 4, 6)
y <- c("apple", "pear")
class(x)
class(y)


## -----------------------------------------------------------------------------
v <- c(1:5)

# access elements via indexing
v[2]
v[c(1,5)]

# access via TRUE/FALSE:
v[c(TRUE,TRUE,FALSE,FALSE,TRUE)]

# access elements via names
v <- c("pear" = "green", "banana" = "yellow", "coconut" = "brown")
v
v["banana"]

#leave out certain elements:
v[c(-2,-3)]


## -----------------------------------------------------------------------------
v1 <- c(1,2,3)
v2 <- c(4,5,6)
# Standard arithmetic
v1 + v2
v1 - v2
v1 * v2


## -----------------------------------------------------------------------------
# Define a short and long vector:
v1 <- c(1, 2, 3, 4, 5)
v2 <- c(1, 2)
v1 + v2
# Note that R does not complain, but rather 'recycles' v2 
# to match the length of v1.


## -----------------------------------------------------------------------------
# Example 1:
v1 <- c(1, -4, 2, 0, pi)
sort(v1)

# Example 2: To make sorting meaningful, all variables need to be of 
# the most complex type:
v1 <- c(1:3, 2 + 2i)
sort(v1)

# Sorting is per increasing numerical or alphabetical order:
v3 <- c("January", "February", "March", "April")
sort(v3)

# Sort order can be reversed:
sort(v3, decreasing = TRUE)


## -----------------------------------------------------------------------------
# Create a matrix.
M = matrix( c(1:6), nrow = 2, ncol = 3, byrow = TRUE)
print(M)
M = matrix( c(1:6), nrow = 2, ncol = 3, byrow = FALSE)
print(M)


## -----------------------------------------------------------------------------
# Unit vector:
matrix (1, 2, 1)

# Zero matrix or vector:
matrix (0, 2, 2)

# Recycling also works for shorter vectors:
matrix (1:2, 4, 4)

# Fortunately, R expects that the vector fits exactly n times in the matrix:
matrix (1:3, 4, 4)
# So, the previous was bound to fail.


## -----------------------------------------------------------------------------
row_names = c("row1", "row2", "row3", "row4")
col_names = c("col1", "col2", "col3")
M <- matrix(c(10:21), nrow = 4, byrow = TRUE, 
            dimnames = list(row_names, col_names))
print(M)



## -----------------------------------------------------------------------------
colnames(M) <- c('C1', 'C2', 'C3')
rownames(M) <- c('R1', 'R2', 'R3', 'R4')
M


## -----------------------------------------------------------------------------
M <- matrix(c(10:21), nrow = 4, byrow = TRUE)
M

# Access one element:
M[1,2]

# The second row:
M[2,]

# The second column:
M[,2]

# Row 1 and 3 only:
M[c(1, 3),]

# Row 2 to 3 with column 3 to 1
M[2:3, 3:1]


## -----------------------------------------------------------------------------
M1 <- matrix(c(10:21), nrow = 4, byrow = TRUE)
M2 <- matrix(c(0:11),  nrow = 4, byrow = TRUE)
M1 + M2
M1 * M2
M1 / M2


## -----------------------------------------------------------------------------
# Example of the dot-product:
a <- c(1:3)
a %*% a
a %*% t(a)
t(a) %*% a

# Define A:
A <- matrix(0:8, nrow = 3, byrow = TRUE)

# Test products:
A %*% a
A %*% t(a) # this is bound to fail!
A %*% A


## -----------------------------------------------------------------------------
A %/% A 


## -----------------------------------------------------------------------------
# Note the difference between the normal product:
A * A 

# and the matrix product %*%:
A %*% A

## -----------------------------------------------------------------------------
# However, there is - of course - only one sum:
A + A

## -----------------------------------------------------------------------------
# Note that the quotients yield almost the same:
A %/% A

A / A


## -----------------------------------------------------------------------------
# This is the matrix A:
A

# The exponential of A:
exp(A)


## -----------------------------------------------------------------------------
# The natural logarithm
log(A)
sin(A)


## -----------------------------------------------------------------------------
# Collapse to a vectore:
colSums(A)
rowSums(A)

# Some functions aggregate the whole matrix to one scalar:
mean(A)
min(A)


## -----------------------------------------------------------------------------
M <- matrix(c(1,1,4,1,2,3,3,2,1), 3, 3)
M

# The diagonal of M:
diag(M)

# Inverse:
solve(M)

# Determinant:
det(M)

# The QR composition:
QR_M <- qr(M)
QR_M$rank

# Number of rows and columns:
nrow(M)
ncol(M)

# Sums of rows and columns:
colSums(M)
rowSums(M)

# Means of rows, columns, and matrix:
colMeans(M)
rowMeans(M)
mean(M)

# Horizontal and vertical concatenation:
rbind(M, M)
cbind(M, M)



## -----------------------------------------------------------------------------
# Create an array:
a <- array(c('A','B'),dim = c(3,3,2))
print(a)

# Access one element:
a[2,2,2]

# Access one layer:
a[,,2]


## ----arrayChunk,eval=TRUE, results='hide'-------------------------------------
# Create two vectors:
v1 <- c(1,1)
v2 <- c(10:13)
col.names <- c("col1","col2", "col3")
row.names <- c("R1","R2")
matrix.names <- c("Matrix1","Matrix2")

# Take these vectors as input to the array.
a <- array(c(v1,v2),dim = c(2,3,2),
     dimnames = list(row.names,col.names, 
     matrix.names))
print(a)

# This allows to address the first row in Matrix 1 as follows:
a['R1',,'Matrix1']


## ----arrayChunk2--------------------------------------------------------------
M1 <- a[,,1]
M2 <- a[,,2]
M2


## -----------------------------------------------------------------------------
x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
apply(x, 2, mean, trim = .2)
col.sums <- apply(x, 2, sum)
row.sums <- apply(x, 1, sum)
rbind(cbind(x, Rtot = row.sums), 
      Ctot = c(col.sums, sum(col.sums)))


## -----------------------------------------------------------------------------
# Re-create the array a (shorter code):
col.names <- c("col1","col2", "col3")
row.names <- c("R1","R2")
matrix.names <- c("Matrix1","Matrix2")
a <- array(c(1,1,10:13),dim = c(2,3,2),
     dimnames = list(row.names,col.names, 
     matrix.names))

# Demonstrate apply:
apply(a, 1, sum)
apply(a, 2, sum)
apply(a, 3, sum)


## -----------------------------------------------------------------------------
  # List is created using list() function.
  myList <- list("Approximation", pi, 3.14, c)
  print(myList)


## -----------------------------------------------------------------------------
# Create the list:
L <- list("Approximation", pi, 3.14, c)

# Assign names to elements:
names(L) <- c("description", "exact", "approx","function")
print(L)

# Addressing elements of the named list:
print(paste("The difference is", L$exact - L$approx))
print(L[3])
print(L$approx)

# However, "function" was a reserved word, so we need to use
# back-ticks in order to address the element:
a <- L$`function`(2,3,pi,5)  # to access the function c(...)
print(a)


## -----------------------------------------------------------------------------
V1 <- c(1,2,3)
L2 <- list(V1, c(2:7))
L3 <- list(L2,V1)
print(L3)
print(L3[[1]][[2]][3])


## -----------------------------------------------------------------------------
# The first object of L2 as a list:
L2[1]
class(L[2])

# The first element of L2 is a numeric vector:
L2[[1]]
class(L2[[2]])
# range
L2[1:2]

# Unexpected result (ranges are not to be used with x[[.]]):
L2[[1:2]] <- 'a'
L2

# Is this what you would expect?
L2[1:2] <- 'a'
L2


## -----------------------------------------------------------------------------
L <- list(1,2)
L[4] <- 4  # position 3 is NULL
L


## -----------------------------------------------------------------------------
L$pi_value <- pi
L


## -----------------------------------------------------------------------------
L[1] <- NULL
L


## -----------------------------------------------------------------------------
L <- L[-2]
L


## -----------------------------------------------------------------------------
L  <- list(c(1:5), c(6:10))
v1 <- unlist(L[1])
v2 <- unlist(L[2])
v2-v1


## -----------------------------------------------------------------------------
# A list of vectors of integers:
L <- list(1L,c(-10L:-8L))
unlist(L)

# Note the named real-valued extra element:
L <- list(c(1:2),c(-10:-8), "pi_value" = pi)
unlist(L)


## -----------------------------------------------------------------------------
# Create a vector containing all your observations:
feedback <- c('Good','Good','Bad','Average','Bad','Good')

# Create a factor object:
factor_feedback <- factor(feedback)

# Print the factor object:
print(factor_feedback)


## ----factors1,fig.cap="The plot-function will result in a bar-chart for a factor-object.\\label{fig:plotfactor}"----
# Plot the histogram -- note the default order is alphabetic
plot(factor_feedback)


## -----------------------------------------------------------------------------
# The nlevels function returns the number of levels:
print(nlevels(factor_feedback))


## ----factors2,fig.cap="The factor objects appear now in a logical order.\\label{fig:orderedfactors}"----
feedback <- c('Good','Good','Bad','Average','Bad','Good')
factor_feedback <- factor(feedback,
                          levels=c("Bad","Average","Good"))
plot(factor_feedback)


## -----------------------------------------------------------------------------
gl(3,2,,c("bad","average","good"),TRUE)


## ----frame1,fig.cap="The standard plot for a data frame in R shows each column printed in function of each other. This is useful to see correlations or how generally the data is structured."----
# Create the data frame.
data_test <- data.frame(
   Name   = c("Piotr", "Pawel","Paula","Lisa","Laura"), 
   Gender = c("Male", "Male","Female", "Female","Female"), 
   Score  = c(78,88,92,89,84), 
   Age    = c(42,38,26,30,35)
   )
print(data_test)

# The standard plot function on a data-frame (Figure 4.3)
# with the pairs() function:
plot(data_test)


## -----------------------------------------------------------------------------
# Get the structure of the data frame:
str(data_test)
# Note that the names became factors (see warning below)

# Get the summary of the data frame:
summary(data_test)

# Get the first rows:
head(data_test)

# Get the last rows:
tail(data_test)

# Extract the column 2 and 4 and keep all rows
data_test.1 <- data_test[,c(2,4)]
print(data_test.1)

# Extract columns by name and keep only selected rows
data_test[c(2:4),c(2,4)]


## -----------------------------------------------------------------------------
d <- data.frame(
   Name   = c("Piotr", "Pawel","Paula","Lisa","Laura"), 
   Gender = c("Male", "Male","Female", "Female","Female"), 
   Score  = c(78,88,92,89,84), 
   Age    = c(42,38,26,30,35),
   stringsAsFactors = FALSE
   )
d$Gender <- factor(d$Gender)  # manually factorize gender
str(d)


## ----eval=FALSE---------------------------------------------------------------
## de(x)                  # fails if x is not defined
## de(x <- c(NA))         # works
## x <- de(x <- c(NA))    # will also save the changes
## data.entry(x)          # de is short for data.entry
## x <- edit(x)           # use the standard editor (vi in *nix)


## ----warning=FALSE------------------------------------------------------------
# The following lines do the same.
data_test$Score[1] <- 80
data_test[3,1]     <- 80


## -----------------------------------------------------------------------------
# Expand the data frame, simply define the additional column:
data_test$End_date <-  as.Date(c("2014-03-01", "2017-02-13",
                   "2014-10-10", "2015-05-10","2010-08-25"))
print(data_test)

# Or use the function cbind() to combine data frames along columns:
Start_date <- as.Date(c("2012-03-01", "2013-02-13",
                   "2012-10-10", "2011-05-10","2001-08-25"))

# Use this vector directly:
df <- cbind(data_test, Start_date)
print(df)

# or first convert to a data frame:
df <- data.frame("Start_date" = t(Start_date))
df <- cbind(data_test, Start_date)
print(df)


## -----------------------------------------------------------------------------
# To add a row, we need the rbind() function:
data_test.to.add <- data.frame(
   Name = c("Ricardo", "Anna"), 
   Gender = c("Male", "Female"), 
   Score = c(66,80), 
   Age = c(70,36),
   End_date = as.Date(c("2016-05-05","2016-07-07"))
   )
data_test.new <- rbind(data_test,data_test.to.add)
print(data_test.new)


## -----------------------------------------------------------------------------
data_test.1 <- data.frame(
   Name = c("Piotr", "Pawel","Paula","Lisa","Laura"), 
   Gender = c("Male", "Male","Female", "Female","Female"), 
   Score = c(78,88,92,89,84), 
   Age = c(42,38,26,30,35)
   )
data_test.2 <- data.frame(
   Name = c("Piotr", "Pawel","notPaula","notLisa","Laura"), 
   Gender = c("Male", "Male","Female", "Female","Female"), 
   Score = c(78,88,92,89,84), 
   Age = c(42,38,26,30,135)
   )
data_test.merged <- merge(x=data_test.1,y=data_test.2,
                          by.x=c("Name","Age"),by.y=c("Name","Age"))

# Only records that match in name and age are in the merged table:
print(data_test.merged)


## -----------------------------------------------------------------------------
data_test$n


## -----------------------------------------------------------------------------
# Get the rownames.
colnames(data_test)
rownames(data_test)
colnames(data_test)[2]
rownames(data_test)[3]

# assign new names
colnames(data_test)[1] <- "first_name"
rownames(data_test) <- LETTERS[1:nrow(data_test)]
print(data_test)


## -----------------------------------------------------------------------------
a <- "Hello"
b <- "world"
paste(a, b, sep = ", ")
c <- "A 'valid' string"


## -----------------------------------------------------------------------------
paste0(12, '%')


## ----eval=FALSE---------------------------------------------------------------
## format(x, trim = FALSE, digits = NULL, nsmall = 0L,
##        justify = c("left", "right", "centre", "none"),
##        width = NULL, na.encode = TRUE, scientific = NA,
##        big.mark   = "",   big.interval = 3L,
##        small.mark = "", small.interval = 5L,
##        decimal.mark = getOption("OutDec"),
##        zero.print = NULL, drop0trailing = FALSE, ...)


## -----------------------------------------------------------------------------
a<-format(100000000,big.mark=" ",
                 nsmall=3,
                 width=20,
                 scientific=FALSE,
                 justify="r")
print(a)


## -----------------------------------------------------------------------------
v1 <- c(2,4,6,8)
v2 <- c(1,2,3,5)
v1 + v2     # addition
v1 - v2     # subtraction
v1 * v2     # multiplication
v1 / v2     # division
v1 %% v2    # remainder of division
v1 %/% v2   # round(v1/v2 -0.5)
v1 ^ v2     # v1 to the power of v2


## -----------------------------------------------------------------------------
v1 %*% v2


## -----------------------------------------------------------------------------
v1 <- c(8,6,3,2)
v2 <- c(1,2,3,5)
v1 > v2     # bigger than
v1 < v2     # smaller than
v1 <= v2    # smaller or equal
v1 >= v2    # bigger or equal
v1 == v2    # equal
v1 != v2    # not equal


## -----------------------------------------------------------------------------
v1 <- c(TRUE, TRUE, FALSE, FALSE)
v2 <- c(TRUE, FALSE, FALSE, TRUE)

v1 & v2     # and
v1 | v2     # or
!v1         # not
v1 && v2    # and applied to the first element
v1 || v2    # or applied to the first element


v1 <- c(TRUE,FALSE,TRUE,FALSE,8,6+3i,-2,0,   NA)
class(v1)  # v1 is a vector or complex numbers
v2 <- c(TRUE)
as.logical(v1)  # coerce to logical (only 0 is FALSE)
v1 & v2
v1 | v2


## -----------------------------------------------------------------------------
FALSE | NA
TRUE  | NA
FALSE & NA
TRUE  & NA
FALSE | NA | TRUE | TRUE
TRUE  & NA & FALSE


## ----eval=FALSE---------------------------------------------------------------
## # left assignment
## x <- 3
## x = 3
## x<<- 3
## 
## # right assignment
## 3 -> x
## 3 ->> x
## 
## #chained assignment
## x <- y <- 4


## -----------------------------------------------------------------------------
mean(v1, na.rm = TRUE)  # works (v1 is defined in previous section)
mean(v1, na.rm <- TRUE) # fails


## -----------------------------------------------------------------------------
# f
# Assigns in the current and superior environment 10 to x,
# then prints it, then makes it 0 only in the function environment
# and prints it again.
# arguments:
#   x   -- numeric
f <- function(x) {x <<- 10; print(x); x <- 0; print(x)}

x <- 3
x

# Run the function f():
f(x)

# Only the value assigned with <<- is available now:
x


## -----------------------------------------------------------------------------
# +-+
# This function is a new operator
# arguments:
#   x -- numeric
#   y -- numeric
# returns:
#   x - y
`+-+` <- function(x, y) x - y
5 +-+ 5
5 +-+ 1

# Remove the new operator:
rm(`+-+`)


## -----------------------------------------------------------------------------
# create a list
x <- c(10:20)
x

# %in% can find an element in a vector
2  %in% x   # FALSE since 2 is not an element of x
11 %in% x   # TRUE since 11 is in x
x[x %in% c(12,13)] # selects elements from x
x[2:4]   # selects the elements with index
            # between 2 and 4


## -----------------------------------------------------------------------------
# %*% the matrix multiplication (or crossproduct)
M = matrix(c(1,2,3,7,8,9,4,5,6), nrow = 3,ncol = 3,
           byrow = TRUE)
M %*% t(M)
M %*% M 
exp(M)


## ----eval=FALSE---------------------------------------------------------------
## if (logical statement) {is equivalent with
##   executed if logical statement is true
## } else {
##   executed if the logical statement if false
## }


## -----------------------------------------------------------------------------
set.seed(1890)
x <- rnorm(1)
if (x < 0) {
  print('x is negative') 
} else if (x > 0) {
  print('x is positive')
} else {
  print('x is zero')
}


## -----------------------------------------------------------------------------
x <- TRUE
y <- pi
y <- if (x) 1 else 2
y  # y is now 1


## -----------------------------------------------------------------------------
z <- 0
y <- if (x) {1; z <- 6} else 2
y  # y is now 6
z  # z is also 6


## -----------------------------------------------------------------------------
x <- 1:6
ifelse(x %% 2 == 0, 'even', 'odd')


## -----------------------------------------------------------------------------
x <- 1:6
y <- LETTERS[1:3]
ifelse(x %% 2 == 0, 'even', y)
# Note that y gets recycled! Use vectors of the same length.


## -----------------------------------------------------------------------------
x <- 'b'
x_info <- switch(x,
    'a' = "One",
    'b' = "Two",
    'c' = "Three",
    stop("Error: invalid `x` value")
  )
# x_info should be 'two' now:
x_info


## -----------------------------------------------------------------------------
x <- 'b'
x_info <- if (x == 'a' ) {
    "One"
  } else if (x == 'b') {
    "Two"
  } else if (x == 'c') {
    "Three" 
  } else {
    stop("Error: invalid `x` value")
  }
# x_info should be 'two' now:
x_info


## ----eval=FALSE---------------------------------------------------------------
## for (value in vector) {
##    statements
## }


## -----------------------------------------------------------------------------
x <- LETTERS[1:5]
for ( j in x) {
   print(j)
}


## -----------------------------------------------------------------------------
x <- c(0, 1, -1, 102, 102)
for ( j in x) {
   print(j)
}


## ----eval=FALSE---------------------------------------------------------------
## repeat {
##    commands
##    if(condition) {break}
## }


## -----------------------------------------------------------------------------
x <- c(1,2)
c <- 2
repeat {
   print(x+c)
   c <- c+1
   if(c > 4) {break}
}


## ----eval=FALSE---------------------------------------------------------------
## while (test_expression) {
##    statement
## }


## -----------------------------------------------------------------------------
x <- c(1,2); c <- 2
while (c < 4) {
   print(x+c)
   c <- c + 1
}


## -----------------------------------------------------------------------------
v <- c(1:5)
for ( j in v) {
   if (j == 3) {
      print("--break--")
      break
   }
   print(j)
}


## -----------------------------------------------------------------------------
v <- c(1:5)
for ( j in v) {
   if (j == 3) {
      print("--skip--")
      next
   }
   print(j)
}


## ----timingLoops--------------------------------------------------------------
n <- 10^7
v1 <- 1:n

# -- using vector arithmetic
t0 <- Sys.time()
v2 <- v1 * 2
t1 <- Sys.time() - t0
print(t1)

# -- using a for loop
rm(v2)
v2 <- c()
t0 <- Sys.time()
for (k in v1) v2[k] <- v1[k] * 2
t2 <- Sys.time() - t0
print(t2)

# t1 and t2 are difftime objects and only differences 
# are defined.
# To get the quotient, we need to coerce them to numeric.
T_rel <- as.numeric(t2, units = "secs") / 
         as.numeric(t1, units = "secs")
T_rel


## ----eval=FALSE---------------------------------------------------------------
## help(c)   # shows help help with the function c
## ?c        # same result
## 
## apropos("cov") #fuzzy search for functions


## ----eval=FALSE---------------------------------------------------------------
## function_name <- function(arg_1, arg_2, ...) {
##    function_body
##    return_value
## }


## -----------------------------------------------------------------------------
# c_surface
# Calculates the surface of a circle
# Arguments:
#   radius -- numeric, the radius of the circle
# Returns
#   the surface of the cicle
c_surface <- function(radius) {
  x <- radius ^ 2 * pi
  return (x)
  }
c_surface(2) + 2


## -----------------------------------------------------------------------------
# c_surface
# Calculates the surface of a circle
# Arguments:
#   radius -- numeric, the radius of the circle
# Returns
#   the surface of the cicle
c_surface <- function(radius) {
  radius ^ 2 * pi
  }

# Test the function:
c_surface(2) + 2


## ----fixEdit,eval=FALSE-------------------------------------------------------
## # Edit the function with vi:
##  fix(c_surface)
## 
## # Or us edit:
##  c_surface <- edit()


## -----------------------------------------------------------------------------
c_surface <- function(radius = 2) {
  radius ^ 2 * pi
  }
c_surface(1)
c_surface()


## ----eval=FALSE---------------------------------------------------------------
## # Download the package (only once):
## install.packages('DiagrammeR')
## 
## # Load it before we can use it (once per session):
## library(DiagrammeR)


## ----eval=FALSE---------------------------------------------------------------
## # See the path where libraries are stored:
## .libPaths()
## 
## # See the list of installed packages:
## library()
## 
## # See the list of currently loaded packages:
## search()


## -----------------------------------------------------------------------------
# available.packages() gets a list:
pkgs <- available.packages(filters = "duplicates")
colnames(pkgs)

# We don't need all, just keep the name:
pkgs <- pkgs[,'Package']

# Show the results:
print(paste('Today, there are', nrow(pkgs), 'packages for R.'))


## -----------------------------------------------------------------------------
# Get the list (only names):
my_pkgs <- library()$results[,1]

# Show the results:
print(paste('I have', length(my_pkgs), 'packages for R.'))


## ----eval=FALSE---------------------------------------------------------------
## # See all installed packages:
## installed.packages()


## ----eval=FALSE---------------------------------------------------------------
## # List all out-dated packages:
## old.packages()


## ----eval=FALSE---------------------------------------------------------------
## # Update all available packages:
## update.packages()


## ----eval=FALSE---------------------------------------------------------------
## # Update all packages in batch mode:
## update.packages(ask = FALSE)


## ----eval=FALSE---------------------------------------------------------------
## # Update one package (example with the TTR package):
## install.packages("TTR")


## ----eval=FALSE---------------------------------------------------------------
## t <- readLines(file.choose())


## ----eval=FALSE---------------------------------------------------------------
## t <- readLines("R.book.txt")


## ----inputCSV,echo=TRUE,results='hide',fig.cap=c("The histogram of the CAD", "A scatter-plot of one variable with another.")----
# To read a CSV-file it needs to be in the current directory
# or we need to supply the full path.
getwd()                           # show actual working directory
setwd("/home/philippe/Downloads") # change working directory
data <- read.csv("eurofxref-hist.csv")
is.data.frame(data)
ncol(data)
nrow(data)
head(data)
hist(data$CAD, col = 'khaki3')
plot(data$USD, data$CAD, col = 'red')


## ----csv1,echo=TRUE,fig.cap="The histogram of the most recent values of the CAD only."----
# get the maximum exchange rate
maxCAD <- max(data$CAD)
# use SQL-like selection
d0 <- subset(data, CAD == maxCAD)
d1 <- subset(data, CAD > maxCAD - 0.1)
d1[,1]
d2<- data.frame(d1$Date,d1$CAD)
d2
hist(d2$d1.CAD, col = 'khaki3')


## -----------------------------------------------------------------------------
write.csv(d2,"output.csv", row.names = FALSE)
new.d2 <- read.csv("output.csv")
print(new.d2)


## ----eval=FALSE---------------------------------------------------------------
## # install the package xlsx if not yet done
## if (!any(grepl("xlsx",installed.packages()))){
##   install.packages("xlsx")}
## library(xlsx)
## data <- read.xlsx("input.xlsx", sheetIndex = 1)


## ----eval=FALSE---------------------------------------------------------------
##   if(!any(grepl("xls", installed.packages()))){
##     install.packages("RMySQL")}
##   library(RMySQL)


## ----eval=FALSE---------------------------------------------------------------
## # the connection is stored in an R object myConnection
## # it needs the database name (db_name), username and password
## myConnection = dbConnect(MySQL(),
##    user = 'root',
##    password = 'xxx',
##    dbname = 'db_name',
##    host = 'localhost')
## 
## # e.g. list the tables available in this database.
## dbListTables(myConnection)


## ----eval=FALSE---------------------------------------------------------------
## # Prepare the query for the database
## result <- dbSendQuery(myConnection,
##   "SELECT * from tbl_students WHERE age > 33")
## 
## # fetch() will get us the results, it takes a parameter n, which
## # is the number of desired records.
## # Fetch all the records(with n = -1) and store it in a data frame.
## data <- fetch(result, n = -1)


## ----eval=FALSE---------------------------------------------------------------
## sSQL = ""
## sSQL[1] <- "UPDATE tbl_students
##             SET score = 'A' WHERE raw_score > 90;"
## sSQL[2] <- "INSERT INTO tbl_students
##             (name, class, score, raw_score)
## 	    VALUES ('Robert', 'Grade 0', 88,NULL);"
## sSQL[3] <- "DROP TABLE IF EXISTS tbl_students;"
## for (k in c(1:3)){
##   dbSendQuery(myConnection, sSQL[k])
##   }


## ----eval=FALSE---------------------------------------------------------------
## dbWriteTable(myConnection, "tbl_name",
##              data_frame_name[, ], overwrite = TRUE)


## -----------------------------------------------------------------------------
environment()  # get the environment
rm(list=ls())  # clear the environment
ls()           # list all objects
a <- "a"
f <- function (x) print(x)
ls()           # note that x is not part of.GlobalEnv


## -----------------------------------------------------------------------------
# f
# Multiple actions and side effects to illustrate environments
# Arguments:
#   x -- single type
f <- function(x){
      # define a local function g() within f()
      g <- function(y){
	      b <- "local variable in g()"
              print(" -- function g() -- ")
              print(environment())
              print(ls())
	      print(paste("b is", b))
	      print(paste("c is", c))
              }
     b <- "LOCAL IN F"    # local variable in f()
     c <- "only defined in f(), not in g()"
     g("parameter to g")
     print(" -- function f() -- ")
     print(environment())
     print(ls())
     print(paste("b is", b))
     }

# Test the function f():
a <- pi
f(1)


## -----------------------------------------------------------------------------
x <- 'Philippe'
rm(x)       # make sure the definition is removed
x           # x is indeed not there (generates an error message)
x <- 2L     # now x is created as a long integer
x <- pi     # coerced to double (real)
x <- c(LETTERS[1:5]) # now it is a vector of strings


## -----------------------------------------------------------------------------
# f
# Demonstrates the scope of variables
f <- function() {
  a <- pi     # define local variable
  print(a)    # print the local variable
  print(b)    # b is not in the scope of the function
}

# Define two variables a and b
a <- 1
b <- 2

# Run the function and note that it knows both a and b.
# For b it cannot find a local definition and hence
# uses the definition of the higher level of scope.
f()

# f() did not change the value of a in the environment that called f():
print(a)


## -----------------------------------------------------------------------------
# Citation from the R documentation:
# Copyright (C) 1997-8 The R Core Team
open.account <- function(total) {
     list(
        deposit = function(amount) {
            if(amount <= 0)
                stop("Deposits must be positive!\n")
            total <<- total + amount
            cat(amount,"deposited. Your balance is", total, "\n\n")
        },
        withdraw = function(amount) {
            if(amount > total)
                stop("You do not have that much money!\n")
            total <<- total - amount
            cat(amount,"withdrawn.  Your balance is", total, "\n\n")
        },
        balance = function() {
            cat("Your balance is", total, "\n\n")
        }
        )
 }

ross <- open.account(100)

robert <- open.account(200)

ross$withdraw(30)

ross$balance()

robert$balance()

ross$deposit(50)

ross$balance()

try(ross$withdraw(500)) # no way..


## -----------------------------------------------------------------------------
rm(list=ls())  # clear the environment


## -----------------------------------------------------------------------------
L <- list(matrix(1:16, nrow=4))
L$data_source <- "mainframe 1"
L$get_data_src <- function(){L$data_source}
print(L$get_data_src())


## -----------------------------------------------------------------------------
# Define a string:
acc <- "Philippe"

# Force an attribute, balance, on it:
acc$balance <- 100

# Inspect the result:
acc


## -----------------------------------------------------------------------------
# a function build in core R
typeof(mean)
is.primitive(mean)
# user defined function are "closures:
add1 <- function(x) {x+1}
typeof(add1)
is.function(add1)
is.object(add1)


## -----------------------------------------------------------------------------
# is.S3
# Determines if an object is S3
# Arguments:
#    x -- an object
# Returns:
#   boolean -- TRUE if x is S3, FALSE otherwise
is.S3 <- function(x){is.object(x) & !isS4(x)}

# Create two test objects:
M  <- matrix(1:16, nrow=4)
df <- data.frame(M)

# Test our new function:
is.S3(M)
is.S3(df)


## -----------------------------------------------------------------------------
library(pryr)
otype(M)
otype(df)
otype(df$X1)            # a vector is not S3
df$fac <-factor(df$X4)
otype(df$fac)           # a factor is S3


## -----------------------------------------------------------------------------
mean
ftype(mean)
sum
ftype(sum)


## -----------------------------------------------------------------------------
apropos("print.")
apropos("mean.")


## -----------------------------------------------------------------------------
methods(methods)
methods(mean)


## ----size='scriptsize'--------------------------------------------------------
getS3method("print","table")


## -----------------------------------------------------------------------------
methods(class = "data.frame")


## -----------------------------------------------------------------------------
my_curr_acc <- list("name" = "Philippe", "balance" <- 100)
class(my_curr_acc) <- "account"  # set the class attribute
otype(my_curr_acc)               # it is an S3 object
class(my_curr_acc)               # the class type is defined above

# Note that the class attribute is not visible in the structure:
str(my_curr_acc)                 


## -----------------------------------------------------------------------------
my_object <- structure(list(), class = "boringClass")


## -----------------------------------------------------------------------------
# print.account
# Print an object of type 'account'
# Arguments:
#   x -- an object of type account
print.account <- function(x){
   print(paste("account holder",x[[1]],sep=": "))
   print(paste("balance       ",x[[2]],sep=": "))
   }
print(my_curr_acc)


## -----------------------------------------------------------------------------
# account
# Constructor function for an object of type account
# Arguments:
#    x -- character (the name of the account holder)
#    y -- numeric (the initial balance of the account
# Returns:
#    Error message in console in case of failure.
account <- function(x,y) {
  if (!is.numeric(y))    stop("Balance must be numeric!")
  if (!is.atomic(x))     stop("Name must be atomic!!")
  if (!is.character(x))  stop("Name must be a string!")
  structure(list("name" = x, "balance" = y), class = "account")
}

# create a new instance for Paul:
paul_account <- account("Paul", 200)

# print the object with print.account():
paul_account


## -----------------------------------------------------------------------------
class(paul_account) <- "data.frame"
print(paul_account)   # R thinks now it is a data.frame
paul_account[[2]]     # the data is still correct
class(paul_account) <- "data.frame"
print(paul_account)   # back to normal: the class is just an attribute


## -----------------------------------------------------------------------------
# add_balance
# Dispatcher function to handle the action of adding a given amount 
# to the balance of an account object.
# Arguments:
#   x      -- account -- the account object 
#   amount -- numeric -- the amount to add to the balance
add_balance <- function(x, amount) UseMethod("add_balance")


## -----------------------------------------------------------------------------
# add_balance.account
# Object specific function for an account for the dispatcher 
# function add_balance()
# Arguments:
#   x      -- account -- the account object 
#   amount -- numeric -- the amount to add to the balance
add_balance.account <- function(x, amount) {
   x[[2]] <- x[[2]] + amount; 
   # Note that much more testing and logic can go here
   # It is not so easy to pass a pointer to a function so we 
   # return the new balance:
   x[[2]]}
my_curr_acc <- add_balance(my_curr_acc, 225) 
print(my_curr_acc)


## -----------------------------------------------------------------------------
# add_balance.default
# The default action for the dispatcher function add_balance
# Arguments:
#   x      -- account -- the account object 
#   amount -- numeric -- the amount to add to the balance
add_balance.default <- function(x, amount) {
  stop("Object provided not of type account.")
}


## -----------------------------------------------------------------------------
# probe
# Dispatcher function
# Arguments:
#    x -- account object
# Returns
#    confirmation of object type
probe <- function(x) UseMethod("probe")

# probe.account
# action for account object for dispatcher function probe()
# Arguments:
#    x -- account object
# Returns
#    confirmation of object "account"
probe.account <- function(x) "This is a bank account"

# probe.default
# action if an incorrect object type is provided to probe()
# Arguments:
#    x -- account object
# Returns
#    error message
probe.default <- function(x) "Sorry. Unknown class"

probe (structure(list(), class = "account"))
# No method for class 'customer', fallback to 'account'
probe(structure(list(), class = c("customer", "account")))
# No method for c class, so falls back to default
probe(structure(list(), class = "customer"))
probe(df)         # fallback to default for data.frame
probe.account(df) # force R to use the account method
my_curr_acc <- account("Philippe", 150) # real account
probe(my_curr_acc)


## -----------------------------------------------------------------------------
# Create the object type Acc to hold bank-accounts:
setClass("Acc", 
  representation(holder       = "character", 
                 branch       = "character",
                 opening_date = "Date"))

# Create the object type Bnk (bank):
setClass("Bnk", 
  representation(name = "character", phone = "numeric"))

# Define investment account as a child of Acc:
setClass("CurrAcc", 
  representation(interest_rate = "numeric",
                 balance  = "numeric"), 
  contains = "Acc")
  
# Define investment account as a child of Acc
setClass("InvAcc", 
  representation(custodian = "Bnk"), contains = "Acc")


## -----------------------------------------------------------------------------
# Create an instance of Bnk:
my_cust_bank <- new("Bnk",
                    name = "HSBC",
                    phone = 123456789)

# Create an instance of Acc:
my_acc <- new("Acc", 
             holder       = "Philippe", 
             branch       = "BXL12",
             opening_date = as.Date("2018-10-02"))


## -----------------------------------------------------------------------------
# Check if it is really an S4 object:
isS4(my_cust_bank)

# Change the phone number and check:
my_cust_bank@phone = 987654321  # change the phone number
print(my_cust_bank@phone)       # check if it changed


## -----------------------------------------------------------------------------
# This will do the same as my_cust_bank@phone:
attr(my_cust_bank, 'phone')

# The function also allows partial matching:
attr(my_cust_bank, which='ph', exact = FALSE)

# attr can also change the value of an attribute.
attr(my_cust_bank, which='phone') <- '123123123'
# Let us verify:
my_cust_bank@phone

# It is even possible to create a new attribute or remove one.
attr(my_cust_bank, 'something') <- 'Philippe'
attr(my_cust_bank, 'something')
attr(my_cust_bank, 'something') <- NULL
attr(my_cust_bank, 'something')
str(my_cust_bank) # the something attribute is totally gone


## -----------------------------------------------------------------------------
x <- 1:9
x        # x is a vector
class(x)

attr(x, "dim") <- c(3,3)
x        # is is now a matrix!
class(x) # but R is not fooled.


## -----------------------------------------------------------------------------
slot(my_acc, "holder")


## -----------------------------------------------------------------------------
my_curr_acc <- new("CurrAcc", 
                  holder = "Philippe", 
                  interest_rate = 0.01, 
                  balance=0, 
                  branch = "LDN12", 
                  opening_date= as.Date("2018-12-01"))

# Note that the following does not work and is bound to fail.
also_an_account <- new("CurrAcc", 
                       holder = "Philippe", 
                       interest_rate = 0.01, 
                       balance=0, Acc=my_acc)


## -----------------------------------------------------------------------------
my_curr_acc@balance <- 500


## -----------------------------------------------------------------------------
my_inv_acc <- new("InvAcc", 
                  custodian = my_cust_bank, 
                  holder="Philippe", 
                  branch="DUB01", 
                  opening_date = as.Date("2019-02-21"))

# note that the first slot is another S4 object:
my_inv_acc


## -----------------------------------------------------------------------------
my_inv_acc@custodian        # our custodian bank is HSBC
my_inv_acc@custodian@name   # note the cascade of @ signs
my_inv_acc@custodian@name <- "DB"  # change it to DB
my_inv_acc@custodian@name   # yes, it is changed
my_cust_bank@name           # but our original bank isn't
my_cust_bank@name <- "HSBC Custody" # try something different
my_inv_acc@custodian@name   # did not affect the account
my_inv_acc@custodian@name <- my_cust_bank@name # change back


## -----------------------------------------------------------------------------
getSlots("Acc")


## -----------------------------------------------------------------------------
# Note the mistake in the following code:
my_curr_acc <- new("CurrAcc", 
                  holder = "Philippe", 
                  interest_rate = 0.01, 
                  balance="0",  # Here is the mistake!
                  branch = "LDN12", 
                  opening_date= as.Date("2018-12-01"))


## -----------------------------------------------------------------------------
x_account <- new("CurrAcc", 
                  holder = "Philippe", 
                  interest_rate = 0.01, 
                  #no balance provided
                  branch = "LDN12", 
                  opening_date= as.Date("2018-12-01"))
x_account@balance  # show what R did with it


## -----------------------------------------------------------------------------
setClass("CurrAcc", 
  representation(interest_rate = "numeric",
                 balance       = "numeric"), 
  contains = "Acc",
  prototype(holder       = NA_character_, 
            interst_rate = NA_real_, 
            balance      = 0))

x_account <- new("CurrAcc", 
                  # no holder
                  # no interest rate
                  # no balance
                  branch = "LDN12", 
                  opening_date= as.Date("2018-12-01"))
x_account         # show what R did:


## -----------------------------------------------------------------------------
# this is the function to open an account
.CurrAcc <- function (holder,
                    interest_rate
                    # branch we know from the user
                    # balance should be 0
                    # opening_date is today
                    ) {

  error_msg = "Invalid input while creating an account\n"
  if (is.atomic(holder) & !is.character(holder)) {
    stop(error_msg, "Invalid holder name.")
    }
  if (!(is.atomic(interest_rate) & is.numeric(interest_rate)
      & (interest_rate >= 0) & (interest_rate < 0.1))) {
    stop(error_msg, "Interest rate invalid.")
    }
  br <- "PAR01"  # pretending to find balance by looking up user
  dt <- as.Date(Sys.Date())
  new("CurrAcc", 
                  holder = holder, 
                  interest_rate = interest_rate, 
                  balance=0, 
                  branch = br, 
                  opening_date= dt)
  }

# Create a new account:
lisa_curr_acc <- .CurrAcc("Lisa", 0.01)
lisa_curr_acc


## -----------------------------------------------------------------------------
# Here is the prototype of a dataset that holds some extra
# information in a structured way.
 setClass("myDataFrame",
          contains = "data.frame",
          slots = list(MySQL_DB   = "character",
                       MySQL_tbl  = "character",
                       data_owner = "character"
                       )
          )
xdf <- new("myDataFrame",
    data.frame(matrix(1:9, nrow=3)),
    MySQL_DB = "myCorporateDB@102.12.12.001",
    MySQL_tbl = "tbl_current_accounts",
    data_owner = "customer relationship team")
   
xdf@.Data   
xdf@data_owner


## -----------------------------------------------------------------------------
str(my_inv_acc)
isS4(my_inv_acc)
pryr::otype(my_inv_acc)


## -----------------------------------------------------------------------------
is(my_inv_acc)
is(my_inv_acc, "Acc")


## -----------------------------------------------------------------------------
is.S3


## -----------------------------------------------------------------------------
# setGeneric needs a function, so we need to create it first.

# credit
# Dispatcher function to credit the ledger of an object of 
# type 'account'.
# Arguments:
#    x -- account object
#    y -- numeric -- the amount to be credited
credit <- function(x,y){useMethod()}

# transform our function credit() to a generic one:
setGeneric("credit")

# Add the credit function to the object CurrAcc
setMethod("credit",
   c("CurrAcc"),
   function (x, y) {
     new_bal <- x@balance + y
     new_bal
     }
   )
   
# Test the function:
my_curr_acc@balance
my_curr_acc@balance <- credit(my_curr_acc, 100)
my_curr_acc@balance


## -----------------------------------------------------------------------------
# debet
# Generic function to debet an account
# Arguments:
#    x -- account object
#    y -- numeric -- the amount to be taken from the account
# Returns
#    confirmation of action or lack thereof
debet <- function(x,y){useMethod()}

# Make it a generic function that verifies the balance
# before the account a debet is booked.
setGeneric("debet")

# Add the debet() function as a method for objects of type CurrAcc
setMethod("debet",
   c("CurrAcc"),
   function (x, y) {
     if(x@balance >= y) {
       new_bal <- x@balance + y} else {
       stop("Not enough balance.")
       }
     new_bal
     }
   )
   
# Test the construction:
my_curr_acc@balance  # for reference
my_curr_acc@balance <- debet(my_curr_acc, 50)
my_curr_acc@balance  # the balance is changed

# We do not have enough balance to debet 5000, so the 
# following should fail:
my_curr_acc@balance <- debet(my_curr_acc, 5000)
my_curr_acc@balance  # the balance is indeed unchanged:


## -----------------------------------------------------------------------------
selectMethod("credit", list("CurrAcc"))


## -----------------------------------------------------------------------------
# Note that we capture the returned value of the setRefClass
# Give this always the same name as the class.
account <- setRefClass("account",
            fields = list(ref_number   = "numeric",
                          holder       = "character",
                          branch       = "character",
                          opening_date = "Date",
                          account_type = "character"
                          ),
            # no method yet.
            )
          
x_acc <- account$new(ref_number   = 321654987,
                    holder        = "Philippe",
                    branch        = "LDN05",
                    opening_date  = as.Date(Sys.Date()),
                    account_type  = "current"
                    )
x_acc


## ----eval=FALSE---------------------------------------------------------------
## setRefClass("account", fields = c("ref_number",
##                                   "holder",
##                                   "branch",
##                                   "opening_date",
##                                   "account_type"
##                                   )
##            )


## ----eval=FALSE---------------------------------------------------------------
## setRefClass("account",
##     fields = list(holder,   # accepts all types
##              branch,        # accepts all types
##              opening_date = "Date" # dates only
##             )
##     )


## -----------------------------------------------------------------------------
isS4(account)
# account is S4 and it has a lot more than we have defined:
account


## -----------------------------------------------------------------------------
custBank    <- setRefClass("custBank",
                   fields = list(name =  "character",
                                 phone = "character"
                                 )
                   )
invAccount  <- setRefClass("invAccount",
                   fields = list(custodian = "custBank"),
                   contains = c("account")
                   # methods go here
                   )


## -----------------------------------------------------------------------------
# Definition of RC object currentAccount
currAccount <- setRefClass("currentAccount",
                   fields = list(interest_rate = "numeric",
                                 balance       = "numeric"),
                   contains = c("account"),
                   methods = list(
                      credit = function(amnt) {
                            balance <<- balance + amnt
                            },
                      debet =  function(amnt) {
                            if (amnt <= balance) {
                               balance <<- balance - amnt
                               } else {
                               stop("Not rich enough!")
                               }
                         }
                      )
                   )
# note how the class reports on itself:
currAccount


## -----------------------------------------------------------------------------
ph_acc <- currAccount$new(ref_number    = 321654987,
                          holder        = "Philippe",
                          branch        = "LDN05",
                          opening_date  = as.Date(Sys.Date()), 
                          account_type  = "current",
                          interest_rate = 0.01,
                          balance       = 0  
                          )


## -----------------------------------------------------------------------------
ph_acc$balance     # after creating balance is 0:
ph_acc$debet(200)  # impossible (not enough balance)
ph_acc$credit(200) # this is bound to fail
ph_acc$balance     # the money arrived in our account
ph_acc$debet(100)  # this is possible
ph_acc$balance     # the money is indeed gone


## -----------------------------------------------------------------------------
alsoCurrAccount <- setRefClass("currentAccount",
                   fields = list(
                             interest_rate = "numeric",
                             balance       = "numeric"),
                   contains = c("account")
                   )
alsoCurrAccount$methods(list(
                      credit = function(amnt) {
                          balance <<- balance + amnt
                          },
                      debet = function(amnt) {
                          if (amnt <= balance) {
                             balance <<- balance - amnt
                             } else {
                             stop("Not rich enough!")
                             }
                            }
                      ))


## -----------------------------------------------------------------------------
# we assume that you installed the package before:
# install.packages("tidyverse")
# so load it:
library(tidyverse)


## ----fig.cap="The sum of sine and cosine illustrated."------------------------
x <- seq(from = 0, to = 2 * pi, length.out = 100)
s <- sin(x)
c <- cos(x)
z <- s + c
plot(x, z, type = "l",col="red", lwd=7)
lines(x, c, col = "blue",  lwd = 1.5)
lines(x, s, col = "darkolivegreen", lwd = 1.5)


## -----------------------------------------------------------------------------
x <- seq(from = 0, to = 2 * pi, length.out = 100)
#df <- as.data.frame((x))
df <- rbind(as.data.frame((x)),cos(x),sin(x), cos(x) + sin(x))
# plot etc.


## -----------------------------------------------------------------------------
library(tidyverse)
x <- seq(from = 0, to = 2 * pi, length.out = 100)
tb <- tibble(x, sin(x), cos(x), cos(x) + sin(x))
# plot and to the rest
# Tibbles are quite similar to data frames, 
# but compare the output:


## ----fig.cap="A tibble plots itself like a data-frame."-----------------------
# not how concise and relevant the output it
print(tb)  

# this does the same as for a data-frame:
plot(tb)

# actually a tibble will still behave as a data frame
is.data.frame(tb)


## -----------------------------------------------------------------------------
tb$x[1]
tb$`sin(x)`[1]


## -----------------------------------------------------------------------------
tb <- tibble(`1` = 1:3, `2` = sin(`1`), `1`*pi, 1*pi)
tb


## -----------------------------------------------------------------------------
# -- data frame -- 
df <- data.frame("value" = pi, "name" = "pi")
df$na        # partial matching of column names

# automatic conversion to factor, plus data frame
# accepts strings:
df[,"name"]   

df[,c("name", "value")] # 

# -- tibble -- 
df <- tibble("value" = pi, "name" = "pi")
df$name       # column name
df$nam        # no partial matching but error msg.
df[,"name"]   # this returns a tibble (no simplification)
df[,c("name", "value")] # no conversion to factor


## -----------------------------------------------------------------------------
tb <- tibble(c("a", "b", "c"), c(1,2,3), 9L,9)
is.data.frame(tb)

# Note also that tibble did no conversion to factors, and
# note that the tibble also recycles the scalars:
tb

# Coerce the tibble to data-frame:
as.data.frame(tb)  

# A tibble does not recycle shorter vectors, so this fails:
fail <- tibble(c("a", "b", "c"), c(1,2))
# That is a major advantage and will save many programming errors.


## -----------------------------------------------------------------------------
t <- tibble("x" = runif(10))                     
t <- within(t, y <- 2 * x + 4 + rnorm(10, mean=0,sd=0.5))


## -----------------------------------------------------------------------------
t <- tibble("x" = runif(10))  %>%
     within(y <- 2 * x + 4 + rnorm(10, mean=0,sd=0.5))


## ----eval=FALSE---------------------------------------------------------------
## # 1. pipe:
## a %>% f()
## # 2. pipe with shortened function:
## a %>% f
## # 3. is equivalent with:
## f(a)


## -----------------------------------------------------------------------------
a <- c(1:10)
a %>% mean()
a %>% mean
mean(a)


## ----eval=FALSE---------------------------------------------------------------
## # The following line
## c <- a    %>%
##      f()
## # is equivalent with:
## c <- f(a)
## 
## # Also, it is easy to see that
## x <- a %>% f(y) %>% g(z)
## # is the same as:
## x <- g(f(a, y), z)


## -----------------------------------------------------------------------------
# f1
# Dummy function that from which only the error throwing part 
# is shown.
f1 <- function() {
    # Here goes the long code that might be doing something risky 
    # (e.g. connecting to a database, uploading file, etc.)
    # and finally, if it goes wrong:
    stop("Early exit from f1!")  # throw error
    }
    
tryCatch(f1(),    # the function to try 
         error   = function(e) {paste("_ERROR_:",e)},
         warning = function(w) {paste("_WARNING_:",w)},
         message = function(m) {paste("_MESSSAGE_:",m)},
         finally="Last command"    # do at the end
         )


## -----------------------------------------------------------------------------
# f1
# Dummy function that from which only the error throwing part 
# is shown.
f1 <- function() {
    # Here goes the long code that might be doing something risky 
    # (e.g. connecting to a database, uploading file, etc.)
    # and finally, if it goes wrong:
    stop("Early exit from f1!")  # something went wrong
    }   %>%
tryCatch(
         error   = function(e) {paste("_ERROR_:",e)},
         warning = function(w) {paste("_WARNING_:",w)},
         message = function(m) {paste("_MESSSAGE_:",m)},
         finally="Last command"    # do at the end
         )
# Note that it fails in silence.


## -----------------------------------------------------------------------------
# This will not work, because lm() is not designed for the pipe.
lm1 <- tibble("x" = runif(10))                            %>%
       within(y <- 2 * x + 4 + rnorm(10, mean=0, sd=0.5)) %>%
       lm(y ~ x)


## -----------------------------------------------------------------------------
# The Tidyverse only makes the %>% pipe available. So, to use the
# special pipes, we need to load magrittr 
library(magrittr) 
lm2 <- tibble("x" = runif(10))                           %>%
       within(y <- 2 * x + 4 + rnorm(10, mean=0,sd=0.5)) %$%
       lm(y ~ x)
summary(lm2)


## -----------------------------------------------------------------------------
coeff <- tibble("x" = runif(10))                           %>%
         within(y <- 2 * x + 4 + rnorm(10, mean=0,sd=0.5)) %$%
         lm(y ~ x)                                         %>%
         summary                                           %>% 
         coefficients
coeff


## ----fig.cap="A linear model fit on generated data to illustrate the piping command."----
library(magrittr)
t <- tibble("x" = runif(100))                           %>%
     within(y <- 2 * x + 4 + rnorm(10, mean=0, sd=0.5)) %T>%
     plot(col="red")   # The function plot does not return anything
                       # so we used the %T>% pipe.
		       
lm3 <-   t                  %$%
         lm(y ~ x)          %T>% # pass on the linear model
         summary            %T>% # further pass on the linear model
         coefficients

tcoef <- lm3 %>% coefficients  # we anyhow need the coefficients

# Add the model (the solid line) to the previous plot:
abline(a = tcoef[1], b=tcoef[2], col="blue", lwd=3)


## -----------------------------------------------------------------------------
x <- c(1,2,3) 

# The following line:
x <- x %>% mean  

# is equivalent with the following:
x %<>% mean

# Show x:
x


## -----------------------------------------------------------------------------
library(pryr)
x <- runif(100)
object_size(x)
y <- x

# x and y together do not take more memory than only x.
object_size(x,y)   

y <- y * 2

# Now, they are different and are stored separately in memory.
object_size(x,y)   


## -----------------------------------------------------------------------------
# The mean of a vector
x <- c(1,2,3,4,5,60)
mean(x)

# Missing values will block the override the result
x <- c(1,2,3,4,5,60,NA)
mean(x)
# Missing values can be ignored with na.rm = TRUE
mean(x, na.rm = TRUE)

# This works also for a matrix
M <- matrix(c(1,2,3,4,5,60), nrow=3)
mean(M)


## -----------------------------------------------------------------------------
v <- c(1,2,3,4,5,6000)
mean(v)
mean(v, trim = 0.2)


## -----------------------------------------------------------------------------
returns <- c(0.5,-0.5,0.5,-0.5)

# arithmetic mean
aritmean <- mean(returns)

# the ln-mean
log_returns <- returns
for(k in 1:length(returns)) {
  log_returns[k] <- log( returns[k] + 1)
  }
logmean <- mean(log_returns)
exp(logmean) - 1

# What is the value of the investment after these returns
V_0 <- 1
V_T <- V_0
for(k in 1:length(returns)) {
  V_T <- V_T * (returns[k] + 1)
  }
V_T

# Compare this to our predictions
## mean of log-returns
V_0 * (exp(logmean) - 1)
## mean of returns
V_0 * (aritmean + 1)


## -----------------------------------------------------------------------------
x <- c(1:5,5e10,NA)
x
median(x)               # no meaningful result with NAs
median(x,na.rm = TRUE)  # ignore the NA
# Note how the median is not impacted by the outlier,
# but the outlier dominates the mean:
mean(x, na.rm = TRUE)


## -----------------------------------------------------------------------------
# my_mode 
# Finds the first mode (only one)
# Arguments:
#   v -- numeric vector or factor
# Returns:
#   the first mode
my_mode <- function(v) {
   uniqv <- unique(v)
   tabv  <- tabulate(match(v, uniqv))
   uniqv[which.max(tabv)]
   }

# now test this function
x <- c(1,2,3,3,4,5,60,NA)
my_mode(x)
x1 <- c("relevant", "N/A", "undesired", "great", "N/A", 
        "undesired", "great", "great")
my_mode(x1)

# text from https://www.r-project.org/about.html
t <- "R is available as Free Software under the terms of the 
Free Software Foundation's GNU General Public License in 
source code form. It compiles and runs on a wide variety of
UNIX platforms and similar systems (including FreeBSD and 
Linux), Windows and MacOS."
v <- unlist(strsplit(t,split=" "))
my_mode(v)


## -----------------------------------------------------------------------------
# my_mode 
# Finds the mode(s) of a vector v
# Arguments:
#   v -- numeric vector or factor
#   return.all -- boolean -- set to true to return all modes
# Returns:
#   the modal elements
my_mode <- function(v, return.all = FALSE) {
  uniqv  <- unique(v)
  tabv   <- tabulate(match(v, uniqv))
  if (return.all) {
    uniqv[tabv == max(tabv)]
  } else {
    uniqv[which.max(tabv)]
  }
}

# example:
x <- c(1,2,2,3,3,4,5)
my_mode(x)
my_mode(x, return.all = TRUE)


## -----------------------------------------------------------------------------
t <- rnorm(100, mean=0, sd=20)
var(t)
sd(t)
sqrt(var(t))
sqrt(sum((t - mean(t))^2)/(length(t) - 1))


## -----------------------------------------------------------------------------
mad(t)
mad(t,constant=1)


## -----------------------------------------------------------------------------
cor(mtcars$hp,mtcars$wt)


## -----------------------------------------------------------------------------
d <- data.frame(mpg = mtcars$mpg, wt = mtcars$wt, hp = mtcars$hp)
# Note that we can feed a whole data-frame in the functions.
var(d)
cov(d)
cor(d)
cov2cor(cov(d))


## -----------------------------------------------------------------------------
x  <- c(-10:10)
df <- data.frame(x=x, x_sq=x^2, x_abs=abs(x), x_exp=exp(x))
cor(df)


## -----------------------------------------------------------------------------
cor(rank(mtcars$hp),rank(mtcars$wt))


## -----------------------------------------------------------------------------
x  <- c(-10:10)
cor(rank(x), rank(x^2))


## -----------------------------------------------------------------------------
# we use the dataset mtcars from MASS
df <- data.frame(mtcars$cyl,mtcars$am)
chisq.test(df)


## ----fig.cap=c("A comparison between a set of random numbers drawn from the normal distribution (khaki) and the theoretical shape of the normal distribution in blue. ")----
 obs <- rnorm(600,10,3)
 hist(obs,col="khaki3",freq=FALSE)
 x <- seq(from=0,to=20,by=0.001)
 lines(x, dnorm(x,10,3),col="blue",lwd=4)


## ----fig.cap=c("The same plot for the returns of the SP500 index seems acceptable, though there are outliers (where the normal distribution converges fast to zero).")----
 library(MASS)
 hist(SP500,col="khaki3",freq=FALSE,border="khaki3")
 x <- seq(from=-5,to=5,by=0.001)
 lines(x, dnorm(x,mean(SP500),sd(SP500)),col="blue",lwd=2)


## ----fig.cap=c("A Q-Q plot is a good way to judge if a set of observations is normally distributed or not.\\label{fig:qq1}")----
 library(MASS)
 qqnorm(SP500,col="red"); qqline(SP500,col="blue")


## ----fig.cap=c("The probability to get maximum x tails when flipping a fair coin, illustrated with the binomial distribution.")----
# Probability of getting 5 or less heads from 10 tosses of 
# a coin.
pbinom(5,10,0.5)

# visualize this for one to 10 numbers of tosses
x <- 1:10
y <- pbinom(x,10,0.5)
plot(x,y,type="b",col="blue", lwd=3,
    xlab="Number of tails",
    ylab="prob of maxium x tails",
    main="Ten tosses of a coin")


## -----------------------------------------------------------------------------
# How many heads should we at least expect (with a probability 
# of 0.25) when a coin is tossed 10 times.
qbinom(0.25,10,1/2)


## -----------------------------------------------------------------------------
# Find 20 random numbers of tails from and event of 10 tosses 
# of a coin
rbinom(20,10,.5)


## -----------------------------------------------------------------------------
N <- 100
t <- data.frame(id = 1:N, result = rnorm(N))
summary(t)


## -----------------------------------------------------------------------------
library(tidyverse) # not only for %>% but also for group_by, etc.
# In mtcars the type of the car is only in the column names,
# so we need to extract it to add it to the data
n <- rownames(mtcars)  

# Now, add a column brand (use the first letters of the type)
t <- mtcars  %>%
     mutate(brand = str_sub(n, 1, 4))    # add column


## -----------------------------------------------------------------------------
# First, we need to find out which are the most abundant brands
# in our dataset (set cutoff at 2: at least 2 cars in database)
top_brands <- count(t, brand) %>% filter(n >= 2)

# top_brands is not simplified to a vector in the tidyverse
print(top_brands)


grouped_cars <- t                      %>% # start with cars
   filter(brand %in% top_brands$brand) %>% # only top-brands
   group_by(brand)                     %>%
   summarise(
       avgDSP = round(mean(disp), 1),
       avgCYL = round(mean(cyl),  1),
       minMPG = min(mpg),
       medMPG = median(mpg),
       avgMPG = round(mean(mpg),2),
       maxMPG = max(mpg),
     )
print(grouped_cars)


## -----------------------------------------------------------------------------
# Each call to summarise() removes a layer of grouping:
by_vs_am <- mtcars %>% group_by(vs, am)
by_vs <- by_vs_am %>% summarise(n = n())
by_vs
by_vs %>% summarise(n = sum(n))

# To removing grouping, use ungroup:
by_vs %>%
  ungroup() %>%
  summarise(n = sum(n))

# You can group by expressions: this is just short-hand for
# a mutate/rename followed by a simple group_by:
mtcars %>% group_by(vsam = vs + am)

# By default, group_by overrides existing grouping:
mtcars             %>%
  group_by(cyl)    %>%
  group_by(vs, am) %>%
  group_vars()

# Use add = TRUE to append grouping levels:
mtcars                         %>%
  group_by(cyl)                %>%
  group_by(vs, am, add = TRUE) %>%
  group_vars()
     


## ----fig.cap=c('The plot-function will generate a scatter-plot for a vector. Note also that the legend is automatically adapted. The \xaxis is the index of the number in the vector and the $y$-axis is the value of the corresponding number in the vector.\\label{fig:plot1}','The plot-function will generate a scatter-plot of each column in function of each other column for a data frame. The main diagonal of the matrix has the names of the columns. Note also that the x and y axis are labelled only on one side and that each row shares the same $y$-axis, while each column shares one \xaxis.\\label{fig:plot2}')----
x <- c(1:20)^2
plot(x)

df <- data.frame('a' = x, 'b' = 1/x, 'c' = log(x), 'd' = sqrt(x))
plot(df)


## ----echo=FALSE,fig.cap="Some plot characters. Most other characters will just plot themselves.\\label{fig:pch}"----
# This sets up an empty plotting field
plot(x = c(0, 4.5),
     y = c(0, 5),
     main = "Some pch arguments",
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     cex.main = 2.6, 
     col = "white"
)

# This will plot all of the standard pch arguments
y = rep(5:0, each=5)
for (i in 0:25) { 
  points(x = i %% 5, y = y[i+1], pch = i,cex = 2, col="blue", bg="khaki3")
  text(0.3 + i %% 5, y = y[i+1], i, cex = 2)
}
for (i in 1:2) {
  ch <- LETTERS[i]
  points(x = i, y = 0, pch = ch,cex = 2, col="red")
  text(0.3 + i, y = 0, ch, cex = 2)
}
for (i in 1:2) {
  ch <- letters[i]
  points(x = i + 2, y = 0, pch = ch,cex = 2, col="red")
  text(0.3 + i + 2, y = 0, ch, cex = 2)
}


## ----fig.cap="A scatter-plot needs an x and a y variable.\\label{fig:scatter1}"----
# prepare the data
library(MASS)

# To make this example more interesting, we convert mpg to l/100km

# mpg2l
# Converts miles per gallon into litres per 100 km
# Arguments:
#    mpg -- numeric -- fuel consumption in MPG
# Returns:
#    Numeric -- fuel consumption in litres per 100 km
mpg2l <- function(mpg = 0) {
100 * 3.785411784 / 1.609344 / mpg}

mtcars$l <- mpg2l(mtcars$mpg)
plot(x = mtcars$hp,y = mtcars$l, xlab = "Horse Power", 
 ylab = "L per 100km", main = "Horse Power vs Milage",
 pch = 22, col="red", bg="yellow")


## ----fig.cap="A line plot of the type b."-------------------------------------
# prepare the data
years <- c(2000,2001,2002,2003,2004,2005)
sales <- c(2000,2101,3002,2803,3500,3450)
plot(x = years,y = sales, type = 'b', 
     xlab = "Years", ylab = "Sales in USD",
     main = "The evolution of our sales")
points(2004,3500,col="red",pch=16) # highlight one point
text(2004,3400,"top sales")        # annotate the highlight


## ----fig.cap="A pie-chart in R."----------------------------------------------
x <- c(10, 20, 12)            # Create data for the graph
labels <- c("good", "average", "bad")
#for saving to a file:
  png(file = "feedback.jpg")  # Give the chart file a name
  pie(x,labels)               # Plot the chart
  dev.off()                   # Save the file
pie(x,labels)                 # Show in the R Graphics screen


## ----fig.cap="A standard bar-chart based on a vector.\\label{fig:bar1}"-------
sales <- c(100,200,150,50,125)
regions <- c("France", "Poland", "UK", "Spain", "Belgium")
barplot(sales, width=1, 
      xlab="Regions", ylab="Sales in EUR", 
	main="Sales 2016", names.arg=regions,
	border="blue", col="brown")


## ----fig.cap="A bar-chart based on a matrix will produce stacked bars. Note how nicely this plot conveys the seasonal trend in the data.\\label{fig:bar2}"----
# Create the input vectors.
colours <- c("orange","green","brown")
regions <- c("Mar","Apr","May","Jun","Jul")
product <- c("License","Maintenance","Consulting")

# Create the matrix of the values.
values <- matrix(c(20,80,0,50,140,10,50,80,20,10,30,
     10,25,60,50), nrow = 3, ncol = 5, byrow = FALSE)

# Create the bar chart.
barplot(values, main = "Sales 2016",
 names.arg = regions, xlab = "Region",
 ylab = "Sales in EUR", col = colours)

# Add the legend to the chart.
legend("topright", product, cex = 1.3, fill = colours)


## ----fig.cap='A boxplot where the total of each bar equals 100\\%. Note how the seasonal trend is obscured in this plot, but how it now tells the story of how the consulting business is growing and the maintenance business, relatively, cannot keep up wih that growth.\\label{fig:boxpl100}'----
# We reuse the matrix 'values' from previous example.

# Add extra space to right of plot area; change clipping to figure
par(mar = c(5, 4, 4, 8) + 0.1, # default margin was c(5, 4, 4, 2) + 0.1
    xpd = TRUE)      # TRUE to restrict all plotting to the plot region

# Create the plot with all totals coerced to 1.0:
barplot(prop.table(values, 2), main = "Sales 2016",
   names.arg = regions, xlab = "Region",
   ylab = "Sales in EUR", col = colours)

# Add the legend, but move it to the right with inset:
legend("topright", product, cex = 1.0, inset=c(-0.3,0), fill = colours,
       title="Business line")


## ----fig.cap=c("Boxplots show information about the central tendency (median) as well as the spread of the data.\\label{fig:boxplot1}")----
library(MASS)
boxplot(mpg ~ cyl,data=mtcars,col="khaki3", 
        main="MPG by number of cylinders")


## ----fig.cap="Violin plot as provided by the function {\tt vioplot} from the package of the same name.\\label{fig:violin1}"----
# install.packages('vioplot') # only do once
library(vioplot)              # load the vioplot library
with(mtcars , vioplot(mpg[cyl==4] , mpg[cyl==6], mpg[cyl==8],
                      col=rgb(0.1,0.4,0.7,0.7) , 
		      names=c("4","6","8") 
		      )
    ) 


## ----fig.cap=c("Violin plot as traced by {\\tt geom\\_violin} provided by the library {\\tt ggplot2}; with the colouring done according to the number of cylinders.\\label{fig:violin2a}", "Violin plot as traced by {\\tt geom\\_violin} provided by the library {\\tt ggplot2}; with the colouring done according to the factoring of the number of cylinders.\\label{fig:violin2b}")----
# Library
library(ggplot2)
 
# First, type of color
ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_violin(aes(fill = cyl))
 
# Second type
ggplot(mtcars, aes(factor(cyl), mpg)) +
  geom_violin(aes(fill = factor(cyl)))


## ----fig.cap=c("A histogram in R is produced by the hist() function.\\label{fig:hist1}","In this histogram, the breaks are changed, and the y-axes is now calibrated as a probability. Note that leaving freq=TRUE would give the wrong impression that there are more observations in the wider brackets.\\label{fig:hist2}")----
library(MASS)
incidents <- ships$incidents
# figure 1: with a rug and fixed breaks
hist(incidents, 
   col=c("red","orange","yellow","green","blue","purple"))
rug(jitter(incidents))  # add the tick-marks

# figure 2: user-defined breaks for the buckets
hist(incidents, 
  col=c("red","orange","yellow","green","blue","purple"),
  ylim=c(0,0.3), breaks=c(0,2,5,10,20,40,80),freq=FALSE)


## ----fig.cap="Two line plots plotted by the function curve()."----------------
fn1 <- function(x) sqrt(1-(abs(x)-1)^2)
fn2 <- function(x) -3*sqrt(1-sqrt(abs(x))/sqrt(2))
curve(fn1,-2,2,ylim=c(-3,1),col="red",lwd=4, 
      ylab = expression(sqrt(1-(abs(x)-1)^2) +++ fn_2))
curve(fn2,-2,2,add=TRUE,lw=4,col="red")
text(0,-1,expression(sqrt(1-(abs(x)-1)^2)))
text(0,-1.25,"++++")
text(0,-1.5,expression(-3*sqrt(1-sqrt(abs(x))/sqrt(2))))


## ----fig.cap='A colour mapping combined with a contour plot provides a nice image of the heights of Auckland\'s Maunga Whau Volcano.\\label{fig:volcano}'----
# A basic plot is possible with the function image()
# image(volcano)
# (try it yourself)

# improved plot as per documentation of image().
x <- 1:nrow(volcano)
y <- 1:ncol(volcano)

# the mapping of colours.
image(x, y, volcano, col = terrain.colors(100), 
      axes = FALSE, xlab='', ylab='')

# add the contour plot
contour(x, y, volcano, levels = seq(90, 200, by = 5),
        add = TRUE, col = "brown")
	
# add axis, a box and a title:
axis(1, at = seq(10, 80, by = 10))
axis(2, at = seq(10, 60, by = 10))
box()
title(main = "Auckland's Maunga Whau Volcano", font.main = 4)


## ----heatmap,fig.cap=c("Heatmap for the ``mtcars'' data.")--------------------
d=as.matrix(mtcars,scale="none")
heatmap(d)


## ----heatmap2,fig.cap=c("Heatmap for the ``mtcars'' data with all columns rescaled")----
heatmap(d,scale="column")


## ----eval=FALSE---------------------------------------------------------------
## heatmap(x, Rowv = NULL, Colv = if(symm) "Rowv" else NULL,
##         distfun = dist, hclustfun = hclust,
##         reorderfun = function(d, w) reorder(d, w),
##         add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
##         scale = c("row", "column", "none"), na.rm = TRUE,
##         margins = c(5, 5), ColSideColors, RowSideColors,
##         cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
##         labRow = NULL, labCol = NULL, main = NULL,
##         xlab = NULL, ylab = NULL,
##         keep.dendro = FALSE, verbose = getOption("verbose"),
## 	...)


## ----eval=FALSE---------------------------------------------------------------
## # if neccesary first download the packages
## install.packages("tm")           # text mining
## install.packages("SnowballC")    # text stemming
## install.packages("RColorBrewer") # colour palettes
## install.packages("wordcloud")    # word-cloud generator


## ----results='hide'-----------------------------------------------------------
# Then load the packages:
library("tm")
library("SnowballC")
library("RColorBrewer")
library("wordcloud")


## -----------------------------------------------------------------------------
# In this example we use a text version of this very book:
t <- readLines("data/r-book.txt")

# Then create a corpus of text 
doc <- Corpus(VectorSource(t))


## ----warning=FALSE------------------------------------------------------------
# The file has still a lot of special characters
# e.g. the following replaces '\', '#', and '|' with space:
toSpace <- content_transformer(function (x , pattern ) 
                               gsub(pattern, " ", x))
doc <- tm_map(doc, toSpace, "\\\\")
doc <- tm_map(doc, toSpace, "#")
doc <- tm_map(doc, toSpace, "\\|")
# Note that the backslash needs to be escaped in R


## ----warning=FALSE,message=FALSE----------------------------------------------
# Convert the text to lower case
doc <- tm_map(doc, content_transformer(tolower))

# Remove numbers
doc <- tm_map(doc, removeNumbers)

# Remove english common stopwords
doc <- tm_map(doc, removeWords, stopwords("english"))

# Remove your own stop words
# specify your stopwords as a character vector:
doc <- tm_map(doc, removeWords, c("can","value","also","price",
              "cost","option","call","need","possible","might",
	      "will","first","etc","one","portfolio", "however",
	      "hence", "want", "simple", "therefore")) 

# Remove punctuations
doc <- tm_map(doc, removePunctuation)

# Eliminate extra white spaces
doc <- tm_map(doc, stripWhitespace)

# Text stemming
#doc <- tm_map(doc, stemDocument)


## -----------------------------------------------------------------------------
library(stringi)

wordToReplace <- c('functions', 'packages')
ReplaceWith   <- c('function',  'package')

doc <- tm_map(doc,  function(x) stri_replace_all_fixed(x, 
           wordToReplace, ReplaceWith, vectorize_all = FALSE))


## -----------------------------------------------------------------------------
dtm <- TermDocumentMatrix(doc)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)


## ----frequBar,fig.cap=c("The frequency of the ten most occurring words in this text. Note tha the name of the software, R, is only one letter, and hence not retained as a word. So R is not not part of this list.")----
barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="khaki3", main ="Most frequent words",
        ylab = "Word frequencies")


## ----wordcloud,fig.cap=c("A word-could for the text of this book. This single image gives a good idea of what this book is about."),out.height="\\textwidth",out.width="\\textwidth",fig.width="\\textwidth",fig.height="\\textwidth"----
set.seed(1879)
wordcloud(words = d$word, freq = d$freq, min.freq = 10,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))


## ----eval=FALSE---------------------------------------------------------------
##   wordcloud(words, freq, scale = c(4,.5), min.freq = 3,
##             max.words=Inf, random.order = TRUE,
## 	    random.color=FALSE, rot.per=.1,
##             colors = "black", ordered.colors = FALSE,
## 	    use.r.layout=FALSE, fixed.asp=TRUE, ...)


## -----------------------------------------------------------------------------
findFreqTerms(dtm, lowfreq = 150)


## -----------------------------------------------------------------------------
# e.g. for the word "function"
findAssocs(dtm, terms = "function", corlimit = 0.15)


## -----------------------------------------------------------------------------
# find colour numbers that contain the word 'khaki'
grep("khaki",colours())

# find the names of those colours
colors()[grep("khaki",colours())]


## -----------------------------------------------------------------------------
# extract the rgb value of a named colour
col2rgb("khaki3")


## ----colSchemes,warning=FALSE,fig.cap='An illustration of six predefined colour schemes in R. This figures uses the dataset {\\tt mtcars} and shows in each sub-plot the number of seconds the car needs to speed up to one fourth of a mile in function of its power ({\\tt hp}). The colouring of each observation (car type) is done in function of the fuel economy {\\tt mpg}. This presentation allows to visualize more than one relationship in one plot.'----
library(ggplot2)
library(gridExtra)
p <- ggplot(mtcars, aes(x = hp, y = qsec, color = mpg)) + 
      geom_point(size=5) 

# no manual colour specification:
p0 <- p + ggtitle('Default colour range')
      
# using named colours:
p1 <- p + scale_color_gradient(low="red", high="green") + 
      ggtitle('Named colour')

# using hexadecimal representation
p2 <- p + scale_color_gradient(low="#0011ee", high="#ff55a0") + 
      ggtitle('Hexadecimal colour')

# RGB definition of colour
p3 <- p + scale_color_gradient(low=rgb(0,0,1), high=rgb(0.8,0.8,0)) + 
      ggtitle('RGB colour')

# rainbow colours
p4 <- p + scale_color_gradientn(colours = rainbow(5)) + 
      ggtitle('Rainbow colours')

# using another colour set
p5 <- p + scale_color_gradientn(colours = terrain.colors(5)) + 
      ggtitle('Terrain colours') + 
      theme_light()  # turn off the grey background

grid.arrange(p0, p1, p2, p3, p4, p5, ncol = 2)


## ----fig.cap='A visualisation of all built in colours in R. Note that the number of the colour can be determined as by taking the y-value minus one times nine plus the x-value.\\label{fig:cols}', out.height="0.9\\textheight",out.width="\\textwidth",fig.width="\\textwidth",fig.height="0.9\\textheight"----
N   <- length(colours())  # this is 657
df  <- data.frame(matrix(1:N, nrow=73, byrow = TRUE))
image(1:(ncol(df)), 1:(nrow(df)), as.matrix(t(df)), 
      col=colours(),
      xlab = "X", ylab = "Y")


## -----------------------------------------------------------------------------
colours()[(3 - 1)  * 9 + 8]  
colours()[(50 - 1) * 9 + 1]


## ----colourSets,fig.cap='Examples of discrete colour sets. The name of the colour-set is the title of each plot.'----
p <- ggplot(mtcars) +
     geom_histogram(aes(cyl, fill=factor(cyl)), bins=3)

# no manual colour specification:
p0 <- p + ggtitle('default colour range')


# using a built-in colour scale
p1 <- p + scale_fill_grey() + 
      ggtitle('Shades of grey')


library(RColorBrewer)
p2 <- p + scale_fill_brewer() + 
      ggtitle('RColorBrewer')
      
p3 <- p + scale_fill_brewer(palette='Set2') + 
      ggtitle('RColorBrewer Set2')
      
p4 <- p + scale_fill_brewer(palette='Accent') + 
      ggtitle('RColorBrewer Accent')
      
grid.arrange(p0, p1, p2, p3, p4, p5, ncol = 2)


## ----eval=FALSE---------------------------------------------------------------
## ts(data = NA, start = 1, end = numeric(), frequency = 1,
##         deltat = 1, ts.eps = getOption("ts.eps"), class = ,
## 	names = )


## -----------------------------------------------------------------------------
library(MASS)
# The SP500 is available as a numeric vector:
str(SP500)


## -----------------------------------------------------------------------------
# Convert it to a time series object.
SP500_ts <- ts(SP500,start = c(1990,1),frequency = 260)


## -----------------------------------------------------------------------------
# Compare the original:
class(SP500)

# with:
class(SP500_ts)


## ----ts,fig.cap=c("The standard plot for a time series object for the returns of the SP500 index in the 1990s.")----
plot(SP500_ts)


## -----------------------------------------------------------------------------
  val = c(339.97)
  for (k in 2:length(SP500)){
    val[k] = val[k-1] * (SP500[k-1] / 100 + 1)
  }
  # Convert both series to a matrix
  M <- matrix(c(SP500,val),nrow=length(SP500))
  
  # Convert the matrix to a time series object
  SP <- ts(M, start=c(1990,1),frequency=260)
  colnames(SP) <- c("Daily Return in Pct","Value")


## ----tsPlot2,fig.cap=c("The standard plot functionality of time series will keep the $z$-axis for both variables the same (even use one common axis).")----
  plot(SP, type = "l", main = "SP500 in the 1990s")


## ----gdpforecast,fig.cap=c("A first plot to show the data before we start. This will allow us to select a suitable method for forecasting.","Using the library forecast with a simple moving average.")----
g <- read.csv('data/gdp/gdp_pol_sel.csv') # get some data
attach(g) # the names of the data are now always available
plot(year, GDP.per.capitia.in.current.USD, type='b', 
     lwd=3, xlab = 'Year', ylab = 'Polish GDP per Capita in USD')


## ----results='hide',warning=FALSE---------------------------------------------
require(forecast)
# make the forecast with the moving average (ma)
g.data <- ts(g$GDP.per.capitia.in.current.USD,start=c(1990))
g.movav = forecast(ma(g.data, order=3), h=5)


## ----forecastMoveAv,fig.cap=c("A forecast based on moving average.")----------
# show the result:
plot(g.movav,col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="Income in current USD")
lines(year,GDP.per.capitia.in.current.USD,col="red",type='b')


## ----warning=FALSE,message=FALSE----------------------------------------------
# testing accuracy of the model by sampling
g.ts.tst <- ts(g.data[1:20],start=c(1990))
g.movav.tst <- forecast(ma(g.ts.tst,order=3),h=5)
accuracy(g.movav.tst, g.data[22:26])


## ----fcmovav1,fig.cap=c("A backtest for our forecast.")-----------------------
plot(g.movav.tst,col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="Income in current USD")
lines(year, GDP.per.capitia.in.current.USD, col="red",type='b')


## -----------------------------------------------------------------------------
train = ts(g.data[1:20],start=c(1990))
test  = ts(g.data[21:26],start=c(2010))
arma_fit <- auto.arima(train)
arma_forecast <- forecast(arma_fit, h = 6)
arma_fit_accuracy <- accuracy(arma_forecast, test)
arma_fit; arma_forecast; arma_fit_accuracy


## ----fig.cap=c("Optimal moving average forecast.")----------------------------
plot(arma_forecast, col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="income in current USD")
lines(year,GDP.per.capitia.in.current.USD,col="red",type='b')


## -----------------------------------------------------------------------------
g.exp <- ses(g.data,5,initial="simple")
g.exp  # simple exponential smoothing uses the last value as
       # the forecast and finds confidence intervals around it


## ----plotExpSm,fig.cap=c("Forecasting with an exponentially smoothed moving average.")----
plot(g.exp,col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="income in current USD")
lines(year,GDP.per.capitia.in.current.USD,col="red",type='b')


## -----------------------------------------------------------------------------
g.exp <- holt(g.data,5,initial="simple")
g.exp  # Holt exponential smoothing


## ----HoltWinters,fig.cap=c("Holt exponentially smoothed moving average.")-----
plot(g.exp,col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="income in current USD")
lines(year,GDP.per.capitia.in.current.USD,col="red",type='b')


## ----stl,fig.cap=c("Using the stl-function to decompose data in a seasonal part and a trend.")----
# we use the data nottem
# Average Monthly Temperatures at Nottingham, 1920-1939
nottem.stl = stl(nottem, s.window="periodic")
plot(nottem.stl)


## ----eval=FALSE---------------------------------------------------------------
## add_ts <- log(exp_ts)


## -----------------------------------------------------------------------------
# Simple exponential: models level
fit <- HoltWinters(g.data, beta=FALSE, gamma=FALSE)

# Double exponential: models level and trend
fit <- HoltWinters(g.data, gamma=FALSE)

# Triple exponential: models level, trend, and seasonal 
# components. This fails on the example, as there is no 
# seasonal trend:
#fit <- HoltWinters(g.data)  

# Predictive accuracy
library(forecast)
accuracy(forecast(fit,5))


## -----------------------------------------------------------------------------
# predict next 5 future values
forecast(fit, 5)


## ----HoltWint,fig.cap=c("The Holt-Winters model fits an exponential trend. Here we plot the double exponential model.")----
plot(forecast(fit, 5),col="blue",lw=4, 
     main="Forecast of GDP per capita of Poland",
     ylab="income in current USD")
lines(year,GDP.per.capitia.in.current.USD,col="red",type='b')


## ----HWnottem,fig.cap=c("The Holt-Winters model applied to the temperatures in Nottingham.")----
# Use the Holt-Winters method for the temperatures
n.hw <- HoltWinters(nottem)
n.hw.fc <- forecast(n.hw,50)
plot(n.hw.fc)


## ----eval=FALSE---------------------------------------------------------------
## # install.packages('RMySQL')
## library(RMySQL)
## # connect to the library
## con <- dbConnect(MySQL(),
##                  user     = "librarian",
##                  password = "librarianPWD",
##                  dbname   = "library",
##                  host     = "localhost"
##                  )
## 
## # in case we would forget to disconnect:
## on.exit(dbDisconnect(con))


## ----eval=FALSE---------------------------------------------------------------
## # show some information
## 
## show(con)
## summary(con, verbose = TRUE)
## # dbGetInfo(con)  # similar as above but in list format
## dbListResults(con)
## dbListTables(con) # check: this might generate too much output
## 
## # get data
## df_books <- dbGetQuery(con, "SELECT COUNT(*) AS nbrBooks
##                      FROM tbl_author_book GROUP BY author;")
## 
## # Now, df_books is a data frame that can be used as usual.
## 
## # close the connection
## dbDisconnect(con)


## -----------------------------------------------------------------------------
# Load the package: 
library(RMySQL)

# db_get_data
# Get data from a MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- the connection info to 
#                                          the MySQL database
#    sSQL     -- character string       -- the SQL statement that 
#                                          selects the records
# Returns
#    data.frame, containing the selected records
db_get_data <- function(con_info, sSQL){
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  df <- dbGetQuery(con, sSQL)
  dbDisconnect(con)
  df
}

# db_run_sql
# Run a query that returns no data in an MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- open connection
#    sSQL     -- character string       -- the SQL statement to run
db_run_sql <-function(con_info, sSQL)
{
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  rs <- dbSendQuery(con,sSQL)
  dbDisconnect(con)
}


## ----fig.cap="Histogram generated with data from the MySQL database.\\label{fig:bookHist}"----
# use the wrapper functions to get data.

# step 1: define the connection info
my_con_info <- list()
my_con_info$user     <- "librarian"
my_con_info$password <- "librarianPWD"
my_con_info$dbname   <- "library"
my_con_info$host     <- "localhost"


# step 2: get the data
my_query <- "SELECT COUNT(*) AS nbrBooks 
                     FROM tbl_author_book GROUP BY author;"
df <- db_get_data(my_con_info, my_query)

# step 3: use this data to produce the histogram:
hist(df$nbrBooks, col='khaki3')


## -----------------------------------------------------------------------------
# -- reset query cache 
sql_reset_query_cache <- function (con_info)
{
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
# remove all cache:
system("sync && echo 3 | sudo tee /proc/sys/vm/drop_caches")
# clear MySQL cache cache and disconnect:
rc <- dbSendQuery(con, "RESET QUERY CACHE;")
dbDisconnect(con)
# once more remove all cache:
system("sync && echo 3 | sudo tee /proc/sys/vm/drop_caches")
}


## -----------------------------------------------------------------------------
# install.packages('sodium') # do only once
                           # fails if you do not have libsodium-dev
library(sodium)

# Create the SHA256 key based on a secret password:
key <- sha256(charToRaw("My sweet secret"))

# Serialize the data to be encrypted:
msg <- serialize("Philippe J.S. De Brouwer", NULL)

# Encrypt:
msg_encr <- data_encrypt(msg, key)


orig <- data_decrypt(msg_encr, key)
stopifnot(identical(msg, orig))

# Tag the message with your key (HMAC):
tag <- data_tag(msg, key)


## -----------------------------------------------------------------------------
# -- 
library(RMySQL)

# -- The functions as mentioned earlier:
# db_get_data
# Get data from a MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- containing the connection 
#                                          info to the MySQL database
#    sSQL     -- character string       -- the SQL statement that selects 
#                                          the records
# Returns
#    data.frame, containing the selected records
db_get_data <- function(con_info, sSQL){
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  df <- dbGetQuery(con, sSQL)
  dbDisconnect(con)
  df
}

# db_run_sql
# Run a query that returns no data in an MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- containing the connection 
#                                          info to the MySQL database
#    sSQL     -- character string       -- the SQL statement to be run
db_run_sql <-function(con_info, sSQL)
{
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  rs <- dbSendQuery(con,sSQL)
  dbDisconnect(con)
}


## ----warning=FALSE,message=FALSE----------------------------------------------
## Load dplyr via tidyverse:
library(tidyverse)

# Define the wrapper functions:

# Step 1: define the connection info.
my_con_info <- list()
my_con_info$user     <- "librarian"
my_con_info$password <- "librarianPWD"
my_con_info$dbname   <- "library"
my_con_info$host     <- "localhost"


# -- The data import was similar to what we had done previously.
# -- However, now we import all tables separately
# Step 2: get the data
my_tables <- c("tbl_authors", "tbl_author_book",
               "tbl_books", "tbl_genres")
my_db_names <- c("authors", "author_book", 
               "books", "genres")

# Loop over the four tables and download their data:
for (n in 1:4) {
  my_sql <- paste("SELECT * FROM `",my_tables[n],"`;", sep="")
  df <- db_get_data(my_con_info, my_sql)
  # the next line uses tibbles are from the tidyverse
  as_tibble(assign(my_db_names[n],df))  
 }

# Step 3: do something with the data
# -- will follow in the remainder of the section


## -----------------------------------------------------------------------------
str(authors)


## -----------------------------------------------------------------------------
# example to illustrate the parser functions

v       <- c("1.0", "2.3", "2.7", ".")
nbrs    <- c("$100,00.00", "12.4%")
s_dte   <- "2018-05-03"

# The parser functions can generate from these strings 
# more specialized types.

# The following will generate an error:
parse_double(v)         # reason: "." is not a number
parse_double(v, na=".") # Tell R what the encoding of NA is
parse_number(nbrs)
parse_date(s_dte)


## -----------------------------------------------------------------------------
parse_guess(v)
parse_guess(v, na=".")
parse_guess(s_dte)
guess_parser(v)
guess_parser(v[1:3])
guess_parser(s_dte)
guess_parser(nbrs)


## -----------------------------------------------------------------------------
s_csv = "'a','b','c'\n001,2.34,.\n2,3.14,55\n3,.,43"
read_csv(s_csv)
read_csv(s_csv, na='.')  # Tell R how to understand the '.'


## -----------------------------------------------------------------------------
# Method 1: before the actual import
spec_csv(s_csv, na='.')

# Method 2: check after facts:
t <- read_csv(s_csv, na='.')
spec(t)


## -----------------------------------------------------------------------------
read_csv(s_csv, na='.',
  col_names = TRUE,
  cols(
    `'a'` = col_character(),
    `'b'` = col_double(),
    `'c'` = col_double()      # coerce to double
    )
  )


## ----eval=FALSE---------------------------------------------------------------
## # Start with:
## t <- read_csv(readr_example("challenge.csv"))
## 
## # Then, to see the issues, do:
## problems(t)
## 
## # Notice that the problems start in row 1001, so
## # the first 1000 rows are special cases. The first improvement
## # can be obtained by increase the guesses
## ## compare
## spec_csv(readr_example("challenge.csv"))
## ## with
## spec_csv(readr_example("challenge.csv"), guess_max = 1001)


## -----------------------------------------------------------------------------
# load readr
library(readr)
# Or load the tidyverse with library(tidyverse), it includes readr.


## -----------------------------------------------------------------------------
# Make a string that looks like a fixed-width table (shortened):
txt <- " book_id  year  title                                           genre                                                                                                                                                                   
       1  1896  Les plaisirs et les jour                        LITmod 
       2  1927  Albertine disparue                              LITmod 
       3  1954  Contre Sainte-Beuve                             LITmod 
       8  1687  PhilosophiÃ¦ Naturalis Principia Mathematica      SCIphy 
       9  -300  Elements (translated )                          SCImat 
      10  2014  Big Data World                                  SCIdat 
      11  2016  Key Business Analytics                          SCIdat 
      12  2011  Maslowian Portfolio Theory                      FINinv 
      13  2016  R for Data Science                              SCIdat"


## -----------------------------------------------------------------------------
fileConn<-file("books.txt")
writeLines(txt, fileConn)
close(fileConn)

my_headers <- c("book_id","year","title","genre")


## -----------------------------------------------------------------------------
# Reading the fixed-width file
# -- > by indicating the widths of the columns
t <- read_fwf(
  file="./books.txt",   
  skip=1,               # skip one line with headers
  fwf_widths(c(8, 6, 48, 8), my_headers)
  )
  
#inspect the input:
print(t)

# -- > same but naming directly
t <- read_fwf(
  file="./books.txt",   
  skip=1,               # skip one line with headers
  fwf_cols(book_id = 8, year = 6, 
           title = 48, genre = 8)
  )

# -- > by selecting columns (by indicating begin and end):
t2 <- read_fwf(
   file = "books.txt", 
   skip = 1,
   fwf_cols(year = c(11, 15), 
            title = c(17, 63))
   )

# -- > by guessing the columns
# The function fwf_empty can help to guess where the columns start 
# based on white space
t3 <- read_fwf(
  file = "books.txt",
  skip = 1,
  fwf_empty("books.txt")
  ) 


## -----------------------------------------------------------------------------
print(t3)


## -----------------------------------------------------------------------------
head(mtcars)
mtcars[1,1]          # mpg is in the first column
rownames(mtcars[1,]) # the name of the car is not a column


## ----warning=FALSE------------------------------------------------------------
## -- 
## -- Load dplyr via tidyverse
library(tidyverse)
library(RMySQL)

## -- The functions as mentioned earlier:

# db_get_data
# Get data from a MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- the connection info
#                                          to the MySQL database
#    sSQL     -- character string       -- the SQL statement that  
#                                          selects the records
# Returns
#    data.frame, containing the selected records
db_get_data <- function(con_info, sSQL){
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  df <- dbGetQuery(con, sSQL)
  dbDisconnect(con)
  df
}

# db_run_sql
# Run a query that returns no data in an MySQL database
# Arguments:
#    con_info -- MySQLConnection object -- the connection info
#                                          to the MySQL database
#    sSQL     -- character string       -- the SQL statement 
#                                          to be run
db_run_sql <-function(con_info, sSQL)
{
  con <- dbConnect(MySQL(), 
                 user     = con_info$user, 
                 password = con_info$password, 
                 dbname   = con_info$dbname, 
                 host     = con_info$host
                 )
  rs <- dbSendQuery(con,sSQL)
  dbDisconnect(con)
}


# use the wrapper functions to get data.

# step 1: define the connection info
my_con_info <- list()
my_con_info$user     <- "librarian"
my_con_info$password <- "librarianPWD"
my_con_info$dbname   <- "library"
my_con_info$host     <- "localhost"


## -- Import 2 tables combined
# step 2: get the data
my_sql <- "SELECT * FROM tbl_authors 
   JOIN tbl_author_book ON author_id = author
   JOIN tbl_books       ON book      = book_id
   JOIN tbl_genres      ON genre     = genre_id;"
t_mix <- db_get_data(my_con_info, my_sql)
t_mix <- as.tibble(t_mix)  %>%  print



## ----warning=FALSE------------------------------------------------------------
# Make a table of how much each author_id occurs:
nbr_auth <- t_mix  %>% count(author_id)

# Do the same and include all fields that are assumed to
# be part of the table authors.
nbr_auth2 <- t_mix    %>% 
  count(author_id, pen_name, full_name, birth_date, death_date, book)

nbr_auth$n - nbr_auth2$n


## -----------------------------------------------------------------------------
# Try without book:
nbr_auth2 <- t_mix    %>% 
  count(author_id, pen_name, full_name, birth_date, death_date)
  
# Now these occurrences are the same:
nbr_auth$n - nbr_auth2$n


## -----------------------------------------------------------------------------
my_authors <- tibble(author_id = t_mix$author_id, 
                    pen_name   = t_mix$pen_name,
                    full_name  = t_mix$full_name,
                    birth_date = t_mix$birth_date,
                    death_date = t_mix$death_date
                    )      %>%
              unique       %>%
              print


## -----------------------------------------------------------------------------
auth <- tibble(
            author_id = as.integer(my_authors$author_id), 
            pen_name   = my_authors$pen_name,
            full_name  = my_authors$full_name,
            birth_date = as.Date(my_authors$birth_date),
            death_date = as.Date(my_authors$death_date)
               )       %>%
        unique         %>%
        print


## -----------------------------------------------------------------------------
auth$full_name <- str_replace(auth$full_name, "\n", "")     %>%
   print


## -----------------------------------------------------------------------------
# First read in some data (using a flat file to remind
# how this works):
 x <- " January   100       102       108
 February  106       105       105
 March     104       104       106
 April     120       122       118
 May       130       100       133
 June      141       139       135
 July      175       176       180
 August    170       188       187
 September 142       148       155
 October   133       137       145
 November  122       128       131
 December  102       108       110"

# Read in the flat file via read_fwf from readr:
t <- read_fwf(x,  fwf_empty(x, col_names = my_headers))

# Set the column names:
colnames(t) <-  c("month", "Sales2017", "Sales2018", "Sales2019")

# Finally, we can show the data as it appeared in the spreadsheet
# from the sales department:
print(t)


## -----------------------------------------------------------------------------
t2 <- gather(t, "year", "sales", 2:4)
t2$year <- str_sub(t2$year,6,9)  # delete the sales word
t2$year <- as.integer(t2$year)   # convert to integer
t2


## ----fig.cap="Finally, we are able to make a plot of the tibble in a way that makes sense and allows to see the trends in the data.\\label{fig:salesdata}"----
library(lubridate)
t2$date <- parse_date_time(paste(t2$year,t$month), orders="ym")
plot(x=t2$date, y=t2$sales, col="red")
lines(t2$date, t2$sales, col="blue")


## -----------------------------------------------------------------------------
library(dplyr)
sales_info <- data.frame(
       time = as.Date('2016-01-01') + 0:9 + rep(c(0,-1), times=5),
       type  = rep(c("bought","sold"),5),
       value = round(runif(10, min = 0, max = 10001))
       )
sales_info

spread(sales_info, type, value)


## -----------------------------------------------------------------------------
sales_info                %>%
  spread(type, value)     %>%
  gather(type, value, 2:3)


## -----------------------------------------------------------------------------
library(tidyr)
turnover <- data.frame(
       what = paste(as.Date('2016-01-01') + 0:9 + rep(c(0,-1), times=5),
                    rep(c("HSBC","JPM"),5), sep="/"),
       value = round(runif(10, min = 0, max = 50))
       )
turnover

separate(turnover, what, into=c("date","counterpart"), sep="/")


## -----------------------------------------------------------------------------
# Use as separator a digit followed by a forward slash 
# and then a capital letter.
separate(turnover, what, into=c("date","counterpart"), 
         sep="[0-9]/[A-Z]")


## -----------------------------------------------------------------------------
library(tidyr)

# Define a data frame:
df <- data.frame(year = 2018, month = 0 + 1:12, day = 5)
print(df)

# Merge the columns to one variable:
unite(df, 'date', 'year', 'month', 'day', sep = '-')


## -----------------------------------------------------------------------------
library(dplyr)


## -----------------------------------------------------------------------------
# Using the example of the library:
dplyr::select(genres,      # the first argument is the tibble
       genre_id, location) # then a list of column names


## -----------------------------------------------------------------------------
a1 <- filter(authors, birth_date > as.Date("1900-01-01"))
paste(a1$pen_name,"--",a1$birth_date)


## -----------------------------------------------------------------------------
authors           %>%
 count(author_id) %>%
 filter(n > 1)    %>%
 nrow()


## -----------------------------------------------------------------------------
author_book     %>%
  count(author) %>%
  filter(n > 1)


## -----------------------------------------------------------------------------
a2 <- books                                         %>%
  inner_join(genres, by = c("genre" = "genre_id"))  
paste(a2$title, "-->", a2$location)


## -----------------------------------------------------------------------------
a3 <- authors                                            %>%
  inner_join(author_book, by = c("author_id"="author"))  %>%
  inner_join(books,       by = c("book"     ="book_id")) %>%
  inner_join(genres,      by = c("genre"    ="genre_id"))%>%
  dplyr::select(pen_name, location)                      %>%
  arrange(location)                                      %>%
  unique()                                               %>%
  print()

# Note the difference with the base-R code below!
b <- merge(authors, author_book, by.x="author_id", 
                                 by.y = "author")
b <- merge(b, books,  by.x="book",  by.y = "book_id")
b <- merge(b, genres, by.x="genre", by.y = "genre_id")
b <- cbind(b$pen_name, b$location)  # colnames disappear 
colnames(b) <- c("pen_name", "location")
b <- as.data.frame(b)
b <- b[order(b$location), ] # sort for data frames is order
b <- unique(b)
print(b)


## ----eval=FALSE---------------------------------------------------------------
## arrange(desc(location))


## -----------------------------------------------------------------------------
a3[!duplicated(a3$location), ]


## ----eval=FALSE---------------------------------------------------------------
## inner_join(A, B, by = c("z", "z")  # ambiguous, but works
## inner_join(A, B, by = "z")         # shorter
## inner_join(A, B)                   # shortest


## -----------------------------------------------------------------------------
ab <- authors                                            %>%
  inner_join(author_book, by = c("author_id"="author"))  %>%
  inner_join(books,       by = c("book"     ="book_id")) %>%
  add_count(author_id)
ab$n


## -----------------------------------------------------------------------------
genres                        %>%
   mutate(genre = genre_id)   %>%  # add column genre
   inner_join(books)          %>%  # leave out the "by="
   dplyr::select(c(title, location))


## -----------------------------------------------------------------------------
t <- authors                                                   %>%
     mutate(short_name = str_sub(pen_name,1,7))                %>%
     mutate(x_name = if_else(str_length(pen_name) > 15,
                            paste(str_sub(pen_name,1,8),
                                  "...",
                                  str_sub(pen_name, 
                                         start = -3),
                                  sep=''),
                            pen_name,
                            "pen_name is NA"
                           )
           )                                                   %>%
     mutate(is_alive = 
       if_else(!is.na(birth_date) & is.na(death_date), 
            "YES",
            if_else(death_date < Sys.Date(), 
                "no", 
                "maybe"),
             "NA")
            )                                                  %>%
    dplyr::select(c(x_name, birth_date, death_date, is_alive)) %>%
    print()


## -----------------------------------------------------------------------------
authors    %>%
  transmute(name = full_name, my_date = as.Date(birth_date) -5)


## -----------------------------------------------------------------------------
authors    %>%
  filter(!is.na(birth_date) & is.na(death_date)) %>%
  transmute(name = full_name, my_date = as.Date(birth_date) -5)  


## ----'setOperators'-----------------------------------------------------------
# Define two sets (with one column):
A <- tibble(col1 = c(1L:4L)) 
B <- tibble(col1 = c(4L,4L,5L))

# Study some of the set-operations:
dplyr::intersect(A,B)
union(A,B)
union_all(A,B)
setdiff(A,B)
setequal(A,B)

# The next example uses a data-frame with two columns:
A <- tibble(col1 = c(1L:4L), 
            col2 = c('a', 'a', 'b', 'b'))
B <- tibble(col1 = c(4L,4L,5L), 
            col2 = c('b', 'b', 'c'))

# Study the same set-operations:
dplyr::intersect(A,B)
union(A,B)
union_all(A,B)
setdiff(A,B)
setequal(A,B)


## ----stringr1-----------------------------------------------------------------
library(tidyverse)
library(stringr)

# define strings
s1 <- "Hello"  # double quotes are fine
s2 <- 'world.' # single quotes are also fine

# Return the length of a string:
str_length(s1)

# Concatenate strings:
str_c(s1, ", ", s2)       # str_c accepts many strings
str_c(s1, s2, sep = ", ") # str_c also has a 


## -----------------------------------------------------------------------------
s <- 'World'
str_c('Hello, ', s, '.')


## -----------------------------------------------------------------------------
library(stringr)                    # or library(tidyverse)
sVector <- c("Hello", ", ", "world", "Philippe")

str_sub (sVector,1,3)               # the first 3 characters
str_sub (sVector,-3,-1)             # the last 3 characters
str_to_lower(sVector[4])            # convert to lowercase
str_to_upper(sVector[4])            # convert to uppercase
str_c(sVector, collapse=" ")        # collapse into one string
str_flatten(sVector, collapse=" ")  # flatten string
str_length(sVector)                 # length of a string

# Nest the functions:
str_c(str_to_upper(str_sub(sVector[4],1,4)),
      str_to_lower(str_sub(sVector[4],5,-1))
     )

# Use pipes:
sVector[4]        %>%
   str_sub(1,4)   %>%
   str_to_upper()


## -----------------------------------------------------------------------------
str <- "abcde"

# Replace from 2nd to 4th character with "+"
str_sub(str, 2, 4) <- "+"  
str


## -----------------------------------------------------------------------------
str <- "F0"
str_dup(str, c(2,3))  # duplicate a string


## -----------------------------------------------------------------------------
str <- c(" 1 ", "  abc", "Philippe De Brouwer   ")
str_pad(str, 5)  # fills with white-space to x characters
# str_pad never makes a string shorter!
# So to make all strings the same length we first truncate:
str      %>%
  str_trunc(10) %>%
  str_pad(10,"right")   %>%
  print
  
# Remove trailing and leading white space:
str_trim(str)
str_trim(str,"left")

# Modify an existing string to fit a line length:
"The quick brown fox jumps over the lazy dog. "  %>%
   str_dup(5)    %>%
   str_c         %>%  # str_flatten also removes existing \n
   str_wrap(50)  %>%  # Make lines of 50 characters long.
   cat                # or writeLines (print shows "\n")


## -----------------------------------------------------------------------------
str <- c("a", "z", "b", "c")

# str_order informs about the order of strings (rank number)
str_order(str)

# Sorting is done with str_sort.
str_sort(str)


## ----stringRegex1-------------------------------------------------------------
library(stringr)   # or library(tidyverse)
sV <- c("philosophy", "physiography", "phis", 
        "Hello world", "Philippe", "Philosophy", 
        "physics", "philology")

# extracting substrings that match a regex pattern
str_extract(sV, regex("Phi"))
str_extract(sV, "Phi")        # the same, regex assumed


## ----stringrRegex2------------------------------------------------------------
str_extract(sV, "(p|P)hi")
# or do it this way
str_extract(sV, "(phi|Phi)")
# but the last form is more limited: for example.
str_extract(sV, "(p|P)h(i|y)")
# is equivalent to
str_extract(sV, "(phi|Phi|phy|Phy)")


## ----stringrRegex3------------------------------------------------------------
str_extract(sV, "(p|P)h(i|y)[^lL]")


## ----stringrRegex4------------------------------------------------------------
str_extract(sV, "(?i)Ph(i|y)[^(?i)L]")


## ----stringrRegexEmail--------------------------------------------------------
emails <- c("god@heaven.org", "philippe@de-brouwer.com",
            "falsemaail@nothingmy", "mistaken.email.@com")
regX <- 
 "^([a-zA-Z0-9_\\-\\.]+)@([a-zA-Z0-9_\\-\\.]+)\\.([a-zA-Z]{2,5})$"
str_extract(emails, regX)


## -----------------------------------------------------------------------------
str_extract("Philippe", "Ph\\w*")  # is greedy
str_extract("Philippe", "Ph\\w*?") # is lazy


## -----------------------------------------------------------------------------
library(rex)
valid_chars <- rex(one_of(regex('a-z0-9\u00a1-\uffff')))

# example: validate an url
expr <- rex(
  start,       # start of the string: ^

  # protocol identifier (optional) + //
  group(list('http', maybe('s')) %or% 'ftp', '://'),

  # user:pass authentication (optional)
  maybe(non_spaces,
    maybe(':', zero_or_more(non_space)),
    '@'),

  #host name
  group(zero_or_more(valid_chars, 
        zero_or_more('-')), 
        one_or_more(valid_chars)),

  #domain name
  zero_or_more('.', 
              zero_or_more(valid_chars, 
              zero_or_more('-')), 
	      one_or_more(valid_chars)),

  #TLD identifier
  group('.', valid_chars %>% at_least(2)),

  # server port number (optional)
  maybe(':', digit %>% between(2, 5)),

  # resource path (optional)
  maybe('/', non_space %>% zero_or_more()),

  end
)

# The rest is only to print it nice (expr can be used)
substring(expr, seq(1,  nchar(expr)-1, 40),
                seq(41, nchar(expr),   40))      %>%
  str_c(sep="\n")                              


## -----------------------------------------------------------------------------
string <- c("one:1", "NO digit", "c5c5c5", "d123d", "123", 6)
pattern <- "\\d"


## -----------------------------------------------------------------------------
# grep() returns the whole string if a match is found
grep(pattern, string, value = TRUE) 

# The default for value is FALSE -> only returns indexes
grep(pattern, string)               

# L for returning a logical variable
grepl(pattern, string)              

# --- stringr ---
# similar to grepl (note order of arguments!)
str_detect(string, pattern)        


## -----------------------------------------------------------------------------
# Locate the first match (the numbers are the position in the string)
regexpr (pattern, string)  

# grepexpr() finds all matches and returns a list
gregexpr(pattern, string)  

# --- stringr ---
# finds the first match and returns a matrix
str_locate(string, pattern) 

# Finds all matches and returns a list (same as grepexpr)
str_locate_all(string, pattern)


## -----------------------------------------------------------------------------
# First, we need additionally a replacement (repl)
repl <- "___"

# sub() replaces the first match:
sub(pattern, repl, string)

# gsub() replaces all matches.
gsub(pattern, repl, string)

# --- stringr ---
# str_replace() replaces the first match.
str_replace(string, pattern, repl)

# str_replace_all() replaces all mathches.
str_replace_all(string, pattern, repl)


## -----------------------------------------------------------------------------
# regmatches() with regexpr() will extract only the first match
regmatches(string, regexpr(pattern, string))

# regmatches() with gregexpr() will extract all matches
regmatches(string, gregexpr(pattern, string)) # all matches

# --- stringr --- 
# Extract the first match
str_extract(string, pattern)

# Similar as str_extract, but returns column instead of row
str_match(string, pattern)

# Extract all matches (list as return)
str_extract_all(string, pattern)

# To get a neat matrix output, add simplify=T
str_extract_all(string, pattern, simplify=TRUE)

# Similar to str_extract_all (but returns column instead of row)
str_match_all(string, pattern)



## -----------------------------------------------------------------------------
# --- base-R ---
strsplit(string, pattern)

# --- stringr ---
str_split(string, pattern)


## ----LoadLubriDate------------------------------------------------------------
# Load the tidyverse for its functionality such as pipes:
library(tidyverse)

# Lubridate is not part of the core-tidyverse, so we need
# to load it separately:
library(lubridate)  


## -----------------------------------------------------------------------------
as.numeric(Sys.time())   # the number of seconds passed 1 January 1970
as.numeric(Sys.time()) / (60 * 60 * 24 * 365.2422)


## -----------------------------------------------------------------------------
# There is a list of functions that convert to a date
mdy("04052018")
mdy("4/5/2018")
mdy("04052018")
mdy("4052018")  # ambiguous formats are refused!
dmy("04052018") # same string, different date


## -----------------------------------------------------------------------------
dt <- ymd(20180505)  %>% print
as.numeric(dt)
ymd(17656)


## -----------------------------------------------------------------------------
yq("201802")


## -----------------------------------------------------------------------------
ymd_hms("2018-11-01T13:59:00") 
dmy_hms("01-11-2018T13:59:00") 

ymd_hm("2018-11-01T13:59") 
ymd_h("2018-11-01T13") 

hms("13:14:15")
hm("13:14")


## -----------------------------------------------------------------------------
as_date("2018-11-12")
as_date(0)
as_date(-365)
as_date(today()) - as_date("1969-02-21")


## -----------------------------------------------------------------------------
today()
now()


## -----------------------------------------------------------------------------
# Note it converts the system time-zone to UTC
as_datetime("2006-07-22T14:00")  

# Force time-zone
as_datetime("2006-07-22T14:00 UTC")

as_datetime("2006-07-22 14:00 Europe/Warsaw") #Fails silently!

dt <- as_datetime("2006-07-22 14:00",tz="Europe/Warsaw") %>% 
      print  

# Get the same date-time numerals in a different time-zone:
force_tz(dt, "Pacific/Tahiti")

# Get the same cosmic moment in a new time-zone
with_tz(dt,  "Pacific/Tahiti")


## -----------------------------------------------------------------------------
today(tzone = "Pacific/Tahiti")
date_decimal(2018.521, tz="UTC")


## -----------------------------------------------------------------------------
dt1 <- make_datetime(year = 1890, month = 12L, day = 29L, 
              hour = 8L, tz = 'MST')
dt1


## -----------------------------------------------------------------------------
# We will use the date from previous hint:
dt1

year(dt)      # extract the year
month(dt)     # extract the month
week(dt)      # extract the week
day(dt)       # extract the day
wday(dt)      # extract the day of the week as number
qday(dt)      # extract the day of the quarter as number
yday(dt)      # extract the day of the year as number
hour(dt)      # extract the hour
minute(dt)    # extract the minutes
second(dt)    # extract the seconds
quarter(dt)   # extract the quarter 
semester(dt)  # extract the  semester
am(dt)        # TRUE if morning
pm(dt)        # TRUE if afternoon
leap_year(dt) # TRUE if leap-year


## -----------------------------------------------------------------------------
# We will use the date from previous example:
dt1

# Experiment changing it:
update(dt, month=5)
update(dt, year=2018)
update(dt, hour=18)


## -----------------------------------------------------------------------------
moment1 <- as_datetime("2018-10-28 01:59:00", tz="Europe/Warsaw")
moment2 <- as_datetime("2018-10-28 02:01:00", tz="Europe/Warsaw")
moment2 - moment1  # could be 2 minutes or 1 hour and 2 minutes
moment3 <- as_datetime("2018-10-28 03:01:00", tz="Europe/Warsaw")

# The clocks were put back in this tz from 3 to 2am.
# So, there is 2 hours difference between 2am and 3am!
moment3 - moment1 


## -----------------------------------------------------------------------------
# Calculate the duration in seconds:
dyears(x = 1/365)
dweeks(x = 1)
ddays(x = 1)
dhours(x = 1)
dminutes(x = 1)
dseconds(x = 1)
dmilliseconds(x = 1) 
dmicroseconds(x = 1) 
dnanoseconds(x = 1) 
dpicoseconds(x = 1) 

# Note that a duration object times a number is again a 
# Duration object:
dpicoseconds(x = 1) * 10^12

# Investigate the object type:
dur <- dnanoseconds(x = 1)
class(dur)
str(dur)
print(dur)


## -----------------------------------------------------------------------------
# useful for automation:
duration(5, unit = "years") 

# coerce and logical
dur <- dyears(x = 10)
as.duration(60 * 60 * 24)
as.duration(dur)
is.duration(dur)
is.difftime(dur)
as.duration(dur)
make_difftime(60, units="minutes")


## -----------------------------------------------------------------------------
years(x = 1)
months(x = 1)
weeks(x = 1)
days(x = 1)
hours(x = 1)
minutes(x = 1)
seconds(x = 1)
milliseconds(x = 1)
microseconds(x = 1)
nanoseconds(x = 1)
picoseconds(x = 1)

# Investigate the object type:
per <- days(x = 1)
class(per)
str(per)
print(per)


## -----------------------------------------------------------------------------
# for automations:
period(5, unit = "years") 

# coerce timespan to period
as.period(5, unit="years") 


as.period(10)
p <- seconds_to_period(10) %>%
     print
period_to_seconds(p)


## -----------------------------------------------------------------------------
years(1) + months(3) + days(13)


## -----------------------------------------------------------------------------
d1 <- ymd_hm("1939-09-01 09:00", tz="Europe/Warsaw")
d2 <- ymd_hm("1945-08-15 12:00", tz="Asia/Tokyo")

interval(d1, d2)  # defines the interval
# or
ww2 <- d1 %--% d2 # defines the same interval

ww2 / days(1)   # the period expressed in days
ww2 / ddays(1)  # duration in terms of days
# The small difference is due to DST and equals one hour:
(ww2 / ddays(1) - ww2 / days(1)) * 24

# Allow the interval to report on its length:
int_length(ww2) / 60 / 60 / 24


## -----------------------------------------------------------------------------
d_date <- ymd("19450430")

# Is a date or ww2erval in another>
d_date %within% ww2
ph <- interval(ymd_hm("1941-12-07 07:48", tz = "US/Hawaii"),
                     ymd_hm("1941-12-07 09:50", tz = "US/Hawaii")
		     )
ph %within% ww2    # is ph in ww2?
int_aligns(ph, ww2) # do ww2 and ph share start or end?

# shift forward or backward
int_shift(ww2, years(1))
int_shift(ww2, years(-1))

# Swap start and end moment
flww2 <- int_flip(ww2)

# coerce all to "positive" (start-date before end-date)
int_standardize(flww2)

# Modify start or end date
int_start(ww2) <- d_date; print(ww2)
int_end(ww2)  <-  d_date; print(ww2)


## -----------------------------------------------------------------------------
dts <- c(ymd("2000-01-10"), ymd("1999-12-28"),
         ymd("1492-01-01"), ymd("2100-10-15")
	 )
round_date(dts, unit="month")
floor_date(dts, unit="month")
ceiling_date(dts, unit="month")
# Change a date to the last day of the previous month or
# to the first day of the month with rollback()
rollback(dts, roll_to_first = FALSE, preserve_hms = TRUE) 


## -----------------------------------------------------------------------------
set.seed(1911)
s <- tibble(reply = runif(n = 1000, min = 0, max = 13))
hml <- function (x = 0) {
  if (x < 0)  return(NA)
  if (x <= 4) return("L")
  if (x <= 8) return("M")
  if (x <= 12) return("H")
  return(NA)
  }
surv <- apply(s, 1, FUN=hml)  # output is a vector
surv <- tibble(reply = surv)  # coerce back to tibble
surv


## -----------------------------------------------------------------------------
# 1. Define the factor-levels in the right order:
f_levels <- c("L", "M", "H")

# 2. Define our data as factors
survey <- parse_factor(surv$reply, levels = f_levels)


## ----fig.cap="The standard plot function on a factored object with some values NA (last block without label).\\label{fig:factorCusSat}"----
summary(survey)
plot(survey, col="khaki3",
     main = "Customer Satisfaction",
     xlab = "Response to the last survey"
     )


## -----------------------------------------------------------------------------
surv2 <- parse_factor(surv$reply, levels = unique(surv$reply))

# Note that the labels are in order of first occurrence
summary(surv2)


## -----------------------------------------------------------------------------
# Count the labels:
fct_count(survey)


## ----'factorCusSat2', fig.cap="Maybe you would prefer to show this plot to the board meeting? This plot takes the two best categories together and creates the impression that more people are happy. Compare this to previous plot."----
# Relabel factors with fct_relabel
HML <- function (x = NULL) {
  x[x == "L"] <- "Low"
  x[x == "M"] <- "Medium/High"
  x[x == "H"] <- "Medium/High"
  x[!(x %in% c("High", "Medium/High", "Low"))] <- NA
  return(x)
  }
f <- fct_relabel(survey, HML)
summary(f)
plot(f, col="khaki3",
     main = "Only one third of customers is not happy",
     xlab = "Response to the expensive survey"
     )


## -----------------------------------------------------------------------------
HMLregex <- function (x = NULL) {
  x[grepl("^L$", x)] <- "Low"
  x[grepl("^M$", x)] <- "Medium/High"
  x[grepl("^H$", x)] <- "Medium/High"
  x[!(x %in% c("High", "Medium/High", "Low"))] <- NA
  return(x)
  }
# This would do exactly the same, but it is a powerful 
# tool with many other possibilities.


## ----fig.cap="A visualisation of how the age of customers impacted the satisfaction in our made-up example.\\label{fig:factorBar}",out.width="\\textwidth"----
num_obs <- 1000  # the number of observations in the survey
# Start from a new survey: srv
srv <- tibble(reply = 1:num_obs)
srv$age <- rnorm(num_obs, mean=50,sd=20)
srv$age[srv$age < 15] <- NA
srv$age[srv$age > 85] <- NA

hml <- function (x = 0) {
  if (x < 0)  return(NA)
  if (x <= 4) return("L")
  if (x <= 8) return("M")
  if (x <= 12) return("H")
  return(NA)
  }

for (n in 1:num_obs) {
  if (!is.na(srv$age[n])) {
     srv$reply[n] <- hml(rnorm(n = 1, mean = srv$age[n] / 7, sd = 2))
   }
   else {
     srv$reply[n] <- hml(runif(n = 1, min = 1, max = 12))
   }
}
f_levels <- c("L", "M", "H")
srv$fct <- parse_factor(srv$reply, levels = f_levels)

# from most frequent to least frequent:
srv$fct                    %>%
fct_infreq(ordered = TRUE) %>%
  levels()
  
# From least frequent to more frequent:
srv$fct      %>%
 fct_infreq  %>%
  fct_rev    %>% 
  levels
  
# Reorder the reply variable in function of median age
fct_reorder(srv$reply, srv$age) %>% 
   levels
   
# Add the function min() to order based on the minimum
# age in each group (instead of default median):
fct_reorder(srv$reply, srv$age, min) %>% 
   levels

# Show the means per class of satisfaction in base-R style:
by(srv$age, srv$fct, mean, na.rm = TRUE)

# Much more accessible result with the dplyr:
satisf <- srv            %>%
          group_by(fct)  %>%
	  summarize(
	     age = median(age, na.rm = TRUE),
	     n = n()
	     )           %>%
         print


# Show the impact of age on satisfaction visually:
par(mfrow=c(1,2))
barplot(satisf$age,  horiz=TRUE, names.arg = satisf$fct,
        col=c("khaki3","khaki3","khaki3","red"), 
	main="Median age per group")
barplot(satisf$n,  horiz=TRUE, names.arg = satisf$fct,
        col=c("khaki3","khaki3","khaki3","red"), 
	main="Frequency per group")


## -----------------------------------------------------------------------------
srv                                 %>%
 mutate("fct_ano" = fct_anon(fct))  %>%
 print


## -----------------------------------------------------------------------------
set.seed(1890)
d1       <- d0  <- iris
i        <- sample(1:nrow(d0), round(0.20 * nrow(d0)))
d1[i,1]  <- NA
i        <- sample(1:nrow(d0), round(0.30 * nrow(d0)))
d1[i,2]  <- NA
head(d1, n=10L)


## ----fig.cap='The visualization of missing data with the function aggr(). The amber colour symbolizes the missing values. Left is a histogram of missing values, and right a view on how correlated missing values are (e.g. missing both sepal width and sepal length only occurs in 8\\% of the values).\\label{fig:iris1}',message=FALSE,warning=FALSE----
# install.packages("VIM") 
library(VIM)
aggr(d1, col=c('darkolivegreen2','goldenrod2'), numbers = TRUE,
     sortVars = TRUE, labels = names(d1), cex.axis = .7, gap = 3,
     ylab=c("Missing data","Pattern"))


## ----fig.cap='The visualization of missing data with the function {\tt md.pattern()}.'----
#install.packages('mice')  # uncomment if necessary
library(mice)              # load the package

# mice provides the improved visualization function md.pattern()
md.pattern(d1)  # function provided by mice


## ----message=FALSE,eval=FALSE-------------------------------------------------
## d2_imp <- mice(d1, m = 5, maxit = 25, method = 'pmm', seed = 1500)


## ----message=FALSE,eval=FALSE-------------------------------------------------
## # e.g. choose set number 3:
## d3_complete <- complete(d2_imp, 3)


## ----message=FALSE------------------------------------------------------------
# install.packages('missForest') # only first time
library(missForest)              # load the library
d_mf <- missForest(d1)           # using the same data as before

# access the imputed data in the ximp attribute:
head(d_mf$ximp)

# normalized MSE of imputation:
d_mf$OOBerror


## -----------------------------------------------------------------------------
# First, time, install the package first via:
# install.packages('Hmisc')
library(Hmisc)

# impute using mean:
SepLImp_mean <- with(d1, impute(Sepal.Length, mean))

# impute a randomly chosen value:
SepLImp_rand <- with(d1, impute(Sepal.Length, 'random'))

# impute the maximum value:
SepLImp_max <- with(d1, impute(Sepal.Length, max))

# impute the minimum value:
SepLImp_min <- with(d1, impute(Sepal.Length, min))

# note the '*' next to the imputed values"
head(SepLImp_min, n = 10L)


## -----------------------------------------------------------------------------
aregImp <- aregImpute(~ Sepal.Length + Sepal.Width 
                        + Petal.Length + Petal.Width + Species, 
			data = d1, n.impute = 4)

print(aregImp)

# n.impute = 4 produced 4 sets of imputed values
# Access the second imputed data set as follows:
head(aregImp$imputed$Sepal.Length[,2], n = 10L)


## ----fig.cap="Two histograms of the same dataset. The histogram with less bins (right) is easier to read and reveals more clearly the underlying model (Gaussian distribution). Using more bins (as in the left plot) can overfit the dataset and obscure the true distribution.\\label{fig:twohists}"----
set.seed(1890)
d <- rnorm(90)
par(mfrow=c(1,2))
hist(d, breaks=70, col="khaki3")
hist(d, breaks=12, col="khaki3")


## -----------------------------------------------------------------------------
# Try a possible cut
c <- cut(d, breaks = c(-3, -1, 0, 1, 2, 3))
table(c)

# This is not good, it will not make solid predictions for the 
# last bind for example
c <- cut(d, breaks = c(-3, -0.5, 0.5, 3))
table(c)

# We have now a similar number of observations in each bin.
# Is that the only thing to think about?


## -----------------------------------------------------------------------------
# install.packages('binr) # do once
library(binr)
b <- bins.quantiles(d, target.bins=5, max.breaks=10)
b$binct


## ----fig.cap="A plot of the fabricated dataset with the spending ratio in function of the age of the customers. The spending ratio is defined as $\\frac{S_{n}}{S_{n-1} + S_n}$, where $S_n$ is the spending in period $n$. If both spends are 0, then the spending ratio is defined as 0."----
set.seed(1890)
age <- rlnorm(1000, meanlog = log(40), sdlog = log(1.3))
y <- rep(NA, length(age))
for(n in 1:length(age)) {
  y[n] <- max(0, 
              dnorm(age[n], mean= 40, sd=10) 
	         + rnorm(1, mean = 0, sd = 10 * dnorm(age[n], 
		   mean= 40, sd=15)) * 0.075)
}
y <- y / max(y)
plot(age, y,
     pch = 21, col = "blue", bg = "red",
     xlab = "age",
     ylab = "spending ratio"
     )

# Assume this data is:
#   age            = age of customer
#   spending_ratio = R : = S_n/ (S_{n-1} + S_n)
#                       (zero if both are zero)
#      with S_n the spending in month n
dt <- tibble (age = age, spending_ratio = y)


## ----fig.cap='A simple aid to select binning borders is plotting a non-parametric fit (left) and the histogram (right). The information from both plots combined can be used to decide on binning.\\label{fig:binvis}'----
# in general we can leave out NAs (in this example redundant)
d1 <- dt[complete.cases(dt),]  

# order() returns sorted indices, so this orders the vector:
d1 <- d1[order(d1$age),]  

# Fit a loess:
d1_loess <- loess(spending_ratio ~ age, d1)

# Add predictions:
d1_pred_loess <- predict(d1_loess)

# Plot the results:
par(mfrow=c(1,2))
plot(d1$age, d1$spending_ratio, pch=16,
     xlab = 'age', ylab = 'spending ratio')
lines(d1$age, d1_pred_loess, lwd = 7, col = 'dodgerblue4')
hist(d1$age, col = 'dodgerblue4', xlab = 'age')
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
# Fit the logistic regression directly on the data without binning:
lReg1 <- glm(formula = spending_ratio ~ age, 
             family = quasibinomial, 
             data = dt)

# Investigate the model:
summary(lReg1)

# Calculate predictions and means square error:
pred1 <- 1 / (1+ exp(-(coef(lReg1)[1] + dt$age * coef(lReg1)[2])))
SE1  <-  (pred1 - dt$spending_ratio)^2
MSE1 <- sum(SE1) / length(SE1)


## -----------------------------------------------------------------------------
# Bin the variable age:
c <- cut(dt$age, breaks = c(15, 30, 55, 90))

# Check the binning:
table(c)
# We have one big bucket and two smaller (with the smallest
# more than 10% of our dataset.

lvls <- unique(c)      # find levels
lvls                   # check levels order

# Create the tibble (a data-frame also works):
dt <- as_tibble(dt)                            %>%
      mutate(is_L = if_else(age <= 30, 1, 0))  %>%
      mutate(is_H = if_else(age > 55 , 1, 0))

# Fit the logistic regression with is_L and is_H:
# (is_M is not used because it is correlated with the previous)
lReg2 <- glm(formula = spending_ratio ~ is_L + is_H, 
             family = quasibinomial, data = dt)

# Investigate the logistic model:
summary(lReg2)

# Calculate predictions for our model and calculate MSE:
pred2 <- 1 / (1+ exp(-(coef(lReg2)[1] + dt$is_L * coef(lReg2)[2] 
                     + dt$is_H * coef(lReg2)[3])))
SE2 <-  (pred2 - dt$spending_ratio)^2
MSE2 <- sum(SE2) / length(SE2)

# Compare the MSE of the two models:
MSE1
MSE2


## ----fig.cap=c("The underlying relation between spending probability for females (left) and males (right) in our fabricated example.\\label{fig:matrixBinning1}")----
# Load libraries and define parameters:
library(tidyverse) # provides tibble (only used in next block)
set.seed(1880)     # to make results reproducible
N <- 500           # number of rows

# Ladies first:
# age will function as our x-value:
age_f   <- rlnorm(N, meanlog = log(40), sdlog = log(1.3))
# x is a temporary variable that will become the propensity to buy:
x_f <- abs(age_f + rnorm(N, 0, 20))    # Add noise & keep positive
x_f <- 1 - (x_f - min(x_f)) / max(x_f) # Scale between 0 and 1
x_f <- 0.5 * x_f / mean(x_f)           # Coerce mean to 0.5
# This last step will produce some outliers above 1
x_f[x_f > 1] <- 1   # Coerce those few that are too big to 1

# Then the gentlemen:
age_m   <- rlnorm(N, meanlog = log(40), sdlog = log(1.3))
x_m <- abs(age_m + rnorm(N, 0, 20))    # Add noise & keep positive
x_m <- 1 - (x_m - min(x_m)) / max(x_m) # Scale between 0 and 1
x_m <- 0.5 * x_m / mean(x_m)           # Coerce mean to 0.5
# This last step will produce some outliers above 1
x_m[x_m > 1] <- 1   # Coerce those few that are too big to 1
x_m <- 1 - x_m                         # relation to be increasing

# Rename (p_x is not the gendered propensity to buy)
p_f <- x_f
p_m <- x_m

# We want a double plot, so change plot params & save old values:
oldparams <- par(mfrow=c(1,2)) 
plot(age_f, p_f,
     pch = 21, col = "blue", bg = "red",
     xlab = "Age",
     ylab = "Spending probability",
     main = "Females"
    )
plot(age_m, p_m,
     pch = 21, col = "blue", bg = "red",
     xlab = "Age",
     ylab = "Spending probability",
     main = "Males"
    )
par(oldparams)   # Reset the plot parameters after plotting


## ----fig.cap=c("The dataset ``as received from the customer service department'' does not show any clear relationship between Age and Sex with the variable that we want to explain: the spending ratio.\\label{fig:MatrixBinning2}")----
# Now, we merge the data and consider this as our input-data
tf <- tibble("age" = age_f, "sex" = "F", "is_good" = p_f) 
tm <- tibble("age" = age_m, "sex" = "M", "is_good" = p_m)
t  <- full_join(tf, tm, by = c("age", "sex", "is_good"))

# Change plot parameters and capture old values:
oldparams <- par(mfrow=c(1,2))  
plot(t$age, t$is_good,
     pch  = 21, col = "black", bg = "khaki3",
     xlab = "Age",
     ylab = "Spending probability",
     main = "Dependence on age"
     )
fct_sex <- factor(t$sex, levels=c("F","M"), labels=c(0,1))
t$sexM  <- as.numeric(fct_sex)    # store for later use
plot(fct_sex, t$is_good,
     col="khaki3",
     main="Dependence on sex",
     xlab="Female       Male")
par(oldparams)   # Reset the plot parameters


## ----fig.cap='The data does not reveal much patterns for any of the variables (Gender and Age).\\label{fig:BinningMatrix3}',warning=FALSE----
d1 <- t[complete.cases(t),]  

d1 <- d1[order(d1$age),]  
d1_age_loess <- loess(is_good ~ age, d1)
d1_age_pred_loess <- predict(d1_age_loess)

d1 <- d1[order(d1$sexM),]  
d1_sex_loess <- loess(is_good ~ sexM, d1)
d1_sex_pred_loess <- predict(d1_sex_loess)

# Plot the results:
par(mfrow=c(2,2))
d1 <- d1[order(d1$age),]  
plot(d1$age, d1$is_good, pch=16,
     xlab = 'Age', ylab = 'Spending probability')
lines(d1$age, d1_age_pred_loess, lwd = 7, col = 'dodgerblue4')
hist(d1$age, col = 'khaki3', xlab = 'age')

d1 <- d1[order(d1$sexM),]  
plot(d1$sexM, d1$is_good, pch=16,
     xlab = 'Gender', ylab = 'Spending probability')
lines(d1$sexM, d1_sex_pred_loess, lwd = 7, col = 'dodgerblue4')
hist(d1$sexM, col = 'khaki3', xlab = 'gender')
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
# Note that we can feed "sex" into the model and it will create
# for us a variable "sexM" (meaning the same as ours)
# To avoid this confusion, we put in our own variable.
regr1 <- glm(formula = is_good ~ age + sexM, 
             family = quasibinomial,
             data = t)

# assess the model:
summary(regr1)

pred1 <- 1 / (1+ exp(-(coef(regr1)[1] + t$age * coef(regr1)[2] 
                     + t$sexM * coef(regr1)[3])))
SE1 <-  (pred1 - t$is_good)^2
MSE1 <- sum(SE1) / length(SE1)


## -----------------------------------------------------------------------------
# 1. Check the potential cut 
c <- cut(t$age, breaks = c(min(t$age), 35, 55, max(t$age)))
table(c)

# 2. Create the matrix variables
t <- as_tibble(t)                                               %>%
    mutate(is_LF = if_else((age <= 35) & (sex == "F"), 1L, 0L)) %>%
    mutate(is_HF = if_else((age >  50) & (sex == "F"), 1L, 0L)) %>%
    mutate(is_LM = if_else((age <= 35) & (sex == "M"), 1L, 0L)) %>%
    mutate(is_HM = if_else((age >  50) & (sex == "M"), 1L, 0L)) %>%
    print

# 3. Check if the final bins aren't too small
t[,5:8] %>% map_int(sum)


## -----------------------------------------------------------------------------
regr2 <- glm(formula = is_good ~ is_LF + is_HF + is_LM + is_HM,
             family = quasibinomial,
             data = t)

# assess the model:
summary(regr2)

pred2 <- 1 / (1+ exp(-(coef(regr2)[1] + 
                     + t$is_LF * coef(regr2)[2] 
		     + t$is_HF * coef(regr2)[3] 
                     + t$is_LM * coef(regr2)[4] 
		     + t$is_HM * coef(regr2)[5] 
                      )))
SE2 <-  (pred2 - t$is_good)^2
MSE2 <- sum(SE2) / length(SE2)


## -----------------------------------------------------------------------------
MSE1
MSE2


## -----------------------------------------------------------------------------
t <- mutate(t, "is_good" = if_else(is_good >= 0.5, 1L, 0L))


## -----------------------------------------------------------------------------
# so we start from this dataset used in previous section:
print(t)


## -----------------------------------------------------------------------------
#install.packages("InformationValue") 
library(InformationValue)

WOETable(X = factor(t$sexM), Y = t$is_good, valueOfGood=1) %>%
   knitr::kable(format.args = list(big.mark = " ", digits=2))

## also functions WOE() and IV(), e.g.
# IV of a categorical variable is the sum of IV of its categories
IV(X = factor(t$sexM), Y = t$is_good, valueOfGood=1)


## -----------------------------------------------------------------------------
WOETable(X = factor(t$is_LF), Y = t$is_good, valueOfGood=1) %>%
   knitr::kable(digits=2)

## also functions WOE() and IV(), e.g.
# IV of a categorical variable is the sum of IV of its categories
IV(X = factor(t$is_LF), Y = t$is_good, valueOfGood=1)


## ----fig.cap=c('A visualization of the loadings of the principal components of the dataset mtcars.\\label{fig:PCA1}','The biplot of the dataset mtcars: all observation and dimensions projected in the plane span by the first two principal components.\\label{fig:PCA2}')----
# PCA: extracting PCs from the correlation matrix
fit <- princomp(mtcars, cor=TRUE)

summary(fit)      # print the variance explained by PC
loadings(fit)     # show PC loadings
head(fit$scores)  # the first principal components 

# plot the loadings (output see figure):
plot(fit,type="b", col='khaki3') 

# show the biplot:
biplot(fit)       


## -----------------------------------------------------------------------------
# Maximum Likelihood Factor Analysis

# Extracting 3 factors with varimax rotation:
fit <- factanal(mtcars, 3, rotation = "varimax")
print(fit, digits = 2, cutoff = .3, sort = TRUE)

# plot factor 1 by factor 2
load <- fit$loadings[,1:2]
plot(load, type = "n")                     # plot the loads
text(load, labels = colnames(mtcars), 
     cex = 1.75, col = 'blue')             # add variable names 


## ----echo=TRUE,results='hide',message=FALSE,warning=FALSE---------------------
# load the library nFactors:
library(nFactors)


## ----fig.cap='Visual aids to select the optimal number of factors.\\label{fig:scree}',warning=FALSE----
# Get the eigenvectors:
eigV <- eigen(cor(mtcars))


# Get mean and selected quantile of the distribution of eigen-
# values of correlation or a covariance matrices of standardized
# normally distributed variables:
aPar <- parallel(subject = nrow(mtcars),var = ncol(mtcars),
                 rep = 100, cent = 0.05)
  
# Get the optimal number of factors analysis:
nScr <- nScree(x = eigV$values, aparallel = aPar$eigen$qevpea)

# See the result
nScr
# and plot it.
plotnScree(nScr) 


## ----eval=FALSE---------------------------------------------------------------
## # PCA Variable Factor Map
## library(FactoMineR)
## 
## # PCA will generate two plots as side effect (not shown here)
## result <- PCA(mtcars)


## ----echo=FALSE,include=FALSE-------------------------------------------------
rm(list=ls())   # reset R at the beginning of a new part


## ----surv,fig.cap=c("A scatter-plot generated by the line ``plot(survey\\$Height, survey\\$Wr.Hnd).''","A lot visualizing the linear regression model (the data in red and the regression in blue.","Using the function abline() and cleaning up the titles.")----
library(MASS)

# Explore the data
plot(survey$Height, survey$Wr.Hnd)


## ----eval=FALSE---------------------------------------------------------------
## < dependent variable> tilde <sum of independent variables>


## -----------------------------------------------------------------------------
# Create the model
lm1 <- lm (formula = Wr.Hnd ~ Height, data = survey)
summary(lm1)


## ----plotlm,fig.cap=c("A plot visualizing the linear regression model (the data in red and the regression in blue).")----
# predictions
h <- data.frame(Height = 150:200)
Wr.lm <- predict(lm1, h)
plot(survey$Height, survey$Wr.Hnd,col="red")
lines(t(h),Wr.lm,col="blue",lwd=3)


## ----abline1,fig.cap=c("Using the function abline() and cleaning up the titles.")----
# Or use the function abline()
plot(survey$Height, survey$Wr.Hnd,col = "red",
     main = "Hand span in function of Height", 
     abline(lm(survey$Wr.Hnd ~ survey$Height ),
            col='blue',lwd=3),
     cex = 1.3,pch = 16,
     xlab = "Height",ylab ="Hand span")


## -----------------------------------------------------------------------------
# We use mtcars from the library MASS
model <- lm(mpg ~ disp + hp + wt, data = mtcars)
print(model)


## -----------------------------------------------------------------------------
# Accessing the coefficients
intercept <- coef(model)[1]
a_disp    <- coef(model)[2]
a_hp      <- coef(model)[3]
a_wt      <- coef(model)[4]

paste('MPG =', intercept, '+', a_disp, 'x disp +', 
     a_hp,'x hp +', a_wt, 'x wt')


## -----------------------------------------------------------------------------
# This allows us to manually predict the fuel consumption
# e.g. for the Mazda Rx4
2.23 + a_disp * 160 + a_hp * 110 + a_wt * 2.62


## ----eval=FALSE---------------------------------------------------------------
## glm(formula, data, family)


## -----------------------------------------------------------------------------
m <- glm(cyl ~ hp + wt, data = mtcars, family = "poisson")
summary(m)


## -----------------------------------------------------------------------------
  m <- glm(cyl ~ hp, data = mtcars, family = "poisson")
  summary(m)


## ----nls1,fig.cap=c("The results of the non-linear regression with nls(). This plot indicates that there is one outlier and you might want to rerun the model without this observation.")----
# Consider observations for dt = d0 + v0 t + 1/2 a t^2
t  <- c(1,2,3,4,5,1.5,2.5,3.5,4.5,1)
dt <- c(8.1,24.9,52,89.2,136.1,15.0,37.0,60.0,111.0,8)

# Plot these values.
plot(t,dt,xlab="time",ylab="distance")

# Take the assumed values and fit into the model.
model <- nls(dt ~ d0 + v0 * t + 1/2 * a * t^2,
             start = list(d0 = 1,v0 = 3,a = 10))

# Plot the model curve
simulation.data <- data.frame(t = seq(min(t),max(t),len = 100))
lines(simulation.data$t,predict(model,
      newdata = simulation.data), col="red", lwd = 3)


## -----------------------------------------------------------------------------
# Learn about the model
summary(model)                # the summary
print(sum(residuals(model)^2))# squared sum of residuals
print(confint(model))         # confidence intervals


## -----------------------------------------------------------------------------
m <- lm(data = mtcars, formula = mpg ~ wt)
summary(m)
summary(m)$r.squared


## -----------------------------------------------------------------------------
# Consider the relation between the hours studied and passing
# an exam (1) or failing it (0):
hours <- c(0,0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 
           1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25,
	   3.50, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50)
pass  <- c(0,0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 
           1, 0, 1, 1, 1, 1, 1, 1)
d <- data.frame(cbind(hours,pass))
m <- glm(formula=pass ~ hours, family = binomial, 
         data = d)


## ----fig.cap=c("The grey diamonds with red border are the data-points (not passed is 0 and passed is 1) and the blue line represents the logistic regression model (or the probability to succeed the exam in function of the hours studied.\\label{fig:logistic.reg}")----
# Visualize the results:
plot(hours, pass, col = "red", pch = 23, bg = "grey",
     xlab = 'Hours studied',
     ylab = 'Passed exam (1) or not (0)')
pred <- 1 / (1+ exp(-(coef(m)[1] + hours * coef(m)[2])))
lines(hours,pred,col="blue",lwd=3)


## -----------------------------------------------------------------------------
# if necessary: install.packages('titanic')
library(titanic)

# This provides a.o. two datasets titanic_train and titanic_test.
# We will work further with the training-dataset.
t <- titanic_train
colnames(t)


## -----------------------------------------------------------------------------
# fit provide a simple model
m <- glm(data    = t, 
         formula = Survived ~ Pclass + Sex + Pclass * Sex + 
	           Age + SibSp, 
	 family  = binomial)
summary(m)


## -----------------------------------------------------------------------------
# building further on the model m
# scores between 0 and 1 (odds)
t2 <- t[complete.cases(t),]
predicScore <- predict(object=m,type="response", newdat = t2)

# Introduce a cut-off level above which we assume survival
predic <- ifelse(predicScore > 0.7, 1, 0)

# The confusion matrix is one line, the headings 2
confusion_matrix <- table(predic, t2$Survived)
rownames(confusion_matrix) <- c("predicted_death",
                                "predicted_survival")
colnames(confusion_matrix) <- c("observed_death",
                                "observed_survival")
print(confusion_matrix)


## ----ROC1,fig.cap='The ROC curve of a logistic regression.',warning=FALSE,message=FALSE----
library(ROCR)
# Re-use the model m and the dataset t2
pred <- prediction(predict(m, type = "response"), t2$Survived)

# Visualize the ROC curve
plot(performance(pred, "tpr", "fpr"), col="blue", lwd=3)
abline(0, 1, lty = 2)


## -----------------------------------------------------------------------------
S4_perf <- performance(pred, "tpr", "fpr")
df <- data.frame(
   x = S4_perf@x.values,
   y = S4_perf@y.values,
   a = S4_perf@alpha.values
   )
colnames(df) <- c(S4_perf@x.name, S4_perf@y.name, S4_perf@alpha.name)
head(df)


## ----ROCgg,fig.cap='The ROC curve plotted with ggplot2.'----------------------
library(ggplot2)
p <- ggplot(data=df, 
            aes(x = `False positive rate`, y = `True positive rate`)) +
            geom_line(lwd=2, col='blue')  + 
            # The next lines add the shading:
            aes(x = `False positive rate`, ymin = 0, 
                ymax = `True positive rate`) + 
            geom_ribbon(, alpha=.5)
p


## ----ROCacc,fig.cap='A plot of the accuracy in function of the cut-off (threshold) level.'----
# Plotting the accuracy (in function of the cut-off)
plot(performance(pred, "acc"), col="blue", lwd=3)


## ----AUCdemo,fig.cap='The area under the curve (AUC) is the area A plus the area B.'----
# Assuming that we have the predictions in the prediction object:
plot(performance(pred, "tpr", "fpr"), col="blue", lwd=3)
abline(0,1,lty=2)
text(0.3,0.5,"A")
text(0.1,0.9,"B")
text(0.8,0.3,"C")


## -----------------------------------------------------------------------------
AUC <- attr(performance(pred, "auc"), "y.values")[[1]]
AUC


## ----fig.cap='The Gini coefficient of the model is area A divided by the sum of A and B.'----
plot(performance(pred, "tpr", "fpr"), col="blue", lwd=3)
abline(0,1,lty=2)
text(0.3,0.5,"A")
text(0.1,0.9,"B")
text(0.8,0.3,"C")


## -----------------------------------------------------------------------------
paste("the Gini is:",round(2 * AUC - 1, 2))


## ----KSvis,fig.cap='The KS as the maximum distance between the cumulative distributions of the positive and negative observations.',echo=FALSE----
# using model m and data frame t2:
predicScore <- predict(object=m,type="response")
d0 <- data.frame(
       score = as.vector(predicScore)[t2$Survived == 0],
       true_result = 'not survived')
d1 <- data.frame(
       score = as.vector(predicScore)[t2$Survived == 1],
       true_result = 'survived')
d <- rbind(d0, d1)
d <- d[complete.cases(d),]

cumDf0 <- ecdf(d0$score)
cumDf1 <- ecdf(d1$score)
x <- sort(d$score)
cumD0 <- cumDf0(x)
cumD1 <- cumDf1(x)
diff <- cumD0 - cumD1
y1  <- gdata::first(cumD0[diff == max(diff)])
y2  <- gdata::first(cumD1[diff == max(diff)])
x1  <- x2 <- quantile(d0$score, probs=y1, na.rm=TRUE)

# plot this with ggplot2
p <- ggplot(d, aes(x = score)) +
     stat_ecdf(geom = "step", aes(col = true_result), lwd=2) +
     ggtitle('Cummulative distributions and KS') +
     geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                  color='navy', lwd=3) + 
     ggplot2::annotate("text", 
              label = paste0("KS=",round((y1-y2)*100,2),"%"), 
              x = x1 + 0.15, y = y2+(y1-y2)/2, color = "navy")
p


## -----------------------------------------------------------------------------
pred <- prediction(predict(m,type="response"), t2$Survived)
ks.test(attr(pred,"predictions")[[1]], 
        t2$Survived,
        alternative = 'greater')


## -----------------------------------------------------------------------------
perf <- performance(pred, "tpr", "fpr")
ks <- max(attr(perf,'y.values')[[1]] - attr(perf,'x.values')[[1]])
ks

# Note: the following line yields the same outcome
ks <- max(perf@y.values[[1]] - perf@x.values[[1]])
ks


## ----KSvisROC,fig.cap='The KS as the maximum distance between the model and a pure random model.'----
pred <- prediction(predict(m,type="response"), t2$Survived)
perf <- performance(pred, "tpr", "fpr")
plot(perf, main = paste0(' KS is',round(ks*100,1),'%'),
     lwd = 4, col = 'red')
lines(x = c(0,1),y=c(0,1), col='blue')

# The KS line:
diff <- perf@y.values[[1]] - perf@x.values[[1]]
xVal <- attr(perf,'x.values')[[1]][diff == max(diff)]
yVal <- attr(perf,'y.values')[[1]][diff == max(diff)]
lines(x = c(xVal, xVal), y = c(xVal, yVal), 
      col = 'khaki4', lwd=8)


## -----------------------------------------------------------------------------
# First, we need a logistic regression. We use the same as before:
library(titanic)
t <- titanic_train
t2 <- t[complete.cases(t),]
m <- glm(data    = t, 
         formula = Survived ~ Pclass + Sex + Pclass*Sex + Age + SibSp, 
         family  = binomial)
library(ROCR)
pred <- prediction(predict(m,type="response", newdat = t2), 
                   t2$Survived)
perf <- performance(pred, "tpr", "fpr")


## -----------------------------------------------------------------------------
# get_best_cutoff
# Finds a cutof for the score so that sensitivity and specificity 
# are optimal.
# Arguments
#   fpr    -- numeric vector -- false positive rate
#   tpr    -- numeric vector -- true positive rate
#   cutoff -- numeric vector -- the associated cutoff values
# Returns:
#   the cutoff value (numeric)
get_best_cutoff <- function(fpr, tpr, cutoff){
        cst <- (fpr - 0)^2 + (tpr - 1)^2
        idx = which(cst == min(cst))
        c(sensitivity = tpr[[idx]], 
          specificity = 1 - fpr[[idx]], 
          cutoff = cutoff[[idx]])
    }
    
# opt_cut_off 
# Wrapper for get_best_cutoff. Finds a cutof for the score so that 
# sensitivity and specificity are optimal.
# Arguments:
#    perf -- performance object (ROCR package)
#    pred -- prediction object (ROCR package)
# Returns:
#   The optimal cutoff value (numeric)
opt_cut_off = function(perf, pred){
    mapply(FUN=get_best_cutoff, 
           perf@x.values, 
           perf@y.values, 
           pred@cutoffs)
   }
opt_cut_off(perf, pred)


## -----------------------------------------------------------------------------
# We introduce cost.fp to be understood as a the cost of a 
# false positive, expressed as a multiple of the cost of a 
# false negative.

# get_best_cutoff
# Finds a cutof for the score so that sensitivity and specificity 
# are optimal.
# Arguments
#   fpr     -- numeric vector -- false positive rate
#   tpr     -- numeric vector -- true positive rate
#   cutoff  -- numeric vector -- the associated cutoff values
#   cost.fp -- numeric        -- cost of false positive divided 
#                                by the cost of a false negative
#				 (default = 1)
# Returns:
#   the cutoff value (numeric)
get_best_cutoff <- function(fpr, tpr, cutoff, cost.fp = 1){
        cst <- (cost.fp * fpr - 0)^2 + (tpr - 1)^2
        idx = which(cst == min(cst))
        c(sensitivity = tpr[[idx]], 
          specificity = 1 - fpr[[idx]], 
          cutoff = cutoff[[idx]])
    }
    
# opt_cut_off 
# Wrapper for get_best_cutoff. Finds a cutof for the score so that 
# sensitivity and specificity are optimal.
# Arguments:
#    perf    -- performance object (ROCR package)
#    pred    -- prediction object (ROCR package)
#    cost.fp -- numeric -- cost of false positive divided by the
#                          cost of a false negative (default = 1)
# Returns:
#   The optimal cutoff value (numeric)
opt_cut_off = function(perf, pred, cost.fp = 1){
    mapply(FUN=get_best_cutoff, 
           perf@x.values, 
           perf@y.values, 
           pred@cutoffs,
	   cost.fp)
   }
opt_cut_off(perf, pred, cost.fp = 5)


## -----------------------------------------------------------------------------
# e.g. cost.fp = 1 x cost.fn
perf_cst1 <- performance(pred, "cost", cost.fp = 1)
str(perf_cst1) # the cost is in the y-values

# the optimal cut-off is then the same as in previous code sample
pred@cutoffs[[1]][which.min(perf_cst1@y.values[[1]])]

# e.g. cost.fp = 5 x cost.fn
perf_cst2 <- performance(pred, "cost", cost.fp = 5)

# the optimal cut-off is now:
pred@cutoffs[[1]][which.min(perf_cst2@y.values[[1]])]


## ----cstCutoff,fig.height='0.9\\textheight',out.height='0.9\\textheight', fig.cap='The cost functions compared different cost structures. In plot (a), we plotted the cost function when the cost of a false positive is equal to the cost of a false negative. In plot (b), a false positive costs five times more than a false negative (valid for a loan in a bank).'----
par(mfrow=c(2,1))
plot(perf_cst1, lwd=2, col='navy', main='(a) cost(FP) = cost(FN)')
plot(perf_cst2, lwd=2, col='navy', main='(b) cost(FP) = 5 x cost(FN)')
par(mfrow=c(1,1))


## ----fig.cap="Three alternatives for the impurity measure in the case of classification problems.\\label{fig:impurity}"----
# Draw Gini, deviance, and misclassification functions

# Define the functions:
gini <- function(x) 2 * x * (1-x)
entr <- function(x) (-x*log(x) - (1-x)*log(1-x))/log(2) / 2
misc <- function(x) {1 - pmax(x,1-x)}

# Plot the curves:
curve(gini, 0, 1, ylim = c(0,0.5), col = "forestgreen", 
      xlab="p", ylab = "Impurity measure", type = "l", lwd = 6)
curve(entr, 0, 1, add = TRUE, lw = 3, col = "black", lwd = 6)
curve(misc, 0, 1, add = TRUE, lw = 3, col = "blue", type = "l", 
      lty = 2, lwd = 6)

# Add the text:
text(0.85, 0.4,   "Gini index",                col = "forestgreen")
text(0.85, 0.485, "Deviance or cross-entropy", col = "black")
text(0.5,  0.3,   "Misclassification index",   col = "blue")


## ----eval=FALSE---------------------------------------------------------------
## rpart(formula, data, weights, subset, na.action = na.rpart,
##       method=c(``class'',``anova''), model = FALSE,
##       x = FALSE, y = TRUE, parms, control, cost, ...)


## ----eval=FALSE---------------------------------------------------------------
## parms = list(prior = c(0.6,0.4))


## ----rpartTree,fig.cap=c("The plot of the complexity parameter (cp) via the function {\\tt plotcp()}.","The decision tree, fitted by rpart. This figure helps to visualize what happens in the decision tree that predicts survival in the Titanic disaster.\\label{fig:tree.rpart}","The same tree as in Figure \\ref{fig:tree.rpart} but now pruned with a complexity parameter $\\rho = 0.01$. Note that the tree is .")----
## example of a regression tree with rpart on the dataset of the Titanic
##
library(rpart)
titanic <- read.csv("data/titanic3.csv")
frm     <- survived ~ pclass + sex + sibsp + parch + embarked + age
t0      <- rpart(frm, data=titanic, na.action = na.rpart, 
  method="class",
  parms = list(prior = c(0.6,0.4)),
  #weights=c(...), # each observation (row) can be weighted
  control = rpart.control(
  minsplit       = 50,  # minimum nbr. of observations required for split
  minbucket      = 20,  # minimum nbr. of observations in a terminal node
  cp             = 0.001,# complexity parameter set to a small value 
                        # this will grow a large (over-fit) tree
  maxcompete     = 4,   # nbr. of competitor splits retained in output
  maxsurrogate   = 5,   # nbr. of surrogate splits retained in output
  usesurrogate   = 2,   # how to use surrogates in the splitting process
  xval           = 7,   # nbr. of cross validations
  surrogatestyle = 0,   # controls the selection of a best surrogate
  maxdepth       = 6)   # maximum depth of any node of the final tree
  )

# Show details about the tree t0:
printcp(t0)             

# Plot the error in function of the complexity parameter
plotcp(t0)              

# print(t0) # to avoid too long output we commented this out
# summary(t0)

# Plot the original decisions tree
plot(t0)


text(t0)

# Prune the tree:
t1 <- prune(t0, cp=0.01)
plot(t1); text(t1)


## ----rpartPlot,fig.cap=c("The decision tree represented by the function {\\tt prp()} from the package {\\tt rpart.plot}. This plot not only looks more elegant, but it is also more informative and less simplified. For example the top node ``sex'' has now two clear options in which descriptions we can recognize the words male and female, and the words are on the branches, so there is no confustion possible which is left and which right.")----
# plot the tree with rpart.plot
library(rpart.plot)
prp(t0, type = 5, extra = 8, box.palette = "auto",
    yesno = 1, yes.text="survived",no.text="dead"
    )


## ----RegTree,fig.cap=c("The plot of the complexity parameter (cp) via the function {\\tt plotcp()}","rpart tree on mpg for the dataset mtcars.\\label{fig:tree.rpartX}","The same tree as in Figure \\ref{fig:tree.rpartX} but now pruned with a complexity parameter $\\rho$ of 0.1. The regression tree is -- in this example -- too simple.")----
# Example of a regression tree with rpart on the dataset mtcars

# The libraries should be loaded by now:
library(rpart); library(MASS); library (rpart.plot)

# Fit the tree:
t <- rpart(mpg ~ cyl + disp + hp + drat + wt + qsec + am + gear,
 data=mtcars, na.action = na.rpart, 
 method     = "anova",
 control    = rpart.control(
   minsplit       = 10,  # minimum nbr. of observations required for split
   minbucket      = 20/3,# minimum nbr. of observations in a terminal node
                         # the default = minsplit/3
   cp             = 0.01,# complexity parameter set to a very small value 
                         # his will grow a large (over-fit) tree
   maxcompete     = 4,   # nbr. of competitor splits retained in output
   maxsurrogate   = 5,   # nbr. of surrogate splits retained in output
   usesurrogate   = 2,   # how to use surrogates in the splitting process
   xval           = 7,   # nbr. of cross validations
   surrogatestyle = 0,   # controls the selection of a best surrogate
   maxdepth       = 30   # maximum depth of any node of the final tree
   )
 )
 
# Investigate the complexity parameter dependence:
printcp(t)
plotcp(t)

# Print the tree:
print(t)
summary(t)

# plot(t) ; text(t)  # This would produce the standard plot from rpart.
# Instead we use:
prp(t, type = 5, extra = 1, box.palette = "Blues", digits = 4,
    shadow.col = 'darkgray', branch = 0.5)

# Prune the tree:
t1 <- prune(t, cp = 0.05)

# Finally, plot the pruned tree:
prp(t1, type = 5, extra = 1, box.palette = "Reds", digits = 4,
    shadow.col = 'darkgray', branch = 0.5)


## -----------------------------------------------------------------------------
# We use the function stats::predict()
predicPerc <- predict(object=t0, newdata=titanic)

# predicPerc is now a matrix with probabilities: in column 1 the 
# probability not to survive, in column 2 the probability to survive:
head(predicPerc)

# This is not what we need. We need to specify that it is a
# classification tree. Here we correct this:
predic <- predict(object=t0, newdata=titanic, type="class")
  # vector with only the fitted class as prediction
head(predic)


## -----------------------------------------------------------------------------
# The confusion matrix:
confusion_matrix <- table(predic, titanic$survived)
rownames(confusion_matrix) <- c("predicted_death",
                                "predicted_survival")
colnames(confusion_matrix) <- c("observed_death",
                                "observed_survival")
confusion_matrix

# As a precentage:
confusion_matrixPerc <- sweep(confusion_matrix, 2, 
                       margin.table(confusion_matrix,2),"/")

# Here is the confusion matrix:
round(confusion_matrixPerc,2)


## ----rocTree,fig.cap=c("ROC curve of the decision tree.", "accuracy")---------
library(ROCR)
pred <- prediction(predict(t0, type = "prob")[,2], 
                                 titanic$survived)

# Visualize the ROC curve:
plot(performance(pred, "tpr", "fpr"), col="blue", lwd=3)
abline(0,1,lty=2)


## ----fig.cap=c("The accuracy for the decision tree on the Titanic data.")-----
plot(performance(pred, "acc"), col="blue", lwd=3)
abline(1,0,lty=2)


## -----------------------------------------------------------------------------
# AUC:
AUC <- attr(performance(pred, "auc"), "y.values")[[1]]
AUC

# GINI:
2 * AUC - 1

# KS:
perf <- performance(pred, "tpr", "fpr")
max(attr(perf,'y.values')[[1]]-attr(perf,'x.values')[[1]])


## ----results='hide'-----------------------------------------------------------
library(randomForest)


## ----forestCars,fig.cap=c("The plot of a randomForest object shows how the model improves in function of the number of trees used.","The importance of each variable in the random-forest model.","Partial dependence on the variables (1 of 3).","Partial dependence on the variables (2 of 3).","Partial dependence on the variables (3 of 3).")----
head(mtcars)
mtcars$l <- NULL  # remove our variable
frm      <- mpg ~ cyl + disp + hp + drat + wt + qsec + am + gear
set.seed(1879)

# Fit the random forest:
forestCars   = randomForest(frm, data = mtcars)

# Show an overview:
print(forestCars)

# Plot the random forest overview:
plot(forestCars)

# Show the summary of fit:
summary(forestCars)

# visualization of the RF:
getTree(forestCars, 1, labelVar=TRUE)

# Show the purity of the nodes:
imp <- importance(forestCars)
imp

# This impurity overview can also be plotted:
plot( imp, lty=2, pch=16)
lines(imp)

# Below we print the partial dependence on each variable.
# We group the plots per 3, to save some space.
impvar = rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op     = par(mfrow=c(1, 3))
for (i in seq_along(impvar)) {
    partialPlot(forestCars, mtcars, impvar[i], xlab=impvar[i],
    main=paste("Partial Dependence on", impvar[i]))
  }


## ----NNlogit,echo=FALSE,fig.cap=c("A logistic regression is actually a neural network with one neuron. Each variable contributes to a sigmoid function in one node, and if that one node gets loadings over a critical threshold, then we predict 1, otherwise 0. The intercept is the ``1'' in a circle and bears a loading of 12.02."),warning=FALSE,message=FALSE----
#install.packages("neuralnet")
library(neuralnet)
nn0 <- neuralnet(mpg ~ wt + qsec + am + hp + disp + cyl + drat + gear + carb,
                 data = mtcars, hidden = c(0),
                 linear.output = TRUE)
plot(nn0, rep = "best", intercept = TRUE, information = FALSE, 
     show.weights = TRUE);


## ----nnDef,eval=FALSE---------------------------------------------------------
## neuralnet(formula, data, hidden = 1, stepmax = 1e+05
##           linear.output = TRUE)


## -----------------------------------------------------------------------------
#install.packages("neuralnet") # Do only once.

# Load the library neuralnet:
library(neuralnet)

# Fit the aNN with 2 hidden layers that have resp. 3 and 2 neurons:
# (neuralnet does not accept a formula wit a dot as in 'y ~ .' )
nn1 <- neuralnet(mpg ~ wt + qsec + am + hp + disp + cyl + drat + 
                       gear + carb,
                 data = mtcars, hidden = c(3,2),
                 linear.output = TRUE)


## ----nnPlot,fig.cap=c("A simple neural net fitted to the dataset of mtcars, predicting the miles per gallon (mpg). In this example we predict the fuel consumption of a car based on some other values in the dataset {\tt mtcars}.")----
plot(nn1, rep = "best", information = FALSE);


## -----------------------------------------------------------------------------
# Get the data about crimes in Boston:
library(MASS)
d <- Boston


## -----------------------------------------------------------------------------
# Inspect if there is missing data:
apply(d,2,function(x) sum(is.na(x)))
# There are no missing values.


## -----------------------------------------------------------------------------
set.seed(1877) # set the seed for the random generator
idx.train <- sample(1:nrow(d),round(0.75*nrow(d)))
d.train   <- d[idx.train,]
d.test    <- d[-idx.train,]


## -----------------------------------------------------------------------------
# Fit the linear model, no default for family, so use 'gaussian':
lm.fit <- glm(medv ~ ., data = d.train)
summary(lm.fit)

# Make predictions:
pr.lm  <- predict(lm.fit,d.test)

# Calculate the MSE:
MSE.lm <- sum((pr.lm - d.test$medv)^2)/nrow(d.test)


## -----------------------------------------------------------------------------
# Store the maxima and minima:
d.maxs <- apply(d, 2, max) 
d.mins <- apply(d, 2, min)

# Rescale the data:
d.sc <- as.data.frame(scale(d, center = d.mins, 
                               scale  = d.maxs - d.mins))

# Split the data in training and testing set:
d.train.sc <- d.sc[idx.train,]
d.test.sc  <- d.sc[-idx.train,]


## -----------------------------------------------------------------------------
library(neuralnet)

# Since the shorthand notation y~. does not work in the 
# neuralnet() function we have to replicate it:
nm  <- names(d.train.sc)
frm <- as.formula(paste("medv ~", paste(nm[!nm %in% "medv"], 
                        collapse = " + ")))

nn2 <- neuralnet(frm, data = d.train.sc, hidden = c(7,5,5),
                 linear.output = T)


## ----fig.cap=c("A visualisation of the ANN. Note that we left out the weights, because there would be too many. With 13 variables, and three layers of respectively 7, 5, and 5 neurons, we have $13\\times 7 + 7 \\times 5 + 1 + 5 \\times 5 + 1 + 5 + 1 = 160$ parameters."),aNNvisualize----
plot(nn2, rep = "best", information = FALSE, 
     show.weights = FALSE)


## -----------------------------------------------------------------------------
# Our independent variable 'medv' is the 14th column, so:
pr.nn2 <- compute(nn2,d.test.sc[,1:13]) 
                                   
# Rescale back to original span:
pr.nn2 <- pr.nn2$net.result*(max(d$medv)-min(d$medv))+min(d$medv)
test.r <- (d.test.sc$medv)*(max(d$medv)-min(d$medv))+min(d$medv)

# Calculate the MSE:
MSE.nn2 <- sum((test.r - pr.nn2)^2)/nrow(d.test.sc)
print(paste(MSE.lm,MSE.nn2))


## ----fig.cap=c("A visualisation of the performance of the ANN (left) compared to the linear regression model (right)."),aNNvisualizePerformance----
par(mfrow=c(1,2))

plot(d.test$medv, pr.nn2, col='red',
     main='Observed vs predicted NN',
     pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright', legend='NN', pch=18, col='red', bty='n')

plot(d.test$medv,pr.lm,col='blue',
     main='Observed vs predicted lm',
     pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright', legend='LM', pch=18,col='blue', bty='n', 
       cex=.95)


## ----fig.cap=c("A visualisation of the performance of the ANN compared to the linear regression model with both models in one plot."),aNNperformance----
plot  (d.test$medv,pr.nn2,col='red',
       main='Observed vs predicted NN',
       pch=18,cex=0.7)
points(d.test$medv,pr.lm,col='blue',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend=c('NN','LM'),pch=18,
       col=c('red','blue'))


## ----xValidationMSE-----------------------------------------------------------
library(boot)
set.seed(1875)
lm.fit <- glm(medv ~ ., data = d)

# The estimate of prediction error is now here:
cv.glm(d, lm.fit, K = 10)$delta[1] 


## ----xValidationANN,message=FALSE---------------------------------------------
# Reminders:
d   <- Boston
nm  <- names(d)
frm <- as.formula(paste("medv ~", paste(nm[!nm %in% "medv"], 
                        collapse = " + ")))
# Store the maxima and minima:
d.maxs <- apply(d, 2, max) 
d.mins <- apply(d, 2, min)

# Rescale the data:
d.sc <- as.data.frame(scale(d, center = d.mins, 
                               scale  = d.maxs - d.mins))

# Set parameters:
set.seed(1873)
cv.error <- NULL  # Initiate to append later
k        <- 10    # The number of repetitions

# This code might be slow, so you can add a progress bar as follows:
#library(plyr) 
#pbar <- create_progress_bar('text')
#pbar$init(k)

# In k-fold cross validation, we must take care to select each
# observation just once in the testing set. This is made easy 
# with modelr:
library(modelr)
kFoldXval <- crossv_kfold(data = d.sc, k = 10, id = '.id')

# Do the k-fold cross validation:
for(i in 1:k){
    # <see note below>
    train.cv  <- kFoldXval$train[i]
    test.cv   <- kFoldXval$test[i]
    test.cv.df <- as.data.frame(test.cv)
    
    # Rebuild the formula (names are changed each run):
    nmKfold <- paste0('X', i, '.', nm)
    medvKfld <- paste0('X', i, '.medv')
    frmKfold <- as.formula(paste(medvKfld, "~", 
                             paste(nmKfold[!nmKfold %in% medvKfld], 
                             collapse = " + ")
		             )
			  )

    # Fit the NN:
    nn2       <- neuralnet(frmKfold, data = train.cv, 
                           hidden = c(7, 5, 5),
                           linear.output=TRUE
			   )

    # The explaining variables are in the first 13 rows, so:
    pr.nn2   <- compute(nn2, test.cv.df[,1:13])   
          
    pr.nn2   <- pr.nn2$net.result * (max(d$medv) - min(d$medv)) + 
                min(d$medv)
    test.cv.df.r <- test.cv.df[[medvKfld]] * 
                    (max(d$medv) - min(d$medv)) + min(d$medv)
    cv.error[i] <- sum((test.cv.df.r - pr.nn2)^2)/nrow(test.cv.df)    
    #pbar$step()  #uncomment to see the progress bar
}


## ----eval=FALSE---------------------------------------------------------------
## index    <- sample(1:nrow(d),round(0.9*nrow(d)))
## train.cv <- d.sc[index,]
## test.cv  <- d.sc[-index,]


## ----eval=FALSE---------------------------------------------------------------
## # <see note below>


## ----fig.cap=c("A boxplot for the MSE cross validation for the ANN."),plotXValidationANN----
# Show the mean of the MSE:
mean(cv.error)
cv.error

# Show the boxplot:
boxplot(cv.error,xlab='MSE',col='gray',
        border='blue',names='CV error (MSE)',
        main='Cross Validation error (MSE) for the ANN',
        horizontal=TRUE)


## ----eval=FALSE---------------------------------------------------------------
## svm(formula, data, subset, na.action = na.omit, scale = TRUE,
##     type = NULL, kernel = 'radial', degree = 3,
##     gamma = if (is.vector(x)) 1 else 1 / ncol(x), coef0 = 0,
##     cost = 1, nu = 0.5,  class.weights = NULL, cachesize = 40,
##     tolerance = 0.001, epsilon = 0.1, shrinking = TRUE,
##     cross = 0, probability = FALSE, fitted = TRUE, ...)


## ----message=FALSE------------------------------------------------------------
library(e1071)
svmCars1 <- svm(cyl ~ ., data = mtcars)
summary(svmCars1)


## -----------------------------------------------------------------------------
# split mtcars in two subsets (not necessary but easier later):
x <- subset(mtcars, select = -cyl)
y <- mtcars$cyl

# fit the model again as a classification model:
svmCars2 <- svm(cyl ~ ., data = mtcars, type = 'C-classification')

# create predictions
pred <- predict(svmCars2, x)

# show the confusion matrix:
table(pred, y)


## -----------------------------------------------------------------------------
svmTune <- tune(svm, train.x=x, train.y=y, kernel = "radial",
                ranges =  list(cost = 10^(-1:2), gamma = c(.5, 1, 2)))

print(svmTune)


## ----kMeansGg,fig.cap='The cars in the dataset mtcars with fuel consumption plotted in function of weight and coloured by the number of cylinders.'----
library(ggplot2)
library(ggrepel)  # provides geom_label_repel()
ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) + 
       geom_point(size = 5) + 
       geom_label_repel(aes(label = rownames(mtcars)),
                  box.padding   = 0.2, 
                  point.padding = 0.25,
                  segment.color = 'grey60')


## -----------------------------------------------------------------------------
ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) + 
       geom_point(size = 5) + 
       geom_text(aes(label = rownames(mtcars)), 
                 hjust = -0.2, vjust = -0.2)


## -----------------------------------------------------------------------------
# normalize weight and mpg
d <- data.frame(matrix(NA, nrow = nrow(mtcars), ncol = 1))
d <- d[,-1]  # d is an empty data frame with 32 rows

rngMpg <- range(mtcars$mpg, na.rm = TRUE)
rngWt  <- range(mtcars$wt,  na.rm = TRUE)

d$mpg_n <- (mtcars$mpg - rngMpg[1]) / rngMpg[2]
d$wt_n  <- (mtcars$wt  - rngWt[1])  / rngWt[2]

# Here is the k-means clustering itself. 
# Note the nstart parameter (the number of random starting sets)
carCluster <- kmeans(d, 3, nstart = 15)

print(carCluster)


## -----------------------------------------------------------------------------
table(carCluster$cluster, mtcars$cyl)
# Note that the rowl are the clusters (1, 2, 3) and the number of 
# cylinders are the columns (4, 6, 8).


## ----ggKMeans1,fig.cap='The result of $k$-means clustering with three clusters on the weight and fuel consumption for the dataset mtcars.'----
# Additionally, we customize colours and legend title. 
# First, we need color names:
my_colors <- if_else(carCluster$cluster == 1, "darkolivegreen3",
                if_else(carCluster$cluster == 2, "coral", "cyan3"))

# We already loaded the libraries as follows:
#library(ggplot2); library(ggrepel)

# Now we can create the plot:
ggplot(mtcars, aes(wt, mpg, fill = factor(carCluster$cluster))) + 
       geom_point(size = 5, colour = my_colors) + 
       scale_fill_manual('Clusters', 
                   values = c("darkolivegreen3","coral", "cyan3"))+
       geom_label_repel(aes(label = rownames(mtcars)),
                  box.padding   = 0.2, 
                  point.padding = 0.25,
                  segment.color = 'grey60') 



## ----fig.cap=c('The variance explained by each principal component in the dataset {\tt mtcars}.', 'A projection of the space of alternatives in the 2D-plane formed by the two largest principal components of mtcars.')----
# Normalize the whole mtcars dataset:
d <- data.frame(matrix(NA, nrow = nrow(mtcars), ncol = 1))
d <- d[,-1]  # d is an empty data frame with 32 rows
for (k in 1:ncol(mtcars)) {
rng <- range(mtcars[, k], na.rm = TRUE)
d[, k]  <- (mtcars[, k]  - rng[1])  / rng[2]
}
colnames(d) <- colnames(mtcars)
rownames(d) <- rownames(mtcars)

# The PCA analysis:
pca1 <- prcomp(d) 
summary(pca1)

# Note also:
class(pca1)


## ----prcompPlot,fig.cap=c("The plot() function applied on a {\\tt prcomp} object visualises the relative importance of the different principal components.", "The custom function {\\tt biplot()} project all data in the plane that is span by the two major PCs.")----
# Plot for the prcomp object shows the variance explained by each PC
plot(pca1, type = 'l')

# biplot shows a projection in the 2D plane (PC1, PC2)
biplot(pca1)


## ----prcompGGPlot, fig.cap='A projection in the plane of the two major principal components via ggplot2. It looks good, but the labels are cluttered.'----
# Same plot with ggplot2:
library(ggplot2)
library(ggfortify)
library(cluster)

autoplot(pca1, data=d, label=TRUE, shape=FALSE, colour='mpg',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3
         )


## ----carCluster---------------------------------------------------------------
carCluster <- kmeans(d, 4, nstart = 10)


## ----fig.cap='The projection of mtcars in the surface formed by the two first principal components and colored by number of cluster.\\label{fig:carsCluster}'----
library(ggplot2)
library(ggrepel)

autoplot(pca1, data=d, label=FALSE, shape=18, size = 5, 
         alpha = 0.6, colour = factor(carCluster$cluster),
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 5) + 
       geom_label_repel(aes(label = rownames(mtcars)),
                  box.padding   = 0.2, 
                  point.padding = 0.25,
                  segment.color = 'grey60')


## -----------------------------------------------------------------------------
my_colors <- if_else(carCluster$cluster == 1,"darkolivegreen3",
                     if_else(carCluster$cluster == 2, "coral", 
		     if_else(carCluster$cluster == 3, "cyan3", 
		             "gray80")))

ggplot(pca1$x,
       aes(x=PC1,y=PC2, fill = factor(carCluster$cluster))) +
       geom_point(size = 5, alpha = 0.7, colour = my_colors)+ 
       scale_fill_manual('Clusters', 
                         values =c("darkolivegreen3","coral", 
			 "cyan3", "gray80")) +
       geom_label_repel(aes(label = rownames(mtcars)),
                  box.padding   = 0.2, 
                  point.padding = 0.25,
                  segment.color = 'grey60')


## -----------------------------------------------------------------------------
pca1 <- prcomp(d)   # unchanged

# Extract the first 3 PCs
PCs <- data.frame(pca1$x[,1:3])

# Then use this in the kmeans function for example.
carCluster <- kmeans(PCs, 4, nstart=10, iter.max=1000)


## ----threePCs,fig.cap='Two dimensional projections of the dependency structure of the data in the first principal components. Note that in this plot, we see different 2D projections of 3D data.'----
# Reminder of previous code:
d <- data.frame(matrix(NA, nrow = nrow(mtcars), ncol = 1))
d <- d[,-1]  # d is an empty data frame with 32 rows
for (k in 1:ncol(mtcars)) {
  rng     <- range(mtcars[, k], na.rm = TRUE)
  d[, k]  <- (mtcars[, k]  - rng[1])  / rng[2]
  }
colnames(d) <- colnames(mtcars)
rownames(d) <- rownames(mtcars)
pca1 <- prcomp(d)               # runs the PCA
PCs <- data.frame(pca1$x[,1:3]) # extracts the first 3 PCs
# Now, PCs holds the three major PCs.

# -- New code below:
# Plot those three PCs:
plot(PCs, col=carCluster$clust, pch=16, cex=3)


## ----plot3DCluster,fig.cap='A three dimensional plot of the cars with on the $z$-axis the first principal component, the second on the $y$-axis and the third along the $z$-axis.'----
library(plot3D)
scatter3D(x = PCs$PC1, y = PCs$PC2, z = PCs$PC3, 
   phi = 45, theta = 45,
   pch = 16, cex = 1.5, bty = "f",
   clab = "cluster", 
   colvar = as.integer(carCluster$cluster), 
   col = c("darkolivegreen3", "coral", "cyan3", "gray"),
   colkey = list(at = c(1, 2, 3, 4),
          addlines = TRUE, length = 0.5, width = 0.5,
          labels = c("1", "2", "3", "4"))
   )
text3D(x = PCs$PC1, y = PCs$PC2, z = PCs$PC3,  labels = rownames(d),
        add = TRUE, colkey = FALSE, cex = 1.2)


## ----eval=FALSE---------------------------------------------------------------
## library(plotly)
## plot_ly(x = PCs$PC1, y = PCs$PC2, z = PCs$PC3,
##         type="scatter3d", mode="markers",
## 	color=factor(carCluster$cluster))


## ----fanny,fig.cap=c('A plot with {\\tt autoplot()}, enhanced with {\\tt ggrepel} of the fuzzy clustering for the dataset mtcars.')----
library(tidyverse)  # provides if_else
library(ggplot2)    # 2D plotting 
library(ggfortify)
library(cluster)    # provides fanny (the fuzzy clustering)
library(ggrepel)    # provides geom_label_repel (de-clutter labels)


carCluster <- fanny(d, 4)
my_colors <- if_else(carCluster$cluster == 1, "coral",
               if_else(carCluster$cluster == 2, "darkolivegreen3", 
	       if_else(carCluster$cluster == 3, "cyan3", 
	         "darkorchid1")))

# Autoplot with visualization of 4 clusters
autoplot(carCluster, label=FALSE, frame=TRUE,  frame.type='norm', 
         shape=16,
         loadings=TRUE,  loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 5,
	 loadings.label.vjust = 1.2, loadings.label.hjust = 1.3) + 
       geom_point(size = 5, alpha = 0.7, colour = my_colors) + 
       geom_label_repel(aes(label = rownames(mtcars)),
                  box.padding   = 0.2, 
                  point.padding = 0.25,
                  segment.color = 'grey40') + 
		  theme_classic()



## ----carsHC,fig.cap='A hierarchical cluster for the dataset mtcars.'----------
# Compute hierarchical clustering
library(tidyverse)
cars_hc <- mtcars                      %>%
           scale                       %>% # scale the data
           dist(method = "euclidean")  %>% # dissimilarity matrix
           hclust(method = "ward.D2")      # hierachical clustering

plot(cars_hc)


## ----eval=FALSE---------------------------------------------------------------
## library(class)
## knn(train, test, cl, k = 1, l = 0, prob = FALSE, use.all = TRUE)


## -----------------------------------------------------------------------------
library(tidyverse)
library(modelr)


## -----------------------------------------------------------------------------
d   <- mtcars
lm1 <- lm(mpg ~ wt + cyl, data = d)


## ----eval=FALSE---------------------------------------------------------------
## add_predictions(data, model, var = "pred", type = NULL)


## -----------------------------------------------------------------------------
library(modelr)

# Use the data defined above:
d1 <- d %>% add_predictions(lm1)

# d1 has now an extra column "pred"
head(d1)


## ----eval=FALSE---------------------------------------------------------------
## add_residuals(data, model, var = "resid")


## -----------------------------------------------------------------------------
d2 <- d1 %>% add_residuals(lm1)

# d2 has now an extra column "resid"
head(d2)


## ----eval=FALSE---------------------------------------------------------------
## bootstrap(data, n, id = ".id")


## ----SCM,fig.cap='The results of the bootstrap exercise: a set of estimates for each coefficient.'----
set.seed(1872)   # make sure that results can be replicated
library(modelr)  # provides bootstrap
library(purrr)   # provides map, map_df, etc.
library(ggplot2) # provides ggplot
d    <- mtcars
boot <- bootstrap(d, 10)

# Now, we can leverage tidyverse functions such as map to create 
# multiple models on the 10 datasets
models <- map(boot$strap, ~ lm(mpg ~ wt + cyl, data = .))

# The function tidy of broom (also tidyverse) allows to create a
# dataset based on the list of models. Broom is not loaded, because
# it also provides a function bootstrap().
tidied <- map_df(models, broom::tidy, .id = "id")

# Visualize the results with ggplot2:
p <- ggplot(tidied, aes(estimate)) + 
     geom_histogram(bins = 5, col = 'red', fill='khaki3', 
                    alpha = 0.5) + 
     ylab('Count') + 
     xlab('Estimate of the coefficient in the plot-title') +
     facet_grid(. ~ term, scales = "free")
p


## -----------------------------------------------------------------------------
# load modelr:
library(modelr)

# Fit a model:
lm1 <- lm(mpg ~ wt + qsec + am, data = mtcars)

# MSE (mean square error):
mse(lm1, mtcars)

# RMSE (root mean square error):
rmse(lm1, mtcars)

# MAD (mean absolute error):
mae(lm1, mtcars)

# Quantiles of absolute error:
qae(lm1, mtcars)

# R-square (variance of predictions divided by the variance of the 
# response variable):
rsquare(lm1, mtcars)


## -----------------------------------------------------------------------------
set.seed(1871)

rs  <- mtcars  %>%
       resample_partition(c(train = 0.6, test = 0.4))
lm2 <- lm(mpg ~ wt + qsec + am, data = rs$train)
rmse(lm2, rs$train); rmse(lm2, rs$test)

# or with the pipe operator:
lm2 %>% rmse(rs$train)
lm2 %>% rmse(rs$test)


## -----------------------------------------------------------------------------
# Fit the model:
lm1 <- lm(mpg ~ wt + qsec + am, data = mtcars)

# Add the predictions and residuals:
df <- mtcars               %>% 
      add_predictions(lm1) %>%
      add_residuals(lm1)

# The predictions are now available in $pred
head(df$pred)

# The residuals are now available in $resid
head(df$resid)

# It is now easy to do something with those predictions and 
# residuals, e.g. the following 3 lines all do the same:
sum((df$pred - df$mpg)^2) / nrow(mtcars)
sum((df$resid)^2) / nrow(mtcars)
mse(lm1, mtcars)


## ----spacingGrid,fig.cap='A spacing grid for the predictions of {\tt mpg}.'----
d <- data_grid(mtcars, wt = seq_range(wt, 10), qsec, am) %>% 
     add_predictions(lm1)
plot(d$wt, d$pred)


## ----bootstrap,fig.cap=c("Bootstrapping the returns of the S\\&P500 index.")----
# Create the sample:
SP500_sample <- sample(SP500,size=100)

# Change plotting to 4 plots in one output:
par(mfrow=c(2,2))

# The histogram of the complete dataset:
hist(SP500,main="(a) Histogram of all data",fr=FALSE,
     breaks=c(-9:5),ylim=c(0,0.4))

# The histogram of the sample:
hist(SP500_sample,main="(b) Histogram of the sample",
     fr=FALSE,breaks=c(-9:5),ylim=c(0,0.4))

# The boxplot of the complete dataset:
boxplot(SP500,main="(c) Boxplot of all data",ylim=c(-9,5))

# The boxplot of the complete sample:
boxplot(SP500_sample,main="(c) Boxplot of the sample",
        ylim=c(-9,5))

# Reset the plot parameters:
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
mean(SP500)
mean(SP500_sample)
sd(SP500)
sd(SP500_sample)


## -----------------------------------------------------------------------------
# Bootstrap generates a number of re-ordered datasets
boot <- bootstrap(mtcars, 3)
# The datasets are now in boot$strap[[n]]
# with n between 1 and 3

# e.g. the 3rd set is addressed as follows:
nrow(boot$strap[[3]])
mean(as.data.frame(boot$strap[[3]])$mpg)

# It is also possible to coerce the selections into a data-frame:
df <- as.data.frame(boot$strap[[3]])
class(df)


## ----MultiHistsBoot,fig.cap="The histograms of the different coefficients of the linear regression model predicting the mpg in the dataset mtcars. We show (a) Estimate for wt., (b) Estimate for hp., (c) Estimate for am:vs., and (d) Estimate for the intercept.", fig.height='0.6\\textheight',out.height='0.6\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
set.seed(1871)
library(purrr)  # to use the function map()
boot <- bootstrap(mtcars, 150)
     
lmodels <- map(boot$strap, ~ lm(mpg ~ wt + hp + am:vs, data = .))

# The function tidy of broom turns a model object in a tibble:
df_mods <- map_df(lmodels, broom::tidy, .id = "id")

# Create the plots of histograms of estimates for the coefficients:
par(mfrow=c(2,2))
hist(subset(df_mods, term == "wt")$estimate, col="khaki3",
     main = '(a) wt', xlab = 'estimate for wt')
hist(subset(df_mods, term == "hp")$estimate, col="khaki3",
     main = '(b) hp', xlab = 'estimate for hp')
hist(subset(df_mods, term == "am:vs")$estimate, col="khaki3",
     main = '(c) am:vs', xlab = 'estimate for am:vs')
hist(subset(df_mods, term == "(Intercept)")$estimate, col="khaki3",
     main = '(d) intercept', xlab = 'estimate for the intercept')
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
d <- mtcars                # get data
set.seed(1871)             # set the seed for the random generator
idx.train <- sample(1:nrow(d),round(0.75*nrow(d)))
d.train <- d[idx.train,]   # positive matches for training set
d.test  <- d[-idx.train,]  # the opposite to the testing set


## -----------------------------------------------------------------------------
set.seed(1870)
sample_cars <- mtcars %>%
               resample(sample(1:nrow(mtcars),5)) # random 5 cars

# This is a resample object (indexes shown, not data):
sample_cars  

#turn it into data:
as.data.frame(sample_cars)

# or into a tibble
as_tibble(sample_cars)

# or use the indices to get to the data:
mtcars[as.integer(sample_cars),]


## -----------------------------------------------------------------------------
library(modelr)
rs <- mtcars  %>%
      resample_partition(c(train = 0.6, test = 0.4))

# address the datasets with: as.data.frame(rs$train)
#                            as.data.frame(rs$test)

# Check execution:
lapply(rs, nrow)


## -----------------------------------------------------------------------------
# 0. Store training and test dataset for further use:
d_train  <- as.data.frame(rs$train)
d_test   <- as.data.frame(rs$test)

# 1. Fit the model on the training dataset:
lm1      <- lm(mpg ~ wt + hp + am:vs, data = rs$train)

# 2. Calculate the desired performance measure (e.g.
# root mean square error (rmse)):
rmse_trn <- lm1 %>% rmse(rs$train)
rmse_tst <- lm1 %>% rmse(rs$test)
print(rmse_trn)
print(rmse_tst)


## -----------------------------------------------------------------------------
# 2. Add predictions and residuals:
x_trn  <- add_predictions(d_train, model = lm1) %>% 
          add_residuals(model = lm1)
x_tst  <- add_predictions(d_test,  model = lm1) %>% 
          add_residuals(model = lm1)


# 3. Calculate the desired risk metrics (via the residuals):
RMSE_trn  <- sqrt(sum(x_trn$resid^2) / nrow(d_train))
RMSE_tst  <- sqrt(sum(x_tst$resid^2) / nrow(d_test))
print(RMSE_trn)
print(RMSE_tst)


## -----------------------------------------------------------------------------
# Monte Carlo cross validation
cv_mc <- crossv_mc(data = mtcars, # the dataset to split
           n = 50,      # n random partitions train and test
	   test = 0.25, # validation set is 25%
	   id = ".id")  # unique identifier for each model

# Example of use:

# Access the 2nd test dataset:
d <- data.frame(cv_mc$test[2])

# Access mpg in that data frame:
data.frame(cv_mc$test[2])$mpg

# More cryptic notations are possible to obtain the same:
mtcars[cv_mc[[2]][[2]][2]$idx,1]


## ----mcCrossVHist,fig.cap="The histogram of the RMSE for a Monte Carlo cross validation on the dataset mtcars."----
set.seed(1868)
library(modelr)     # sample functions
library(purrr)      # to use the function map()

cv_mc <- crossv_mc(mtcars, n = 50, test = 0.40)
mods  <- map(cv_mc$train, ~ lm(mpg ~ wt + hp + am:vs, data = .))
RMSE  <- map2_dbl(mods, cv_mc$test, rmse)
hist(RMSE, col="khaki3")


## ----eval=FALSE---------------------------------------------------------------
## library(magrittr)  # to access the %T>% pipe
## crossv <- mtcars                                          %>%
##           crossv_mc(n = 50, test = 0.40)
## RMSE <- crossv                                            %$%
##         map(train, ~ lm(mpg ~ wt + hp + am:vs, data = .)) %>%
##         map2_dbl(crossv$test, rmse)                      %T>%
##         hist(col = "khaki3", main ="Histogram of RMSE",
## 	     xlab = "RMSE")


## -----------------------------------------------------------------------------
library(modelr)
# k-fold cross validation
cv_k  <- crossv_kfold(data = mtcars, 
           k = 5,      # number of folds
           id = ".id") # unique identifier for each


## -----------------------------------------------------------------------------
cv_k$test


## ----fig.cap='Histogram of the RMSE based on a 5-fold cross validation. The histogram indeed shows that there were 5 observations. Note the significant spread of RMSE: the largest one is about four times the smallest.'----
set.seed(1868)
library(modelr)
library(magrittr)  # to access the %T>% pipe
crossv <- mtcars                                           %>% 
          crossv_kfold(k = 5)                             
RMSE <- crossv                                             %$%
        map(train, ~ lm(mpg ~ wt + hp + am:vs, data = .))  %>%
        map2_dbl(crossv$test, rmse)                        %T>%
        hist(col = "khaki3", main ="Histogram of RMSE", 
	     xlab = "RMSE")


## ----initQuantmod,warning=FALSE,message=FALSE---------------------------------
# Install quantmod:
if(!any(grepl("quantmod", installed.packages()))){
    install.packages("quantmod")}

# Load the library:
library(quantmod)

## ----echo=FALSE,include=FALSE-------------------------------------------------
options("getSymbols.yahoo.warning"=FALSE)


## ----warning=FALSE,message=FALSE----------------------------------------------
# Download historic data of the Google share price:
getSymbols("GOOG",src="yahoo")         # get Google's history
getSymbols(c("GS","GOOG"),src="yahoo") # to load more than one


## -----------------------------------------------------------------------------
setSymbolLookup(HSBC='yahoo',GOOG='yahoo')
setSymbolLookup(DEXJPUS='FRED')
setSymbolLookup(XPTUSD=list(name="XPT/USD",src="oanda"))

# Save the settings in a file:
saveSymbolLookup(file="qmdata.rda")  
# Use this in new sessions calling:
loadSymbolLookup(file="qmdata.rda")

# We can also download a list of symbols as follows:
getSymbols(c("HSBC","GOOG","DEXJPUS","XPTUSD")) 


## ----GetSymbols---------------------------------------------------------------
stockList <- stockSymbols()  # get all symbols
nrow(stockList)     # number of symbols
colnames(stockList) # information in this list


## ----getFX,warning=FALSE------------------------------------------------------
getFX("EUR/PLN",from="2019-01-01")


## ----getHSBC------------------------------------------------------------------
getSymbols("HSBC",src="yahoo") #get HSBC's data from Yahoo


## ----demoQuantmod,fig.cap=c("Demonstration of the barChart() function of the package quantmod.","Demonstration of the lineChart() function of the package quandmod.","Demonstration of the candleChart() function of the package quantmod.")----

# 1. The bar chart:
barChart(HSBC)     

# 2. The line chart:
lineChart(HSBC)
# Note: the lineChart is also the default that yields the same
#       result as plot(HSBC)

# 3. The candle chart:
candleChart(HSBC, subset='last 1 years',theme="white",
            multi.col=TRUE)


## ----eval=FALSE---------------------------------------------------------------
## candleChart(HSBC,subset='2018::2018-01')
## candleChart(HSBC,subset='last 1 months')


## ----bollingerBands,fig.keep='last',fig.cap=c("Bollinger bands with the package quandmod.")----
getSymbols(c("HSBC"))
chartSeries(HSBC, subset='last 4 months')
addBBands(n = 20, sd = 2, maType = "SMA", draw = 'bands', 
          on = -1)


## ----eval=FALSE---------------------------------------------------------------
## myxtsdata["2008-01-01/2010-12-31"]  # between 2 date-stamps
## 
## # All data before or after a certain time-stamp:
## xtsdata["/2007"]  # from start of data until end of 2007
## xtsdata["2009/"]  # from 2009 until the end of the data
## 
## # Select the data between different hours:
## xtsdata["T07:15/T09:45"]


## ----eval=FALSE---------------------------------------------------------------
## HSBC['2017']    #returns HSBC's OHLC data for 2017
## HSBC['2017-08'] #returns HSBC's OHLC data for August 2017
## HSBC['2017-06::2018-01-15'] # from June 2017 to Jan 15 2018
## 
## HSBC['::']     # returns all data
## HSBC['2017::'] # returns all data in HSBC, from 2017 onward
## my.selection <- c('2017-01','2017-03','2017-11')
## HSBC[my.selection]


## ----eval=FALSE---------------------------------------------------------------
## last(HSBC)               #returns the last quotes
## last(HSBC,5)             #returns the last 5 quotes
## last(HSBC, '6 weeks')    # the last 6 weeks
## last(HSBC, '-1 weeks')   # all but the last week
## last(HSBC, '6 months')   # the last 6 months
## last(HSBC, '3 years')    # the last 3 years
## 
## # these functions can also be combined:
## last(first(HSBC, '3 weeks'), '5 days')
## #


## ----eval=FALSE---------------------------------------------------------------
## periodicity(HSBC)
## unclass(periodicity(HSBC))
## to.weekly(HSBC)
## to.monthly(HSBC)
## periodicity(to.monthly(HSBC))
## ndays(HSBC); nweeks(HSBC); nyears(HSBC)


## -----------------------------------------------------------------------------
getFX("USD/EUR")
periodicity(USDEUR)
to.monthly(USDEUR)
periodicity(to.monthly(USDEUR))


## -----------------------------------------------------------------------------
endpoints(HSBC,on="years") 

# find the maximum closing price each year
apply.yearly(HSBC,FUN=function(x) {max(Cl(x)) } )

# the same thing - only more general
subHSBC <- HSBC['2012::']
period.apply(subHSBC,endpoints(subHSBC,on='years'),
             FUN=function(x) {max(Cl(x))} )

# the following line does the same but is faster:
as.numeric(period.max(Cl(subHSBC),endpoints(subHSBC,
           on='years')))


## -----------------------------------------------------------------------------
seriesHi(HSBC)
has.Cl(HSBC)
tail(Cl(HSBC))


## ----eval=FALSE---------------------------------------------------------------
## Lag(Cl(HSBC))
## Lag(Cl(HSBC),c(1,5,10)) #One, five and ten period lags
## Next(OpCl(HSBC))
## 
## # Open to close one, two and three-day lags:
## Delt(Op(HSBC),Cl(HSBC),k=1:3)


## ----eval=FALSE---------------------------------------------------------------
## dailyReturn(HSBC)
## weeklyReturn(HSBC)
## monthlyReturn(HSBC)
## quarterlyReturn(HSBC)
## yearlyReturn(HSBC)
## allReturns(HSBC)     # all previous returns


## ----quantmodModel1,warning=FALSE---------------------------------------------
# First, we create a quantmod object. 
# At this point, we do not need to load data.

setSymbolLookup(SPY='yahoo',
   VXN=list(name='^VIX',src='yahoo'))

qmModel <- specifyModel(Next(OpCl(SPY)) ~ OpCl(SPY) + Cl(VIX))
head(modelData(qmModel))


## ----qmodMod,fig.cap=c("The evolution of the HSBC share for the last ten years.")----
getSymbols('HSBC',src='yahoo') #google doesn't carry the adjusted price
lineChart(HSBC)


## -----------------------------------------------------------------------------
HSBC.tmp   <- HSBC["2010/"]          #see: subsetting for xts objects    


## -----------------------------------------------------------------------------
# use 70% of the data to train the model:
n          <- floor(nrow(HSBC.tmp) * 0.7)
HSBC.train <- HSBC.tmp[1:n]               # training data
HSBC.test  <- HSBC[(n+1):nrow(HSBC.tmp)]  # test-data
# head(HSBC.train)


## -----------------------------------------------------------------------------
# making sure that whenever we re-run this the latest data 
# is pulled in:
m.qm.tr <- specifyModel(Next(Op(HSBC.train)) ~ Ad(HSBC.train)
          + Hi(HSBC.train) - Lo(HSBC.train) + Vo(HSBC.train))

D <- modelData(m.qm.tr)


## -----------------------------------------------------------------------------
# Add the additional column:
D$diff.HSBC <- D$Hi.HSBC.train - D$Lo.HSBC.train 

# Note that the last value is NA:
tail(D, n = 3L)

# Since the last value is NA, let us remove it:
D <- D[-nrow(D),]                           


## -----------------------------------------------------------------------------
colnames(D) <- c("Next.Op","Ad","Hi","Lo","Vo","Diff")


## -----------------------------------------------------------------------------
m1 <- lm(D$Next.Op ~ D$Ad + D$Diff + D$Vo)
summary(m1)


## -----------------------------------------------------------------------------
m2 <- lm(D$Next.Op ~ D$Ad + D$Diff)
summary(m2)


## ----fig.cap=c("The Q-Q plot of our naive model to forecast the next opening price of the HSBC stock. The results seem to be reasonable.\\label{fig:qq:m2}")----
qqnorm(m2$residuals)
qqline(m2$residuals,col='blue',lwd=2)


## -----------------------------------------------------------------------------
m.qm.tst <- specifyModel(Next(Op(HSBC.test)) ~  Ad(HSBC.test) 
             + Hi(HSBC.test) - Lo(HSBC.test) + Vo(HSBC.test))

D.tst <- modelData(m.qm.tst)
D.tst$diff.HSBC.test <- D.tst$Hi.HSBC.test-D.tst$Lo.HSBC.test  
#tail(D.tst)                           # the last value is NA
D.tst <- D[-nrow(D.tst),]  # remove the last value that is NA

colnames(D.tst) <- c("Next.Op","Ad","Hi","Lo","Vo","Diff")


## -----------------------------------------------------------------------------
a   <- coef(m2)['(Intercept)']
bAd <- coef(m2)['D$Ad']
bD  <- coef(m2)['D$Diff']
est <- a + bAd * D.tst$Ad + bD * D.tst$Diff


## -----------------------------------------------------------------------------
# -- Mean squared prediction error (MSPE)
#sqrt(mean(((predict(m2,newdata = D.tst) - D.tst$Next.Op)^2)))
sqrt(mean(((est - D.tst$Next.Op)^2)))

# -- Mean absolute errors (MAE)
mean((abs(est - D.tst$Next.Op)))

# -- Mean absolute percentage error (MAPE)
mean((abs(est - D.tst$Next.Op))/D.tst$Next.Op)

# -- squared sum of residuals
print(sum(residuals(m2)^2))  

# -- confidence intervals for the model
print(confint(m2)) 


## -----------------------------------------------------------------------------
# Compare the coefficients in a refit:
m3 <- lm(D.tst$Next.Op ~ D.tst$Ad + D.tst$Diff)
summary(m3)


## -----------------------------------------------------------------------------
M0 <- matrix(c(
 1.6 , -0.83 , 1.4 , 4.7 , 1 , 0.9 , 1.1 ,
 1.8 , -0.83 , 1.0 , 4.7 , 1 , 0.9 , 0.8 ,
 1.8 , -0.83 , 1.2 , 4.7 , 1 , 0.9 , 0.6 ,
 1.6 , -1.24 , 1.4 , 2.8 , 1 , 0.9 , 0.8 ,
 0.9 , -0.83 , 1.4 , 4.7 , 1 , 0.7 , 0.8 ,
 0.9 , -0.83 , 0.8 , 4.7 , 1 , 0.7 , 0.6 ,
 0.7 ,  1.02 , 0.2 , 2.0 , 3 , 1.1 , 1.3 ,
 1.1 ,  0.52 , 1.0 , 1.3 , 3 , 0.6 , 0.9 ,
 1.2 , -0.83 , 1.3 , 4.7 , 1 , 0.8 , 0.5 ,
 0.9,  0.18 , 0.9 , 7.3 ,  1 , 0.8 , 0.6 ),
 byrow = TRUE,
 ncol = 7)
colnames(M0) <- c("tlnt","stab","cost","infl","trvl","infr","life")
# We use the IATA code of a nearby airport as abbreviation,
# so, instead of:
#rownames(M0) <- c("Bangalore", "Mumbai", "Delhi", "Manilla", 
#                 "Hyderabad", "Sao Polo", "Dublin", "Krakow", 
#                 "Chennai", "Buenos Aires")
# ... we use this:
rownames(M0) <- c("BLR", "BOM", "DEL", "MNL", "HYD", "GRU", 
                  "DUB", "KRK", "MAA", "EZE")

M0


## -----------------------------------------------------------------------------
# Political stability is a number between -2.5 and 2.5
# So, we make it all positive by adding 2.5:
M0[,2] <- M0[,2] + 2.5

# Lower wage inflation is better, so invert the data:
M0[,4] <- 1 / M0[,4]

# then we define a function:

# mcda_rescale_dm 
# Rescales a decision matrix M
# Arguments:
#    M -- decision matrix
#         criteria in columns and higher numbers are better.
# Returns
#    M -- normalised decision matrix
mcda_rescale_dm <- function (M) {
  colMaxs <- function(M) apply(M, 2, max, na.rm = TRUE)
  colMins <- function(M) apply(M, 2, min, na.rm = TRUE)
  M <- sweep(M, 2, colMins(M), FUN="-")
  M <- sweep(M, 2, colMaxs(M) - colMins(M), FUN="/")
  M
}

# Use this function:
M <- mcda_rescale_dm(M0)

# Show the new decision matrix elegantly:
knitr::kable(round(M,2))


## -----------------------------------------------------------------------------
# mcda_get_dominated
# Finds the alternatives that are dominated by others
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
# Returns
#    Dom -- prefM -- a preference matrix with 1 in position ij 
#                    if alternative i is dominated by alternative j.
mcda_get_dominated <- function(M) {
  Dom  <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  dominatedOnes <- c()
  for (i in 1:nrow(M)) {
    for (j in 1:nrow(M)) {
      isDom <- TRUE
      for (k in 1:ncol(M)) {
        isDom <- isDom && (M[i,k] >= M[j,k])
      }
      if(isDom && (i != j)) {
        Dom[j,i] <- 1
        dominatedOnes <- c(dominatedOnes,j)
      }
    }
  }
  colnames(Dom) <- rownames(Dom) <- rownames(M)
  class(Dom) <- "prefM"
  Dom
}


## -----------------------------------------------------------------------------
# mcda_get_dominants
# Finds the alternatives that dominate others
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
# Returns
#    Dom -- prefM -- a preference matrix with 1 in position ij 
#                    if alternative i dominates alternative j.
mcda_get_dominants <- function (M) t(mcda_get_dominated(M))


## -----------------------------------------------------------------------------
Dom <- mcda_get_dominants(M)
print(Dom)


## -----------------------------------------------------------------------------
# mcda_del_dominated
# Removes the dominated alternatives from a decision matrix
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
# Returns
#    A decision matrix without the dominated alternatives
mcda_del_dominated <- function(M) {
  Dom <- mcda_get_dominated(M)
  M[rowSums(Dom) == 0,]
}


## -----------------------------------------------------------------------------
M1 <- mcda_del_dominated(M)
round(M1,2)


## -----------------------------------------------------------------------------
M1 <- mcda_rescale_dm(M1)


## -----------------------------------------------------------------------------
# Firs, we load diagram:
require(diagram)

# plot.prefM
# Specific function to handle objects of class prefM for the
# generic function plot()
# Arguments:
#    PM  -- prefM -- preference matrix
#    ... -- additional arguments passed to plotmat()
#           of the package diagram.
plot.prefM <- function(PM, ...)
{
  X <- t(PM) # We want arrows to mean '... is better than ...'
             # plotmat uses the opposite convention.
  plotmat(X,  
          box.size    = 0.1, 
	  cex.txt     = 0, 
	  lwd         = 5 * X,  # lwd proportional to preference
	  self.lwd    = 3,
	  lcol        = 'blue',
	  self.shiftx = c(0.06, -0.06, -0.06, 0.06),
	  box.lcol    = 'blue',
	  box.col     = 'khaki3',
	  box.lwd     = 2,
	  relsize     = 0.9,
	  box.prop    = 0.5,
	  endhead     = FALSE,
	  main        = "",
	  ...)
}


## ----fig.cap='A visualization of the dominance relationship.'-----------------
# We pass the argument 'curve = 0' to the function plotmat,
# since otherwise in this case the arrow from BLR to MAA 
# would be hidden after the box of EZE.
plot(Dom, curve = 0)


## -----------------------------------------------------------------------------
# mcda_get_dominated
# Finds the alternatives that are dominated by others
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
#    w -- numeric vector of weights for the criteria
# Returns
#    a vector with a score for each alternative
mcda_wsm <- function(M, w) {
  X <- M %*% w
  colnames(X) <- 'pref'
  X
}


## -----------------------------------------------------------------------------
# the critia: "tlnt" "stab" "cost" "infl" "trvl" "infr" "life"
w <- c(        0.125, 0.2,   0.2,   0.2,  0.175,  0.05,  0.05)
w <- w / sum(w)  # the sum was 1 already, but just to be sure.

# now we can execute our function mcda_wsm()
mcda_wsm(M1, w)


## -----------------------------------------------------------------------------
# mcda_wsm_score
# Returns the scores for each of the alternative for each of 
# the criteria weighted by their weights.
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
#    w -- numeric vector of weights for the criteria
# Returns
#    a score-matrix of class scoreM
mcda_wsm_score <- function(M, w) {
  X <- sweep(M1,MARGIN=2,w,`*`)
  class(X) <- 'scoreM'
  X
}

# plot.scoreM
# Specific function for an object of class scoreM for the 
# generic function plot().
# Arguments:
#    M -- scoreM -- score matrix
# Returns:
#    plot
plot.scoreM <- function (M) {
  # 1. order the rows according to rowSums
  M <- M[order(rowSums(M),decreasing=T),]
  
  # 2. use a bar-plot on the transposed matrix
  barplot(t(M), 
     legend = colnames(M),
     xlab   = 'Score',
     col    = rainbow(ncol(M))
     )
}


## ----fig.cap='The scores of different cities according to the WSM.\\label{fig:wsm_score}'----
sM <- mcda_wsm_score(M1, w)
plot(sM)


## -----------------------------------------------------------------------------
# mcda_electre Type 2
# Push the preference matrixes PI.plus, PI.min and 
# PI.indif in the environment that calls this function.
# Arguments:
#    M -- normalized decision matrix with alternatives in rows,
#         criteria in columns and higher numbers are better.
#    w -- numeric vector of weights for the criteria
# Returns nothing but leaves as side effect:
#    PI.plus   -- the matrix of preference 
#    PI.min    -- the matrix of non-preference
#    PI.indif  -- the indifference matrix
mcda_electre <- function(M,  w) {
  # initializations
  PI.plus  <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  PI.min   <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  PI.indif <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))

  # calculate the preference matrix
  for (i in 1:nrow(M)){
    for (j in 1:nrow(M)) {
      for (k in 1:ncol(M)) {
        if (M[i,k] > M[j,k]) {
          PI.plus[i,j] <<- PI.plus[i,j] + w[k]
        }
        if (M[j,k] > M[i,k]) {
          PI.min[i,j] <<- PI.min[i,j] + w[k]
        }
        if (M[j,k] == M[i,k]) {
          PI.indif[j,i] <<- PI.indif[j,i] + w[k]
        }
      }
    }
  }
}


## -----------------------------------------------------------------------------
# mcda_electre1
# Calculates the preference matrix for the ELECTRE method
# Arguments:
#    M -- decision matrix (colnames are criteria, rownames are alternatives)
#    w -- vector of weights
#    Lambda -- the cutoff for the levels of preference
#    r -- the vector of maximum inverse preferences allowed
#    index -- one of ['C1', 'C2']
# Returns:
#    object of class prefM (preference matrix)
mcda_electre1 <- function(M,  w, Lambda, r, index='C1') {
  # get PI.plus, PI.min and PI.indif 
  mcda_electre(M,w)
  
  # initializations
  CM <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  PM <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  colnames(PM) <- rownames(PM) <- rownames(M)
  
  # calcualte the preference matrix
  if (index == 'C1') {
    # for similarity index C1
    for (i in 1:nrow(M)){
      for (j in 1:nrow(M)) {
        CM[i,j] <- (PI.plus[i,j] + PI.indif[i,j]) / (PI.plus[i,j] + 
	            PI.indif[i,j] + PI.min[i,j])
        if((CM[i,j] > Lambda) && ((M[j,] - M[i,]) <= r) && 
	  (PI.plus[i,j] > PI.min[i,j])) PM[i,j] = 1
      }
    }
  } else {
    # for similarity index C2
    for (i in 1:nrow(M)){
      for (j in 1:nrow(M)) {
        if (PI.min[i,j] != 0) 
        {CM[i,j] <- (PI.plus[i,j]) / (PI.min[i,j])}
        else
        {CM[i,j] = 1000 * PI.plus[i,j]} # to avoid dividing by 0
        if((CM[i,j] > Lambda) && ((M[j,] - M[i,]) <= r) && 
	   (PI.plus[i,j] > PI.min[i,j])) {PM[i,j] = 1}
      }
    }
  }
  for (i in 1:nrow(PM)) PM[i,i] = 0
  class(PM) <- 'prefM'
  PM
}


## ----eval=FALSE---------------------------------------------------------------
## list(PI.plus = PI.plus, PI.min = PI.min, PI.indif = PI.indif)


## ----eval=FALSE---------------------------------------------------------------
## # If we did not push the values PI.plus, PI.min, and
## # PI.indif into the environment of this functions, we would write:
## X <- mcda_electre(M,w)
## # and then address X as in the following code as follows:
## X$PI.min[i,j]


## ----fig.cap='The preference structure as found by the ELECTRE I method given all parameters in the code.\\label{fig:electre1}'----
# the critia: "tlnt" "stab" "cost" "infl" "trvl" "infr" "life"
w <- c(        0.125, 0.2,   0.2,   0.2,  0.175,  0.05,  0.05)
w <- w / sum(w)  # the sum was 1 already, but just to be sure.
r  <- c(0.3,    0.5,  0.5,   0.5,   1,     0.9,   0.5)

eM <- mcda_electre1(M1, w, Lambda=0.6, r=r)
print(eM)
plot(eM)


## ----electre1C2,fig.cap='The results of ELECTRE I with comparability index C2 and parameters as in the aforementioned code.'----
# the critia: "tlnt" "stab" "cost" "infl" "trvl" "infr" "life"
w <- c(        0.125, 0.2,   0.2,   0.2,  0.175,  0.05,  0.05)
w <- w / sum(w)  # the sum was 1 already, but just to be sure.
r  <- c(0.3,    0.5,  0.5,   0.5,   1,     0.9,   0.5)

eM <- mcda_electre1(M1, w, Lambda=1.25, r=r, index='C2')
plot(eM)


## ----fig.cap='The preference structure as found by the ELECTRE II method given all parameters in the code.\\label{fig:electre2}'----
# the critia: "tlnt" "stab" "cost" "infl" "trvl" "infr" "life"
w <- c(        0.125, 0.2,   0.2,   0.2,  0.175,  0.05,  0.05)
w <- w / sum(w)  # the sum was 1 already, but just to be sure.
r  <- c(1,    1,  1,   1,   1,     1,   1)

eM <- mcda_electre1(M1, w, Lambda=0.0, r=r)
print(eM)
plot(eM)


## -----------------------------------------------------------------------------
mcda_electre2 <- function (M1, w) {
  r <- rep(1L, ncol(M))
  mcda_electre1(M1, w, Lambda=0.0, r=r)
  }


## ----eval=FALSE---------------------------------------------------------------
## sum(rowSums(prefM)) == A * A - A
## # with prefM the preference matrix and
## #      A the number of alternatives.


## ----fig.cap="Examples of smooth transition schemes for preference functions $\\pi(d)$. The ``d'' is to be understood as the difference in score for a given criterion.\\label{fig:pref2}", echo=FALSE----
library(ggplot2)
library(latex2exp)
d <- seq(from = -3, to = +3, length.out = 100)

## error function family:
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
# (see Abramowitz and Stegun 29.2.29)
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)
 
## Gudermannian function
gd <- function(x) asin(tanh(x))

f1 <- function(x) erf( sqrt(pi) / 2 * x)
f2 <- function(x) tanh(x)
f3 <- function(x) 2 / pi * gd(pi / 2 * x)
f4 <- function(x) x / sqrt(1 + x^2)
f5 <- function(x) 2 / pi * atan(pi / 2 * x)
f6 <- function(x) x / (1 + abs(x))

df <- data.frame(d = d, y = f1(d), Function = "erf( sqrt(pi) / 2 * d)")
df <- rbind(df, data.frame(d = d, y = f2(d), Function = "tanh(d)"))
df <- rbind(df, data.frame(d = d, y = f3(d), Function = "2 / pi * gd(pi / 2 * d)"))
df <- rbind(df, data.frame(d = d, y = f4(d), Function = "d / (1 + d^2)"))
df <- rbind(df, data.frame(d = d, y = f5(d), Function = "2 / pi * atan(pi / 2 * d)"))
df <- rbind(df, data.frame(d = d, y = f6(d), Function = "x / (1 + abs(d))"))

fn <- ""
fn[1] <- "erf \\left(\\frac{\\sqrt{\\pi} d}{2}\\right)"
fn[2] <- "tanh(x)"
fn[3] <- "\\frac{2}{\\pi} gd\\left( \\frac{\\pi d}{2} \\right)"
fn[4] <- "\\frac{d}{1 + d^2}"
fn[5] <- "\\frac{2}{\\pi} atan\\left(\\frac{\\pi d}{2}\\right)"
fn[6] <- "\\frac{x}{1+ |x|}"


ggplot(data = df, aes(x = d, y = y, color = Function)) +
   geom_line(aes(col = Function), lwd=2) +
   guides(color=guide_legend(title=NULL)) +
   scale_color_discrete(labels=lapply(sprintf('$\\pi(d) = %s$', fn), TeX)) + 
   theme(legend.justification = c(1, 0), 
        legend.position = c(1, 0),   # south east
        legend.box.margin=ggplot2::margin(rep(20, times=4)),
        legend.key.size = unit(1.5, "cm")  # increase vertical space between legend items
	) +
   ylab(TeX('Preference --- $\\pi$'))


## ----fig.cap="Examples of practically applicable preferences functions $P(d)$. The ``x'' is to be understood as the difference in score for a given criterion. Note that all scaling factors are optimized for a normalized decision matrix. The function {\\tt gaus()} refers to {\\tt $exp (-(x-0)^2 / 0.5)$}.\\label{fig:pref}", fig.height='0.9\\textheight', out.height='0.9\\textheight', fig.width='0.99\\textwidth', out.width='0.99\\textwidth', echo=FALSE----
f_curve <- function(f) {
  g <- Vectorize(f)
  s <- deparse(f)[2]
   curve(g, xlab = '', ylab = '', col = 'red', lwd = 3, 
         from = -1, to = +1,
         main = bquote(bold(.(s)))
#          main = s
	  )
   }

gaus <- function(x) exp (-(x-0)^2 / (0.5)^2)
f1 <- function(x) ifelse(x<0, 0, - 3/2 * x^5 + 5/2 * x^3)
f2 <- function(x) ifelse(x<0, 0, sin(pi * x / 2))
f3 <- function(x) min(1, max(1.5*x-0.2, 0))
f4 <- function(x) ifelse(x<0, 0, x)
f5 <- function(x) ifelse(x < 0, 0, 1 - gaus(x))
f6 <- function(x) 1+tanh(6*(x-0.6))

par(mfrow=c(3,2))
f_curve(f1)
f_curve(f2)
f_curve(f3)
f_curve(f4)
f_curve(f5)
f_curve(f6)
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
# mcda_promethee
# delivers the preference flow matrices for the Promethee method
# Arguments:
#    M      -- decision matrix
#    w      -- weights
#    piFUNs -- a list of preference functions, 
#              if not provided min(1,max(0,d)) is assumed.
# Returns (as side effect)
# phi_plus <<- rowSums(PI.plus)
# phi_min  <<- rowSums(PI.min)
# phi_     <<- phi_plus - phi_min
#
mcda_promethee <- function(M, w, piFUNs='x')
{
  if (piFUNs == 'x') {
       # create a factory function:
       makeFUN <- function(x) {x; function(x) max(0,x) }
       P <- list()
       for (k in 1:ncol(M)) P[[k]] <- makeFUN(k)
       } # in all other cases we assume a vector of functions
# initializations
PI.plus  <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
PI.min   <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
# calculate the preference matrix
for (i in 1:nrow(M)){
  for (j in 1:nrow(M)) {
    for (k in 1:ncol(M)) {
      if (M[i,k] > M[j,k]) {
        PI.plus[i,j] = PI.plus[i,j] + w[k] * P[[k]](M[i,k] - M[j,k])
      }
      if (M[j,k] > M[i,k]) {
        PI.min[i,j] = PI.min[i,j] + w[k] * P[[k]](M[j,k] - M[i,k])
      }
    }
  }
}
# note the <<- which pushes the results to the upwards environment
phi_plus <<- rowSums(PI.plus)
phi_min  <<- rowSums(PI.min)
phi_     <<- phi_plus - phi_min
}


## -----------------------------------------------------------------------------
# mcda_promethee1
# Calculates the preference matrix for the Promethee1 method
# Arguments:
#    M      -- decision matrix
#    w      -- weights
#    piFUNs -- a list of preference functions, 
#              if not provided min(1,max(0,d)) is assumed.
# Returns:
#    prefM object -- the preference matrix
#
mcda_promethee1 <- function(M, w, piFUNs='x') {
  # mcda_promethee adds phi_min, phi_plus & phi_ to this environment:
  mcda_promethee(M, w, piFUNs='x') 
  
  # Now, calculate the preference relations:
  pref     <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
    for (i in 1:nrow(M)){
      for (j in 1:nrow(M)) {
        if (phi_plus[i] == phi_plus[j] && phi_min[i]==phi_min[j]) {
	    pref[i,j] <- 0
	  }
	  else if ((phi_plus[i] > phi_plus[j] && 
	            phi_min[i] < phi_min[j] ) || 
	          (phi_plus[i] >= phi_plus[j] && 
		    phi_min[i] < phi_min[j] )) {
	    pref[i,j] <- 1
	  }
	  else {
	    pref[i,j] = NA
	  }
	}
      }
  rownames(pref) <- colnames(pref) <- rownames(M)
  class(pref)    <- 'prefM'
  pref
}


## -----------------------------------------------------------------------------
# We reuse the decision matrix M1 and weights w as defined above.
m <- mcda_promethee1(M1, w)


## ----plotProm1,fig.cap='The hierarchy between alternatives as found by PROMethEE I.'----
# We reuse the decision matrix M1 and weights w as defined above.
m <- mcda_promethee1(M1, w)
plot(m)


## ----promethe1FUN,fig.cap='The result for PROMethEE I with different preference functions provided.'----
# Make shortcuts for some of the functions that we will use:
gauss_val <- function(d) 1 - exp(-(d - 0.1)^2 / (2 * 0.5^2))
x         <- function(d) max(0,d)
minmax    <- function(d) min(1, max(0,2*(d-0.5)))
step      <- function(d) ifelse(d > 0.5, 1,0)

# Create a list of 7 functions (one per criterion):
f <- list()
f[[1]] <- gauss_val
f[[2]] <- x
f[[3]] <- x
f[[4]] <- gauss_val
f[[5]] <- step
f[[6]] <- x
f[[7]] <- minmax

# Use the functions in mcda_promethee1:
m <- mcda_promethee1(M1, w, f)

# Plot the results:
plot(m)


## ----fig.cap="Promethee II can also be seen as using a richer preference structure that can become negative. Here are some examples of practically applicable preferences functions $P(d)$. The function {\\tt gaus()} refers to {\\tt $exp (-(x-0)^2 / 0.5)$}.\\label{fig:pref:MCDA}", fig.height='0.9\\textheight', out.height='0.9\\textheight', fig.width='0.99\\textwidth', out.width='0.99\\textwidth', echo=FALSE----
library(latex2exp)
f_curve <- function(f) {
  g <- Vectorize(f)
   curve(g, xlab = '', ylab = '', col = 'red', lwd = 3, 
         from = -1, to = +1,
         main = TeX(toString(deparse(f)[2])))
   }

gaus <- function(x) exp (-(x-0)^2 / 0.5)
f1 <- function(x) - 3/2 * x^5 + 5/2 * x^3
f2 <- function(x) sin(pi * x / 2)
f3 <- function(x) min(1, max(2*x, -1))
f4 <- function(x) x
f5 <- function(x) ifelse(x < 0 , gaus(x) - 1, 1 - gaus(x))
f6 <- function(x) tanh(3 * x)

par(mfrow=c(3,2))
f_curve(f1)
f_curve(f2)
f_curve(f3)
f_curve(f4)
f_curve(f5)
f_curve(f6)
par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
# mcda_promethee2
# Calculates the Promethee2 preference matrix
# Arguments:
#    M      -- decision matrix
#    w      -- weights
#    piFUNs -- a list of preference functions, 
#              if not provided min(1,max(0,d)) is assumed.
# Returns:
#    prefM object -- the preference matrix
#
mcda_promethee2 <- function(M, w, piFUNs='x')
 { # promethee II
  mcda_promethee(M, w, piFUNs='x')
  pref     <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
    for (i in 1:nrow(M)){
      for (j in 1:nrow(M)) {
          pref[i,j] <- max(phi_[i] - phi_[j],0)
        }
      }
rownames(pref) <- colnames(pref) <- rownames(M)
class(pref) <- 'prefM'
pref
}


## ----fig.cap='The hierarchy between alternatives as found by PROMethEE II. The thickness of the lines corresponds to the strength of the preference.\\label{fig:prometheeII}'----
m <- mcda_promethee2(M1, w)
plot(m)


## ----prom2Bar,fig.cap='PROMethEE II provides a full ranking. Here we show how much each alternative is preferable over its competitors. The size of the blocks is relative to the amount of preference over the other alternative.'----
# We can consider the rowSums as a "score".
rowSums(m)

# So, consider the prefM as a score-matrix (scoreM):
plot.scoreM(m)


## ----gaia,fig.cap=c('The variance explained by each principal component.', 'A projection of the space of alternatives in the 2D-plane formed by the two most dominating principal components.')----
pca1 <- prcomp(M1) 
summary(pca1)

# plot for the prcomp object shows the variance explained by each PC
plot(pca1, type = 'l')

# biplot shows a projection in the 2D plane (PC1, PC2)
biplot(pca1)


## ----gaiaGG,fig.cap=c('A standard plot with {\\tt autoplot()} with labels coloured','Autoplot with visualization of two clusters')----
library(ggplot2)
library(ggfortify)
library(cluster)

# Autoplot with labels colored
autoplot(pca1, data=M1, label=TRUE, shape=FALSE, colour='cost',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3
         )

# Autoplot with visualization of 2 clusters
autoplot(fanny(M1,2), label=TRUE, frame=TRUE, shape=FALSE,
         loadings=TRUE,  loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


## ----fig.cap="Clustering with elliptoid borders, labels of alternative, projections of the criteria and a ``decision vector'' -- the projection of the weights -- constitute a ``Gaia-plot.''\\label{fig:gaia}"----
# Use the weights as defined above:
w

# Calculate coordinates
dv1 <- sum( w * pca1$rotation[,1])  # decision vector PC1 component
dv2 <- sum( w * pca1$rotation[,2])  # decision vector PC2 component

p <- autoplot(pam(M1,2), frame=TRUE, frame.type='norm', label=TRUE, 
         shape=FALSE,
         label.colour='blue',label.face='bold', label.size=6,
         loadings=TRUE,  loadings.colour = 'dodgerblue4',
         loadings.label = TRUE, loadings.label.size = 6, 
	 loadings.label.colour='dodgerblue4',
         loadings.label.vjust = 1.2, loadings.label.hjust = 1.3
         )
p <- p + scale_y_continuous(breaks = 
                        round(seq(from = -1, to = +1, by = 0.2), 2))
p <- p + scale_x_continuous(breaks = 
                        round(seq(from = -1, to = +1, by = 0.2), 2))
p <- p + geom_segment(aes(x=0, y=0, xend=dv1, yend=dv2), size = 2,
                      arrow = arrow(length = unit(0.5, "cm")))
p <- p + ggplot2::annotate("text", x = dv1+0.2, y = dv2-0.01, 
                  label = "decision vector",
                  colour = "black", fontface =2)
p


## -----------------------------------------------------------------------------
### Outrank
# M is the decision matrix (formulated for a maximum problem)
# w the weights to be used for each rank
# method = [order | invorder | median | average ]
outrank <- function (M, w, method='order')
{
  order      <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  order.inv  <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  order.pref <- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
  
  for (i in 1:nrow(M)){
    for (j in 1:nrow(M)) {
      for (k in 1:ncol(M)) {
        if (M[i,k] > M[j,k]) { order[i,j] = order[i,j] + w[k] }
        if (M[j,k] > M[i,k]) { order.inv[i,j] = order.inv[i,j] + w[k] }
      }
    }
  }
  for (i in 1:nrow(M)){
    for (j in 1:nrow(M)) {
      if (order[i,j] > order[j,i]){
        order.pref[i,j] = 1
        order.pref[j,i] = 0
      }
      else if (order[i,j] < order[j,i]) {
        order.pref[i,j] = 0
        order.pref[j,i] = 1
      }
      else {
        order.pref[i,j] = 0
        order.pref[j,i] = 0
      }
    }
  }
 class(order.pref) <- 'prefM'
 order.pref
}


## -----------------------------------------------------------------------------
r_w <- 0.05              # the interest rate per week 
r_y <- (1 + r_w)^52 - 1  # assume 52 weeks per year
paste0("the APR is: ", round(r_y * 100, 2),"%!")


## -----------------------------------------------------------------------------
r   <- 0.1
CFs <- c(-100, 100, 100)
t   <- c(   0,   5,   7)
NPV <- sum(CFs / (1 + r)^t)
print(round(NPV, 2))


## -----------------------------------------------------------------------------
# bond_value
# Calculates the fair value of a bond
# Arguments:
#    time_to_mat -- time to maturity in years
#    coupon      -- annual coupon in $
#    disc_rate   -- discount rate (risk free + risk premium)
#    nominal     -- face value of the bond in $
# Returns:
#    the value of the bond in $
bond_value <- function(time_to_mat, coupon, disc_rate, nominal){
  value  <- 0
  # 1/ all coupons
  for (t in 1:time_to_mat) {
     value <- value + coupon * (1 + disc_rate)^(-t)
     }
  # 2/ end payment of face value
  value <- value + nominal * (1 + disc_rate)^(-time_to_mat)
  value
}

# We assume that the required interest rate is the 
# risk free interest rate of 3%.

# The fair value of the bond is then:
bond_value(time_to_mat = 5, coupon = 5, disc_rate = 0.03, 
           nominal = 100)


## -----------------------------------------------------------------------------
bond_value(time_to_mat = 5, coupon = 5, disc_rate = 0.035, 
           nominal = 100)


## -----------------------------------------------------------------------------
V <- bond_value(time_to_mat = 5, coupon = 5, disc_rate = 0.03, 
           nominal = 100)
CFs <- c(seq(5, 5, length.out=4), 105)
t   <- c(1:5)
r   <- 0.03
MacD <- 1/V * sum(t * CFs / (1 + r)^t)
print(MacD)


## -----------------------------------------------------------------------------
R_A <- 0.02 + 1.25 * (0.10 - 0.02)
print(paste0('The RR for company A is: ', 
             round(R_A, 2) * 100, '%'))


## -----------------------------------------------------------------------------
R_B <- 0.02 + 0.75 * (0.10 - 0.02)
print(paste0('The RR for B is: ', round(R_B, 2) * 100, '%'))
print(paste0('The beta changed by ',
             round((0.75 / 1.25 - 1) * 100, 2), 
             '% and the RR by ', 
	     round((R_B / R_A - 1) * 100, 2), '%.'))


## -----------------------------------------------------------------------------
V_0 <- 10 * (1 + 0.00) / (0.01 + 0.05 - 0.00)
print(round(V_0,2))


## -----------------------------------------------------------------------------
V_0 <- 10 * (1 + 0.02) / (0.01 + 0.05 -0.02)
print(round(V_0,2))


## -----------------------------------------------------------------------------
V_0 <- 10 * (1 + 0.02) / (0.01 + 1.5 * 0.05 - 0.02)
print(round(V_0,2))


## -----------------------------------------------------------------------------
V_0 <- 10 * (1 + 0.02) / (0.01 + 1.5 * 0.05 - 0.10)
print(round(V_0,2))


## ----longCall, fig.cap='The intrinsic value of a long call illustrated with its payoff and profit. The profit is lower, since it takes into account that the option buyer has paid a fixed premium for the option.'----
# Let us plot the value of the call in function of the strike
FS <- seq(80, 120, length.out=150) # future spot price
X  <- 100                          # strike
P  <- 5                            # option premium
T  <- 3                            # time to maturity
r  <- 0.03                         # discount rate
payoff <- mapply(max, FS-X, 0)
profit <- payoff - P * (1 + r)^T

# Plot the results:
plot(FS, payoff,
     col='red', lwd=3, type='l',
     main='LONG CALL value at maturity',
     xlab='Future strike price',
     ylab='$',
     ylim=c(-10,20)
     )
lines(FS, profit,
      col='blue', lwd=2)
text(105,8, 'Payoff', col='red')
text(115,5, 'Profit', col='blue')


## ----shortCall, fig.cap='The intrinsic value of a short call illustrated with its payoff and profit. The profit is higher, since this the position of the option-wirter, so this party has got the premium at the start of the contract. Note that the loss is unlimited.'----
FS <- seq(80, 120, length.out=150) # future spot price
X  <- 100                          # strike
P  <- 5                            # option premium
T  <- 3                            # time to maturity
r  <- 0.03                         # discount rate
payoff <- - mapply(max, FS-X, 0)
profit <- P * (1 + r)^T + payoff

# Plot the results:
plot(FS, payoff,
     col='red', lwd=3, type='l',
     main='SHORT CALL value at maturity',
     xlab='Future spot price',
     ylab='$',
     ylim=c(-20,10)
     )
lines(FS, profit,
      col='blue', lwd=2)
text(90,1.5, 'Payoff', col='red')
text(90,7, 'Profit', col='blue')


## ----longShortPut, fig.cap='The payoff and profit for a long put (left) and a short put (right).'----
FS <- seq(80, 120, length.out=150) # future spot price
X  <- 100                          # strike
P  <- 5                            # option premium
T  <- 3                            # time to maturity
r  <- 0.03                         # discount rate

# the long put:
payoff <- mapply(max, X - FS, 0)
profit <- payoff - P * (1 + r)^T

par(mfrow=c(1,2))

plot(FS, payoff,
     col='red', lwd=3, type='l',
     main='LONG PUT at maturity',
     xlab='Future spot price',
     ylab='$',
     ylim=c(-20,20)
     )
lines(FS, profit,
      col='blue', lwd=2)
text(110,1,  'Payoff', col='red')
text(110,-4, 'Profit', col='blue')

# the short put:
payoff <- - mapply(max, X - FS, 0)
profit <- payoff + P * (1 + r)^T

plot(FS, payoff,
     col='red', lwd=3, type='l',
     main='SHORT PUT at maturity',
     xlab='Future spot price',
     ylab='',
     ylim=c(-20,20)
     )
lines(FS, profit,
      col='blue', lwd=2)
text(110,1, 'Payoff', col='red')
text(110,6, 'Profit', col='blue')

par(mfrow=c(1,1))  # reset the plot interface


## -----------------------------------------------------------------------------
# call_intrinsicVal
# Calculates the intrinsic value for a call option
# Arguments:
#    Spot   -- numeric -- spot price
#    Strike -- numeric -- the strike price of the option
# Returns
#    numeric -- intrinsic value of the call option.
call_intrinsicVal <- function(Spot, Strike) {max(Spot - Strike, 0)}



# put_intrinsicVal
# Calculates the intrinsic value for a put option
# Arguments:
#    Spot   -- numeric -- spot price
#    Strike -- numeric -- the strike price of the option
# Returns
#    numeric -- intrinsic value of the put option.
put_intrinsicVal <- function(Spot, Strike) {max(-Spot + Strike, 0)}



# call_price
# The B&S price of a call option before maturity
# Arguments:
#    Spot   -- numeric -- spot price in $ or %
#    Strike -- numeric -- the strike price of the option  in $ or %
#    T      -- numeric -- time to maturity in years
#    r      -- numeric -- interest rates (e.g. 0.02 = 2%)
#    vol    -- numeric -- standard deviation of underlying in $ or %
# Returns
#    numeric -- value of the call option in $ or %
#
call_price <- function (Spot, Strike, T, r, vol)
 {
  d1 <- (log(Spot / Strike) + (r + vol ^ 2/2) * T) / (vol * sqrt(T))
  d2 <- (log(Spot / Strike) + (r - vol ^ 2/2) * T) / (vol * sqrt(T))
  pnorm(d1) * Spot - pnorm(d2) * Strike * exp(-r * T)
  }

# put_price
# The B&S price of a put option before maturity
# Arguments:
#    Spot   -- numeric -- spot price in $ or %
#    Strike -- numeric -- the strike price of the option  in $ or %
#    T      -- numeric -- time to maturity in years
#    r      -- numeric -- interest rates (e.g. 0.02 = 2%)
#    vol    -- numeric -- standard deviation of underlying in $ or %
# Returns
#    numeric -- value of the put option in $ or %
#
put_price <- function(Spot, Strike, T, r, vol)
 {
 Strike * exp(-r * T) - Spot + call_price(Spot, Strike, T, r, vol)
 }


## -----------------------------------------------------------------------------
# Examples:
call_price (Spot = 100, Strike = 100, T = 1, r = 0.02, vol = 0.2)
put_price  (Spot = 100, Strike = 100, T = 1, r = 0.02, vol = 0.2)


## ----fig.cap='The price of a long call compared to its intrinsic value. The market value is always positive.\\label{fig:call}'----
# Long call
spot <- seq(50,150, length.out=150)
intrinsic_value_call <- apply(as.data.frame(spot), 
                              MARGIN=1, 
			      FUN=call_intrinsicVal, 
			      Strike=100)
market_value_call    <- call_price(Spot = spot, Strike = 100, 
                                   T = 3, r = 0.03, vol = 0.2)
plot(spot, market_value_call,
     type = 'l', col= 'red', lwd = 4,
     main = 'European Call option',
     xlab = 'Spot price',
     ylab = 'Option value')
text(115, 40, 'Market value', col='red')
lines(spot, intrinsic_value_call,
      col= 'forestgreen', lwd = 4)
text(130,15, 'Intrinsic value', col='forestgreen')


## ----fig.cap='The price of a long put compared to its intrinsic value. Note that the market price of a put can be lower than its intrinsic value. This is because of the cost of carry.\\label{fig:put}'----
# Long put
spot <- seq(50,150, length.out=150)
intrinsic_value_put <- apply(as.data.frame(spot), 
                             MARGIN=1, 
			     FUN=put_intrinsicVal, 
			     Strike=100)
market_value_put    <- put_price(Spot = spot, Strike = 100, 
                                 T = 3, r = 0.03, vol = 0.2)
plot(spot, market_value_put,
     type = 'l', col= 'red', lwd = 4,
     main = 'European Put option',
     xlab = 'Spot price',
     ylab = 'Option value')
text(120, 8, 'market value', col='red')
lines(spot, intrinsic_value_put,
      col= 'forestgreen', lwd = 4)
text(75,10, 'intrinsic value', col='forestgreen')


## -----------------------------------------------------------------------------
# CRR_price
# Calculates the CRR binomial model for an option
# Arguments:
#         S0 -- numeric -- spot price today (start value)
#         SX -- numeric -- strike, e.g. 100
#      sigma -- numeric -- the volatility over the maturity period, 
#                          e.g. ca. 0.2 for shares on 1 yr
#        Rrf -- numeric -- the risk free interest rate (log-return)
# optionType -- character -- 'lookback' for lookback option, 
#                            otherwise vanilla call is assumed
#    maxIter -- numeric -- number of iterations 
# Returns:
#    numeric -- the value of the option given parameters above
CRR_price <- function(S0, SX, sigma, Rrf, optionType, maxIter)
  {
  Svals <- mat.or.vec(2^(maxIter), maxIter+1)
  probs <- mat.or.vec(2^(maxIter), maxIter+1)
  Smax  <- mat.or.vec(2^(maxIter), maxIter+1)
  Svals[1,1] <- S0
  probs[1,1] <- 1
  Smax[1,1]  <- S0
  dt <-  1 / maxIter
  u <- exp(sigma * sqrt(dt))
  d <- exp(-sigma * sqrt(dt))
  p = (exp(Rrf * dt) - d) / (u - d)
  for (n in 1:(maxIter))
    {
    for (m in 1:2^(n-1))
     {
     Svals[2*m-1,n+1] <- Svals[m,n] * u
     Svals[2*m,n+1]   <- Svals[m,n] * d
     probs[2*m-1,n+1] <- probs[m,n] * p
     probs[2*m,n+1]   <- probs[m,n] * (1 - p)
     Smax[2*m-1,n+1]  <- max(Smax[m,n], Svals[2*m-1,n+1])
     Smax[2*m,n+1]    <- max(Smax[m,n], Svals[2*m,n+1])
     }
    }
  if (optionType == 'lookback')
   {
     exp.payoff <- (Smax - SX)[,maxIter + 1]  * probs[,maxIter + 1]
   }  # lookback call option
   else
    {
     optVal <- sapply(Svals[,maxIter + 1] - SX,max,0)
     exp.payoff <- optVal * probs[,maxIter + 1]
    }  # vanilla call option
  sum(exp.payoff) / (1 + Rrf)
  }


## -----------------------------------------------------------------------------
# plot_CRR
# This function will call the CRR function iteratively for 
# number of iterations increasing from 1 to maxIter and 
# plot the results on screen (or if desired uncomment the 
# relevant lines to save to disk).
# Arguments:
# optionType -- character -- 'lookback' for lookback option, 
#                            otherwise vanilla call is assumed
#    maxIter -- numeric -- maximal number of iterations 
#   saveFile -- boolean -- TRUE to save the plot as pdf
# Returns:
#    numeric -- the value of the option given parameters above

plot_CRR <- function(optionType, maxIter, saveFile = FALSE) 
  {
  x <- seq(1,maxIter)
  y <- mat.or.vec(maxIter,1)
  for (k in 1:maxIter) {y[k] <- CRR_price(100,
                                          100,
                                          sigma      = 0.2,
					  Rrf        = 0.02,
					  optionType = optionType,
					  maxIter    = k)}
  d <- data.frame(x, y)
  colnames(d) <- c('maxIter', 'value')
  p <- qplot(maxIter, value, data=d, geom = "point",size=I(3) )
  p <- p + geom_line() + ggtitle(optionType)
  p <- p + xlab('Number of iterations' ) + 
       ylab('Estimated value of the option')

  if(saveFile) {
    p
    ggsave(paste('img/binomial_CRR_',optionType,'.pdf',sep=''))
    }
  # Return the plot:
  p
  } 


## ----CRRcall,fig.cap='The Cox--Ross--Rubinstein model for the binomial model applied to a call option. Note how the process converges smooth and quick.'----
library(ggplot2)
# Plot the convergence of the CRR algorithm for a call option.
plot_CRR("Call", maxIter = 20)


## ----CRRrus,fig.cap='The Cox--Ross--Rubinstein model for the binomial model applied to an unlimited look-back option (aka Russian Option). For the Russian option, the convergence process is not obvious. In fact the more steps we allow, the more expensive the option becomes, because the more there are moments that can have a lower price.'----
# Plot the convergence of the CRR algorithm for a call option.
plot_CRR("lookback", maxIter = 15)


## -----------------------------------------------------------------------------
# We still use ggplot2
library(ggplot2)


## -----------------------------------------------------------------------------
# plot_price_evol
# Plots the evolution of Call price in function of a given variable
# Arguments:
#    var       -- numeric   -- vector of values of the variable
#    varName   -- character -- name of the variable to be studied
#    price     -- numeric   -- vector of prices of the option
#    priceName -- character -- the name of the option
#    reverseX  -- boolean   -- TRUE to plot x-axis from high to low
# Returns
#    ggplot2 plot
#
plot_price_evol <- function(var, varName, price, priceName, 
                            reverseX = FALSE)
{
  d <- data.frame(var, price)
  colnames(d) <- c('x','y')
  p <- qplot(x, y, data=d, geom = "line",size=I(2) )
  p <- p + geom_line()
  if (reverseX) {p <- p + xlim(max(var),min(var))}  # reverse axis
  p <- p + xlab(varName ) + ylab(priceName)
  p   # return the plot
}


## ----callDependencies,fig.cap='The value of a call option depends on many variables. Some are illustrated in these plots.',fig.height='0.9\\textheight',out.height='0.9\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# Define the default values:
t      <- 1
Spot   <- 100
Strike <- 100
r      <- log(1 + 0.03)
vol    <- 0.2


## ... time
T <- seq(5, 0.0001, -0.01)
Call <- c(call_price (Spot, Strike, T, r, vol))
p1 <- plot_price_evol(T, "Time to maturity (years)", Call, "Call", 
                      TRUE)

## ... interest
R <- seq(0.001, 0.3, 0.001)
Call <- c(call_price (Spot, Strike, t, R, vol))
p2 <- plot_price_evol(R, "Interest rate", Call, "Call")

## ... volatility
vol <- seq(0.00, 0.2, 0.001)
Call <- c(call_price (Spot, Strike, t, r, vol))
p3 <- plot_price_evol(vol, "Volatility", Call, "Call")

## ... strike
X <- seq(0, 200, 1)
Call <- c(call_price (Spot, X, t, r, vol))
p4 <- plot_price_evol(X, "Strike", Call, "Call")

## ... Spot
spot <- seq(0, 200, 1)
Call <- c(call_price (spot, Strike, t, r, vol))
p5 <- plot_price_evol(spot, "Spot price", Call, "Call")

# In the next line we use the function grid.arrange()
# from the gridExtra package
library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 3)


## ----putDependencies,fig.cap='The value of a put option depends on many variables. Some are illustrated in these plots.',fig.height='0.9\\textheight',out.height='0.9\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# Define the default values:
t      <- 1
Spot   <- 100
Strike <- 100
r      <- log(1 + 0.03)
vol    <- 0.2


## ... time
T <- seq(5, 0.0001, -0.01)
Call <- c(put_price (Spot, Strike, T, r, vol))
p1 <- plot_price_evol(T, "Time to maturity (years)", 
                      Call, "Call", TRUE)

## ... interest
R <- seq(0.001, 0.3, 0.001)
Call <- c(put_price (Spot, Strike, t, R, vol))
p2 <- plot_price_evol(R, "Interest rate", Call, "Call")

## ... volatility
vol <- seq(0.00, 0.2, 0.001)
Call <- c(put_price (Spot, Strike, t, r, vol))
p3 <- plot_price_evol(vol, "Volatility", Call, "Call")

## ... strike
X <- seq(0, 200, 1)
Call <- c(put_price (Spot, X, t, r, vol))
p4 <- plot_price_evol(X, "Strike", Call, "Call")

## ... Spot
spot <- seq(0, 200, 1)
Call <- c(put_price (spot, Strike, t, r, vol))
p5 <- plot_price_evol(spot, "Spot price", Call, "Call")

# In the next line we use the function grid.arrange()
# from the gridExtra package
library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 3)


## ----deltaCP,fig.cap='An illustration of how the delta of a call and put compare in function of the spot price. Note the difference in scale on the $y$-axis.',fig.height='0.9\\textheight',out.height='0.9\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# Define the functions to calculate the price of the delta

# call_delta
# Calculates the delta of a call option
# Arguments:
#    S      -- numeric -- spot price
#    Strike -- numeric -- strike price
#    T      -- numeric -- time to maturity
#    r      -- numeric -- interest rate
#    vol    -- numeric -- standard deviation of underlying
call_delta  <- function (S, Strike, T, r, vol)
 {
  d1 <- (log (S / Strike)+(r + vol ^2 / 2) * T) / (vol * sqrt(T))
  pnorm(d1)
  }


# put_delta
# Calculates the delta of a put option
# Arguments:
#    S      -- numeric -- spot price
#    Strike -- numeric -- strike price
#    T      -- numeric -- time to maturity
#    r      -- numeric -- interest rate
#    vol    -- numeric -- standard deviation of underlying
put_delta  <- function (S, Strike, T, r, vol)
 {
  d1 <- (log (S / Strike)+(r + vol ^2 / 2) * T) / (vol * sqrt(T))
  pnorm(d1) - 1
  }


## DELTA CALL
spot <- seq(0,200, 1)
delta <- c(call_delta(spot, Strike, t, r, vol))
p1 <- plot_price_evol(spot, "Spot price", delta, "Call delta")

## DELTA PUT
spot <- seq(0,200, 1)
delta <- c(put_delta(spot, Strike, t, r, vol))
p2 <- plot_price_evol(spot, "Spot price", delta, "Put delta")

# plot the two visualizations:
grid.arrange(p1, p2, nrow = 2)


## -----------------------------------------------------------------------------
# load the plotting library
library(ggplot2)


## -----------------------------------------------------------------------------
# portfolio_plot
# Produces a plot of a portfolio of the value in function of the 
# spot price of the underlying asset.
# Arguments:
#   portf            - data.frame - composition of the portfolio
#                   with one row per option, structured as follows:
#                       - ['long', 'short'] - position
#                       - ['call', 'put']   - option type
#                       - numeric           - strike
#                       - numeric           - gearing (1 = 100%)
#   structureName="" - character - label of the portfolio
#   T = 1            - numeric   - time to maturity (in years)
#   r = log(1 + 0.02) - numeric   - interest rate (per year) 
#                                    as log-return
#   vol = 0.2        - numeric   - annual volatility of underlying
#   spot.min = NULL  - NULL for automatic scaling x-axis, value for min
#   spot.max = NULL  - NULL for automatic scaling x-axis, value for max
#   legendPos=c(.25,0.6) - numeric vector - set to 'none' to turn off
#   yLims = NULL     - numeric vector - limits y-axis, e.g. c(80, 120)
#   fileName = NULL  - character - filename, NULL for no saving
#   xlab = "default" - character - x axis label, NULL to turn off
#   ylab = "default" - character - y axis label, NULL to turn off
# Returns (as side effect)
#   ggplot plot
#   pdf file of this plot (in subdirectory ./img/)

portfolio_plot <- function(portf, 
	   structureName="", # name of the option strategy
	   T = 1,            # time to maturity (in years)
	   r = log(1 + 0.02),# interest rate (per year)
	   vol = 0.2,        # annual volatility of the underlying
	   spot.min = NULL,  # NULL for automatic scaling x-axis
	   spot.max = NULL,  # NULL for automatic scaling x-axis
	   legendPos=c(.25,0.6), # set to 'none' to turn off
	   yLims = NULL,     # limits of y-axis, e.g. c(80, 120)
	   fileName = NULL,  # NULL for no saving plot into file
	   xlab = "default", # set to NULL to turn off (or string)
	   ylab = "default"  # set to NULL to turn off (or string)
	   ) {
# portf = data frame with: long/short, call/put, strike, gearing

the_S = 100       # The spot price today is always 100
# 
  strikes <- as.numeric(portf[,3])
  strike.min <- min(strikes)
  strike.max <- max(strikes)
  if (is.null(spot.min)) {
    spot.min <- min(0.8*strike.min, max(0,2*strike.min - strike.max))
    }
  if (is.null(spot.max)) {
    spot.max <- max(1.2 * strike.max, 2 * strike.max - strike.min)}
  if (structureName == ""){
    structureName<- paste(deparse(substitute(fileName)), 
                          collapse = "", sep="")
    }
  nbrObs   <- 200
  spot     <- seq(spot.min,spot.max,len=nbrObs)
  val.now  <- seq(0,0,len=nbrObs)
  val.end  <- seq(0,0,len=nbrObs)
  for (k in 1:nrow(portf))
    {
     Strike  <- as.numeric(portf[k,3])
     gearing <- as.numeric(portf[k,4])
     if (portf[k,1] == 'long'){theSign <- 1}else{theSign = -1}
     if (portf[k,2] == 'call')
       {
        purchasePrice <- call_price(the_S, Strike, T, r, vol)
        callVal  <- sapply(spot, call_price, Strike=Strike, T=T, 
	                   r=r, vol=vol)
        val.now.incr <- callVal - purchasePrice 
        val.end.incr <- sapply(spot, call_intrinsicVal, 
	                Strike = Strike) - purchasePrice
       }
       else
       {
        if (portf[k,2] == 'put')
          {
          purchasePrice <- put_price(the_S, Strike, T, r, vol)
          callVal  <- sapply(spot, put_price, Strike=Strike, T=T, 
	                     r=r, vol=vol)
          val.now.incr <- callVal - purchasePrice 
          val.end.incr <- sapply(spot, put_intrinsicVal, 
	                  Strike = Strike) - purchasePrice
          }
          else # then it is 'underlying'
          {
          val.now.incr <- spot - Strike
          val.end.incr <- spot - Strike
          }
       }
     val.now <- val.now + val.now.incr * gearing * theSign
     val.end <- val.end + val.end.incr * gearing * theSign
     }
  d1 <- data.frame(spot, val.end, 
          paste('intrinsic value',structureName,sep=" "), 3)
  d2 <- data.frame(spot, val.now,  
          paste('value 1 year to maturity',structureName,sep=" "), 2)
  colnames(d1) <- c('spot', 'value', 'legend')
  colnames(d2) <- c('spot', 'value', 'legend')
  dd <- rbind(d1,d2)
  p <- qplot(spot, value, data=dd, color = legend, 
             geom = "line",size=I(2) )
  if(is.null(xlab)) {
      p <- p + theme(axis.title.x = element_blank())
      } else { 
	if(xlab == "default") {p <- p + xlab('spot price')
	   } else {p <- p + xlab(xlab)}}
  if(is.null(ylab)) {
      p <- p + theme(axis.title.y = element_blank())
      } else { 
	if(ylab == "default") {p <- p + ylab('Value')
	   } else {p <- p + ylab(ylab)}}
  p <- p + ylab('Value')
  p <- p + theme(legend.position=legendPos)
  if(legendPos == "none") {p <- p + ggtitle(structureName)}
  if (!is.null(yLims)) {p <- p + scale_y_continuous(limits=yLims)}
  if(!is.null(fileName)) {
    # remove punctuation:
    fileName <- str_replace_all(fileName, "[[:punct:]]", "")
    # remove spaces:
    fileName <- str_replace_all(fileName, " ", "")            
    # save file in sub-directory img
    ggsave(paste('img/',fileName,'.pdf',sep=''), 
           width = 6, height = 3)
    }
  # return the plot:
  p
}


## ----fig.cap='Linear option strategies illustrated. The red line is the intrinsic value and the green line is the value of today if the spot price would move away from 100. Part 1 (basic strategies).\\label{fig:lin_opt_strat1}',fig.height='0.9\\textheight',out.height='0.9\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# long call
portfolio <- rbind(c('long','call',100,1))
p1 <- portfolio_plot(portfolio, 'Long call', 
                     legendPos="none", xlab = NULL)

# short call
portfolio <- rbind(c('short','call',100,1))
p2 <- portfolio_plot(portfolio, 'Short call', legendPos="none", 
                     xlab = NULL, ylab = NULL)

# long put
portfolio <- rbind(c('long','put',100,1))
p3 <- portfolio_plot(portfolio, 'Long put', legendPos="none", 
                     xlab=NULL)

# short put
portfolio <- rbind(c('short','put',100,1))
p4 <- portfolio_plot(portfolio, 'Short put', legendPos="none", 
                     xlab = NULL, , ylab = NULL)

# -- long call and short put
portfolio <- rbind(c('long','call',100,1))
portfolio <- rbind(portfolio, c('short','put',100,1))
p5 <- portfolio_plot(portfolio, 'Long call + short put', 
                     legendPos="none", xlab = NULL)

# -- call
portfolio <- rbind(c('long','call',100,1))
p6 <- portfolio_plot(portfolio, 'Call', legendPos="none", 
                     xlab = NULL, ylab = NULL)

# -- put
portfolio <- rbind(c('long','put',100,1))
p7 <- portfolio_plot(portfolio, 'Put', legendPos="none")

# -- callput
portfolio <- rbind(c('short','put',100,1))
portfolio <- rbind(portfolio, c('long','call',100,1))
p8 <- portfolio_plot(portfolio, 'Call + Put', legendPos="none",
                     ylab = NULL)


# show all visualizations:
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)


## ----fig.cap='Linear option strategies illustrated. Part 2 (basic composite structures).\\label{fig:lin_opt_strat2}',fig.height='0.9\\textheight',out.height='0.9\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# -- callspread
portfolio <- rbind(c('short','call',120,1))
portfolio <- rbind(portfolio, c('long','call',100,1))
p1 <- portfolio_plot(portfolio, 'CallSpread',
                     legendPos="none", xlab = NULL)

# -- short callspread
portfolio <- rbind(c('long','call',120,1))
portfolio <- rbind(portfolio, c('short','call',100,1))
p2 <- portfolio_plot(portfolio, 'Short allSpread', 
                     legendPos="none", xlab = NULL, ylab = NULL)

# -- callspread differently
portfolio <- rbind(c('short','put',120,1))
portfolio <- rbind(portfolio, c('long','put',100,1))
p3 <- portfolio_plot(portfolio, 'Short putSpread',
                     legendPos="none", xlab = NULL)

# -- putspread
portfolio <- rbind(c('short','put',80,1))
portfolio <- rbind(portfolio, c('long','put',100,1))
p4 <- portfolio_plot(portfolio, 'PutSpread',
                     legendPos="none", xlab = NULL, ylab = NULL)

# -- straddle
portfolio <- rbind(c('long','call',100,1))
portfolio <- rbind(portfolio, c('long','put',100,1))
p5 <- portfolio_plot(portfolio, 'Straddle', spot.min = 50, 
                     spot.max = 150,legendPos="none", xlab = NULL)
# Note that our default choices for x-axis range are not suitable 
# for this structure. Hence, we add spot.min and spot.max

# -- short straddle
portfolio <- rbind(c('short','call',100,1))
portfolio <- rbind(portfolio, c('short','put',100,1))
p6 <- portfolio_plot(portfolio, 'Short straddle',spot.min = 50, 
                     spot.max = 150, legendPos="none", 
		     xlab = NULL, ylab = NULL)

# -- strangle
portfolio <- rbind(c('long','call',110,1))
portfolio <- rbind(portfolio, c('long','put',90,1))
p7 <- portfolio_plot(portfolio, 'Strangle', 
                     spot.min = 50, spot.max = 150,
                     legendPos="none", xlab = NULL)

# -- butterfly
portfolio <- rbind(c('long','call',120,1))
portfolio <- rbind(portfolio, c('short','call',100,1))
portfolio <- rbind(portfolio, c('long','put',80,1))
portfolio <- rbind(portfolio, c('short','put',100,1))
p8 <- portfolio_plot(portfolio, 'Butterfly',
                     spot.min = 50, spot.max = 150,
                     legendPos="none", xlab = NULL, ylab = NULL)


# show all visualizations:
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)


## ----fig.cap='Linear option strategies illustrated. Part 3 (some more complex structures and extra one structure that is entirely made up).\\label{fig:lin_opt_strat3}',fig.height='0.75\\textheight',out.height='0.75\\textheight',fig.width='0.99\\textwidth',out.width='0.99\\textwidth'----
# -- condor
portfolio <- rbind(c('long','call',140,1))
portfolio <- rbind(portfolio, c('short','call',120,1))
portfolio <- rbind(portfolio, c('long','put',60,1))
portfolio <- rbind(portfolio, c('short','put',80,1))
p1 <- portfolio_plot(portfolio, 'Condor',spot.min = 40, 
                     spot.max = 160, legendPos="none", 
		     xlab = NULL)

# -- short condor
portfolio <- rbind(c('short','call',140,1))
portfolio <- rbind(portfolio, c('long','call',120,1))
portfolio <- rbind(portfolio, c('short','put',60,1))
portfolio <- rbind(portfolio, c('long','put',80,1))
p2 <- portfolio_plot(portfolio, 'Short Condor',spot.min = 40, 
                     spot.max = 160, legendPos="none", 
		     xlab = NULL, ylab = NULL)

# -- geared call
portfolio <- rbind(c('long','call',100.0,2))
p3 <- portfolio_plot(portfolio, 
                     structureName="Call with a gearing of 2",
		     legendPos="none", xlab = NULL)

# -- nearDigital (approximate a digital option with a geared call)
portfolio <- rbind(c('short','call',100.1,10))
portfolio <- rbind(portfolio, c('long','call',100,10))
p4 <- portfolio_plot(portfolio, 'Near digital',
                     legendPos="none", xlab = NULL, ylab = NULL)

# -- a complex structure:
portfolio <- rbind(c('long','call',110,1))
portfolio <- rbind(portfolio, c('short','call',105,1))
portfolio <- rbind(portfolio, c('short','put',95,1))
portfolio <- rbind(portfolio, c('long','put',90,1))
portfolio <- rbind(portfolio, c('long','put',80,1))
portfolio <- rbind(portfolio, c('long','call',120,1))
portfolio <- rbind(portfolio, c('short','call',125,1))
portfolio <- rbind(portfolio, c('short','put',70,10))
portfolio <- rbind(portfolio, c('short','put',75,1))
portfolio <- rbind(portfolio, c('short','call',130,10))
portfolio <- rbind(portfolio, c('long','call',99,10))
portfolio <- rbind(portfolio, c('short','call',100,10))
portfolio <- rbind(portfolio, c('short','put',100,10))
portfolio <- rbind(portfolio, c('long','put',101,10))
p5 <- portfolio_plot(portfolio, 'Fun',legendPos='none', 
                     spot.min=60, spot.max=140,
		     yLims=c(-0,25))

# show all visualizations:
# Pasing a layout_matrix to the function grid.arrange()
# allows to make the lower plot bigger:
layoutM <- rbind(c(1,2),
                 c(3,4),
                 c(5,5))
grid.arrange(p1, p2, p3, p4, p5, nrow = 3, layout_matrix = layoutM)


## ----coveredCall,fig.cap='A covered call is a short call where the losses are protected by having the underlying asset in the same portfolio.'----
## --- covered call ----

nbrObs     <- 100
the.S      <- 100
the.Strike <- 100
the.r      <- log (1 + 0.03)
the.T      <- 1
the.vol    <- 0.2
Spot.min   = 80
Spot.max   = 120
LegendPos  = c(.5,0.2)
Spot   <- seq(Spot.min,Spot.max,len=nbrObs)
val.end.call <-  - sapply(Spot, call_intrinsicVal, Strike = 100)
call.value <- call_price(the.S, the.Strike, the.T, the.r, the.vol)

d.underlying <- data.frame(Spot, Spot - 100,  'Underlying',   1)
d.shortcall  <- data.frame(Spot, val.end.call,  'Short call', 1)
d.portfolio  <- data.frame(Spot, 
                           Spot + val.end.call + call.value - 100,
			   'portfolio', 1.1)
colnames(d.underlying) <- c('Spot', 'value', 'Legend','size')
colnames(d.shortcall)  <- c('Spot', 'value', 'Legend','size')
colnames(d.portfolio)  <- c('Spot', 'value', 'Legend','size')
dd <- rbind(d.underlying,d.shortcall,d.portfolio)
p <- qplot(Spot, value, data = dd, color = Legend, geom = "line",
           size=size )
p <- p + xlab('Value of the underlying' ) + ylab('Profit at maturity')
p <- p + theme(legend.position = LegendPos)
p <- p + scale_size(guide = 'none')
print(p) 


## ----marriedPut,fig.cap='A married put is a put option combined with the underlying asset.'----
## --- married put ----

LegendPos = c(.8,0.2)
Spot        <- seq(Spot.min,Spot.max,len=nbrObs)
val.end.put <-  sapply(Spot, put_intrinsicVal, Strike = 100)
put.value   <- - put_price(the.S, the.Strike, the.T, the.r, the.vol)

d.underlying <- data.frame(Spot, Spot - 100,  'Underlying', 1)
d.shortput   <- data.frame(Spot, val.end.put,  'Long put',  1)
d.portfolio  <- data.frame(Spot, 
                           Spot + val.end.put + put.value - 100,
			   'portfolio',
			   1.1)
colnames(d.underlying) <- c('Spot', 'value', 'Legend','size')
colnames(d.shortput)   <- c('Spot', 'value', 'Legend','size')
colnames(d.portfolio)  <- c('Spot', 'value', 'Legend','size')
dd <- rbind(d.underlying,d.shortput,d.portfolio)
p  <- qplot(Spot, value, data = dd, color = Legend, geom = "line",
            size = size )
p <- p + xlab('Value of the underlying' ) + ylab('Profit at maturity')
p <- p + theme(legend.position = LegendPos)
p <- p + scale_size(guide = 'none')
print(p) 


## ----collar,fig.cap='A collar is a structure that protects us from strong downwards movements at the cost of a limited upside potential.'----
## --- collar ----

# using the same default values as for the previous code block.
LegendPos = c(.6,0.25)
Spot         <-   seq(Spot.min,Spot.max,len=nbrObs)
val.end.call <- - sapply(Spot, call_intrinsicVal, Strike = 110)
val.end.put  <- + sapply(Spot, put_intrinsicVal, Strike = 95)
call.value   <- call_price(the.S, the.Strike, the.T, the.r, the.vol)
put.value    <- put_price(the.S, the.Strike, the.T, the.r, the.vol)

d.underlying <- data.frame(Spot, Spot - 100,   'Underlying', 1)
d.shortcall  <- data.frame(Spot, val.end.call, 'Short call', 1)
d.longput    <- data.frame(Spot, val.end.put,  'Long call',  1)
d.portfolio  <- data.frame(Spot, Spot + val.end.call + call.value + 
                     val.end.put - put.value - 100, 'portfolio',1.1)
colnames(d.underlying) <- c('Spot', 'value', 'Legend','size')
colnames(d.shortcall)  <- c('Spot', 'value', 'Legend','size')
colnames(d.longput)    <- c('Spot', 'value', 'Legend','size')
colnames(d.portfolio)  <- c('Spot', 'value', 'Legend','size')
dd <- rbind(d.underlying,d.shortcall,d.longput,d.portfolio)
p  <- qplot(Spot, value, data=dd, color = Legend, geom = "line",
            size=size )
p <- p + xlab('Value of the underlying' ) + ylab('Profit at maturity')
p <- p + theme(legend.position = LegendPos)
p <- p + scale_size(guide = 'none')
print(p) 


## -----------------------------------------------------------------------------
##--------- the example of the capital protected structure
N           <- 5
nominal     <- 1000
inDeposit   <- nominal * (1.02)^(-N)
cst         <- 0.01  *0        # in PERCENT
pvCosts     <- N * cst * nominal
rest4option <- 1000 - inDeposit - pvCosts
callPrice   <- call_price (100, Strike=100, T=5, r=0.02, vol=0.02)
# reformulate this price as a percentage and then adjust to nominal
callPrice   <- callPrice / 100 * 1000   
gearing     <- rest4option / callPrice
paste('The gearing is:', round(gearing * 100,2))


## -----------------------------------------------------------------------------
# install once: install.packages('ggplot2')
library(ggplot2)


## ----ggSimplest,fig.cap='A basic and simple scatter-plot generated with {\\tt ggplot2}.'----
p <- ggplot(mtcars, aes(x=wt, y=mpg))
# So far printing p would result in an empty plot.
# We need to add a geom to tell ggplot how to plot.
p <- p + geom_point()

# Now, print the plot:
p


## ----ggSimpleTitle, fig.cap='The same plot as in previous figure, but now enhanced with Loess estimates complete with their confidence interval, custom title, axis labels, and improved font properties.'----
# Note that we can reuse the object p
p <- p + geom_smooth(method = "loess", span = 0.9, alpha = 0.3)
p <- p + xlab('Weight') + ggtitle('Fuel consumption explained')
p <- p + ylab('Fuel consumption in mpg')
p <- p + theme(axis.text  = element_text(size=14),
               axis.title = element_text(size=16),
	       plot.title = element_text(size=22, face = 'bold'))
p


## ----eval=FALSE---------------------------------------------------------------
## # Try *one* of the following:
## p <- p + theme_minimal(base_size = 16)  # change also font size
## p <- p + theme_light()                  # light theme
## p <- p + theme_dark()                   # dark theme
## p <- p + theme_classic()
## p <- p + theme_minimal()
## # Much more themes are available right after loading ggplot2.
## 
## # Change the theme for the entire session:
## theme_set(theme_classic(base_size = 12))


## ----ggAddLayers,fig.cap='The same plot as in previous Figure, but now enhanced with different colours that depend on the gearbox (the parameter {\\tt am}, and the size of the dot corresponds to the time the car needs to get from standstill to the distance of 0.2~Miles (ca. 322~m).'----
# We start from the previously generated plot p.
p <- p + aes(colour = factor(am))
p <- p + aes(size   = qsec)
p


## ----ggFacet,message=FALSE,warning=FALSE,fig.cap='A facet plot will create sub-plots per discrete value of one or more variables. In this case we used the number of cylinders, so each sub-plot is for the number of cylinders that is in its title.'----
library(tidyverse)  # provides the pipe: %>%
mtcars %>%
  ggplot(aes(x = wt, y = mpg)) + 
        geom_point() +
        geom_smooth(method = "loess", span = 0.99, alpha=0.3) + 
	xlab('Weight') + 
	ggtitle('MPG in function of weight per nbr cylinders') + 
	facet_wrap( ~ cyl, nrow=2) + 
	theme_classic(base_size = 14) # helps the grey to stand out


## ----bigScatterPlot,fig.cap='The standard functionality for scatterplots is not optimal for large datasets. In this case, it is not clear what the relation is between LTI and DPD. ggplot2 will provide us some more tools to handle this situation.'----
# Set the seed to allow results to be replicated:
set.seed(1868)

# Generate the data:
LTI <- runif(10000, min = 0, max = 1)
DPD <- abs(rnorm(10000, 
                 mean = 70 * ifelse(LTI < 0.5, LTI, 0.5), 
		 sd = 30 * sqrt(LTI) + 5))

# Plot the newly generated data and try to make it not cluttered:
plot(LTI, DPD, 
     pch=19,          # small dot
     ylim =c(0, 100)) # not show outliers, hence zooming in


## ----ggContour,fig.cap='The contour plot is able to show where the density of points is highest by a visually attractive gradient in colour.'----
# We add also the colour schemes of viridisLite.
library(viridisLite)

d <- data.frame(LTI = LTI, DPD = DPD)
p <- ggplot(d, aes(x = LTI, y = DPD)) + 
     stat_density_2d(geom = "raster", aes(fill = ..density..), 
                     contour = FALSE) + 
     geom_density_2d() + 
     scale_fill_gradientn(colours = viridis(256, option = "D")) + 
     ylim(0,100)
p
# Note that ggplot will warn us about the data that is not shown
# due to the cut-off on the y-axis -- via the function ylim(0, 100).


## ----geomSmooth,fig.cap='Adding a Loess estimate is a good idea to visualize the general trend in the data. The algorithm, however needs a time that increases as the square of the number of observations and will be too slow for larger datasets.'----
# Now plot with ggplot. 
# This will generate warnings, because the observations that are not 
# plotted will not be used for our smoothing algorithm.
ggplot(d, aes(x = LTI, y = DPD)) + geom_point(alpha=0.25) +
       geom_smooth(method='loess') + 
       ylim(0,100)


## ----viridisLite1,fig.cap='This plot shows a facet plot of a contour plot with customised colour scheme.'----
library(ggplot2)
library(viridisLite)

# take a subset of the dataset diamonds
set.seed(1867)
d <- diamonds[sample(nrow(diamonds), 1500),]
p <- ggplot(d, aes(x, price)) +
     stat_density_2d(geom = "raster", aes(fill = ..density..), 
                     contour = FALSE) + 
#     geom_density_2d() + 
     facet_grid(. ~ cut) +
     scale_fill_gradientn(colours = viridis(256, option = "D")) +
     ggtitle('Diamonds per cut type')
p


## ----eval=FALSE---------------------------------------------------------------
## ---
## title: "R Markdown"
## author: "Philippe De Brouwer"
## date: "January 1, 2020"
## output: beamer_presentation
## ---
## 
## ```{r setup, include=FALSE}
## knitr::opts_chunk$set(echo = FALSE)
## ```
## 
## ## R Markdown
## 
## This is an R Markdown presentation. Markdown is a simple formatting
## syntax for authoring HTML, PDF, and MS Word documents. For more
## details on using R Markdown see <http://rmarkdown.rstudio.com>.
## 
## When you click the **Knit** button a document will be generated
## that includes both content as well as the output of any embedded R
## code chunks within the document.
## 
## ## Slide with Bullets
## 
## - Bullet 1
## - Bullet 2
## - Bullet 3
## 
## ## Slide with R Output
## 
## ```{r cars, echo = TRUE}
## summary(cars)
## ```
## 
## ## Slide with Plot
## 
## ```{r pressure}
## plot(pressure)
## ```


## ----eval=FALSE---------------------------------------------------------------
## # Level 1: title slide
## 
## ## level 2: title or new slide
## This line goes on slide the slide with title 'level 2: title
## or new slide'
## 
## ### level 3: box on slide
## This line goes in the body of the box that will have the title
## 'level 3: box on slide'


## ----eval=FALSE,include=FALSE,echo=FALSE--------------------------------------
## $

## ----eval=FALSE---------------------------------------------------------------
## ## 50 random numbers
## ```{r showPlot}
## hist(runif(50))
## ```


## -----------------------------------------------------------------------------


## ----shiny1, eval=FALSE-------------------------------------------------------
## library(shiny)
## runExample("01_hello")


## ----eval=FALSE---------------------------------------------------------------
## # The filename must be app.R to execute it on a server.
## 
## # The name of the server function must match the argument in the
## # shinyApp function also it must take the arguments input and output.
## server <- function(input, output) {
##   output$distPlot <- renderPlot({
##     # any plot must be passed through renderPlot()
##     hist(rnorm(input$nbr_obs), col = 'khaki3', border = 'white',
##          breaks=input$breaks,
## 	 main = input$title,
## 	 xlab = "random observations")
##     })
##  output$nbr_txt <- renderText({
##    paste("You have selected", input$nbr_obs, "observations.")
##  })  # note the use of brackets ({ })
## }
## 
## # The name of the ui object must match the argument in the shinyApp
## # function and we must provide an object ui (that holds the html
## # code for our page).
## ui <- fluidPage(
##   titlePanel("Our random simulator"),
##   sidebarLayout(
##     sidebarPanel(
##       sliderInput("nbr_obs", "Number of observations:",
##                   min = 10, max = 500, value = 100),
##       sliderInput("breaks", "Number of bins:",
##                   min = 4, max = 50, value = 10),
##       textInput("title", "Title of the plot",
##                 value = "title goes here")
##     ),
##     mainPanel(plotOutput("distPlot"),
##               h4("Conclusion"),
## 	      p("Small sample sizes combined with a high number of
## 	        bins might provide a visual image of the
## 		distribution that does not resemble the underlying
## 		dynamics."),
## 	       "Note that we can provide text, but not
## 	       <b>html code</b> directly.",
## 	       textOutput("nbr_txt")  # object name in quotes
## 	      )
##   )
## )
## 
## # finally we call the shinyApp function
## shinyApp(ui = ui, server = server)


## ----eval=FALSE---------------------------------------------------------------
## # Load the library:
## library(rsconnect)
## 
## # Upload the app (default filename is app.R)
## rsconnect::deployApp('path/to/your/app',
##                      server="shinyapps.io", account="xxx")


## ----eval=FALSE---------------------------------------------------------------
## library(leaflet)
## content <- paste(sep = "<br/>",
##   "<b><a href='http://www.de-brouwer.com/honcon/'>
##       Honorary Consulate of Belgium</a></b>",
##   "ul. Marii Grzegorzewskiej 33A",
##   "30-394 Krakow"
##   )
## map <- leaflet()                                    %>%
##   addProviderTiles(providers$OpenStreetMap)         %>%
##   addMarkers(lng = 19.870188, lat = 50.009159)      %>%
##   addPopups(lng = 19.870188, lat = 50.009159, content,
##     options = popupOptions(closeButton = TRUE))     %>%
##   setView(lat = 50.009159,lng = 19.870188, zoom = 12)
## map


## ----eval=FALSE---------------------------------------------------------------
## library(titanic)    # for the data
## library(tidyverse)  # for the tibble
## library(ggvis)      # for the plot
## 
## titanic_train$Age     %>%
##     as_tibble         %>%
##     na.omit           %>%
##     ggvis(x = ~value) %>%
##     layer_densities(
##       adjust = input_slider(.1, 2, value = 1, step = .1,
##                             label = "Bandwidth"),
##       kernel = input_select(
##         c("Gaussian"     = "gaussian",
##           "Epanechnikov" = "epanechnikov",
##           "Rectangular"  = "rectangular",
##           "Triangular"   = "triangular",
##           "Biweight"     = "biweight",
##           "Cosine"       = "cosine",
##           "Optcosine"    = "optcosine"),
##         label = "Kernel")
##     )


## ----eval=FALSE---------------------------------------------------------------
## library(ggvis)
## 
## function(input, output, session) {
##   # A reactive subset of mtcars
##   reacCars <- reactive({ mtcars[1:input$n, ] })
## 
##   # register observers and place the controls
##   reacCars %>%
##      ggvis(~wt, ~mpg, fill=(~cyl)) %>%
##      layer_points() %>%
##      layer_smooths(span = input_slider(0.5, 1, value = 1,
##         label = 'smoothing span:')) %>%
##      bind_shiny("plot1", "plot_ui_div")
## 
##   output$carsData <- renderTable({ reacCars()[, c("wt", "mpg")] })
## }


## ----eval=FALSE---------------------------------------------------------------
## library(ggvis)
## 
## fluidPage(sidebarLayout(
##   sidebarPanel(
##     # explicit code for a slider-bar:
##     sliderInput("n", "Number of points", min = 1, max = nrow(mtcars),
##                 value = 10, step = 1),
##     # no code needed for the smoothing span, ggvis does this
##     uiOutput("plot_ui_div") # produces a <div> with id corresponding
##                             # to argument in bind_shiny
##   ),
##   mainPanel(
##     # place the plot "plot1" here:
##     ggvisOutput("plot1"),    # matches argument to bind_shiny()
##     # under this the table of selected card models:
##     tableOutput("carsData")  # parses the result of renderTable()
##   )
## ))


## ----eval=FALSE---------------------------------------------------------------
## shiny::runApp("/path/to/my/app")


## ----googleVis1, eval=FALSE---------------------------------------------------
## library(googleVis)
## demo(package='googleVis')


## -----------------------------------------------------------------------------
# diversity
# Calculates the entropy of a system with equiprobable states.
# Arguments:
#    x -- numeric vector -- observed probabilities of classes
# Returns:
#    numeric -- the entropy / diversity measure
diversity <- function(x) {
  f <- function(x) x * log(x)
  x1 <- mapply(FUN = f, x)
  - sum(x1) / log(length(x))
  }


## -----------------------------------------------------------------------------
# diversity
# Calculates the entropy of a system with discrete states.
# Arguments:
#   x     -- numeric vector -- observed probabilities of classes
#   prior -- numeric vector -- prior probabilities of the classes
# Returns:
#    numeric -- the entropy / diversity measure
diversity <- function(x, prior = NULL) {
  if (min(x) <= 0) {return(0);} # the log will fail for 0
  # if the numbers are higher than 1, then not probabilities but 
  # populations are given, so we rescale to probabilities:
  if (sum(x) != 1) {x <- x / sum(x)}
  N <- length(x)
  if(!is.null(prior)) {
    for (i in (1:N)) {
      a <- (1 - 1 / (N * prior[i])) / (1 - prior[i])
      b <- (1 - N * prior[i]^2) / (N * prior[i] * (1 - prior[i]))
      x[i] <- a * x[i]^2 + b * x[i]
    }
   }
  f <- function(x) x * log(x)
  x1 <- mapply(FUN = f, x)
  - sum(x1) / log(N)
  }


## -----------------------------------------------------------------------------
# Consider the following prior probabilities:
pri <- c(0.1,0.5,0.4)

# No prior priorities supplied, so 1/N is most diverse:
diversity(c(0.1,0.5,0.4))      

# The population matches prior probabilities, so index should be 1
diversity(c(0.1,0.5,0.4), prior = pri) 

# Very non-diverse population:
diversity(c(0.999,0.0005,0.0005), prior = pri) 

# Only one sub-group is represented (no diversity):
diversity(c(1,0,0), prior = pri)         

# Numbers instead of probabilities provided, also this works:
diversity(c(100,150,200)) 


## ----divFunction,fig.cap="The evolution of the gender-diversity-index in function of one of the representation of one of the genders in our population."----
females <- seq(from=0,to=1, length.out = 100)
div <- numeric(0)
for (i in (1:length(females))) {
  div[i] <- diversity (c(females[i], 1 - females[i]))
  }
  
d <- as.data.frame(cbind(females, div)  )
colnames(d) <- c('percentage females', 'diversity index')
library(ggplot2)
p <- ggplot(data = d, 
        aes(x = `percentage females`, y = `diversity index`)) +
     geom_line(color = 'red', lwd = 3) + 
     ggtitle('Diversity Index') + 
     xlab('percentage females') + ylab('diversity index')
p


## -----------------------------------------------------------------------------
library(tidyverse)
N <- 200
set.seed(1866)

d0 <- data.frame("ID"      = 1:N,

          # Log-normal age distribution
          "age"         = round(rlnorm(N, log(30), log(1.25))),
	  # A significant bias towards the old continent:
          "continent"   = ifelse(runif(N) < 0.3, "America", 
	                   ifelse(runif(N) < 0.7,"Europe","Other")),
          # A mild bias towards males:
          "gender"      = ifelse(runif(N) < 0.45, "F", "M"),
	  # Grade will be filled in later:
          "grade"       = 0,
	  # Three teams of different sizes:
          "team"        = ifelse(runif(N) < 0.6, "bigTeam", 
	                    ifelse(runif(N) < 0.6, 
                            "mediumTeam", 
		            ifelse(runif(N) < 0.8, "smallTeam", 
			           "XsmallTeam"))),
          # Most people have little people depending on them:
          "dependents"  = round(rlnorm(N,log(0.75),log(2.5))),
	  # Random performance (no bias linked, but different group sizes):
          "performance" = ifelse(runif(N) < 0.1, "L", 
	                    ifelse(runif(N) < 0.6, "M", 
                            ifelse(runif(N) < 0.7, "H", "XH"))),
          # Salary will be filled in later:
          "salary"      = 0,
	  # We make just a snapshot dashboard, so we do not need this now, 
	  # but we could use this later to show evolution:
          "timestamp"   = as.Date("2020-01-01")
          )

# Now we clean up age and fill in grade, salary and lastPromoted without 
# any bias for gender, origin -- but with a bias for age.
d1 <- d0                                                  %>%
  mutate(age    = ifelse((age < 18), age + 10, age))      %>%
  mutate(grade  = ifelse(runif(N) * age < 20, 0, 
                    ifelse(runif(N) * age < 25, 1, 
                    ifelse(runif(N) * age < 30, 2, 3))))  %>%
  mutate(salary = round(exp(0.75*grade)*4000 + 
                    rnorm(N,0,1500)))                     %>%
  mutate(lastPromoted = round(exp(0.05*(3-grade))*1 + 
                    abs(rnorm(N,0,5))) -1)


## -----------------------------------------------------------------------------
# If not done yet, install the package:
# install.packages('flexdashboard')

# Then load the package:
library(flexdashboard)


## ----flexdash,eval=FALSE,echo=FALSE-------------------------------------------
## ---
## title: "Divsersity in Action"
## output:
##   flexdashboard::flex_dashboard:
##     theme: cosmo
##     orientation: rows
##     vertical_layout: fill
##     #storyboard: true
##     social: menu
##     source: embed
## ---
## 
## ```{r}
## # (C) Philippe J.S. De Brouwer -- 2019
## # demo: http://rpubs.com/phdb/diversity_dash01
## ```
## 
## ```{r setup, include=FALSE}
## library(flexdashboard)
## library(tidyverse)
## library(ggplot2)
## library(knitr)
## library(gridExtra)
## #install.packages('plotly')
## library(plotly)
## N <- 150
## set.seed(1865)
## 
## d0 <- data.frame("ID"         = 1:N,
##                 "age"         = round(rlnorm(N, log(30), log(1.25))),
##                 "continent"   = ifelse(runif(N) < 0.3, "America", ifelse(runif(N) < 0.7, "Europe","Other")),
##                 "gender"      = ifelse(runif(N) < 0.4, "F", "M"),
##                 "grade"       = 0,
##                 "team"        = ifelse(runif(N) < 0.6, "bigTeam", ifelse(runif(N) < 0.6,
##                                "mediumTeam", ifelse(runif(N) < 0.8, "smallTeam", "XsmallTeam"))),
##                 "dependents"  = round(rlnorm(N,log(0.65),log(1.5))),
##                 "performance" = ifelse(runif(N) < 0.1, "L", ifelse(runif(N) < 0.6, "M",
##                                 ifelse(runif(N) < 0.7, "H", "XH"))),
##                 "salary"      = 0,
##                 "timestamp"   = as.Date("2020-01-01")
##                 )
## 
## d1 <- d0 %>%
##   mutate(age    = ifelse((age < 18), age + 10, age)) %>%
##   mutate(grade  = ifelse(runif(N) * age < 20, 0, ifelse(runif(N) * age < 25, 1, ifelse(runif(N) * age < 30, 2, 3)))) %>%
##   mutate(salary = round(exp(0.75*grade)*4000 + rnorm(N,0,2500)))  %>%
##   mutate(lastPromoted = round(exp(0.05*(3-grade))*1 + abs(rnorm(N,0,5))) -1)
## ```
## 
## Overview
## ========
## 
## Row
## -------------------------------------
## 
## ```{r}
## # our diversity function
## diversity <- function(x, prior = NULL) {
##   if (min(x) <= 0) {return(0);} # the log will fail for 0
##   # if the numbers are higher than 1, then not probabilities but
##   # populations are given, so we rescale to probabilities:
##   if (sum(x) != 1) {x <- x / sum(x)}
##   N <- length(x)
##   if(!is.null(prior)) {
##     for (i in (1:N)) {
##       a <- (1 - 1 / (N * prior[i])) / (1 - prior[i])
##       b <- (1 - N * prior[i]^2) / (N * prior[i] * (1 - prior[i]))
##       x[i] <- a * x[i]^2 + b * x[i]
##     }
##    }
##   f <- function(x) x * log(x)
##   x1 <- mapply(FUN = f, x)
##   - sum(x1) / log(N)
## }
## # the gauges for the different dimensions
## ```
## 
## ### Gender
## ```{r}
## # ranges:
## rGreen <- c(0.900001, 1)
## rAmber <- c(0.800001, 0.9)
## rRed   <- c(0, 0.8)
## iGender <- round(diversity(table(d1$gender)),3)
## gauge(iGender, min = 0, max = 1, gaugeSectors(
##   success = rGreen, warning = rAmber, danger = rRed
##   ))
## kable(table(d1$gender))
## ```
## 
## ### Age
## ```{r}
## # consider each band of ten years as a group
## iAge <- round(diversity(table(round(d1$age/10))),3)
## gauge(iAge, min = 0, max = 1, gaugeSectors(
##   success = rGreen, warning = rAmber, danger = rRed
##   ))
## kable(table(round(d1$age/10)*10))
## ```
## 
## 
## ### Roots
## ```{r}
## iRoots <- round(diversity(table(d1$continent)),3)
## gauge(iRoots, min = 0, max = 1, gaugeSectors(
##   success = rGreen, warning = rAmber, danger = rRed
##   ))
## kable(table(d1$continent))
## ```
## 
## ### Dependents
## ```{r}
## # we only monitor if someone has dependents
## xdep <- ifelse(d1$dependents >= 1, 1, 0)
## iDep <- round(diversity(table(xdep)),3)
## gauge(iDep, min = 0, max = 1, gaugeSectors(
##   success = rGreen, warning = rAmber, danger = rRed
##   ))When the aforementioned code is executed,
## kable(table(d1$dependents))
## ```
## 
## Row
## -------------------------------------
## ### Info
## The diversity indeces show how diverse our workforce is. They are calculated similar to entropy: $I = -\frac{1}{\log(N)} \sum_i^N {p_i \log p_i}$, where there are $N$ possible and mutually exclusive states $i$. They range from $0$ to $1$.
## 
## ### Average Diversity Index
## ```{r}
## x  <- mean(c(iGender, iAge, iDep, iRoots))
## valueBox(x,
##          icon  = ifelse(x > 0.9, "fa-smile" , ifelse(x > 0.8, "fa-meh", "fa-sad-tear")),
##          color = ifelse(x > 0.9, "success" , ifelse(x > 0.8, "warning", "danger"))
##          )
## ```
## 
## 
## Gender
## ======================================
## 
## Row {.tabset}
## -------------------------------------
## 
## ### Composition
## ```{r}
## p2 <- ggplot(data = d1, aes(x=gender, fill=gender)) +
##   geom_bar(stat="count", width=0.7) +
##   facet_grid(rows=d1$grade) +
##   ggtitle('workforce composition i.f.o. salary grade (level in the company)')
## ggplotly(p2)
## ```
## 
## ### Salary
## ```{r}
## p1 <- ggplot(data = d1, aes(x=gender, y=salary, fill=gender)) +
##   geom_boxplot() +
##   facet_grid(rows=d1$grade) +
##   ggtitle('The salary gap per salary grade (level in the company)')
## ggplotly(p1)
## ```
## 
## ### Promotions
## ```{r}
## d1$promoted = ifelse(d1$lastPromoted <= 2,1,0)
## p <- ggplot(data = d1, aes(x=gender, fill=gender, y=promoted)) +
##   stat_summary(fun.y="mean", geom="bar") +
##   facet_grid(rows=d1$grade) +
##   ggtitle('promotion propensity per grae')
## p
## ggplotly(p)
## ```
## 
## Age
## ===
## Row {.tabset}
## -------------------------------------
## 
## ### Histogram
## 
## ```{r}
## qplot(d1$age, geom="histogram", binwidth=5,
##         fill=I("steelblue"), col=I("black")
##         ) +
##   xlab('Age') +
##   ggtitle('Histogram for Age')
## ```
## 
## Roots {.tabset}
## ===============
## Row {.tabset}
## -------------------------------------
## 
## ### rworldmap
## 
## ```{r warning=FALSE}When the aforementioned code is executed,
## # the default R-approach:
## install.packages('rworldmap')
## library(rworldmap)
## 
## nbrPerCountry = read.table(text="
## country value
## Poland 100
## Ukraine 65
## UK 2
## USA 1
## China 3
## Germany 0
## France 1
## Italy 20
## Greece 25
## Spain 13
## Portugal 7
## Mexico 55
## Belarus 5
## Russia 7
## Vietnam 1
## India 25
## Belgium 2
## Chile 6
## Congo 1
## ", header=T)
## 
## x <- joinCountryData2Map(nbrPerCountry, joinCode="NAME", nameJoinColumn="country")
## 
## mapCountryData(x, nameColumnToPlot="value", catMethod="fixedWidth")
## ```
## 
## ### leaflet
## ```{r}
## #install.packages('maps')
## library(maps)
## #install.packages('sp')
## library(sp)
## #install.packages('leaflet')
## library(leaflet)
## map <- leaflet() %>%
##   addProviderTiles(providers$OpenStreetMap) %>%
##   setView(lng = 0, lat = 0, zoom = 2)
## map
## ```
## 
## Dependents
## ==========
## Row {.tabset}
## -------------------------------------
## 
## ### Histogram
## 
## ```{r}
## qplot(d1$dependents, geom="histogram", binwidth=1,
##         fill=I("steelblue"), col=I("black")
##         ) +
##   xlab('Dependents') +
##   ggtitle('Histogram for number of dependents')
## ```


## ----flexdashStructure,eval=FALSE,echo=TRUE-----------------------------------
## ---
## title: "Diversity in Action"
## output:
##   flexdashboard::flex_dashboard:
##     theme: cosmo
##     orientation: rows
##     vertical_layout: fill
##     social: menu
##     source: embed
## ---
## 
## ```{r setup, include=FALSE}
## # load packages and functions
## ```
## Overview
## ========
## Row
## -------------------------------------
## ### Gender
## ```{r}
## code to generate the first gauge
## ```
## 
## Row
## -------------------------------------
## 
## Gender
## ======================================
## 
## Row {.tabset}
## -------------------------------------
## 
## Age
## ======================================
## Row {.tabset}
## -------------------------------------
## ### first tab
## ### second tab
## 
## Roots
## ======================================
## Row {.tabset}
## -------------------------------------
## 
## Dependants
## ======================================
## Row {.tabset}
## -------------------------------------


## ----flexdashOverview,eval=FALSE,echo=TRUE------------------------------------
## Overview
## ========
## 
## Row
## -------------------------------------
## 
## ```{r}
## # here goes the R-code to get data and calculate diversity indices.
## ```
## 
## ### Gender
## ```{r genderGauge}
## # ranges:
## rGreen   <- c(0.900001, 1)
## rAmber   <- c(0.800001, 0.9)
## rRed     <- c(0, 0.8)
## iGender  <- round(diversity(table(d1$gender)),3)
## gauge(iGender, min = 0, max = 1, gaugeSectors(
##       success = rGreen, warning = rAmber, danger = rRed
##       ))
## kable(table(d1$gender))
## ```
## ### This is only the first gauge, the next one can be described below.


## ----eval=FALSE---------------------------------------------------------------
## gauge <-  gvisGauge(iGender,
##              options=list(min = 0, max = 1,
## 	             greenFrom = 0.9,greenTo = 1,
## 		     yellowFrom = 0.8, yellowTo = 0.8,
## 		     redFrom = 0, redTo = 0.8,
## 		     width = 400, height = 300))


## ----flexdashGender,eval=FALSE,echo=TRUE--------------------------------------
## Gender
## ======================================
## 
## Row {.tabset}
## -------------------------------------
## 
## ### Composition
## ```{r}
## p2 <- ggplot(data = d1, aes(x=gender, fill=gender)) +
##   geom_bar(stat="count", width=0.7) +
##   facet_grid(rows=d1$grade) +
##   ggtitle('workforce composition i.f.o. salary grade')
## ggplotly(p2)
## ```
## 
## ### Tab2
## ```{r RCodeForTab2}
## # etc.
## ```


## ----shinyDash1,eval=FALSE,echo=TRUE------------------------------------------
## # this file should be called app.R
## library(shinydashboard)
## 
## # general code (not part of server or ui function)
## my_seed <- 1865
## set.seed(my_seed)
## N <- 150
## 
## # user interface ui
## ui <- dashboardPage(
##   dashboardHeader(title = "ShinyDashBoard"),
##   dashboardSidebar(
##       title = "Choose ...",
##       numericInput('N', 'Number of data points:', N),
##       numericInput('my_seed', 'seed:', my_seed),
##       sliderInput("bins", "Number of bins:",
##                   min = 1, max = 50, value = 30)
##   ),
##   dashboardBody(
##     fluidRow(
##       valueBoxOutput("box1"),
##       valueBoxOutput("box2")
##       ),
##     plotOutput("p1", height = 250)
##     )
##   )
## 
## # server function
## server <- function(input, output) {
##   d <- reactive({
##     set.seed(input$my_seed)
##     rnorm(input$N)
##     })
##   output$p1 <- renderPlot({
##     x    <- d()
##     bins <- seq(min(x), max(x), length.out = input$bins + 1)
##     hist(x, breaks = bins, col = 'deepskyblue3', border = 'gray')
##     })
##   output$box1 <- renderValueBox({
##     valueBox(
##       value    = formatC(mean(d()), digits = 2, format = "f"),
##       subtitle = "mean", icon     = icon("globe"),
##       color    = "light-blue")
##     })
##   output$box2 <- renderValueBox({
##       valueBox(
##         value    = formatC(sd(d()), digits = 2, format = "f"),
##         subtitle = "standard deviation",  icon     = icon("table"),
##         color    = "light-blue")
##       })
## 
## }
## 
## # load the app
## shinyApp(ui, server)


## -----------------------------------------------------------------------------
# Load the library:
library(parallel)

# The function detectCores finds the number of cores available
numCores <- detectCores() 
numCores 


## ----mclapply1----------------------------------------------------------------
# Set the sample size:
N <- 50

# Set the starting points for the k-means algorithm:
starts <- rep(100, N)   # 50 times the value 100

# Prepare data as numeric only:
t <- as.data.frame(titanic_train)
t <- t[complete.cases(t),]
t <- t[,c(2,3,5,6,7,8)]
t$Sex <- ifelse(t$Sex == 'male', 1, 0)

# Prepare the functions to be executed:
f <- function(nstart) kmeans(t, 4, nstart = nstart)

# Now, we are ready to go, we try two approaches:

# 1. With the standard function base::lapply
system.time(results <- lapply(starts, f))

# 2. With parallel::mclapply
system.time(results <- mclapply(starts, f, mc.cores = numCores))


## -----------------------------------------------------------------------------
# 1. The regular for loop:
for (n in 1:4) print(gamma(n))

# 2. The foreach loop:
library(foreach)
foreach (n = 1:4) %do% print(gamma(n))


## -----------------------------------------------------------------------------
class(foreach (n = 1:4))


## -----------------------------------------------------------------------------
# Once installed, the package must be loaded once per session:
library(doParallel)

# Find the number of available cores (as in previous section):
numCores <- parallel::detectCores() 

# Register doParallel:
registerDoParallel(numCores)  

# Now, we are ready to put it in action:
foreach (n = 1:4) %dopar% print(gamma(n))


## -----------------------------------------------------------------------------
# Collapse to wide data.frame:
foreach (n = 1:4, .combine = cbind) %dopar% print(gamma(n))

# Collapse to long data.frame:
foreach (n = 1:4, .combine = rbind) %dopar% print(gamma(n))

# Collapse to vector:
foreach (n = 1:4, .combine = c)     %dopar% print(gamma(n))


## -----------------------------------------------------------------------------
# We reuse parameters and data as defined in the previous section:
system.time(
foreach (n = 1:N) %dopar% {
    kmeans(t, 4, nstart=100)
    }
  )


## ----echo=FALSE---------------------------------------------------------------
detach("package:doParallel", unload=TRUE)
#detach("package:parallel", unload=TRUE)


## -----------------------------------------------------------------------------
library(snow)


## ----eval=FALSE---------------------------------------------------------------
## detach("package:parallel", unload = TRUE)
## library(snow)


## -----------------------------------------------------------------------------
cl <- makeCluster(c("localhost", "localhost"), type = "SOCK")


## ----snowTitanic--------------------------------------------------------------
library(titanic)  # provides the dataset: titanic_train
N <- 50
starts <- rep(100, N)   # 50 times the value 100

# Prepare data as numeric only:
t <- as.data.frame(titanic_train)
t <- t[complete.cases(t),]
t <- t[,c(2,3,5,6,7,8)]
t$Sex <- ifelse(t$Sex == 'male', 1, 0)

# Prepare the functions to be executed:
f <- function(nstart) kmeans(t, 4, nstart=nstart)

# 1. with the standard function 
system.time(results <- lapply(starts, f))

# 2. with the cluster
# First, we must export the object t,  so that it can 
# be used by the cluster:
clusterExport(cl, "t")
#clusterExport(cl, "starts") # Not needed since it is in the function f
system.time(
  result2 <- parLapply(cl, starts, f)
  )


## -----------------------------------------------------------------------------
f <- function (x, y) x + y + rnorm(1)
clusterCall(cl, function (x, y) {x + y + rnorm(1)}, 0, pi)
clusterCall(cl, f, 0, pi)

# Both forms are semantically similar to:
clusterCall(cl, eval, f(0,pi))
# Note that the random numbers are the same on both clusters.


## -----------------------------------------------------------------------------
# Note that ...
clusterEvalQ(cl, rnorm(1))

# ... is almost the same as
clusterCall(cl, evalq, rnorm(1))


## -----------------------------------------------------------------------------
clusterApply(cl, c(0, 10), sum, pi)


## -----------------------------------------------------------------------------
stopCluster(cl)


## ----echo=FALSE---------------------------------------------------------------
detach("package:snow", unload=TRUE)


## ----eval=FALSE---------------------------------------------------------------
## # Remove a list of selected variables:
## rm (list = c('A,, 'gpuA', 'gpuB', 'vlcA', 'vlcB'))
## 
## # Remove *all* variables previously defined (be careful!):
## rm(list = ls())
## 
## # Run the garbage collector:
## gc()
## 
## # If not longer needed, unload gpuR:
## detach('package:gpuR', unload = TRUE)


## ----eval=FALSE---------------------------------------------------------------
## # install from github
## devtools::install_github(repo = "rstudio/spark-install",
##                          subdir = "R")
## library(sparkinstall)
## 
## # lists the versions available to install
## spark_available_versions()
## 
## # installs an specific version
## spark_install(version = "2.4.3")
## 
## # uninstalls an specific version
## spark_uninstall(version = "2.4.3", hadoop_version = "2.6")


## -----------------------------------------------------------------------------
system('start-master.sh')


## ----sparkrLoad,message=FALSE,warning=FALSE-----------------------------------
library(tidyverse)
library(dplyr)
library(SparkR)
# Note that loading SparkR will generate many warning messages,
# because it overrrides many functions such as summary, first,
# last, corr, ceil, rbind, expr, cov, sd and many more.

sc <- sparkR.session(master = "local", appName = 'first test',
                     sparkConfig = list(spark.driver.memory = '2g'))

# Show the session:
sc

DF <- as.DataFrame(mtcars)

# The DataFrame is for big data, so the attempt to print all data,
# might surprise us a little:
DF
# R assumes that the data-frame is big data and does not even
# start printing all data.

# head() will collapse the first lines to a data.frame:
head(DF)


## ----eval=FALSE---------------------------------------------------------------
## # If DDF is a distributed DataFrame defined by SparkR,
## # we can add a checkpoint as follows:
## DDF <- SparkR::checkpoint(DDF)


## -----------------------------------------------------------------------------
library(titanic)
library(tidyverse)

# This provides a.o. two datasets titanic_train and titanic_test.
# We will work further with the training-dataset.
T <- as.DataFrame(titanic_train)

# The SparkDataFrame inherits from data.frame, so most function 
# work as expected on a DataFrame:
colnames(T)
str(T)
summary(T)
class(T)

# But has more functions.
printSchema(T)

# Truncated information collapses to data.frame:
T %>% head %>% class


## -----------------------------------------------------------------------------
X <- T %>% SparkR::select(T$Age)         %>% head
Y <- T %>% SparkR::select(column('Age')) %>% head
Z <- T %>% SparkR::select(expr('Age'))   %>% head
cbind(X, Y, Z)


## -----------------------------------------------------------------------------
T %>% 
  SparkR::filter("Age < 20 AND Sex == 'male' AND Survived == 1") %>%
  SparkR::select(expr('PassengerId'), expr('Pclass'), expr('Age'), 
                 expr('Survived'), expr('Embarked'))             %>%
  head
      
# The following is another approach. The end-result is the same, 
# however, we bring the data to the R's working memory and then 
# use dplyr. Note the subtle differences in syntax.
SparkR::collect(T)                                              %>%
      dplyr::filter(Age < 20 & Sex == 'male' & Survived == 1)   %>%
      dplyr::select(PassengerId, Pclass, Age, Survived,
                   Embarked)                                    %>%
      head


## -----------------------------------------------------------------------------
# Most functions work as expected
TMP <- T                                              %>% 
       SparkR::group_by(expr('Pclass'), expr('Sex'))  %>%
             summarize(count = n(expr('Survived')))
N   <- nrow(T)
TMP %>% mutate(Percentage = expr('count') / N  * 100) %>% 
        arrange('Pclass', 'Sex') %>% SparkR::collect()


## ----eval=FALSE---------------------------------------------------------------
## library(tidyverse)
## library(SparkR)
## library(titanic)
## sc <- sparkR.session(master = "local", appName = 'first test',
## sparkConfig = list(spark.driver.memory = '2g'))
## T <- as.DataFrame(titanic_train)


## ----eval=FALSE---------------------------------------------------------------
## # The data:
## T <- as.DataFrame(titanic_train)
## 
## # The schema can be a structType:
## schema <- SparkR::structType(SparkR::structField("Age", "double"),
##                         SparkR::structField("ageGroup", "string"))
## 
## # Or (since Spark 2.3) it can also be a DDL-formatted string
## schema <- "age DOUBLE, ageGroup STRING"
## 
## # The function to be applied:
## f <- function(x) {
##   data.frame(x$Age, if_else(x$Age < 30, "youth", "mature"))
##   }
## 
## # Run the function f on the Spark cluster:
## T2 <- SparkR::dapply(T, f, schema)
## 
## # Inspect the results:
## head(SparkR::collect(T2))
## 
## ##     age           ageGroup
## ## 1    22              youth
## ## 2    38             mature
## ## 3    26              youth
## ## 4    35             mature
## ## 5    35             mature
## ## 6    NA               <NA>


## -----------------------------------------------------------------------------
DFcars  <- createDataFrame(mtcars)
DFcars2 <- dapply(DFcars, function(x) {x}, schema(DFcars))
head(collect(DFcars2))


## ----echo=FALSE,results='hide'------------------------------------------------
# clean up
rm(T)
gc()

## -----------------------------------------------------------------------------
# The data:
T <- as.DataFrame(titanic_train)

# The function to be applied:
f <- function(x) { 
  y <- data.frame(x$Age, ifelse(x$Age < 30, "youth", "mature")) 
  
  # We specify now column names in the data.frame to be returned:
    colnames(y) <- c("age", "ageGroup")
  
  # and we return the data.frame (base R type):
  y
  }

# Run the function f on the Spark cluster:
T2_DF <- dapplyCollect(T, f)

# Inspect the results (T2_DF is now a data.frame, no collect needed):
head(T2_DF)


## -----------------------------------------------------------------------------
# define the function to be used:
f <- function (key, x) {
  data.frame(key, min(x$Age,  na.rm = TRUE), 
                  mean(x$Age, na.rm = TRUE), 
		  max(x$Age,  na.rm = TRUE))
  }

# The schema also can be specified via a DDL-formatted string
schema <- "class INT, min_age DOUBLE, avg_age DOUBLE, max_age DOUBLE"

maxAge <- gapply(T, "Pclass", f, schema)

head(collect(arrange(maxAge, "class", decreasing = TRUE)))


## -----------------------------------------------------------------------------
# define the function to be used:
f <- function (key, x) {
  y <- data.frame(key, min(x$Age,  na.rm = TRUE), 
                       mean(x$Age, na.rm = TRUE), 
	               max(x$Age,  na.rm = TRUE))
  colnames(y) <- c("class", "min_age", "avg_age", "max_age")
  y
  }

maxAge <- gapplyCollect(T, "Pclass", f)

head(maxAge[order(maxAge$class, decreasing = TRUE), ])


## -----------------------------------------------------------------------------
# first a trivial example to show how spark.lapply works:
surfaces <- spark.lapply(1:3, function(x){pi * x^2})
print(surfaces)


## ----eval=FALSE---------------------------------------------------------------
## mFamilies   <- c("binomial", "gaussian")
## trainModels <- function(fam) {
##   m <- SparkR::glm(Survived ~ Age + Pclass,
##                    data = T,
## 		   family = fam)
##   summary(m)
##   }
## mSummaries <- spark.lapply(mFamilies, trainModels)


## ----eval=FALSE---------------------------------------------------------------
## # df is a data.frame (R-data.frame)
## # DF is a DataFrame (distributed Spark data-frame)
## 
## # We already saw how to get data from Spark to R:
## df <- collect(DF)
## 
## # From R to Spark:
## DF <- createDataFrame(df)


## ----eval=FALSE---------------------------------------------------------------
## loadDF(fileName,
##        source = "csv",
##        header = "true",
##        sep = ",")


## -----------------------------------------------------------------------------
T %>% withColumn("AgeGroup", column("Age") / lit(10))        %>%
      SparkR::select(expr('PassengerId'), expr('Age'), 
                      expr('AgeGroup'))                      %>%
      head


## -----------------------------------------------------------------------------
T %>% cube("Pclass", "Sex")  %>%
      agg(avg(T$Age))        %>%
      collect()


## ----warning=FALSE,error=FALSE,message=FALSE----------------------------------
# Prepare training and testing data:
T_split <- randomSplit(T, c(8,2), 2)
T_train <- T_split[[1]]
T_test  <- T_split[[2]]

# Fit the model:
M1 <- spark.glm(T_train, 
                Survived ~ Pclass + Sex, 
		family = "binomial")

# Save the model:
path1 <- tempfile(pattern = 'ml', fileext = '.tmp')
write.ml(M1, path1)

# Retrieve the model
M2 <- read.ml(path1)

# Do something with M2
summary(M2)

# Add predictions to the model for the test data:
PRED1 <- predict(M2, T_test)

# Show the results:
x <- head(collect(PRED1))
head(cbind(x$Survived, x$prediction))

# Close the connection:
unlink(path1)


## ----eval=FALSE---------------------------------------------------------------
## # We unload the package SparkR.
## detach("package:SparkR", unload=TRUE)
## detach("package:tidyverse", unload=TRUE)
## detach("package:dplyr", unload=TRUE)


## ----eval=FALSE---------------------------------------------------------------
## install.packages('sparklyr')


## ----eval=FALSE---------------------------------------------------------------
## # First, load the tidyverse packages:
## library(tidyverse)
## library(dplyr)
## 
## # Load the package sparklyr:
## library(sparklyr)
## 
## # Load the data:
## library(titanic)
## 
## # Our spark-master is already running, so no need for:
## # system('start-master.sh')
## 
## # Connect to the local Spark master
## sc <- spark_connect(master = "local")


## ----eval=FALSE---------------------------------------------------------------
## Titanic_tbl <- copy_to(sc, titanic_train)
## Titanic_tbl
## 
## ## # Source: spark<titanic_train> [?? x 12]
## ##    PassengerId Survived Pclass Name    Sex     Age SibSp Parch ...
## ##          <int>    <int>  <int> <chr>   <chr> <dbl> <int> <int> ...
## ##  1           1        0      3 Brau... male     22     1     0 ...
## ##  2           2        1      1 Cumi... female   38     1     0 ...
## ##  3           3        1      3 Heik... female   26     0     0 ...
## ##  4           4        1      1 Futr... female   35     1     0 ...
## ##  5           5        0      3 Alle... male     35     0     0 ...
## ##  6           6        0      3 Mora... male    NaN     0     0 ...
## ##  7           7        0      1 McCa... male     54     0     0 ...
## ##  8           8        0      3 Pals... male      2     3     1 ...
## ##  9           9        1      3 John... female   27     0     2 ...
## ## 10          10        1      2 Nass... female   14     1     0 ...
## ## # ... with more rows, and 1 more variable: Embarked <chr>
## 
## 
## # More datasets can be stored in the same connection:
## cars_tbl <- copy_to(sc, mtcars)
## 
## # List the available tables:
## src_tbls(sc)
## ## [1] "mtcars"        "titanic_train"


## ----eval=FALSE---------------------------------------------------------------
## Titanic_tbl %>% summarise(n = n())
## ## # Source: spark<?> [?? x 1]
## ##       n
## ##   <dbl>
## ## 1   891
## 
## # Alternatively:
## Titanic_tbl %>% spark_dataframe() %>% invoke("count")
## ## [1] 891
## 
## Titanic_tbl                                   %>%
##   dplyr::group_by(Sex, Embarked)              %>%
##   summarise(count = n(), AgeMean = mean(Age)) %>%
##   collect
## ## # A tibble: 7 x 4
## ##   Sex    Embarked count AgeMean
## ##   <chr>  <chr>    <dbl>   <dbl>
## ## 1 male   C           95    33.0
## ## 2 male   S          441    30.3
## ## 3 female C           73    28.3
## ## 4 female ""           2    50
## ## 5 female S          203    27.8
## ## 6 male   Q           41    30.9
## ## 7 female Q           36    24.3


## ----eval=FALSE---------------------------------------------------------------
## library(DBI)
## sSQL <- "SELECT Name, Age, Sex, Embarked FROM titanic_train
##          WHERE Embarked = 'Q' LIMIT 10"
## dbGetQuery(sc, sSQL)
## ##                            Name Age    Sex Embarked
## ## 1              Moran, Mr. James NaN   male        Q
## ## 2          Rice, Master. Eugene   2   male        Q
## ## 3   McGowan, Miss. Anna "Annie"  15 female        Q
## ## 4 O'Dwyer, Miss. Ellen "Nellie" NaN female        Q
## ## 5      Glynn, Miss. Mary Agatha NaN female        Q


## ----eval=FALSE---------------------------------------------------------------
## # sdf_len creates a DataFrame of a given length (5 in the example)
## x <- sdf_len(sc, 5, repartition = 1) %>%
##      spark_apply(function(x) pi * x^2)
## print(x)
## ## # Source: spark<?> [?? x 1]
## ##      id
## ##   <dbl>
## ## 1  3.14
## ## 2 12.6
## ## 3 28.3
## ## 4 50.3
## ## 5 78.5


## ----eval=FALSE---------------------------------------------------------------
## # Transform data and create buckets for Age:
## t <- Titanic_tbl %>%
##   ft_bucketizer(input_col  = "Age",
##                 output_col = "AgeBucket",
##                 splits     = c(0, 10, 30, 90))
## 
## # Split data in training and testing set:
## partitions <- t %>%
##   sdf_random_split(training = 0.7, test = 0.3, seed = 1942)
## t_training <- partitions$training
## t_test     <- partitions$test
## 
## # Fit a logistic regression model
## M <- t_training %>%
##     ml_logistic_regression(Survived ~ AgeBucket + Pclass + Sex)
## 
## # Add predictions:
## pred <- sdf_predict(t_test, M)
## 
## # Show details of the quality of the fit of the model:
## ml_binary_classification_evaluator(pred)


## ----echo=FALSE---------------------------------------------------------------
system('stop-master.sh')
#detach("package:sparklyr", unload=TRUE)
detach("package:SparkR", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
#detach("package:dplyr", unload=TRUE)


## -----------------------------------------------------------------------------
x <- 1:1e4

# tmiming the function mean:
system.time(mean(x))


## -----------------------------------------------------------------------------
N <- 2500

# repeating it to gather longer time
system.time({for (i in 1:N) mean(x)})


## -----------------------------------------------------------------------------
system.time({for (i in 1:N) mean(x)})
system.time({for (i in 1:N) mean(x)})


## -----------------------------------------------------------------------------
# timing a challenger model
system.time({for (i in 1:N) sum(x) / length(x)})


## -----------------------------------------------------------------------------
N <- 1500

# Load microbenchmark:
library(microbenchmark)

# Create a microbenchmark object:
comp <- microbenchmark(mean(x),              # 1st code block
                       {sum(x) / length(x)}, # 2nd code block
		       times = N)            # number times to run


## -----------------------------------------------------------------------------
summary(comp)


## ----microbenchAuto,fig.cap='The package {\\tt microbenchmark} also defines a suible visualisation bia {\\tt autoplot}. This shows us the violin plots of the different competing methods. This also makes clear how heavy the right tails of the distribution of run-times are.',warning=FALSE,message=FALSE----
# Load ggplot2:
library(ggplot2)

# Use autoplot():
autoplot(comp)


## -----------------------------------------------------------------------------
x <- 1:1e+4
y <- 0

# Here we recalculate the sum multiple times.
system.time(for(i in x) {y <- y + sum(x)})

# Here we calculate it once and store it.
system.time({sum_x <- sum(x); for(i in x) {y <- y + sum_x}})


## -----------------------------------------------------------------------------
# Define some numbers.
x1 <- 3.00e8
x2 <- 663e-34
x3 <- 2.718282
x4 <- 6.64e-11

y1 <- pi
y2 <- 2.718282
y3 <- 9.869604
y4 <- 1.772454

N <- 1e5

# 1. Adding some of them directly:
f1 <- function () {
  for (i in 1:N) {
    x1 + y1
    x2 + y2
    x3 + y3
    x4 + y4
    }
  }
system.time(f1())

# 2. Converting first to a vector and then adding the vectors:
f2 <- function () {
  x <- c(x1, x2, x3, x4)
  y <- c(y1, y2, y3, y4)
  for (i in 1:N) x + y
}
system.time(f2())

# 3. Working with the elements of the vectors:
f3 <- function () {
  x <- c(x1, x2, x3, x4)
  y <- c(y1, y2, y3, y4)
  for (i in 1:N) {
    x[1] + y[1]
    x[2] + y[2]
    x[3] + y[3]
    x[4] + y[4]
  }
  }
system.time(f3())

# 4. Working with the elements of the vectors and code shorter:
f4 <- function () {
  x <- c(x1, x2, x3, x4)
  y <- c(y1, y2, y3, y4)
  for (i in 1:N) for(n in 1:4) x[n] + y[n]
  }
system.time(f4())


## ----ffff,fig.cap='This figure shows that working with vectors is fastest (f2()).',warning=FALSE,message=FALSE----
# Remind the packages:
library(microbenchmark)
library(ggplot2)

# Compare all four functions:
comp <- microbenchmark(f1(), f2(), f3(), f4(), times = 15)

# Create the violin plots:
autoplot(comp)


## ----appendTime---------------------------------------------------------------
N <- 2e4
   
# Method 1: using append():
system.time ({
    lst <- list()
    for(i in 1:N) {lst <- append(lst, pi)}
   })
    
# Method 2: increasing length while counting with length():
system.time ({
    lst <- list()
    for(i in 1:N) {
        lst[[length(lst) + 1]] <- pi}
   })
   
# Method 3: increasing length using a counter:
system.time ({
    lst <- list()
    for(i in 1:N) {lst[[i]] <- pi}
   })
   
# Method 4: pre-allocate memory:
system.time({
    lst <- vector("list", N)
    for(i in 1:N) {lst[[i]] <- pi}
    })


## -----------------------------------------------------------------------------
N <- 500
# simple operations on a matrix
M <- matrix(1:36, nrow = 6)
system.time({for (i in 1:N) {x1 <- t(M); x2 <- M + pi}})

# simple operations on a data-frame
D <- as.data.frame(M)
system.time({for (i in 1:N) {x1 <- t(D); x2 <- D + pi}})


## -----------------------------------------------------------------------------
x <- 1:1e4
N <- 1000
system.time({for (i in 1:N) mean(x)})
system.time({for (i in 1:N) sum(x) / length(x)})


## -----------------------------------------------------------------------------
N <- 732
# Use ts from stats to create the a time series object:
t1 <- stats::ts(rnorm(N), start=c(2000,1), end=c(2060,12), 
                frequency=12)
t2 <- stats::ts(rnorm(N), start=c(2010,1), end=c(2050,12), 
                frequency=12)

# Create matching zoo and xts objects:
zoo_t1 <- zoo::zoo(t1)
zoo_t2 <- zoo::zoo(t2)
xts_t1 <- xts::as.xts(t1)
xts_t2 <- xts::as.xts(t2)

# Run a merge on them.
# Note that base::merge() is a dispatcher function.
system.time({zoo_t <- merge(zoo_t1, zoo_t2)})
system.time({xts_t <- merge(xts_t1, xts_t2)})

# Calculate the lags.
system.time({for(i in 1:100) lag(zoo_t1)})
system.time({for(i in 1:100) lag(xts_t1)})


## ----eval=FALSE---------------------------------------------------------------
## # standard function:
## f1 <- function(n, x = pi) for(i in 1:n) x = 1 / (1+x)
## 
## # using curly brackets:
## f2 <- function(n, x = pi) for(i in 1:n) x = 1 / {1+x}
## 
## # adding unnecessary round brackets:
## f3 <- function(n, x = pi) for(i in 1:n) x = (1 / (1+x))
## 
## # adding unnecessary curly brackets:
## f4 <- function(n, x = pi) for(i in 1:n) x = {1 / {1+x}}
## 
## # achieving the same result by raising to a power
## f5 <- function(n, x = pi) for(i in 1:n) x = (1+x)^(-1)
## 
## # performing the power with curly brackets
## f6 <- function(n, x = pi) for(i in 1:n) x = {1+x}^{-1}
## 
## N <- 1e6
## library(microbenchmark)
## comp <- microbenchmark(f1(N), f2(N), f3(N), f4(N), f5(N), f6(N),
##                        times = 150)
## 
## comp
## ## Unit: milliseconds
## ##   expr      min       lq     mean   median       uq      max  ...
## ##   f1(N) 37.37476 37.49228 37.76950 37.57212 37.79876 39.99120 ...
## ##   f2(N) 37.29297 37.50435 37.79612 37.63191 37.81497 41.09414 ...
## ##   f3(N) 37.96886 38.18751 38.59619 38.28713 38.68162 47.66612 ...
## ##   f4(N) 37.88111 38.06787 38.41134 38.16297 38.36706 42.53103 ...
## ##   f5(N) 45.12742 45.31632 45.67364 45.45465 45.69882 49.65297 ...
## ##   f6(N) 45.93406 46.03159 46.51151 46.15287 46.64509 52.95426 ...
## 
## # Plot the results:
## library(ggplot2)
## autoplot(comp)


## ----eval=FALSE---------------------------------------------------------------
## library(compiler)
## N <- as.double(1:1e7)
## 
## # Create a *bad* function to calculate mean:
## f <- function(x) {
##   xx = 0
##   l = length(x)
##   for(i in 1:l)
##     xx = xx + x[i]/l
##   xx
## }
## 
## # Time the function:
## system.time(f(N))
## ##   user  system elapsed
## ##  0.61    0.000    0.61
## 
## # Compile the function:
## cmp_f <- cmpfun(f)
## 
## # Time the compiled version
## system.time(cmp_f(N))
## ##  user  system elapsed
## ##  0.596  0.00    0.596
## 
## # The difference is small, so we use microbenchmark
## library(microbenchmark)
## comp <- microbenchmark(f(N), cmp_f(N), times = 150)
## 
## # See the results:
## comp
## ## Unit: milliseconds
## ##      expr      min       lq     mean   median       uq      ...
## ##      f(N) 552.2785 553.9911 559.6025 556.1511 562.7207 601.5...
## ##  cmp_f(N) 552.2152 553.9453 558.1029 555.8457 560.0771 588.4...
## 
## 
## # Plot the results.
## library(ggplot2)
## autoplot(comp)


## -----------------------------------------------------------------------------
options(R_COMPILE_PKGS = 1)


## -----------------------------------------------------------------------------
options(R_ENABLE_JIT = 0)


## ----eval=FALSE---------------------------------------------------------------
## # Naive implementation of the Fibonacci numbers in R:
## Fib_R <- function (n) {
##   if ((n == 0) | (n == 1)) return(1)
##   return (Fib_R(n - 1) + Fib_R(n - 2))
##   }
## 
## # The R-function compiled via cmpfun():
## library(compiler)
## Fib_R_cmp <- cmpfun(Fib_R)
## 
## # Using native C++ via cppFunction():
## Rcpp::cppFunction('int Fib_Cpp(int n) {
##   if ((n == 0) || (n == 1)) return(1);
##   return (Fib_Cpp(n - 1) + Fib_Cpp(n - 2));
## }')
## 
## library(microbenchmark)
## N <- 30
## comp <- microbenchmark(Fib_R(N), Fib_R_cmp(N),
##                        Fib_Cpp(N), times = 25)
## 
## comp
## ## Unit: milliseconds
## ##          expr         min          lq        mean      median          uq
## ##      Fib_R(N) 1449.755022 1453.320560 1474.679341 1456.202559 1472.447928
## ##  Fib_R_cmp(N) 1444.145773 1454.127022 1489.742750 1459.170600 1554.450501
## ##    Fib_Cpp(N)    2.678766    2.694425    2.729571    2.711567    2.749208
## ##          max neval cld
## ##  1596.226483    25   b
## ##  1569.764246    25   b
## ##     2.858784    25  a
## 
## library(ggplot2)
## autoplot(comp)


## ----eval=FALSE---------------------------------------------------------------
## Fib_R2 <- function (n) {
##  x = 1
##  x_prev = 1
##  for (i in 2:n) {
##    x <- x + x_prev
##    x_prev = x
##    }
##  x
##  }
## 
## Rcpp::cppFunction('int Fib_Cpp2(int n) {
##   int x = 1, x_prev = 1, i;
##   for (i = 2; i <= n; i++) {
##     x += x_prev;
##     x_prev = x;
##     }
##   return x;
##   }')
## 
## N <- 30
## comp <- microbenchmark(Fib_R(N), Fib_R2(N),
##                        Fib_Cpp(N), Fib_Cpp2(N),
## 		       times = 20)
## 
## comp
## ## Unit: microseconds
## ##         expr         min           lq         mean     median ...
## ##     Fib_R(N) 1453850.637 1460021.5865 1.495407e+06 1471455.852...
## ##    Fib_R2(N)       2.057       2.4185 1.508404e+02       4.792...
## ##   Fib_Cpp(N)    2677.347    2691.5255 2.757781e+03    2697.519...
## ##  Fib_Cpp2(N)       1.067       1.4405 5.209175e+01       2.622...
## ##          max neval cld
## ##  1603991.462    20   b
## ##     2925.070    20  a
## ##     3322.077    20  a
## ##      964.378    20  a
## 
## library(ggplot2)
## autoplot(comp)


## ----eval=FALSE---------------------------------------------------------------
## sourceCpp("/path/cppSource.cpp")


## ----eval=FALSE---------------------------------------------------------------
## Rprof("/path/to/my/logfile")
## ... code goes here
## Rprof(NULL)


## -----------------------------------------------------------------------------
f0 <- function() x <- pi^2
f1 <- function() x <- pi^2 + exp(pi)
f2 <- function(n) for (i in 1:n) {f0(); f1()}
f3 <- function(n) for (i in 1:(2*n)) {f0(); f1()}
f4 <- function(n) for (i in 1:n) {f2(n); f3(n)}


## ----eval=FALSE---------------------------------------------------------------
## # Start the profiling:
## Rprof("prof_f4.txt")
## 
## # Run our functions:
## N <- 500
## f4(N)
## 
## # Stop the profiling process:
## Rprof(NULL)


## -----------------------------------------------------------------------------
# show the summary:
summaryRprof("prof_f4.txt")


## ----eval=FALSE---------------------------------------------------------------
## N <- 1000
## require(profr)
## pr <- profr({f4(N)}, 0.01)
## plot(pr)


## ----flameCallee,fig.cap=c('A flame-plot produced by the function {\\tt flameGraph()} from the package {\\tt proftools}.','A Callee Tee Map produced by the function {\\tt calleeTreeMap()} from the package {\\tt proftools}. This visualization shows boxes with surfaces that are relative to the time that the function takes.')----
library(proftools)

# Read in the existing profile data from Rprof:
pd <- readProfileData("prof_f4.txt")

# Print the hot-path:
hotPaths(pd, total.pct = 10.0)

# A flame-graph (stacked time-bars: first following plot)
flameGraph(pd)

# Callee tree-map (intersecting boxes with area related to time
# spent in a function: see plot below)
calleeTreeMap(pd)


## ----eval=FALSE---------------------------------------------------------------
## unlink("prof_f4.txt")


## ----eval=FALSE---------------------------------------------------------------
## install.packages("devtools")
## install.packages("roxygen2")


## -----------------------------------------------------------------------------
library(devtools)

## ----eval=FALSE---------------------------------------------------------------
## library(roxygen2)


## ----eval=FALSE---------------------------------------------------------------
## setwd(".")  # replace this with your desired directory
## 
## # Step 1: define the function(s)
## # -- below we repeat the define the function as defined earlier.
## # diversity
## # Calculates the entropy of a system with discrete states.
## # Arguments:
## #   x     -- numeric vector -- observed probabilities of classes
## #   prior -- numeric vector -- prior probabilities of the classes
## # Returns:
## #    numeric -- the entropy / diversity measure
## diversity <- function(x, prior = NULL) {
##   if (min(x) <= 0) {return(0);} # the log will fail for 0
##   # if the numbers are higher than 1, then not probabilities but
##   # populations are given, so we rescale to probabilities:
##   if (sum(x) != 1) {x <- x / sum(x)}
##   N <- length(x)
##   if(!is.null(prior)) {
##     for (i in (1:N)) {
##       a <- (1 - 1 / (N * prior[i])) / (1 - prior[i])
##       b <- (1 - N * prior[i]^2) / (N * prior[i] * (1 - prior[i]))
##       x[i] <- a * x[i]^2 + b * x[i]
##     }
##    }
##   f <- function(x) x * log(x)
##   x1 <- mapply(FUN = f, x)
##   - sum(x1) / log(N)
##   }


## ----eval=FALSE---------------------------------------------------------------
## # Step 2: create the empty package directory:
## package.skeleton(list = c("diversity"), # list all the functions
##                                         # from the current
## 					# environment that we want
## 					# in the package.
##                 name = "div"            # name of the package
## 		)


## -----------------------------------------------------------------------------
list.files(path = "div", all.files = TRUE)


## ----eval=FALSE---------------------------------------------------------------
## ?diversity


## -----------------------------------------------------------------------------
#' Function to calculate a diversity index based on entropy.
#'
#' This function returns entropy of a system with discrete states
#' @param x numeric vector, observed probabilities
#' @param prior numeric vector, the prior probabilities
#' @keywords diversity, entropy
#' @return the entropy or diversity measure
#' @examples
#' x <- c(0.4, 0.6)
#' diversity(x)


## ----eval=FALSE---------------------------------------------------------------
## setwd("./div")  # or your choice of package directory
## document()


## ----eval=FALSE---------------------------------------------------------------
## Warning: The existing 'NAMESPACE' file was not generated by roxygen2,
## and will not be overwritten.
## Warning: The existing 'diversity.Rd' file was not generated by
## roxygen2, and will not be overwritten.


## ----echo=FALSE---------------------------------------------------------------
setwd("~/Documents/science/books/R-book")


## ----echo=TRUE,results='hide',message=FALSE,warning=FALSE---------------------
setwd(".") # the directory in which ./div is a subdirectory
install("div")


## -----------------------------------------------------------------------------
x  <- c(0.3, 0.2, 0.5)
pr <- c(0.33, 0.33, 0.34)
diversity(x, prior = pr)


## ----eval=FALSE---------------------------------------------------------------
## ?diversity


## ----eval=FALSE---------------------------------------------------------------
## install_github('div','your_github_username')


## -----------------------------------------------------------------------------
# When you use a R-package, you should cite it in your work.
# The information can be found as follows:
citation('base')

# You can even extract the BibTex information:
toBibtex(citation('base'))


## ----snipBin,fig.cap='A visual aid to select binning borders is plotting a non-parametric fit and the histogram.'----
# -- remind the data fabrication --- 
set.seed(1865)
age <- rlnorm(1000, meanlog = log(40), sdlog = log(1.3))
y <- rep(NA, length(age))
for(n in 1:length(age)) {
  y[n] <- max(0, 
              dnorm(age[n], mean= 40, sd=10) 
	         + rnorm(1, mean = 0, sd = 10 * dnorm(age[n], 
		   mean= 40, sd=15)) * 0.075)
}
y <- y / max(y)
dt <- tibble (age = age, spending_ratio = y)

# -- plotting loesss and histogram with ggplot2 ---
library(ggplot2)
library(gridExtra)
p1 <- ggplot(dt, mapping = aes(x = age, y = spending_ratio)) +
     geom_point(alpha=0.2, size = 2) +
     geom_smooth(method = "loess")
p2 <- qplot(dt$age, geom="histogram", binwidth=5) 
grid.arrange(p1, p2, ncol=2)


## ----eval=FALSE,fig.cap="Some plot characters. Most other characters will just plot themselves."----
## # This sets up an empty plotting field
## plot(x = c(0, 4.5),
##      y = c(0, 5),
##      main = "Some pch arguments",
##      xaxt = "n",
##      yaxt = "n",
##      xlab = "",
##      ylab = "",
##      cex.main = 2.6,
##      col = "white"
## )
## 
## # This will plot all of the standard pch arguments
## y = rep(5:0, each=5)
## for (i in 0:25) {
##   points(x = i %% 5, y = y[i+1], pch = i,cex = 2, col="blue", bg="khaki3")
##   text(0.3 + i %% 5, y = y[i+1], i, cex = 2)
## }
## for (i in 1:2) {
##   ch <- LETTERS[i]
##   points(x = i, y = 0, pch = ch,cex = 2, col="red")
##   text(0.3 + i, y = 0, ch, cex = 2)
## }
## for (i in 1:2) {
##   ch <- letters[i]
##   points(x = i + 2, y = 0, pch = ch,cex = 2, col="red")
##   text(0.3 + i + 2, y = 0, ch, cex = 2)
## }


## ----KSvisCodeOnly,fig.cap='The KS as the maximum distance between the cumulative distributions of the positive and negative observations.',eval=FALSE----
## # Using model m and data frame t2:
## predicScore <- predict(object=m,type="response")
## d0 <- data.frame(
##        score = as.vector(predicScore)[t2$Survived == 0],
##        true_result = 'not survived')
## d1 <- data.frame(
##        score = as.vector(predicScore)[t2$Survived == 1],
##        true_result = 'survived')
## d  <- rbind(d0, d1)
## d  <- d[complete.cases(d),]
## 
## cumDf0 <- ecdf(d0$score)
## cumDf1 <- ecdf(d1$score)
## x      <- sort(d$score)
## cumD0  <- cumDf0(x)
## cumD1  <- cumDf1(x)
## diff   <- cumD0 - cumD1
## y1     <- gdata::first(cumD0[diff == max(diff)])
## y2     <- gdata::first(cumD1[diff == max(diff)])
## x1     <- x2 <- quantile(d0$score, probs=y1, na.rm=TRUE)
## 
## # Plot this with ggplot2:
## p <- ggplot(d, aes(x = score)) +
##      stat_ecdf(geom = "step", aes(col = true_result), lwd=2) +
##      ggtitle('Cummulative distributions and KS') +
##      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
##                   color='navy', lwd=3) +
##      ggplot2::annotate("text",
##               label = paste0("KS=",round((y1-y2)*100,2),"%"),
##               x = x1 + 0.15, y = y2+(y1-y2)/2, color = "navy")
## p


## ----eval=FALSE---------------------------------------------------------------
## library(ggplot2)
## library(latex2exp)
## d <- seq(from = -3, to = +3, length.out = 100)
## 
## ## error function family:
## erf     <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## # (see Abramowitz and Stegun 29.2.29)
## erfc    <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## erfinv  <- function (x) qnorm((1 + x)/2)/sqrt(2)
## erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)
## 
## ## Gudermannian function
## gd <- function(x) asin(tanh(x))
## 
## f1 <- function(x) erf( sqrt(pi) / 2 * x)
## f2 <- function(x) tanh(x)
## f3 <- function(x) 2 / pi * gd(pi / 2 * x)
## f4 <- function(x) x / sqrt(1 + x^2)
## f5 <- function(x) 2 / pi * atan(pi / 2 * x)
## f6 <- function(x) x / (1 + abs(x))
## 
## df <- data.frame(d = d, y = f1(d),
##                  Function = "erf( sqrt(pi) / 2 * d)")
## df <- rbind(df, data.frame(d = d, y = f2(d), Function = "tanh(d)"))
## df <- rbind(df, data.frame(d = d, y = f3(d),
##                            Function = "2 / pi * gd(pi / 2 * d)"))
## df <- rbind(df, data.frame(d = d, y = f4(d),
##                            Function = "d / (1 + d^2)"))
## df <- rbind(df, data.frame(d = d, y = f5(d),
##                            Function = "2 / pi * atan(pi / 2 * d)"))
## df <- rbind(df, data.frame(d = d, y = f6(d),
##                            Function = "x / (1 + abs(d))"))
## 
## fn <- ""
## fn[1] <- "erf \\left(\\frac{\\sqrt{\\pi} d}{2}\\right)"
## fn[2] <- "tanh(x)"
## fn[3] <- "\\frac{2}{\\pi} gd\\left( \\frac{\\pi d}{2} \\right)"
## fn[4] <- "\\frac{d}{1 + d^2}"
## fn[5] <- "\\frac{2}{\\pi} atan\\left(\\frac{\\pi d}{2}\\right)"
## fn[6] <- "\\frac{x}{1+ |x|}"
## 
## 
## ggplot(data = df, aes(x = d, y = y, color = Function)) +
##    geom_line(aes(col = Function), lwd=2) +
##    guides(color=guide_legend(title=NULL)) +
##    scale_color_discrete(labels=lapply(sprintf('$\\pi(d) = %s$', fn),
##                         TeX)) +
##    theme(legend.justification = c(1, 0),
##         legend.position = c(1, 0),   # south east
##         legend.box.margin=ggplot2::margin(rep(20, times=4)),
## 	# increase vertical space between legend items:
##         legend.key.size = unit(1.5, "cm")
## 	) +
##    ylab(TeX('Preference --- $\\pi$'))


## -----------------------------------------------------------------------------
f <- function(x) x^2 + 1
f
Vectorize(f)


## ----eval=FALSE---------------------------------------------------------------
## f_curve <- function(f) {
##   g <- Vectorize(f)
##   s <- deparse(f)[2]
##    curve(g, xlab = '', ylab = '', col = 'red', lwd = 3,
##          from = -1, to = +1,
## #-1-         main = bquote(bold(.(s))))
##           main = s
## 	 )
##    }
## 
## gaus <- function(x) exp (-(x-0)^2 / 0.5)
## f1 <- function(x) - 3/2 * x^5 + 5/2 * x^3
## f2 <- function(x) sin(pi * x / 2)
## f3 <- function(x) min(1, max(2*x, -1))
## f4 <- function(x) x
## f5 <- function(x) ifelse(x < 0 , gaus(x) - 1, 1 - gaus(x))
## f6 <- function(x) tanh(3 * x)
## 
## par(mfrow=c(3,2))
## f_curve(f1)
## f_curve(f2)
## f_curve(f3)
## f_curve(f4)
## f_curve(f5)
## f_curve(f6)
## par(mfrow=c(1,1))


## ----eval=FALSE---------------------------------------------------------------
## nottemC <- 5/9 * nottem - 32


## -----------------------------------------------------------------------------
dotproduct <- function (m1, m2) {
  if(!is.matrix(m1) || !is.matrix(m2)) {
    print("ERROR 1: m1 and m2 must be matrices.");
    return(-1);
  }
  if (ncol(m1) != nrow(m2)) {
    print("ERROR 2: ncol(m1) must match nrow(m2).")
    return(-2);
  }
  # Note that we do not check for values NA or Inf. In R they work
  # as one would expect.
  m <- matrix(nrow=nrow(m1), ncol=ncol(m2))
  for (k in 1:nrow(m1)) {
    for (l in 1:ncol(m1)){
      m[k,l] <- 0
      for (n in 1:ncol(m1)) {
        m[k,l] <- m[k,l] + m1[k,n] * m2[n,l]
        }
      }
  }
  return(m);
}
  
M1 <- matrix(1:6, ncol=2)
M2 <- matrix(5:8, nrow=2)

# Try our function
dotproduct(M1, M2)

# Compare with
M1 %*% M2  

# The following must fail:
dotproduct("apple", M2)

# This must fail too:
dotproduct(as.matrix(pi),M1)


## -----------------------------------------------------------------------------
# Loading mtcars  
library(MASS)     # However, this is probably already loaded

# Exploring the number of gears
summary(mtcars$gear)

## show the histogram
# since the only valuse are 3,4 and 5 we format the 
# breaks around those values.
hist(mtcars$gear, col="khaki3", 
     breaks=c(2.5, 2.9,3.1,3.9,4.1,4.9,5.1, 5.5))

# Study the correlation between gears and transmission
# (see section on correlations)
cor(mtcars$gear, mtcars$am)
cor(rank(mtcars$gear), rank(mtcars$am))

# We can also show the correlation between gears and transmission
# Since the data is mainly overlapping we use the function 
# jitter() to add some noise so we can see individual points.
plot(jitter(mtcars$am,0.125), jitter(mtcars$gear,0.3), 
            pch=21, col="blue", bg="red",
    xlab = "Transmission (0 = automatic, 1 = manual)",
    ylab = "Number of forward gears",
    main = "Jittered view of gears and transmission"
    )
# Add a blue line (linear fit -- see linear regression)
abline(lm(mtcars$gear ~ mtcars$am),
       col='blue',lwd=3)


## -----------------------------------------------------------------------------
# Create the factors, taking care that we know what will
# be manual and what will be automatic.
f <- factor(mtcars$am, levels=c(0,1),labels=c("automatic", "manual"))

# To manually inspect the result:
head(cbind(f,mtcars$am))

# Show off the result
plot(f, col="khaki3")


## -----------------------------------------------------------------------------
# show the distribution
hist(mtcars$hp, col="khaki3")

# Try a possible cut
table(cut(mtcars$hp, breaks=c(50,100,175,350)))

# Assume we like the cut and use this one:
c <- cut(mtcars$hp, breaks=c(50,100,175,350)) # add the cut
l <- unique(c)                                # find levels
l                                             # check levels order

# Note that we provide the labels in the order of levels
f <- factor(c, levels=l, labels=c("M", "L", "H"))
plot(f, col="khaki3", 
     main="Horsepower as a factor",
     xlab="bins: L=low, M=Medium, H=high",
     ylab="number of cars in that bin")


## -----------------------------------------------------------------------------
M <- matrix(c(1:9),nrow=3); 
D <- data.frame(M); 
rownames(D) <- c("2016","2017","2018"); 
colnames(D) <- c("Belgium", "France", "Poland"); 
cbind(D,rowSums(D)); 
D <- D[,-2]


## ----warning=FALSE------------------------------------------------------------
# Start a standard input
t <- read_csv(readr_example("challenge.csv"))

# Then, to see the issues, do:
problems(t)

# Notice that the problems start in row 1001, so
# the first 1000 rows are special cases. The first improvement
# can be obtained by increase the guesses
## compare 
spec_csv(readr_example("challenge.csv"))
## with 
spec_csv(readr_example("challenge.csv"), guess_max = 1001)

# second step:
t <- read_csv(readr_example("challenge.csv"), guess_max = 1001)

# Let us see:
head(t)
tail(t)


## -----------------------------------------------------------------------------
# We use the same seed so our results will be comparable.
set.seed(1865)

# We use the function mutate() from dplyr:
library(dplyr)

# For completeness sake, we generate the same data again.
N <- 500
age_f   <- rlnorm(N, meanlog = log(40), sdlog = log(1.3))
x_f <- abs(age_f + rnorm(N, 0, 20))    # Add noise & keep positive
x_f <- 1 - (x_f - min(x_f)) / max(x_f) # Scale between 0 and 1
x_f <- 0.5 * x_f / mean(x_f)           # Coerce mean to 0.5
# This last step will produce some outliers above 1
x_f[x_f > 1] <- 1                      # Coerce those > 1 to 1

age_m   <- rlnorm(N, meanlog = log(40), sdlog = log(1.3))
x_m <- abs(age_m + rnorm(N, 0, 20))    # Add noise & keep positive
x_m <- 1 - (x_m - min(x_m)) / max(x_m) # Scale between 0 and 1
x_m <- 0.5 * x_m / mean(x_m)           # Coerce mean to 0.5
# This last step will produce some outliers above 1
x_m[x_m > 1] <- 1                      # Coerce those > 1 to 1
x_m <- 1 - x_m                         # Make the relation increasing

p_f <- x_f
p_m <- x_m

tf <- tibble("age" = age_f, "sex" = "F", "is_good" = p_f) 
tm <- tibble("age" = age_m, "sex" = "M", "is_good" = p_m)
t  <- full_join(tf, tm, by = c("age", "sex", "is_good")) %>%
      mutate("is_good" = if_else(is_good >= 0.5, 1L, 0L))%>%  
      mutate("sexM"    = if_else(sex == "F", 0, 1)) 

###########
# Model 1 #
###########
regr1 <- glm(formula = is_good ~ age + sexM, 
             family = binomial,
             data = t)

# assess the model:
summary(regr1)

pred1 <- 1 / (1+ exp(-(coef(regr1)[1] + t$age * coef(regr1)[2] 
                     + t$sexM * coef(regr1)[3])))
SE1 <-  (pred1 - t$is_good)^2
MSE1 <- sum(SE1) / length(SE1)

###########
# Model 2 #
###########
# make the same cut
t <- as_tibble(t)                                               %>%
    mutate(is_LF = if_else((age <= 35) & (sex == "F"), 1L, 0L)) %>%
    mutate(is_HF = if_else((age >  50) & (sex == "F"), 1L, 0L)) %>%
    mutate(is_LM = if_else((age <= 35) & (sex == "M"), 1L, 0L)) %>%
    mutate(is_HM = if_else((age >  50) & (sex == "M"), 1L, 0L))

regr2 <- glm(formula = is_good ~ is_LF + is_HF + is_LM + is_HM,
             family = binomial,
             data = t)

# assess the model:
summary(regr2)

pred2 <- 1 / (1 + exp(-(coef(regr2)[1] + 
                     + t$is_LF * coef(regr2)[2] 
		     + t$is_HF * coef(regr2)[3] 
                     + t$is_LM * coef(regr2)[4] 
		     + t$is_HM * coef(regr2)[5] 
                      )))
SE2 <-  (pred2 - t$is_good)^2
MSE2 <- sum(SE2) / length(SE2)

# Finally, we also note that the MSE has improved too.
MSE1
MSE2


## -----------------------------------------------------------------------------
# mcda_promethee_list
# delivers the preference flow matrices for the Promethee method
# Arguments:
#    M      -- decision matrix
#    w      -- weights
#    piFUNs -- a list of preference functions, 
#              if not provided min(1,max(0,d)) is assumed.
# Returns (as side effect)
# phi_plus <<- rowSums(PI.plus)
# phi_min  <<- rowSums(PI.min)
# phi_     <<- phi_plus - phi_min
#
mcda_promethee_list <- function(M, w, piFUNs='x')
{
  if (piFUNs == 'x') {
       # create a factory function:
       makeFUN <- function(x) {x; function(x) max(0,x) }
       P <- list()
       for (k in 1:ncol(M)) P[[k]] <- makeFUN(k)
       } # else, we assume a vector of functions is provided
# initializations
PI.plus  <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
PI.min   <<- matrix(data=0, nrow=nrow(M), ncol=nrow(M))
# calculate the preference matrix
for (i in 1:nrow(M)){
  for (j in 1:nrow(M)) {
    for (k in 1:ncol(M)) {
      if (M[i,k] > M[j,k]) {
        PI.plus[i,j] = PI.plus[i,j] + w[k] * P[[k]](M[i,k] - M[j,k])
      }
      if (M[j,k] > M[i,k]) {
        PI.min[i,j] = PI.min[i,j] + w[k] * P[[k]](M[j,k] - M[i,k])
      }
    }
  }
}
# Till this line the code was exactly the same.
# In the following line we create the list to return.
L <- list("phi_plus" <- rowSums(PI.plus), 
          "phi_min"  <- rowSums(PI.min),
          "phi_"     <- phi_plus - phi_min
	  )
return(L)
}


## -----------------------------------------------------------------------------
set.seed(1492)
M <- matrix (runif(9), nrow = 3)
w <- c(runif(3))
L <- mcda_promethee_list(M, w)


## -----------------------------------------------------------------------------
IV(X = factor(mtcars$vs), Y = factor(mtcars$am))
WOETable(X = factor(mtcars$vs), Y = factor(mtcars$am))
