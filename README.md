# Econ565
Here I provide the codes that I developed for Econ 565 by Sam Hwang at University of British Columbia.

## Exercise 1
Econ565_exercise1.ipynb: I code a monte carlo simulation to produce data for a two-alternative discrete choice model. Then I use the data to estimate the parameters of the model.

## Exercise 2
data.jl: It generates the data based on the model provided in section 8 of Berry (1994) RAND. Then, I use this data to estimate the parameters of the model

main.do: It uses data from "data.jl" to replicate the Table 1 of Berry (1994). The estimates for column 1 and 2 are very similar to estimates reported in the paper, however, the estimates in column 3 and 4 are not.
