/*******************************************************/
/*******************************************************/
/********************SET UP ****************************/
clear
clear matrix
cd "C:\Users\amirh\Dropbox\PhD\3\Econ 565\Exercise 2"
cap log close
log using "Berry_1994RAND_discrete.log", replace

set more off
pause on 

/********************DATA WITH SIGMA_D = 1 ****************************/
import delimited data_market.csv, encoding(Big5) clear

rename v1 y
rename v3 x
rename v4 p
rename v2 market
rename v5 w
rename v6 c

local iter 	= 100
local M 	= 500
local iter1 = 2
mat r_ols	= J(6,`iter',.)
mat r_iv 	= J(6,`iter',.)


// Drop invalid draws in the sample
drop if market == 0			//Invalid draw in which the monte carlo didn't have an answer
drop if y == "-Inf" | y == "Inf" //Share of one of the products is zero
destring, replace
duplicates tag market, g(tag)
drop if tag == 0
drop tag
unique market
local all = 2*`M'*`iter'
dis `all'
keep in 1/`all'			//Keep only 2*500*100 = M*iter observations
gen obs = _n

gen temp1 = x[_n - 1]
gen temp2 = x[_n + 1]
replace temp1 = . if market ~= market[_n-1]
replace temp2 = . if market ~= market[_n+1]

gen x_rival = temp1
replace x_rival = temp2 if x_rival == .
drop temp1 temp2

/******************** REGRESSIOINS ****************************/
forvalues i = 1(1)`iter'{
	local l = 2*`M'*(`i' - 1) + 1
	local h = 2*`M'*`i'
	qui reg y x p if obs >= `l' & obs <= `h'
	mat r_ols[1, `i'] = -_b[p]
	mat r_ols[2, `i'] = _se[p]
	mat r_ols[3, `i'] = _b[x]
	mat r_ols[4, `i'] = _se[x]
	mat r_ols[5, `i'] = _b[_cons]
	mat r_ols[6, `i'] = _se[_cons]

	qui ivreg y x (p = w x_rival) if obs >= `l' & obs <= `h'
	mat r_iv[1, `i'] = -_b[p]
	mat r_iv[2, `i'] = _se[p]
	mat r_iv[3, `i'] = _b[x]
	mat r_iv[4, `i'] = _se[x]
	mat r_iv[5, `i'] = _b[_cons]
	mat r_iv[6, `i'] = _se[_cons]
}

mata : st_matrix("mean_ols", rowsum(st_matrix("r_ols")))
mat mean_ols = mean_ols/`iter'
mata : st_matrix("mean_iv", rowsum(st_matrix("r_iv")))
mat mean_iv = mean_iv/`iter'

matlist mean_ols
matlist mean_iv


/********************DATA WITH SIGMA_D = 3 ****************************/
import delimited data_market3.csv, encoding(Big5) clear

rename v1 y
rename v3 x
rename v4 p
rename v2 market
rename v5 w
rename v6 c

local iter 	= 100
local M 	= 500
local iter1 = 2
mat r_ols	= J(6,`iter',.)
mat r_iv 	= J(6,`iter',.)


// Drop invalid draws in the sample
drop if market == 0			//Invalid draw in which the monte carlo didn't have an answer
drop if y == "-Inf" | y == "Inf" //Share of one of the products is zero
destring, replace
duplicates tag market, g(tag)
drop if tag == 0
drop tag
unique market
local all = 2*`M'*`iter'
dis `all'
keep in 1/`all'			//Keep only 2*500*100 = M*iter observations
gen obs = _n

gen temp1 = x[_n - 1]
gen temp2 = x[_n + 1]
replace temp1 = . if market ~= market[_n-1]
replace temp2 = . if market ~= market[_n+1]

gen x_rival = temp1
replace x_rival = temp2 if x_rival == .
drop temp1 temp2

/******************** REGRESSIOINS ****************************/
forvalues i = 1(1)`iter'{
	local l = 2*`M'*(`i' - 1) + 1
	local h = 2*`M'*`i'
	qui reg y x p if obs >= `l' & obs <= `h'
	mat r_ols[1, `i'] = _b[p]
	mat r_ols[2, `i'] = _se[p]
	mat r_ols[3, `i'] = _b[x]
	mat r_ols[4, `i'] = _se[x]
	mat r_ols[5, `i'] = _b[_cons]
	mat r_ols[6, `i'] = _se[_cons]

	qui ivreg y x (p = w x_rival) if obs >= `l' & obs <= `h'
	mat r_iv[1, `i'] = _b[p]
	mat r_iv[2, `i'] = _se[p]
	mat r_iv[3, `i'] = _b[x]
	mat r_iv[4, `i'] = _se[x]
	mat r_iv[5, `i'] = _b[_cons]
	mat r_iv[6, `i'] = _se[_cons]
}

mata : st_matrix("mean_ols", rowsum(st_matrix("r_ols")))
mat mean_ols = mean_ols/`iter'
mata : st_matrix("mean_iv", rowsum(st_matrix("r_iv")))
mat mean_iv = mean_iv/`iter'

matlist mean_ols
matlist mean_iv

log close
