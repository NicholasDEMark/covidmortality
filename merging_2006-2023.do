
*********SEER Pops***********
clear
infix year 1-4 st_fips 7-8 ct_fips 9-11  race 14 sex 16 age 17-18 population 19-29 ///
using /Users/nicholasmark/Dropbox/data/populations/us.1969_2023.singleages.through89.90plus.adjusted.txt"


label define race 1 "White" 2 "Black" 3 "Other"
label values race race

bys year sex age: egen tot_pop = total(population) 
keep year sex age tot_pop
duplicates drop
save "/Users/nicholasmark/Dropbox/covid/2025/data/sex_age_1969-2023"

use "/Users/nicholasmark/Dropbox/covid/2025/data/sex_age_1969-2023", clear
sort sex year age 
bys sex: g jan_pop = (tot_pop+tot_pop[_n-1])/2

expand 12
g cohort = year-age
bys sex year cohort: g month = _n

gen date = mdy(month, 15, year)  // 15th of each month
format date %td

replace jan_pop = . if month >1

bys sex cohort: ipolate jan_pop date, g(newpop)
save "/Users/nicholasmark/Dropbox/covid/2025/data/sex_cohort_1969-2023"
**************************

use "/Users/nicholasmark/Dropbox/covid/2025/data/sex_cohort_1969-2023", clear
g fem = (sex ==2)
keep month fem year cohort newpop
ren month monthdth
merge 1:1 cohort year monthdth fem using "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_cohort_sex_2006-2023", keep(3)
gen cumu_month = (year - 2006) * 12 + month
g mort_rate = 100000*(tot_deaths/newpop)
g age = year-cohort
keep cohort cumu_month year age month tot_deaths mort_rate fem newpop
order cohort cumu_month year age newpop month tot_deaths mort_rate fem 


preserve
keep if fem ==0 
drop fem
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023.csv", replace
restore

preserve
keep if fem ==1
drop fem
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023.csv", replace
restore

twoway line mort_rate cumu_month if cohort == 1950, by(fem)

/* OLD ESTIMATES USING ACS
*************
*Overall
************

*The population file does have numbers for col and non col, it just has %
*So tot deaths/tot pop is mortality rate
use "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_age_sex_2006-2023", clear
merge m:1 year age fem using "/Users/nicholasmark/Dropbox/covid/2025/data/ACS_educ_sex_age_2006-2023"
g cohort = year-age
g mort_rate = 100000*(tot_deaths/tot_pop)
gen cumu_month = (year - 2006) * 12 + month
sort age cumu_month
twoway line mort_rate cumu_month if age == 80, by(sex)


preserve
keep if fem ==0
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023.csv", replace
restore

preserve
keep if fem ==1
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023.csv", replace
restore

*************
*By Education
*************
use "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_age_sex_educ_2006-2023", clear
merge m:1 year age fem using "/Users/nicholasmark/Dropbox/covid/2025/data/ACS_educ_sex_age_2006-2023"
g cohort = year-age
gen cumu_month = (year - 2006) * 12 + month

g pop_col1 = tot_pop*pct_col
g pop_col0 = tot_pop*(1-pct_col)

g mort_rate = 100000*(tot_deaths/pop_col1) if col ==1
replace mort_rate = 100000*(tot_deaths/pop_col0) if col ==0

drop if col ==.

sort age cumu_month
twoway line mort_rate cumu_month if age == 80, by(sex col)


foreach col in 0 1 {
preserve
keep if fem ==1 & col ==`col'
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_col`col'_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_col`col'_2006-2023.csv", replace
restore
}

foreach col in 0 1 {
preserve
keep if fem ==0 & col ==`col'
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_col`col'_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_col`col'_2006-2023.csv", replace
restore
}

/* The Race denominators need some editing - with the 2020 Census, people could select multiple races, so need to go back and recode things for other races. 
*************
*By Race
************
use "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_age_sex_race_2006-2023", clear
merge m:1 year age fem race using "/Users/nicholasmark/Dropbox/covid/2025/data/ACS_race_sex_age_2006-2023"
g cohort = year-age
g mort_rate = 100000*(tot_deaths/tot_pop)
gen cumu_month = (year - 2006) * 12 + month
sort age cumu_month
twoway line mort_rate year if age == 80, by(sex race)


preserve
keep if fem ==0
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_2006-2023.csv", replace
restore

preserve
keep if fem ==1
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_female_2006-2023.csv", replace
restore



foreach race in 0 1 {
preserve
keep if fem ==0 & race ==`race'
keep age year cumu_month month mort_rate cohort
order cohort cumu_month year month age mort_rate
save "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_race`race'_2006-2023", replace
export delimited "/Users/nicholasmark/Dropbox/covid/2025/data/mort_rates_male_race`race'_2006-2023.csv", replace
restore
}

