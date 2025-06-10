* ==============================================================================
*  COHORT APPORTIONMENT COMPARISON
use "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_age_sex_2006-2023", clear

* Prepare the data
* ---------------
* Ensure variables are properly formatted
destring monthdth fem age year tot_deaths, replace ignore("NA" ".")

* Create birth cohort identifiers
gen cohort_young = year - age - 1  // Born later in birth year (younger cohort)
gen cohort_old = year - age        // Born earlier in birth year (older cohort)

* ==============================================================================
* UNIFORM BIRTH DISTRIBUTION ASSUMPTION  
* ==============================================================================
* Uniform birth distribution: proportion depends on how much of the year remains
* For someone dying in month M at age A, what's the probability they were born in cohort (year-age-1)?
* This equals the fraction of their birth year that falls AFTER month M of the previous year
gen prop_uniform_young = (2*monthdth - 1) / 24
gen prop_uniform_old = 1 - prop_uniform_young

gen deaths_uniform_young = tot_deaths * prop_uniform_young
gen deaths_uniform_old = tot_deaths * prop_uniform_old

* ==============================================================================
* ORGANIZE AND LABEL VARIABLES
* ==============================================================================

* Label variables
label var year "Year of death"
label var monthdth "Month of death"
label var age "Age at death (years)"
label var fem "Female indicator"
label var tot_deaths "Total deaths"

label var cohort_young "Younger birth cohort (year-age+1)"
label var cohort_old "Older birth cohort (year-age)"

* ==============================================================================
* SAVE YOUNGER COHORT AS SEPARATE FILE
* ==============================================================================

preserve
keep monthdth fem year cohort_young deaths_uniform_young
ren cohort_young cohort
tempfile young
save `young'
restore

* ==============================================================================
* MERGE WITH OLDER COHORT
* ==============================================================================

keep monthdth fem year cohort_old deaths_uniform_old
ren cohort_old cohort

merge 1:1 monthdth fem year cohort using `young'

* ==============================================================================
* SUM FOR COHORT TOTALS
* ==============================================================================

egen tot_deaths_app = rsum(deaths_uniform_*)

*twoway line tot_deaths year if cohort == 1950 & monthdth ==1 & fem ==1
keep monthdth fem year cohort tot_deaths
save "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_cohort_sex_2006-2023", replace

* ==============================================================================
* Check difference with age data
* ==============================================================================

use "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_age_sex_2006-2023", clear
g cohort = year-age
merge 1:1 monthdth cohort year fem using "/Users/nicholasmark/Dropbox/covid/2025/data/deaths_cohort_sex_2006-2023"
sort fem age year month

g diff = tot_deaths-tot_deaths_app
twoway line diff year if cohort == 1932 & monthdth ==1 & fem ==1
