******************************************************************************************************************************************************************************************
* Non-Manufacturing Oil Exposure
******************************************************************************************************************************************************************************************

*** Convert 6-digit NAICS to 3-digit NAICS ***
use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\czone_industry1990.dta" 
merge m:m ussic87 using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\naicstosic87.dta"
drop _merge
drop if ussic87 == .
drop if usnaics == .
drop if imp_emp1990 == .
gen emp1990 = weight*imp_emp1990
sort czone usnaics
gen double naics = floor(usnaics/1000)
collapse (sum) emp1990, by(naics czone)
sort naics czone
drop if naics > 310 & naics < 340
gen match = 2

*** Calculate l_{j,i,t}/l_{i,t} ***
egen emp_czone = sum(emp1990), by(czone)
gen cz_share = emp1990/emp_czone

sort czone naics

*** Matching dataset with change in labor due to oil shocks ***
merge m:m match using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor Due to Oil 1990-2000.dta"
drop if yr == .
drop if czone == .
drop _merge
sort naics czone

*** Create oil exposure variable ***
egen emp_ind = sum(emp1990), by(naics)

egen total_emp_ind = sum(emp1990), by(match)

gen ind_share = emp_ind/total_emp_ind
gen dl_nonman = cz_share * dl * ind_share

*** Sum up across all industries for each CZ ***
collapse (sum) dl_nonman, by(czone)
gen yr = 1990
save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1990-2000.dta", replace

*** Repeat the same exercise for 2000 ***

*** Convert 6-digit NAICS to 3-digit NAICS ***
use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\czone_industry2000.dta" 
merge m:m ussic87 using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\naicstosic87.dta"
drop _merge
drop if ussic87 == .
drop if usnaics == .
drop if imp_emp2000 == .
gen emp2000 = weight*imp_emp2000
sort czone usnaics
gen double naics = floor(usnaics/1000)
collapse (sum) emp2000, by(naics czone)
sort naics czone
drop if naics > 310 & naics < 340
gen match = 2

*** Calculate l_{j,i,t}/l_{i,t} ***
egen emp_czone = sum(emp2000), by(czone)
gen cz_share = emp2000/emp_czone

*** Matching dataset with change in labor due to oil shocks ***
merge m:m match using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor Due to Oil 2000-2007.dta"
drop if yr == .
drop if czone == .
drop _merge
sort naics czone

*** Create oil exposure variable ***
egen emp_ind = sum(emp2000), by(naics)

egen total_emp_ind = sum(emp2000), by(match)

gen ind_share = emp_ind/total_emp_ind
gen dl_nonman = cz_share * dl * ind_share

*** Sum up across all industries for each CZ ***
collapse (sum) dl_nonman, by(czone)
gen yr = 2000
save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 2000-2007.dta", replace

append using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1990-2000.dta"
sort czone yr
save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1990-2007.dta", replace

use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1990-2007.dta" 
merge m:m czone yr using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset.dta"
drop if statefip == .

save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset.dta", replace
