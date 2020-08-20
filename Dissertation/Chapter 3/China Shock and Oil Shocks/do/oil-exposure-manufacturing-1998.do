******************************************************************************************************************************************************************************************
* Manufacturing Oil Exposure
******************************************************************************************************************************************************************************************

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

*** Calculate l_{j,i,t}/l_{i,t} ***
egen emp_czone = sum(emp2000), by(czone)
gen share = emp2000/emp_czone

*** Matching dataset with change in labor due to oil shocks ***
merge m:m naics using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor Due to Oil 1998-2007_Man.dta"
drop if yr == .
drop if czone == .
drop _merge
sort czone naics

*** Create oil exposure variable ***
gen dl_man = share * dl

*** Sum up across all industries for each CZ ***
collapse (sum) dl_man, by(czone)
gen yr = 2000
save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1998-2007_Man.dta", replace

use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\Change in Labor due to Oil Shocks Across CZ 1998-2007_Man.dta" 
merge m:m czone yr using "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\workfile_china.dta"
drop if statefip == .
drop _merge
drop if dl_man == .
save "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset 1998-2007.dta", replace
