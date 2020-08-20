******************************************************************************************************************************************************************************************
* Tables 8-10
******************************************************************************************************************************************************************************************

use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset.dta", replace

use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset 1998-2007.dta", replace

use "C:\Users\naafe\China\Codes\China Shock and Oil Shocks\dta\final dataset 1993-1998.dta", replace

******************************************************************************************************************************************************************************************
* Manufacturing Employment
******************************************************************************************************************************************************************************************

ivregress 2sls d_sh_empl_mfg (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man l_shind_manuf_cbp reg* l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_mfg dl_man l_shind_manuf_cbp reg* l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_empl_mfg_edu_c (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_mfg_edu_c dl_man l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_empl_mfg_edu_nc (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_mfg_edu_nc dl_man l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

******************************************************************************************************************************************************************************************
* Non-Manufacturing Employment
******************************************************************************************************************************************************************************************

ivregress 2sls d_sh_empl_nmfg (d_tradeusch_pw=d_tradeotch_pw_lag) dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_nmfg dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_empl_nmfg_edu_c (d_tradeusch_pw=d_tradeotch_pw_lag) dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_nmfg_edu_c dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_empl_nmfg_edu_nc (d_tradeusch_pw=d_tradeotch_pw_lag) dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_empl_nmfg_edu_nc dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

******************************************************************************************************************************************************************************************
* Unemployment Rate
******************************************************************************************************************************************************************************************

ivregress 2sls d_sh_unempl (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_unempl dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_unempl_edu_c (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_unempl_edu_c dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_unempl_edu_nc (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_unempl_edu_nc dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

******************************************************************************************************************************************************************************************
* NILF
******************************************************************************************************************************************************************************************

ivregress 2sls d_sh_nilf (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_nilf dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_nilf_edu_c (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_nilf_edu_c dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)

ivregress 2sls d_sh_nilf_edu_nc (d_tradeusch_pw=d_tradeotch_pw_lag) dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
ivregress 2sls d_sh_nilf_edu_nc dl_man dl_nonman l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource reg* t2 [aw=timepwt48], cluster(statefip)
