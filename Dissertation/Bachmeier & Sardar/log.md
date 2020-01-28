[Project Roadmap](https://lancebachmeier.com/active/projects/oilthreshold/doc/trunk/project-roadmap.md)
based on our discussion May 11, 2018.

[Discussion of GC lag length selection](https://lancebachmeier.com/active/projects/oilthreshold/doc/trunk/discussion-gc-lags.md)

Please add comments to that file so that everything is in one place when
we return to it in the future.

*2018/10/04 21:12:22*

NS: I have added impulse responses for aggregate industrial production
(manufacturing and overall).

Looked at the test using nominal price of oil as well. There is no evidence
of Granger Causality in that case.

The next thing I am going to do is use Kilian's oil shocks for our analysis.
I recall you mentioning that we need to control for the precautionary demand
shocks while using NOPI as the price shock. Why are we using demand shock as 
a control variable instead of using it as an alternative to NOPI?

*2018/04/04 20:35:12*

NS: I have added aggregate industrial production (manufacturing and overall)
to the analysis. 

_AIC v. SIC_
I did the comparison of lag lengths obtained from AIC and SIC. Using SIC
fewer lags were chosen, as opposed to AIC in which almost 20 lags were
chosen. I played around with the code and noticed that SIC was more stable
than AIC. For example, in our function we are testing for a max of 24 lags.
If we test for a larger no. (let's say 48 instead of 24), the lags chosen using 
AIC keep on increasing. Whereas for SIC, the no. of lags chosen stay the
same even if we test for more lags.That explains the motivation to use SIC
instead of AIC.

_Heteroskedasticity_
Conducted White test and concluded there is heteroskedasticity. Hence 
Newey West standard errors were used while conducting Wald Test. I
concluded Granger causality for overall IP, manfacturing IP. At the 
sectoral level, apparel, electrical, electronics, furniture, machinery,
metals, motor, paper also exhibit non-linearity.



*2018/04/02 11:29:24*

LB: I played with the code some more to figure out what is going on. Inside
the `gctest` function, you can choose whether to use the AIC or SIC for
model selection. You can also choose whether or not to use NW corrected
standard errors.

The results depend heavily on the lag selection criteria and whether or
not you use NW corrected standard errors. That's going to be a problem
when we try to publish. We have to provide some justification for that.

You will need to dig into this. The first thing to do is to see if there
is heteroskedasticity using the White test. You should be able to modify
the `gctest` function to do that. Please do not use copy and paste for all
sectors. That is not a sustainable approach. Get a comparison of the lag
lengths chosen using the AIC and SIC. If we don't figure this out, we won't
be able to publish the paper, regardless of what else we do elsewhere in
the paper.

Also, add the overall industrial production index to the analysis. We
can't report just sector-level results without the aggregate results.

*2018/03/27 16:58:25*

I made changes to the file `threshold effects of oil shocks.Rmd`. I
removed code duplication and generally made it easier for me to work with.

*Dec 5 10:33 AM*

You can add to the top of this file in order to save a message to the project.
