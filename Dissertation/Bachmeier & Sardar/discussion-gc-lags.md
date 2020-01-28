**NS started discussion 2018/26/04 17:06:12**

I have considered predictions of other variables. I have looked at data from
consumption, investment, unemployment rate and personal savings. You can find
these variables in the two new datasets I have uploaded. 

Unfortunately the disaggregated data on consumption is not as detailed as 
industrial production. For IP, we had 17 industries whereas for consumption I was
only able to find disaggregated data for 6 sectors. 

Monthly Data: Variables for which the data was available at the monthly level. 
Includes aggregate pce, durable goods (further divided into motors, furnishings, 
recreational), nondurable goods (further divided into food, clothing, gasoline), 
personal savings and unemployment rate. Data for consumption was taken from the 
BEA website. I had to construct the data myself by deflating (used price indexes
by major type of product) nominal pce for each sector since it wasn't available 
at such a high frequency. 

I got a mixed responses. Evidence of Granger Causality for aggregate pce, motor 
vehicles & parts, aggreagate non-durables, food, clothing and personal savings.

On the other hand durables, furnishings, gasoline and unemployment rate do 
not respond asymmetrically to oil shocks across business cycles.

Quarterly Data: Investment data is only available at a quarterly frequency.
These included pvt. residential fixed investment, private domestic investment
& private nonresidential fixed investment. There was no evidence in Granger
Causality for any of these varaibles. We can argue that quarterly frequency 
does not capture the response to oil shocks given the high volaitlity of 
oil prices.


**LB started discussion 2018/04/23 09:07:46**

I don't know that the argument that the SIC is "more stable" than the AIC
will be considered acceptable to a referee. I think a better strategy
would be to do the following:

- Use the SIC for everything
- Fix the lag length at a longer value like 12 and redo the test

If we are using the Newey-West corrected covariance matrix, we should be
okay, because omitting lags introduces serial correlation. From a prediction
perspective, we can argue that the threshold model does better. Publishing
in a top journal requires that we deal carefully with issues of specification
at each step of the analysis. The goal is to make it "feel" reliable, so
minor changes in the analysis make a big difference, and in the case of
a Granger causality test, that means taking the choice of lag length
seriously.

There are also multi-step causality tests. It might take several periods
before the nonlinearity has an effect on forecasts. We can consider those
tests after we have something that will look good to the referees using
the current tests that we're doing.

The most important next step in the analysis will be to consider predictions
of other variables (consumption data, for example). After that, we need
to worry about the exact form of the oil shock and the controls we put into
the regression.

Essentially, we need to focus on the Granger causality testing to the point
that we could publish a paper that looks at just that issue. Only then should
we move on to the IRF analysis. I think there are too many loose ends right
now that have to be tied up before we can move on. Looking at what has been
done, it feels like this is something we just threw together in an attempt
to add length to the paper, rather than the result of careful analysis being
done to satisfy our curiosity about predictions using the nonlinear model.
That may or may not be the case, but if that is how our paper reads, we will
not be able to publish in a good journal.
