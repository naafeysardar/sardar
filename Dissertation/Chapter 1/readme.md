# Changes made by Naafey

I have changed the naming schemes. Data in the full-dataset.csv file contains 
consumption in nominal terms. I have made sure that each transformed variable 
has its own unique name so I do not run into the issue of over writing and 
ending up measuring something else. The new naming scheme is as follows: 

n.durable <= Nominal durables PCE 

durable <= Percentage change in real durables PCE

All other consumption measures have been defined similarly. I have noticed 
a difference in the values computed using '100*pctChange(n.durable/price)'
and 'durable/stats::lag(durable,-1) - 1'. 

When I use the (durable/stats::lag(durable,-1) - 1) command it calculates the 
percentage change in durables PCE between current and previous quarter. I am 
certain of that because I redid the calculations manually to see if I get the
same answer. In that case, I am calculating the following:

$$ Durables_t = 100 \times \frac{Real.Durables_t - Real.Durables_{t-1}}{Real.Durables_{t-1}} $$

I am not sure why the pctChange function is giving me different values. I 
couldn't find any documentation on it either. 

The following link gives a description of the method I used earlier to find
the percentage differece,

https://stackoverflow.com/questions/14614710/calculate-percentage-change-in-an-r-data-frame


# Changes made by Lance Bachmeier

*This file created 2019/11/01*

I tried to use the programs that you had, but it was really inconvenient
and hard to follow due to unclear dependencies and such. Ideally, everything
should go into an R package and be accessed by calling functions. If that's
not possible, you need at a minimum to have one file that does everything
related to that part of the analysis. That's a recipe for disaster as your
program grows and you spend time away from this paper.

I created the file `linear.R` to hold all of the linear VAR model analysis.
That includes sourcing other files first as needed.

Another disaster waiting to happen is renaming variables. It may or may
not do what you think it's doing. Those errors show up at bad times and
in unexpected ways.
