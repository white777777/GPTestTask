# GPTestTask
regression analysis app
In input we have csv file with:
1) calendar month
2) borehole name
3) working hours of borehole at this month
4) quantity of oil from borehole
5) quantity of water from this borehole

This program trys to approximate quantity of oil or water with some special function q_i(t) = q_{0i} f(t) wher f(t) with parameters and find this parameters:
1) f(t) = exp(-Dt), D>0
2) f(t) = (1+t)^(-a),  a>0
3) f(t) = (1+bDt)^(-1/b), b>0, D>0
4) f(t) = (1+t)^(-a) if (t<tau)
f(t) = (1+bD(t-tau))^(-1/b) if t>tau
a>0, b>0, tau>0

regression model is ln q_i (t_j) = ln q_{0i} +ln f(t_j) +\Ksi_{ij}
