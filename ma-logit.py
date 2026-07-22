import itertools
import time
from multiprocessing import Pool

import sympy as sy
import numpy as np

from scipy.integrate import quad
from scipy.special import expit

class LogisticSigmoid(sy.Function):
    def fdiff(self, argindex=1):
        return self * (1 - self)
    
    # def _eval_evalf(self, prec):
    #     # 1. Get the argument as a SymPy float
    #     arg = self.args[0]._eval_evalf(prec)

m, x, theta, theta2 = sy.symbols('m x theta theta2')
beta = sy.IndexedBase('beta')
# n = 10

pdf_normal = 1 / sy.sqrt(2 * sy.pi) * sy.exp(-theta**2 / 2)
pr = LogisticSigmoid(beta[2] + sy.sqrt(beta[1]) * theta)

# integrand = sy.binomial(m, x) * (pr**x) * ((1 - pr)**(m - x)) * pdf_normal
integrand =  (pr**x) * ((1 - pr)**(m - x)) * pdf_normal
L_unit = sy.Integral(integrand, (theta,-sy.oo,sy.oo))
l_unit = sy.log(L_unit)

# def k_term(args):
#     order, mu_null, tau2_null, m_obs = args
#     print(order)
#     inner_order, outer_order = order
#     inner_deriv = l_unit.diff(*[beta[i] for i in inner_order])
#     inner_deriv = inner_deriv.replace(sy.Integral, lambda expr, var: sy.Integral(sy.factor(expr), var))
#     binomial_pmf =  sy.binomial(m_obs,x) * (pr**x * (1-pr)**(m_obs-x) * pdf_normal).subs(theta,theta2)
#     expected = sy.Integral(
#         sy.summation(inner_deriv * binomial_pmf, (x,0,m_obs)),
#                         (theta2,-sy.oo,sy.oo)
#                         )
#     if len(outer_order)==0:
#         out = eval_with_numeric_integrals(expected,values={beta[1]:tau2_null,beta[2]:mu_null,m:m_obs})
#     else:
#         outer_deriv = expected.diff(*[beta[i] for i in outer_order])
#         out = eval_with_numeric_integrals(outer_deriv,values={beta[1]:tau2_null,beta[2]:mu_null,m:m_obs})
#     print(out)
#     return(out)


