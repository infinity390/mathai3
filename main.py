import mathai
import parser
from base import *

# Expansion example
expr_expansion = parser.take_input("(x+y+z)^3")
expanded_expr = mathai.solve(mathai.expand2(expr_expansion))
print(expanded_expr)

# Integration example
expr_integration = parser.take_input("x*(1+sin(x)/x)")
integrated_expr = mathai.integratex(expr_integration, 4, "v_0", [], True, False, True)  # integrate wrt x
integrated_expr = mathai.solve(mathai.expand2(integrated_expr))
print(integrated_expr)

# Differentiate what we integrated
differentiated_expr = mathai.diffx(integrated_expr)

# Verify that differentiation of the integral returns the original expression
is_correct = mathai.solve(mathai.expand2(differentiated_expr - expr_integration)) == 0
print(is_correct)  # True if d/dx ∫F(x)dx = F(x)
