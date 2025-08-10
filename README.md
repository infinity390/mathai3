## Comparison: mathai vs. sympy expansion

This project demonstrates symbolic expansion using both the custom `mathai` library and the well-known `sympy` library.

### Using mathai
```python
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
```

### Using sympy
```python
from sympy import symbols, sin, expand, integrate, diff, simplify

# Define variable
x, y, z = symbols('x y z')

# Expansion example
expr_expansion = (x + y + z) ** 3
expanded_expr = expand(expr_expansion)
print(expanded_expr)

# Integration example
expr_integration = x * (1 + sin(x) / x)
integrated_expr = integrate(expr_integration, x)
print(simplify(integrated_expr))

# Differentiate what we integrated
differentiated_expr = diff(integrated_expr, x)

# Verify that differentiation of the integral returns the original expression
is_correct = simplify(differentiated_expr - expr_integration) == 0
print(is_correct)  # True if d/dx ∫F(x)dx = F(x)
```

Sympy output: ``` `x**3 + 3*x**2*y + 3*x**2*z + 3*x*y**2 + 6*x*y*z + 3*x*z**2 + y**3 + 3*y**2*z + 3*y*z**2 + z**3`

x**2/2 - cos(x)

True``` <br>
My Math Ai output: ```((3*(x^2)*y)+(3*(x^2)*z)+(3*(y^2)*x)+(3*(y^2)*z)+(3*(z^2)*x)+(3*(z^2)*y)+(6*x*y*z)+(x^3)+(y^3)+(z^3))

((-1*cos(x))+((2^-1)*(x^2)))

True```

Run the main.py file for getting the Math Ai's output
