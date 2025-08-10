## Comparison: mathai vs. sympy Expansion

This project demonstrates symbolic expansion using both the custom `mathai` library and the well-known `sympy` library.

### Using mathai
```python
import mathai
import parser
from base import *

equation = parser.take_input("(x+y+z)^3")
equation = mathai.solve(mathai.expand2(equation))

print(equation)
```

### Using sympy
```python
from sympy import symbols, expand

x, y, z = symbols('x y z')
equation = expand((x + y + z) ** 3)

print(equation)
```

Sympy expanded: ``x**3 + y**3 + z**3 + 3*x**2*y + 3*x*y**2 + 3*x**2*z + 3*x*z**2 + 3*y**2*z + 3*y*z**2 + 6*x*y*z`` <br>
My Math Ai output: ``((3*(x^2)*y)+(3*(x^2)*z)+(3*(y^2)*x)+(3*(y^2)*z)+(3*(z^2)*x)+(3*(z^2)*y)+(6*x*y*z)+(x^3)+(y^3)+(z^3))``

Run the main.py file for getting the Math Ai's output
