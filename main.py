import mathai
import parser
from base import *

equation = parser.take_input("(x+y+z)^3")
equation = mathai.solve(mathai.expand2(equation))

print(equation)
