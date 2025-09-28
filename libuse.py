import sys
from mathai import *
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Increase recursion limit for deep symbolic tasks ---
sys.setrecursionlimit(10000)

# --- Mathai functions ---
def integration_byparts(item):
    return simplify(fraction(simplify(byparts(simplify(parse(item)))[0])))

def integration_apart(item):
    return simplify(fraction(integrate(apart(factor2(simplify(parse(item)))))[0]))

def integration_direct(item):
    return simplify(fraction(simplify(integrate(simplify(parse(item)))[0])))

def integration_trig(item):
    return simplify(trig0(integrate(trig1(simplify(parse(item))))[0]))

def algebra(item):
    return logic0(simplify(expand(simplify(parse(item)))))

def trig_basic(item):
    return logic0(simplify(expand(trig3(simplify(parse(item))))))

def trig_advanced(item):
    return logic0(simplify(trig0(trig1(trig4(simplify(fraction(trig0(simplify(parse(item))))))))))

# --- All tasks flattened into a single list ---
all_tasks = [
    # Advanced trigonometry first
    *[(item, trig_advanced) for item in [
        "cos(x)/(1+sin(x)) + (1+sin(x))/cos(x) = 2*sec(x)",
        "(1+sec(x))/sec(x) = sin(x)^2/(1-cos(x))"
    ]],
    
    # integration by parts
    *[(item, integration_byparts) for item in ["sin(x)*x", "x*sin(3*x)", "x*log(abs(x))", "arctan(x)"]],
    
    # partial fractions
    *[(item, integration_apart) for item in ["x/((x+1)*(x+2))", "1/(x^2-9)"]],
    
    # direct integration
    *[(item, integration_direct) for item in [
        "x*sqrt(x+2)", "sin(cos(x))*sin(x)", "2*x/(1+x^2)",
        "sqrt(a*x+b)", "cos(sqrt(x))/sqrt(x)",
        "e^(arctan(x))/(1+x^2)", "sqrt(sin(2*x))*cos(2*x)"
    ]],
    
    # trig integration
    *[(item, integration_trig) for item in ["sin(2*x+5)^2", "sin(x)^4", "cos(2*x)^4"]],
    
    # algebra
    *[(item, algebra) for item in ["(x+1)^2 = x^2+2*x+1", "(x+1)*(x-1) = x^2-1"]],
    
    # trig basic
    *[(item, trig_basic) for item in ["2*sin(x)*cos(x)=sin(2*x)"]],
]

# --- Worker wrapper ---
def run_task(task):
    item, func = task
    try:
        result = func(item)
    except Exception as e:
        result = str(e)
    return item, result

# --- Parallel execution ---
if __name__ == "__main__":
    print("Running tasks in parallel with recursion limit =", sys.getrecursionlimit(), "\n")
    results = []
    with ProcessPoolExecutor() as executor:
        future_to_task = {executor.submit(run_task, task): task for task in all_tasks}
        for future in as_completed(future_to_task):
            item, result = future.result()
            results.append((item, result))
            print(f"{item}  =>  {result}")
