import sys
import os
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

# --- All tasks ---
all_tasks = [
    *[(item, trig_advanced) for item in [
        "cos(x)/(1+sin(x)) + (1+sin(x))/cos(x) = 2*sec(x)",
        "(1+sec(x))/sec(x) = sin(x)^2/(1-cos(x))"
    ]],
    *[(item, integration_byparts) for item in ["sin(x)*x", "x*sin(3*x)", "x*log(abs(x))", "arctan(x)"]],
    *[(item, integration_apart) for item in ["x/((x+1)*(x+2))", "1/(x^2-9)"]],
    *[(item, integration_direct) for item in [
        "x*sqrt(x+2)", "sin(cos(x))*sin(x)", "2*x/(1+x^2)",
        "sqrt(a*x+b)", "cos(sqrt(x))/sqrt(x)",
        "e^(arctan(x))/(1+x^2)", "sqrt(sin(2*x))*cos(2*x)"
    ]],
    *[(item, integration_trig) for item in ["sin(2*x+5)^2", "sin(x)^4", "cos(2*x)^4"]],
    *[(item, algebra) for item in ["(x+1)^2 = x^2+2*x+1", "(x+1)*(x-1) = x^2-1"]],
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

# --- Parallel execution using all cores ---
if __name__ == "__main__":
    num_cores = os.cpu_count()  # detect available cores
    print(f"Running tasks in parallel on {num_cores} cores with recursion limit = {sys.getrecursionlimit()}\n")
    
    # Store results in original order
    results = [None] * len(all_tasks)

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        # Map each task to a future
        future_to_index = {executor.submit(run_task, task): idx for idx, task in enumerate(all_tasks)}
        for future in as_completed(future_to_index):
            idx = future_to_index[future]
            item, result = future.result()
            results[idx] = (item, result)

    # Print results in original task order
    for item, result in results:
        print(f"{item}  =>  {result}")
