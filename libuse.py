import sys
import os
import time
from mathai import *
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Increase recursion limit ---
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

# --- Batch execution ---
if __name__ == "__main__":
    num_cores = os.cpu_count()
    print(f"Running tasks in batches of {num_cores} (one task per core)...\n")
    
    start_time = time.time()

    # Process tasks in batches
    for i in range(0, len(all_tasks), num_cores):
        batch = all_tasks[i:i+num_cores]  # next batch
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            futures = {executor.submit(run_task, task): task for task in batch}
            for future in as_completed(futures):
                item, result = future.result()
                print(f"{item}  =>  {result}\n")

    total_time = time.time() - start_time
    print(f"All tasks completed in {total_time:.2f} seconds")
