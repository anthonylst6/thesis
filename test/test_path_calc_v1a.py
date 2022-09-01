from pathlib import Path
ver = Path(__file__).stem[-3:]

def func_from_calc():
    print("version from calc: " + ver)
    
func_from_calc()