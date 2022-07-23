#!/bin/env python3
from datetime import datetime
import math
import platform

print("Writing the date and square root of 42 to file out.txt...")
fp = open('out.txt', 'w')
today = datetime.now()

fp.write("This is a test script for python on Katana\n")
fp.write(str(math.sqrt(42))+ '\n')
fp.write(str(today)+'\n')

fp.close()

print("\n")
print("This is Python version {}".format(platform.python_version()) )
