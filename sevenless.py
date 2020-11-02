#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)

for i in range(start, end + 1):
    if i%7 != 0 or i == 0:
        print(i)