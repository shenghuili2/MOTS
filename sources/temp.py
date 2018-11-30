import math
x = 0
a = 1389416494410643649
for i in range(2,math.ceil(1389416494410643649**0.5)):
    if a % i == 0:
        x += i
        x += (a // i)
        print ("i: ",i)
        print ("a/i: ", a//i)
    if i % 10000000 == 0:
        print ("check: ", i)
print ("x: ", x)
x:  18998163238110