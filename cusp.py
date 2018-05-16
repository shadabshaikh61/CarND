# Hello World program in Python
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L

x=[]
y=[]
    
for i in frange(-1,1,0.0001):
    for j in frange(-1,1,0.0001):
        if(i*i*i+j*j*j ==i*j):
            x.append(i)
            y.append(j)
print(x)
print(y)
