import math as m
def isfloat(Value):
  try:
    float(Value)
    return True
  except ValueError:
    return False

inpt = [ '0.0' ]
n = 0
while isfloat(inpt[n]):
 n += 1
 inpt.append(input(str(n).zfill(2)+' '))
mv = 0.0
for i in range(1,n):
 mv += float(inpt[i])
mv /= n-1
se = 0.0
for i in range(1,n):
 se += pow((mv - float(inpt[i])),2)
se /= (n-2)
se = m.sqrt(se)
sm = se/m.sqrt(n-1)
print('')
print('mv '+str(mv))
print('se '+str(se))
print('sm '+str(sm))
