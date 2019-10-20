import numpy as np
import matplotlib.pyplot as plt

#D Exponential function
def f(x):
    f = 0.5*(x + abs(x))
    return f

domain = np.linspace(-1, 1, 100)


### Chebshev nodes with 20
n = 20
i = np.arange(n, dtype=np.float64)
xn= []

for i in range(n):
   x = np.cos((2 *(i+ 1) - 1) / (2 * n) * np.pi)
   xn.append(x)
  
x= np.asarray(xn)

# Cubic polynomial
p_cubic=np.polyfit(x,f(x),3)
y_cubic=np.polyval(p_cubic,domain) 

# Monomial of order 5
p_order5=np.polyfit(x,f(x),5)
y_order5=np.polyval(p_order5,domain)

#3 Monomial of order 10
p_order10=np.polyfit(x,f(x),10)
y_order10=np.polyval(p_order10,domain)

#Calculation of approximation errors 
error_cubic=abs(f(domain)-y_cubic)
error_order5=abs(f(domain)-y_order5)
error_order10=abs(f(domain)-y_order10)

#Plot 
plt.figure(1)

plt.subplot(121)
plt.plot(domain, f(domain), color='pink', label='ramp function')
plt.plot(domain, y_cubic,'r', label='cubic')
plt.plot(domain, y_order5, color='orange' , label='monomial 5')
plt.plot(domain, y_order10, color='blue', label='monomial 10')
plt.legend(['ramp', 'cubic polynomial', 'monomial 5', 'monomial 10'], loc='best')
plt.title('Ramp function Chebychev nodes')

plt.subplot(122)
plt.plot(domain, error_cubic,'r', label='cubic')
plt.plot(domain, error_order5, color='orange', label='monomial 5')
plt.plot(domain, error_order10, color='blue', label='monomial 10')
plt.legend(loc='best')
plt.title('Ramp function errors Chebychev nodes')

plt.subplots_adjust(top=1.0, bottom=0.1, left=0, right=1.8, hspace=0.7, wspace=0.7)

# Ramp function
def f(x):
    f = 1/(1+25*(x**2))
    return f

domain = np.linspace(-1, 1, 100)

# Chebychev interpolation nodes
n = 20
i = np.arange(n, dtype=np.float64)
xn= []

for i in range(n):
   x = np.cos((2 *(i+ 1) - 1) / (2 * n) * np.pi) 
   xn.append(x)
  
x= np.asarray(xn)

# Cubic polynomial
p_cubic=np.polyfit(x,f(x),3) 
y_cubic=np.polyval(p_cubic,domain) 

# Monomial of order 5
p_order5=np.polyfit(x,f(x),5)
y_order5=np.polyval(p_order5,domain)

# Monomial of order 10
p_order10=np.polyfit(x,f(x),10)
y_order10=np.polyval(p_order10,domain)

# Calculation of approximation errors
error_cubic=abs(f(domain)-y_cubic)
error_order5=abs(f(domain)-y_order5)
error_order10=abs(f(domain)-y_order10)

# Plot 
plt.figure(2)

plt.subplot(221)
plt.plot(domain, f(domain), color='pink', label='runge function')
plt.plot(domain, y_cubic,'r', label='cubic')
plt.plot(domain, y_order5, color='orange' , label='monomial 5')
plt.plot(domain, y_order10, color='blue', label='monomial 10')
plt.legend(['runge', 'cubic polynomial', 'monomial 5', 'monomial 10'], loc='best')
plt.title('Runge function Chebychev nodes')

plt.subplot(222)
plt.plot(domain, error_cubic,'r', label='cubic')
plt.plot(domain, error_order5, color='orange', label='monomial 5')
plt.plot(domain, error_order10, color='blue', label='monomial 10')
plt.legend(loc='best')
plt.title('Runge function errors Chebychev nodes')

plt.subplots_adjust(top=2.5, bottom=0.1, left=0, right=1.8, hspace=0.7, wspace=0.7)

#Definition of the function we want to approximate 
def f(x):
    f = np.exp(1/x)
    return f
domain = np.linspace(-1, 1, 100)

# Chebychev interpolation nodes
n = 20
i = np.arange(n, dtype=np.float64)
xn= []

for i in range(n):
   x = np.cos((2 *(i+ 1) - 1) / (2 * n) * np.pi)
   xn.append(x)
  
x= np.asarray(xn)

# Cubic polynomial
p_cubic=np.polyfit(x,f(x),3) 
y_cubic=np.polyval(p_cubic,domain) 

#2 Monomial of order 5
p_order5=np.polyfit(x,f(x),5)
y_order5=np.polyval(p_order5,domain)

#3 CASE: Monomial of order 10
p_order10=np.polyfit(x,f(x),10)
y_order10=np.polyval(p_order10,domain)

#Calculation of approximation errors 
error_cubic=abs(f(domain)-y_cubic)
error_order5=abs(f(domain)-y_order5)
error_order10=abs(f(domain)-y_order10)

#Plot 
plt.figure(3)

plt.subplot(321)
plt.plot(domain, f(domain), color='pink', label='exponential function')
plt.plot(domain, y_cubic,'r', label='cubic')
plt.plot(domain, y_order5, color='orange' , label='monomial 5')
plt.plot(domain, y_order10, color='blue', label='monomial 10')
plt.ylim([-20000,200000])
plt.legend(['exponential f', 'cubic polynomial', 'monomial 5', 'monomial 10'], loc='best')
plt.title('Exponential function approximation Chebychev nodes')

plt.subplot(322)
plt.plot(domain, error_cubic,'r', label='cubic')
plt.plot(domain, error_order5, color='orange', label='monomial 5')
plt.plot(domain, error_order10, color='blue', label='monomial 10')
plt.legend(loc='best')
plt.ylim([0,200000])
plt.title('Exponential function errors Chebychev nodes')

plt.subplots_adjust(top=3.7, bottom=0.1, left=0, right=1.8, hspace=0.7, wspace=0.7)


