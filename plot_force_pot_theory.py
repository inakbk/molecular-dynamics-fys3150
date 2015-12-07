#plotting the force and the potential 
from pylab import *

rCut = 2.5 #sigma = 1
r_marked = linspace(1,rCut+0.5,1000) 

def F_marked(r_marked):
	r_marked6 = 1/r_marked**(6)
	r_marked8 = r_marked6/r_marked**2
	return r_marked8*(2*r_marked6 - 1)

def U_marked(r_marked):
	r_marked6 = 1/r_marked**(6)
	return r_marked6*(r_marked6 - 1)

figure(1)
subplot(2,1,1)
plot(r_marked,U_marked(r_marked), 'r')
plot(r_marked,U_marked(r_marked) - U_marked(rCut), '--g')
hold('on')
plot(rCut,U_marked(rCut), 'r-o')
plot(r_marked,zeros(len(r_marked)), 'k')
axis([0.99,rCut+0.5,-0.26,0.01])
ylabel('$U/4\epsilon$', fontsize=18)

subplot(2,1,2)
plot(r_marked,F_marked(r_marked), 'b')
plot(rCut,F_marked(rCut), 'b-o')
plot(r_marked,zeros(len(r_marked)), 'k')

axis([0.99,rCut+0.5,-0.1,0.01])
xlabel('$\sigma/r_{ij}$', fontsize=18)
ylabel('$F\sigma^2/24\epsilon$', fontsize=18)

print U_marked(rCut), F_marked(rCut)

show()


