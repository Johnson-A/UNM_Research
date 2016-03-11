from time import clock
from dolfin import *
from math import *
from shutil import copyfile
from os import makedirs
import numpy

"""This version of the code runs a swarm of simulations of various viscosities and temperatures per viscosity. Since 5-29, it also includes an adiabatic temperature variance at the LAB."""

def RunJob(Tb,mu_value,path):
	set_log_level(ERROR)

	runtimeInit = clock() #to time the calculation
	tfile = File(path+'/t6t.pvd')
	mufile = File(path+"/mu.pvd")
	ufile = File(path+'/velocity.pvd')
	gradpfile = File(path+'/gradp.pvd')
	pfile = File(path+'/pstar.pvd')
	parameters = open(path+'/parameters','w',0)

	temp_values = [27.+273, Tb+273, 1300.+273, 1305.+273]

	rho_0 = 3300.
	alpha = 2.5e-5
	g = 9.81
	kappa = 1.E-6 #not actually used.. (in cian's notes he is using the non dimensional kappa in the equations, which here i have defined as 1 (as kappa/kappa_0) so i think i can get away with this)

	dTemp = temp_values[3] - temp_values[0]
	b    = 12.7
	cc   = log(128)
	#Ep  = 0.014057
	Ep   = b/dTemp 
	theta = 0.5
	h = 1000000.
	mu_a = mu_value #this was taken from the blankenbach paper, can change..
	kappa_0 = 1.E-6
	
	#above - 'freely chosen' variables 
	temp_values = [x/dTemp for x in temp_values] #non dimensionalising temp

	mu_bot = exp(-Ep*(temp_values[3]*dTemp - 1573) + cc )*mu_a
	print mu_a
	print mu_bot
	Ra = rho_0*alpha*g*dTemp*h**3/(kappa_0*mu_a)
	w0 = rho_0*alpha*g*dTemp*h**2/mu_a 
	tau = h/w0
	p0 = mu_a*w0/h
	print Ra
	print w0
	print p0
	vslipx = 1.6e-09/w0
	dt = 3.E11/tau
	t = dt
	tEnd = 3.E15/tau #non dimensionalising times.  Note these are not in million years,
	nx = 50
	ny = 50
	MeshHeight = 400000./h #nondimensional (nd)
	MeshWidth = 1000000./h #nd
	LABHeight = 0.75*MeshHeight
	for name in dir():
		ev = str(eval(name))
		if name[0] != '_' and ev[0] != '<':
			parameters.write(name+' = '+ev+'\n')
	vslip = Constant((1.6e-09/w0, 0.0)) #nondimensional
	noslip = Constant((0.0, 0.0))

	class PeriodicBoundary(SubDomain):
	    def inside(self, x, on_boundary):
		return abs(x[0]) < DOLFIN_EPS and on_boundary

	    def map(self, x, y):
		y[0] = x[0] - MeshWidth
		y[1] = x[1]

	pbc = PeriodicBoundary()

	class LithosExp(Expression):
	    def eval(self, values, x):
	#        values[0] = LABHeight - 0.1*MeshHeight + cos(10*pi*x[0]/MeshWidth)*0.1*MeshHeight
        #		disc = (0.2 * MeshHeight) ** 2 - (x[0] - 0.5 * MeshWidth) ** 2
                height = 0.05
		width = 0.2
		scale = 0.05
		ridge1 = height*(1-tanh((x[0] - (0.5 + width)* MeshWidth)/scale))
		ridge2 = height*(1-tanh((x[0] - (0.5 - width)* MeshWidth)/scale))

	#	hump = 0
	#	if disc > 0:
	#	    hump = sqrt(disc)
		
		hump = ridge1 - ridge2
		values[0] = LABHeight - hump

	LAB = LithosExp()

	class TempExp(Expression):
		def eval(self, value, x):
				if x[1] >= LAB(x):
				    value[0] = temp_values[0] + (temp_values[1]-temp_values[0])*(MeshHeight - x[1])/(MeshHeight - LAB(x))
				else:
				    	value[0] = temp_values[3] - (temp_values[3]-temp_values[2])*(x[1])/(LAB(x))

	class muExp(Expression):
	    def __init__(self, muTemp):
		self.muTemp = muTemp
	    def eval(self, value, x):
		t1 = numpy.zeros(1, dtype='d')
		self.muTemp.eval(t1,x)
		value[0] = exp(-Ep*(t1*dTemp - 1573) + cc*x[1]/MeshHeight)*mu_a
		#if value[0] > 1e+26:
		#    value[0] = 1e+26 ###########################################
		value[0] = value[0] / mu_a

	    
	    


	def ZerothStokes(zerothMu, zerothTemp):

	    solver = KrylovSolver('tfqmr', 'amg')
	    R = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
	    Llll = VectorFunctionSpace(mesh, 'CG', 2)
	    Q = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
	    W = R * Q

	    def top(x, on_boundary):return abs(x[1] - MeshHeight) <= DOLFIN_EPS
	    def bottom(x, on_boundary):return abs(x[1]) <= DOLFIN_EPS
	    def left(x, on_boundary): return ((x[1] <= LAB(x)) and (x[0] <= DOLFIN_EPS))
	    def right(x, on_boundary): return (((x[1])<= LAB(x)) and ((x[0] - MeshWidth) <= DOLFIN_EPS))
	    noslip = Constant((0.0, 0.0))
	    bc0 = DirichletBC(W.sub(0), noslip, top)
	    bc1 = DirichletBC(W.sub(0), vslip, bottom)
	    bcp1 = DirichletBC(W.sub(1), Constant(0.0), bottom)
	    bcs = [bc0, bc1, bcp1]

	    u = Function(W)
	    v,p = split(u)
	    u_t = TestFunction(W) #not actually used here
	    v_t, p_t = TestFunctions(W)
	    u_a = TrialFunction(W)
		
	    ################# Form!!! #########

	    r_v = (inner(sym(grad(v_t)), 2.*zerothMu*sym(grad(v))) \
		- div(v_t)*p  \
		- zerothTemp*v_t[1])*dx

	    r_p = p_t*div(v)*dx

	    ###############
	    #not both this and below forms use cian's scaling.
		
	    r = r_v + r_p
	 
	    solve(r == 0, u, bcs)
	    vel, pre = u.split()
	    return vel




	mesh = RectangleMesh(Point(0.0,0.0),Point(MeshWidth,MeshHeight),nx,ny,'left/right')


	Svel = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
	Spre = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
	Stemp = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
	Smu = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
	Sgradp = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)

	S0 = MixedFunctionSpace([Svel,Spre,Stemp])

	u = Function(S0)
	v,p,T = split(u)

	v_t, p_t, T_t = TestFunctions(S0)
	du = TrialFunction(S0)

	T0 = Function(Stemp)
	T0.interpolate(TempExp())


	mu = Function(Smu)
	mu.interpolate(muExp(T0))

	v0 = Function(Svel)
	v0 = ZerothStokes(mu, T0)

	v_theta = (1. - theta)*v0 + theta*v

	T_theta = (1. - theta)*T + theta*T0

	r_v = (inner(sym(grad(v_t)), 2.*mu*sym(grad(v))) \
	    - div(v_t)*p \
	    - T*v_t[1] )*dx

	r_p = p_t*div(v)*dx

	r_T = (T_t*((T - T0) \
	    + dt*inner(v_theta, grad(T_theta))) \
	    + (dt/Ra)*inner(grad(T_t), grad(T_theta)) )*dx

	r = r_v + r_p + r_T

	def top(x, on_boundary): return x[1] >= MeshHeight - DOLFIN_EPS 
	def bottom(x, on_boundary): return x[1] <= DOLFIN_EPS
	def left(x, on_boundary): return ((x[1] <= LAB(x)) and (x[0] <= DOLFIN_EPS))
	def right(x, on_boundary): return (((x[1])<= LAB(x)) and ((x[0] - MeshWidth) <= DOLFIN_EPS))

	bcv0 = DirichletBC(S0.sub(0), noslip, top)
	bcv1 = DirichletBC(S0.sub(0), vslip, bottom)
	bcp0 = DirichletBC(S0.sub(1), Constant(0.0), bottom)
	bct0 = DirichletBC(S0.sub(2), Constant(temp_values[0]), top)
	bct1 = DirichletBC(S0.sub(2), Constant(temp_values[3]), bottom)

	bcs = [bcv0, bcv1, bcp0, bct0, bct1]

	t = 0
        count = 0
	while (t < tEnd):
		solve(r == 0, u, bcs)
		t += dt
		nV, nP, nT = u.split()
		if (count % 100 == 0):
			pfile << nP
			ufile << nV
			tfile << nT
			mufile << mu
			gradpfile << project(grad(nP),Sgradp)
		count +=1
		assign(T0, nT)
		assign(v0, nV)
		mu.interpolate(muExp(T0))

	print 'Case mu=%g, Tb=%g complete.'%(mu_a, Tb), ' Run time =', clock() - runtimeInit, 's'

if __name__ == "__main__":
	Tbs = [800,1000,1300]
	#Mus = numpy.linspace(1e19,1e21,5) 
	Mus = [1e19]
	for mu in Mus:
		mufolder='mu='+str(mu)
		try:
			makedirs(mufolder)
		except:
			pass
		for temp in Tbs:
			tempfolder='Tb='+str(temp)
			workpath=mufolder+'/'+tempfolder	
			RunJob(temp,mu,workpath)

