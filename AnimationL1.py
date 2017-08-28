import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi
from matplotlib import animation

G=6.67e-11
M_sun=2e30               
M_earth=6e24
M_jupiter=1.898e27
R_peri_jupiter=741e9        
R_peri_earth=147e9
R_api_jupiter=817e9
R_api_earth=152e9
R_sun=6.9e8
year_earth=3.14e7
year_jupiter=11.86*year_earth

RJ=R_peri_jupiter
M_planet=M_jupiter
R_peri=R_peri_jupiter
R_api=R_api_jupiter
year=year_jupiter

M_total=M_planet+M_sun
alpha=M_planet/M_total


initial_x_planet= R_peri
initial_y_planet=0
initial_vx_planet=0
initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]

T=year
div=10000
t=sp.linspace(0.,30*T,div)
t2=t

initial_x_asteroidref=R_peri*(1-((alpha/3)**(1./3.)))
initial_y_asteroidref=0
initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5                 
initial_vx_asteroidref=0
initial_vy_asteroidref=initial_speed_asteroidref
initial_conditions_asteroidref=[initial_x_asteroidref,initial_vx_asteroidref,initial_y_asteroidref,initial_vy_asteroidref]
    
def planet(planetpos,t):
    x_planet=planetpos[0]
    y_planet=planetpos[2]
    vx_planet=planetpos[1]
    vy_planet=planetpos[3]
    ax_planet=-((G*M_sun)/(((x_planet**2)+(y_planet**2))**1.5))*x_planet
    ay_planet=-((G*M_sun)/(((x_planet**2)+(y_planet**2))**1.5))*y_planet
    
    return [vx_planet, ax_planet,vy_planet, ay_planet]
    
solnplanet=spi.odeint(planet,initial_conditions_planet,t)
x_planet_old=solnplanet[:,0]
y_planet_old=solnplanet[:,2]

solnasteroidref=spi.odeint(planet,initial_conditions_asteroidref,t)
x_asteroidref=solnasteroidref[:,0]
y_asteroidref=solnasteroidref[:,2]

def asteroid(asteroidpos,t2,x_planet_old,y_planet_old):
        
    x_planet_new=np.interp(t2,t,x_planet_old)
    y_planet_new=np.interp(t2,t,y_planet_old)
    r_plan=np.sqrt((x_planet_new**2)+(y_planet_new**2))
        
    x_asteroid=asteroidpos[0]
    y_asteroid=asteroidpos[2]
    r_ast=np.sqrt((x_asteroid**2)+(y_asteroid**2))
        
    if r_ast==r_plan or r_ast==0:
        vx_asteroid=0
        vy_asteroid=0
        ax_asteroid=0
        ay_asteroid=0
    else:
        vx_asteroid=asteroidpos[1]
        vy_asteroid=asteroidpos[3]
        
        ax_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*x_asteroid -(((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(x_asteroid-x_planet_new))
        ay_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*y_asteroid- (((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(y_asteroid-y_planet_new))
        
            
    return [vx_asteroid,ax_asteroid,vy_asteroid,ay_asteroid]
        
i=1.01*initial_x_asteroidref
    
initial_x_asteroid=i
initial_y_asteroid=0
initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
x_asteroid=solnasteroid[:,0]
y_asteroid=solnasteroid[:,2]
    
x_asteroidframe=x_asteroid-x_asteroidref
y_asteroidframe=y_asteroid-y_asteroidref
        
r_ast=np.sqrt((x_asteroidframe**2)+(y_asteroidframe**2))
        
fig = pl.figure()

ax = pl.axes(xlim=(-5e12, 5e12), ylim=(-5e12, 5e12))
ax.patch.set_facecolor('black')

line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    
    return line,

def animate(f):
    f=5*f
    x1 = x_asteroid[:f]
    y1 = y_asteroid[:f]
    x2 = x_planet_old[:f]
    y2 = y_planet_old[:f]
    line.set_data(x1, y1)
    
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=10000000, interval=0.001, blit=True)
circ=pl.Circle((0,0), radius=25*R_sun, color='orange')
ax.add_patch(circ)
pl.plot(x_planet_old, y_planet_old)
pl.show()
        