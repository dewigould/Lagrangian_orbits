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
R_sun=6.9e8
year_earth=3.14e7
year_jupiter=11.86*year_earth

year=year_jupiter
R_peri=R_peri_jupiter
Initial=0.005*M_sun
Final=0.06*M_sun
Number=10
M_planet=np.arange(Initial,Final,(Final-Initial)/Number)
Max_displacements=[]


#initial conditions for planet
initial_x_planet= R_peri
initial_y_planet=0
initial_vx_planet=0
initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]

#Time Period for Integration
T=1*year
div=10000
t=sp.linspace(0.,50*T,div)
t2=t

#Initial Condition for Reference Asteroid at L4
initial_x_asteroidref=R_peri*0.5
initial_y_asteroidref=R_peri*(3**0.5)*0.5
initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5                 
initial_vx_asteroidref=-1*initial_speed_asteroidref*(3**0.5)*0.5
initial_vy_asteroidref=initial_speed_asteroidref*0.5
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
 
i=0.05*M_sun
M_total=i+M_sun  
          
def asteroid(asteroidpos,t2,x_planet_old,y_planet_old):
            
        x_planet_new=np.interp(t2,t,x_planet_old)
        y_planet_new=np.interp(t2,t,y_planet_old)
        
        x_asteroid=asteroidpos[0]
        y_asteroid=asteroidpos[2]

        vx_asteroid=asteroidpos[1]
        vy_asteroid=asteroidpos[3]

        ax_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*x_asteroid -(((G*i)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(x_asteroid-x_planet_new))
        ay_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*y_asteroid- (((G*i)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(y_asteroid-y_planet_new))

        
        return [vx_asteroid,ax_asteroid,vy_asteroid,ay_asteroid]

solnasteroidref=spi.odeint(asteroid,initial_conditions_asteroidref,t2,args=(x_planet_old,y_planet_old))
x_asteroidref=solnasteroidref[:,0]
y_asteroidref=solnasteroidref[:,2]


initial_x_asteroid=initial_x_asteroidref*1.01
initial_y_asteroid=initial_y_asteroidref*1.01
initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5

initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))

x_asteroid=solnasteroid[:,0]
y_asteroid=solnasteroid[:,2]

j=0.07*M_sun
M_total=j+M_sun  
          
def asteroid2(asteroidpos,t2,x_planet_old,y_planet_old):
            
        x_planet_new=np.interp(t2,t,x_planet_old)
        y_planet_new=np.interp(t2,t,y_planet_old)
        
        x_asteroid2=asteroidpos[0]
        y_asteroid2=asteroidpos[2]

        vx_asteroid2=asteroidpos[1]
        vy_asteroid2=asteroidpos[3]

        ax_asteroid2=-((G*M_sun)/(((x_asteroid2**2)+(y_asteroid2**2))**1.5))*x_asteroid2 -(((G*j)/((((x_asteroid2-x_planet_new)**2)+((y_asteroid2-y_planet_new)**2))**1.5))*(x_asteroid2-x_planet_new))
        ay_asteroid2=-((G*M_sun)/(((x_asteroid2**2)+(y_asteroid2**2))**1.5))*y_asteroid2- (((G*j)/((((x_asteroid2-x_planet_new)**2)+((y_asteroid2-y_planet_new)**2))**1.5))*(y_asteroid2-y_planet_new))

        
        return [vx_asteroid2,ax_asteroid2,vy_asteroid2,ay_asteroid2]

initial_x_asteroid=initial_x_asteroidref*1.01
initial_y_asteroid=initial_y_asteroidref*1.01
initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5

initial_conditions_asteroid2=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
solnasteroid2=spi.odeint(asteroid2,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))

x_asteroid2=solnasteroid2[:,0]
y_asteroid2=solnasteroid2[:,2]

k=0.06*M_sun
M_total=k+M_sun  
          
def asteroid3(asteroidpos,t2,x_planet_old,y_planet_old):
            
        x_planet_new=np.interp(t2,t,x_planet_old)
        y_planet_new=np.interp(t2,t,y_planet_old)
        
        x_asteroid3=asteroidpos[0]
        y_asteroid3=asteroidpos[2]

        vx_asteroid3=asteroidpos[1]
        vy_asteroid3=asteroidpos[3]

        ax_asteroid3=-((G*M_sun)/(((x_asteroid3**2)+(y_asteroid3**2))**1.5))*x_asteroid3 -(((G*k)/((((x_asteroid3-x_planet_new)**2)+((y_asteroid3-y_planet_new)**2))**1.5))*(x_asteroid3-x_planet_new))
        ay_asteroid3=-((G*M_sun)/(((x_asteroid3**2)+(y_asteroid3**2))**1.5))*y_asteroid3- (((G*k)/((((x_asteroid3-x_planet_new)**2)+((y_asteroid3-y_planet_new)**2))**1.5))*(y_asteroid3-y_planet_new))

        
        return [vx_asteroid3,ax_asteroid3,vy_asteroid3,ay_asteroid3]

initial_x_asteroid=initial_x_asteroidref*1.01
initial_y_asteroid=initial_y_asteroidref*1.01
initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5

initial_conditions_asteroid3=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
solnasteroid3=spi.odeint(asteroid3,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))

x_asteroid3=solnasteroid3[:,0]
y_asteroid3=solnasteroid3[:,2]
 
fig = pl.figure()

ax = pl.axes(xlim=(-4e12, 4e12), ylim=(-4e12, 4e12))
ax.patch.set_facecolor('black')
line, = ax.plot([], [], lw=2)
line2,=ax.plot([],[],lw=2)
line3,=ax.plot([],[],lw=2)
line4,=ax.plot([],[],lw=2)
def init():
    line.set_data([], [])
    line2.set_data([],[])
    line3.set_data([],[])
    line4.set_data([],[])
    return line,line2, line3, line4

def animate(f):
    f=8*f
    x1 = x_asteroid[:f]
    y1 = y_asteroid[:f]
    x2 = x_planet_old[:f]
    y2 = y_planet_old[:f]
    x3 = x_asteroid2[:f]
    y3 = y_asteroid2[:f]
    x4 = x_asteroid3[:f]
    y4 = y_asteroid3[:f]
    line.set_data(x1, y1)
    line2.set_data(x2,y2)
    line3.set_data(x3,y3)
    line4.set_data(x4,y4)
    return line, line2, line3, line4

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1000000000, interval=0.001, blit=True)
circ=pl.Circle((0,0), radius=25*R_sun, color='orange')
ax.add_patch(circ)

pl.show()
        