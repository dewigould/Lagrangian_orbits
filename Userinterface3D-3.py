import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi
from mpl_toolkits.mplot3d import Axes3D

G=6.67e-11
M_sun=2e30               
M_earth=6e24
M_jupiter=1.898e27
R_peri_jupiter=741e9        
R_peri_earth=147e9
R_api_jupiter=817e9
R_api_earth=152e9
R_sun=69
year_earth=3.14e7
year_jupiter=11.86*year_earth

M_planet=M_jupiter
R_peri=R_peri_jupiter
R_api=R_api_jupiter
year=year_jupiter

M_total=M_planet+M_sun
alpha=M_planet/M_total

Asteroid=str(raw_input('Would you like to explore asteroids at L1, L2, L3, L4 or L5? '))
if Asteroid in ['l4', 'L4', 'L5', 'l5']:

    time=float(raw_input('For how many Jupiter years would you like to run the simulation? '))
    T=time*year
    Number=float(raw_input('How many asteroids would you like to inlcude in the simulation? '))

    #initial conditions for planet
    initial_x_planet= R_peri
    initial_y_planet=0
    initial_vx_planet=0
    initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
    initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]

    #Time Period for Integration
    div=10000
    t=sp.linspace(0.,T,div)
    t2=t

    if Asteroid in ['L4', 'l4']:
    #Initial Condition for Reference Asteroid at L4
        initial_x_asteroidref=R_peri*0.5
        initial_y_asteroidref=R_peri*(3**0.5)*0.5
        initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
        initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5                 
        initial_vx_asteroidref=-1*initial_speed_asteroidref*(3**0.5)*0.5
        initial_vy_asteroidref=initial_speed_asteroidref*0.5
        initial_conditions_asteroidref=[initial_x_asteroidref,initial_vx_asteroidref,initial_y_asteroidref,initial_vy_asteroidref]

    if Asteroid in ['L5', 'l5']:
        #Initial Condition for Reference Asteroid at L5
        initial_x_asteroidref=R_peri*0.5
        initial_y_asteroidref=-1*R_peri*(3**0.5)*0.5
        initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
        initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5                 
        initial_vx_asteroidref=initial_speed_asteroidref*(3**0.5)*0.5
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

    x_sunCOM=-1*(M_planet/M_total)*x_planet_old
    y_sunCOM=-1*(M_planet/M_total)*y_planet_old

    x_planetCOM=(1-(M_planet/M_total))*x_planet_old
    y_planetCOM=(1-(M_planet/M_total))*y_planet_old

    def asteroid(asteroidpos,t2,x_planet_old,y_planet_old):
        
            x_planet_new=np.interp(t2,t,x_planet_old)
            y_planet_new=np.interp(t2,t,y_planet_old)
        
            x_asteroid=asteroidpos[0]
            y_asteroid=asteroidpos[2]
    
            vx_asteroid=asteroidpos[1]
            vy_asteroid=asteroidpos[3]
    
            ax_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*x_asteroid -(((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(x_asteroid-x_planet_new))
            ay_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*y_asteroid- (((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(y_asteroid-y_planet_new))
    
            return [vx_asteroid,ax_asteroid,vy_asteroid,ay_asteroid]
    
    solnasteroidref=spi.odeint(asteroid,initial_conditions_asteroidref,t2,args=(x_planet_old,y_planet_old))
    x_asteroidref=solnasteroidref[:,0]
    y_asteroidref=solnasteroidref[:,2]
    if Asteroid in ['L4', 'l4', 'L5', 'l5']:
        Initial=0.9*R_peri
        Final=1.1*R_peri
    loopy=np.arange(Initial,Final,(Final-Initial)/Number)
            
    if Asteroid in ['L4', 'l4']:
        for i in loopy:
            initial_x_asteroid=i*0.5
            initial_y_asteroid=i*(3**0.5)*0.5
            initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
            initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
            solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
            x_asteroid=solnasteroid[:,0]
            y_asteroid=solnasteroid[:,2]
    
            x_asteroidframe=x_asteroid-x_asteroidref
            y_asteroidframe=y_asteroid-y_asteroidref
        
            fig1=pl.figure(1)
            ax1 = fig1.gca(projection='3d')
            
            ax1.set_xlim((-10e2, 10e2))
            ax1.set_ylim((-10e2, 10e2))
            ax1.set_zlim((-10e2, 10e2))
    
            ax1.patch.set_facecolor('black')
            ax1.axis('off')
            pl.title('Planet and Asteroids in SS Frame')
            #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
            pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
            #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
            ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
            zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
            ax1.plot_surface(xs, ys, zs, color='r')

            pl.show()
        
            #pl.figure(2)
            #pl.title('Sun and Planet about their COM')
            #pl.plot(x_sunCOM/1e9,y_sunCOM/1e9)
            #pl.plot(x_planetCOM/1e9,y_planetCOM/1e9)
            #pl.show()
        
            #fig2=pl.figure(3)
           
            #ax2 = fig2.add_subplot(1,1,1)
            #ax2.patch.set_facecolor('black')
            
            #pl.title('Asteroid in L4 Frame')
            #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
            #pl.show()
    
    if Asteroid in ['L5', 'l5']:
        for i in loopy:
            initial_x_asteroid=i*0.5
            initial_y_asteroid=-1*i*(3**0.5)*0.5
            initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
            initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
            solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
            x_asteroid=solnasteroid[:,0]
            y_asteroid=solnasteroid[:,2]
    
            x_asteroidframe=x_asteroid-x_asteroidref
            y_asteroidframe=y_asteroid-y_asteroidref
        
            fig3=pl.figure(1)
            ax3 = fig3.gca(projection='3d')
            ax3.set_xlim((-10e2, 10e2))
            ax3.set_ylim((-10e2, 10e2))
            ax3.set_zlim((-10e2, 10e2))
    
            ax3.patch.set_facecolor('black')
            ax3.axis('off')
            pl.title('Planet and Asteroids in SS Frame')
            #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
            pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
            #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
            ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
            zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
            ax3.plot_surface(xs, ys, zs, color='r')
            pl.show()
        
            #pl.figure(2)
            #pl.title('Sun and Planet about their COM')
            #pl.plot(x_sunCOM/1e9,y_sunCOM/1e9)
            #pl.plot(x_planetCOM/1e9,y_planetCOM/1e9)
            #pl.show()
        
            #fig4=pl.figure(3)
            
            #ax4 = fig4.add_subplot(1,1,1)
            #ax4.patch.set_facecolor('black')
            
            #pl.title('Asteroid in L5 Frame')
            #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
            #pl.show()

if Asteroid in ['l1', 'L1']:
    #initial conditions for planet
    initial_x_planet= R_peri
    initial_y_planet=0
    initial_vx_planet=0
    initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
    initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]

    #Time Period for Integration
    T=1*year
    div=10000
    t=sp.linspace(0.,T,div)
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

    x_sunCOM=-1*(M_planet/M_total)*x_planet_old
    y_sunCOM=-1*(M_planet/M_total)*y_planet_old

    x_planetCOM=(1-(M_planet/M_total))*x_planet_old
    y_planetCOM=(1-(M_planet/M_total))*y_planet_old

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
        
    loopy=np.arange(1.01*initial_x_asteroidref,1.4*initial_x_asteroidref,0.04*initial_x_asteroidref)
    
    Max_dispnew=[]
    for i in loopy:
        initial_x_asteroid=i
        initial_y_asteroid=0
        initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
        initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
        solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
        x_asteroid=solnasteroid[:,0]
        y_asteroid=solnasteroid[:,2]
    
        x_asteroidframe=x_asteroid-x_asteroidref
        y_asteroidframe=y_asteroid-y_asteroidref
        
        fig5=pl.figure(1)
        ax5 = fig5.gca(projection='3d')
        ax5.set_xlim((-10e2, 10e2))
        ax5.set_ylim((-10e2, 10e2))
        ax5.set_zlim((-10e2, 10e2))

        ax5.patch.set_facecolor('black')
        ax5.axis('off')
        pl.title('Planet and Asteroids in SS Frame')
        #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
        pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
        #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
        ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
        zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
        ax5.plot_surface(xs, ys, zs, color='r')
        pl.show()
        
        #pl.figure(2)
        #pl.title('Sun and Planet about their COM')
        #pl.plot(x_sunCOM/1e9,y_sunCOM/1e9)
        #pl.plot(x_planetCOM/1e9,y_planetCOM/1e9)
        #pl.show()
        
        #fig6=pl.figure(3)
        #ax6 = fig6.add_subplot(1,1,1)
        #ax6.patch.set_facecolor('black')
        
        #pl.title('Asteroid in L1 Frame')
        #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
        #pl.show()

if Asteroid in ['l2', 'L2']:
    initial_x_planet= R_peri
    initial_y_planet=0
    initial_vx_planet=0
    initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
    initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]

    #Time Period for Integration
    T=1*year
    div=30000
    t=sp.linspace(0.,3*T,div)
    t2=t
    
    initial_x_asteroidref=R_peri*(1+((alpha/3)**(1./3.)))
    initial_y_asteroidref=0
    initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
    initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5 
    initial_speed_asteroidrefnew=(2*3.14*initial_x_asteroidref)/(1*year)             
    initial_vx_asteroidref=0
    initial_vy_asteroidref=initial_speed_asteroidref
    initial_conditions_asteroidref=[initial_x_asteroidref,initial_vx_asteroidref,initial_y_asteroidref,initial_speed_asteroidrefnew]
    
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
    
    x_sunCOM=-1*(M_planet/M_total)*x_planet_old
    y_sunCOM=-1*(M_planet/M_total)*y_planet_old
    
    x_planetCOM=(1-(M_planet/M_total))*x_planet_old
    y_planetCOM=(1-(M_planet/M_total))*y_planet_old
    
    def asteroid(asteroidpos,t2,x_planet_old,y_planet_old):
            
            x_planet_new=np.interp(t2,t,x_planet_old)
            y_planet_new=np.interp(t2,t,y_planet_old)
            
            r_plan=np.sqrt((x_planet_new**2)+(y_planet_new**2))
            
            x_asteroid=asteroidpos[0]
            y_asteroid=asteroidpos[2]
            r_ast=np.sqrt((x_asteroid**2)+(y_asteroid**2))
        
            if r_ast==r_plan:
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
    
    loopy=np.arange(1.01*initial_x_asteroidref,1.4*initial_x_asteroidref,0.04*initial_x_asteroidref)
    Max_disp=[]
    for i in loopy:
        initial_x_asteroid=i
        initial_y_asteroid=0
        initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
        initial_conditions_asteroid=[initial_x_asteroid, 0, initial_y_asteroid, initial_speed_asteroidrefnew]
        solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
        x_asteroid=solnasteroid[:,0]
        y_asteroid=solnasteroid[:,2]
        
        x_asteroidframe=x_asteroid-x_asteroidref
        y_asteroidframe=y_asteroid-y_asteroidref
        
        fig7=pl.figure(1)
        ax7 = fig7.gca(projection='3d')
        ax7.set_xlim((-10e2, 10e2))
        ax7.set_ylim((-10e2, 10e2))
        ax7.set_zlim((-10e2, 10e2))

        ax7.patch.set_facecolor('black')
        ax7.axis('off')
        pl.title('Planet and Asteroids in SS Frame')
        #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
        pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
        #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
        ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
        zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
        ax7.plot_surface(xs, ys, zs, color='r')
        pl.show()
        
        #pl.figure(2)
        #pl.title('Sun and Planet about their COM')
        #pl.plot(x_sunCOM/1e9,y_sunCOM/1e9)
        #pl.plot(x_planetCOM/1e9,y_planetCOM/1e9)
        #pl.show()
        
        #fig8=pl.figure(3)
        #ax8 = fig8.add_subplot(1,1,1)
        #ax8.patch.set_facecolor('black')
        #pl.title('Asteroid in L2 Frame')
        #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
        #pl.show()

if Asteroid in ['l3', 'L3']:
    #initial conditions for planet
    initial_x_planet= R_peri
    initial_y_planet=0
    initial_vx_planet=0
    initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
    initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]
    
    #Time Period for Integration
    T=1*year
    div=10000
    t=sp.linspace(0.,2*T,div)
    t2=t
    
    #Initial Condition for Reference Asteroid at L1
    initial_x_asteroidref=-1*R_peri*(1+((5./12.)*alpha))
    initial_y_asteroidref=0
    initial_R_asteroidref=((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5
    initial_speed_asteroidref=(((G*M_sun)/(((initial_x_asteroidref**2)+(initial_y_asteroidref**2))**0.5)))**0.5                 
    initial_vx_asteroidref=0
    initial_vy_asteroidref=-1*initial_speed_asteroidref
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
    
    x_sunCOM=-1*(M_planet/M_total)*x_planet_old
    y_sunCOM=-1*(M_planet/M_total)*y_planet_old
    
    x_planetCOM=(1-(M_planet/M_total))*x_planet_old
    y_planetCOM=(1-(M_planet/M_total))*y_planet_old
    
    def asteroid(asteroidpos,t2,x_planet_old,y_planet_old):
            
            x_planet_new=np.interp(t2,t,x_planet_old)
            y_planet_new=np.interp(t2,t,y_planet_old)
            
            x_asteroid=asteroidpos[0]
            y_asteroid=asteroidpos[2]
        
            vx_asteroid=asteroidpos[1]
            vy_asteroid=asteroidpos[3]
        
            ax_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*x_asteroid -(((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(x_asteroid-x_planet_new))
            ay_asteroid=-((G*M_sun)/(((x_asteroid**2)+(y_asteroid**2))**1.5))*y_asteroid- (((G*M_planet)/((((x_asteroid-x_planet_new)**2)+((y_asteroid-y_planet_new)**2))**1.5))*(y_asteroid-y_planet_new))
        
            
            return [vx_asteroid,ax_asteroid,vy_asteroid,ay_asteroid]
        
    
    solnasteroidref=spi.odeint(asteroid,initial_conditions_asteroidref,t2,args=(x_planet_old,y_planet_old))
    x_asteroidref=solnasteroidref[:,0]
    y_asteroidref=solnasteroidref[:,2]
    
    loopy=np.arange(0.98*initial_x_asteroidref,1.1*initial_x_asteroidref,0.01*initial_x_asteroidref)
    Max_dispnew=[]
    for i in loopy:
        initial_x_asteroid=i
        initial_y_asteroid=0
        initial_R_asteroid=((initial_x_asteroid**2)+(initial_y_asteroid**2))**0.5
    
        initial_conditions_asteroid=[initial_x_asteroid, initial_vx_asteroidref, initial_y_asteroid, initial_vy_asteroidref]
        solnasteroid=spi.odeint(asteroid,initial_conditions_asteroid,t2,args=(x_planet_old,y_planet_old))
    
        x_asteroid=solnasteroid[:,0]
        y_asteroid=solnasteroid[:,2]
    
        x_asteroidframe=x_asteroid-x_asteroidref
        y_asteroidframe=y_asteroid-y_asteroidref
        
        fig9=pl.figure(1)
        ax9 = fig9.gca(projection='3d')
        ax9.set_xlim((-10e2, 10e2))
        ax9.set_ylim((-10e2, 10e2))
        ax9.set_zlim((-10e2, 10e2))

        ax9.patch.set_facecolor('black')
        ax9.axis('off')
        pl.title('Planet and Asteroids in SS Frame')
        #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
        pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
        #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
        ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
        zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
        ax9.plot_surface(xs, ys, zs, color='r')
        pl.show()
        
        #pl.figure(2)
        #pl.title('Sun and Planet about their COM')
        #pl.plot(x_sunCOM/1e9,y_sunCOM/1e9)
        #pl.plot(x_planetCOM/1e9,y_planetCOM/1e9)
        #pl.show()
        
        #fig0=pl.figure(3)
        #ax0 = fig0.add_subplot(1,1,1)
        #ax0.patch.set_facecolor('black')
        #pl.title('Asteroid in L3 Frame')
        #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
        #pl.show()

answer=str(raw_input('Would you like to see what changing the mass of the planet does? (yes or no) '))  
if answer in ['yes', 'Yes', 'YES']:
    
    year=year_jupiter
    R_peri=R_peri_jupiter
    Initial=0.005*M_sun
    Final=0.06*M_sun
    Number=10
    M_planet=np.arange(Initial,Final,(Final-Initial)/Number)
    Max_displacements=[]
    
    for i in M_planet:
    
        M_total=i+M_sun
    
        #initial conditions for planet
        initial_x_planet= R_peri
        initial_y_planet=0
        initial_vx_planet=0
        initial_vy_planet= (((G*M_sun)/(((initial_x_planet**2)+(initial_y_planet**2))**0.5)))**0.5
        initial_conditions_planet = [initial_x_planet, initial_vx_planet, initial_y_planet, initial_vy_planet]
    
        #Time Period for Integration
        T=1*year
        div=10000
        t=sp.linspace(0.,5*T,div)
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
    
        x_asteroidframe=x_asteroid-x_asteroidref
        y_asteroidframe=y_asteroid-y_asteroidref
    
        r_asteroidFrm=np.sqrt((x_asteroidframe**2)+(y_asteroidframe**2))
        Max_displacements.append(np.amax(r_asteroidFrm))
        
        fig01=pl.figure(6)
        ax01 = fig01.gca(projection='3d')
        ax01.set_xlim((-10e2, 10e2))
        ax01.set_ylim((-10e2, 10e2))
        ax01.set_zlim((-10e2, 10e2))

        ax01.patch.set_facecolor('black')
        ax01.axis('off')
        
        #pl.plot(x_planet_old/1e9,y_planet_old/1e9)
        pl.plot(x_asteroid/1e9,y_asteroid/1e9,0)
        #pl.plot(x_asteroidref/1e9,y_asteroidref/1e9)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        xs = (R_sun) * np.outer(np.cos(u), np.sin(v))
        ys = (R_sun) * np.outer(np.sin(u), np.sin(v))
        zs = (R_sun) * np.outer(np.ones(np.size(u)), np.cos(v))
        ax01.plot_surface(xs, ys, zs, color='r')
        pl.show()
        
        #fig02=pl.figure(7)
        #ax02 = fig02.add_subplot(1,1,1)
        #ax02.patch.set_facecolor('black') 
        #pl.title('Asteroid in L4 Frame (CHANGING MASS)')
        #pl.plot(x_asteroidframe/1e9,y_asteroidframe/1e9)
        #pl.show()
    
    #pl.figure(8)
    #pl.title('Max Displacement in AsteroidRef Frame vs Mass of Planet (CHANGING MASS)')
    #pl.xlabel('Mass of planet in units of Msun')
    #pl.ylabel('Maximum Displacement in units of Million kms')
    #pl.plot(M_planet/M_sun,Max_displacements)
    #pl.show()




    

        

    
    
    
    
    
    

    