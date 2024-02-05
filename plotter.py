import matplotlib.pyplot as plt
import numpy as np
import csv



NAME="/home/chris/Desktop/dynamical/data"

time_steps=20
l=50
sims=40;
particles=200;
X=1000;

def extrapolate(step,file):
    name=NAME+str(file)+".txt"
    b_val=np.zeros(X)
    u_val=np.zeros(X)
    
    f = open(name, "r")
    with f as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for i,row in enumerate(plots):
           
            if particles*step<=i<particles*(step+1):
                if int(row[1])==1 :
                    x=int(row[0])
                    b_val[x-l:x+l]+=0.5/l               
                else:
                    x=int(row[0])
                    u_val[x-l:x+l]+=0.5/l
                
    return b_val,u_val

def create_densities(file):
    
    b_sum= np.zeros(X)
    u_sum= np.zeros(X)
    
    b_var= np.zeros(X)
    u_var =np.zeros(X)
    for i in range(0,sims):
        b,u =np.array(extrapolate(i,file))
        
        b_sum+=b
        u_sum+=u
        
        b_var+=np.square(b)
        u_var+=np.square(u)
    
    b_var=b_var/sims
    u_var=u_var/sims
   
    
    b_avg=b_sum/sims
    u_avg=u_sum/sims
    
    b_var-=np.square(b_avg)
    u_var-=np.square(u_avg)
    
    return b_avg, b_var ,u_avg, u_var  

def plot_stats_for_file(file):
    b,B,u,U=create_densities(file)
   
    g=np.sqrt(B/sims)
    v=np.sqrt(U/sims)
    
    x=np.linspace(0,10,X)
    shower=np.full((X),0.618)
    shower2=np.full((X),0.38)
    fig,(ax2,ax4)= plt.subplots(2)
   
    ax2.plot(x,b-g,'b')
    ax2.plot(x,b,'r')
    ax2.plot(x,b+g,'b')
    ax2.plot(x,shower2,'k')
    ax2.yaxis.set_ticks(np.arange(0,1.2,0.2))
    
    ax2.xaxis.set_ticks(np.arange(0,10,2))
   
    
       
    
    ax4.plot(x,u-v,'b')
    ax4.plot(x,u,'r')
    ax4.plot(x,u+v,'b')
    ax4.plot(x,shower,'k')
    ax4.yaxis.set_ticks(np.arange(0,1.2,0.2))   
    
    ax4.xaxis.set_ticks(np.arange(0,10,5))
    
   
    fig.savefig(str(file)+".png")

def plot_stats ():    
    for i in range(0,time_steps):
        plot_stats_for_file(i)


plot_stats()


