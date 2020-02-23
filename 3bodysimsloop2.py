
import matplotlib.pyplot as plt
'''functions'''
def sqrt(n):   #Calculate squareroot
    return n**0.5


m2s=[]#the masses to run in Mearth

for x in [100,200,300,400,500,600,700,800,900,1000]:
    m2s.append(x)



        
for n in range(len(m2s)):

    '''constants'''
    pi=3.1415927
    G=1.967e-44 #AU^3/kg s^2

    msun=2e30 #kg
    mearth=5.97e24 #kg
    Mjupiter=318 #Mearth

    spd=86400 #sec per Earth day
    dpy=365 #earth day per Earth year
    spy=dpy*spd #sec per Earth year


    '''inputs'''
    mstaro=1 #Msun
    m1o=1   #Mearth
    m2o=m2s[n] #Mearth
    p1o=100 #days          ##KOI872b=33.6
    p2o=200 #days            ##KOI872c=57
    tstep=864 #time steps in sec
    tlimitd=100000 #days


    '''convert input units'''
    mstar=mstaro*msun #mass of star in kg
    m1=m1o*mearth #mass of planet 1 in kg
    m2=m2o*mearth #mass of planet 2 in kg
    p1=p1o/dpy  #period of p1 in year
    p2=p2o/dpy  #period of p2 in year

    p1s=p1o*spd #period of p1 in second
    p2s=p2o*spd #period of p2 in second


    '''init velo&distance in AU'''
    r1=(p1**2*mstaro)**0.333  #semimajor axis of planet 1 in AU :P^2=r^3/m
    v1x=0
    v1y=2*pi*r1/p1
    v1y=v1y/spy  #velocity of p1 in AU/s; P=2pi R/v

    r2=(p2**2*mstaro)**0.333
    v2x=0
    v2y=2*pi*r2/p2
    v2y=v2y/spy


    '''init position'''
    x1=r1
    y1=0

    x2=r2
    y2=0


    '''time controls'''
    tlimit=tlimitd*spd #time limit in sec


    '''arrays to store postions'''
    arax1=[]
    aray1=[]

    arax2=[]
    aray2=[]




    '''initialize values for calculating periods'''
    yold1=y1
    yold2=y2
    told1=0
    told2=0     #the last time a period is finished

    arap1=[p1s]  #periods

    arap2=[p2s]  #periods

    f=0 #flag
    magnitude=0  #The difference between longest and shortest period of p1


    '''tracing'''
    for t in range(0,tlimit,tstep):

        x12=x1-x2
        y12=y1-y2
        r21=sqrt(x12**2+y12**2) #distance between 1 and 2
        
        #planet 1
        r1=sqrt(x1**2+y1**2) #distance between planet and star
        
        a1=G*mstar/r1**2
        a12=G*m2/r21**2
        
        a1x=-a1*x1/r1
        a1y=-a1*y1/r1
        a12x=-a12*x12/r21  #acceration due to gravity of second star
        a12y=-a12*y12/r21
        
        v1x=v1x+a1x*tstep+a12x*tstep #v=u+at
        v1y=v1y+a1y*tstep+a12y*tstep

        x1=x1+v1x*tstep #s=s+vt
        y1=y1+v1y*tstep


        arax1.append(x1)
        aray1.append(y1)

        #planet 2
        r2=sqrt(x2**2+y2**2) #distance between planet and star
        
        a2=G*mstar/r2**2
        a21=G*m1/r21**2
        
        a2x=-a2*x2/r2
        a2y=-a2*y2/r2
        a21x=a21*x12/r21
        a21y=a21*y12/r21
        
        v2x=v2x+a2x*tstep+a21x*tstep #v=u+at
        v2y=v2y+a2y*tstep+a21y*tstep

        x2=x2+v2x*tstep #s=s+vt
        y2=y2+v2y*tstep

        arax2.append(x2)
        aray2.append(y2)

        #periods
        if y1>=0 and yold1<=0 and f==1: #change in sign of y axis indicates a transit
            temp=t-told1 #the period
            arap1.append(temp)
            told1=t
        if y2>=0 and yold2<=0 and f==1:
            temp=t-told2
            arap2.append(temp)
            told2=t
            
        yold1=y1
        yold2=y2
        f=1

    m=max(arap1)-min(arap1)
    print(m)
    po=plt.figure()
    plt.suptitle("p1:blue, period:"+str(p1o)+"days"
                 "\np2:red, period:"+str(p2o)+"days"
                 "\nmass1:"+str(m1o)+"Mearth:"
                 "\nmass2:"+str(m2o)+"Mearth:")
    plt.xlabel("orbits")
    plt.ylabel("period(second)")
    plt.plot(arap1,'b')
    plt.plot(arap2,'r')
    plt.grid()
    plt.savefig(str(m2o)+'.pdf')



