# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:58:11 2020

@author: Angad Bajwa,Aditya Pethkar,Mandar Burande
"""
#Libraries required
import cmath
import math
import matplotlib.pyplot as plt
from datetime import datetime
from decimal import Decimal

#Function to convert floating point numbers into scientific notation
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

#Function to input variables
def user_input():
    variables = {}

#extracts information from every line of .txt file by splitting the line from the '=' symbol
    with open("input.txt") as f:
        for line in f:
            name, value = line.split("=")
            variables[name] = float(value)

    system = variables["Type of the system( 1-Symmetrical ,2-Unsymmetrical )"]
    phase_spacing_sym = variables["Spacing between the phase conductors(m)(Symmetrical)"]
    a1 = variables["Spacing between the phase conductors(enter distance)(m) a(Unsymmetrical)"]
    b1 = variables["Spacing between the phase conductors(enter distance)(m) b(Unsymmetrical)"]
    c1 = variables["Spacing between the phase conductors(enter distance)(m) c(Unsymmetrical)"]
    no_sub_cond = variables["Number of sub-conductors per bundle"]
    sub_cond_spacing = variables["Spacing between the sub-conductors(m)"]
    no_strands = variables["Number of strands in each sub-conductor"]
    d = variables["Diameter of each strand(m)"]
    length_of_line = variables["Length of line(km)"]
    model = variables["Model of Line(1-Short , 2-Nominal Pi ,3-Long )"]
    resistance = variables["Resistance of line per km(Ohms)"]
    freq = variables["Power Frequency(Hz)"]
    nom_volt = variables["Nominal Voltage(V)"]
    power_rec = variables["Receiving end load(MW)"]
    pf = variables["Power factor of receiving end load"]



#dictionary to pass the input parameters to calculations function
    parameters={
                "system":system,
                "phase_spacing_sym":phase_spacing_sym,
                "a1":a1,
                "b1":b1,
                "c1":c1,
                "no_subconductors":no_sub_cond,
                "sub_cond_spacing":sub_cond_spacing,
                "no_strands":no_strands,
                "diameter":d,
                "length_of_line":length_of_line,
                "model":model,
                "resistance":resistance,
                "freq":freq,
                "nom_volt":nom_volt,
                "load":power_rec,
                "pf":pf,
        }
    return parameters

#Function to perform calculations
def calculations(parameters):

#Receving the input parameters from the user_input function

    system=parameters["system"]
    phase_spacing_sym=parameters['phase_spacing_sym']
    a1=parameters["a1"]
    b1=parameters["b1"]
    c1=parameters["c1"]
    no_sub_conductors=parameters["no_subconductors"]
    sub_cond_spacing=parameters["sub_cond_spacing"]
    no_strands=parameters["no_strands"]
    d=parameters["diameter"]
    length=parameters["length_of_line"]
    model=parameters["model"]
    resistance=parameters["resistance"]
    freq=parameters["freq"]
    nom_volt=parameters["nom_volt"]
    power_rec=parameters["load"]
    pf=parameters["pf"]


#Receiving End Voltage

    VR=nom_volt/math.sqrt(3)


#Calculating radius---r
    a = 3
    b = -3
    c = 1-no_strands
 # calculate the discriminant
    di = (b ** 2) - (4 * a * c)

 # find two solutions
    sol1 = (-b - cmath.sqrt(di)) / (2 * a)
    sol2 = (-b + cmath.sqrt(di)) / (2 * a)

    if (sol1.real >= sol2.real):
        n = sol1.real
    else:
        n= sol2.real


    D=(2*n-1)*d
    r=D/2

# Inductance---L

    if(system==1):
        if(no_sub_conductors==2):
            h=math.pow(no_sub_conductors,2)
            i=math.pow(phase_spacing_sym,6)/(math.pow(r,1)*math.exp(-1/4)*math.pow(sub_cond_spacing,1))
            q=math.pow(i,1/h)
            print(q)

        elif(no_sub_conductors==3):
            h = math.pow(no_sub_conductors, 2)
            i = math.pow(phase_spacing_sym, 6) / math.pow(r*math.pow(sub_cond_spacing,2)*math.exp(-1/4),3)
            q = math.pow(i, 1 / h)
            print(q)

        elif(no_sub_conductors==4):
            h = math.pow(no_sub_conductors, 2)
            i = math.pow(phase_spacing_sym, 6) / math.pow(r*math.pow(sub_cond_spacing,3)*math.exp(-1/4)*math.sqrt(2),4)
            q = math.pow(i, 1 / h)

        else:
            h = math.pow(no_sub_conductors, 2)
            i = math.pow(phase_spacing_sym, 6) / math.pow(r*math.exp(-1 / 4) * math.pow(sub_cond_spacing, 5) * 6,6)
            q = math.pow(i, 1 / h)


    else:
        if(no_sub_conductors==2):

            h=math.pow(a1*b1*c1,1/3)
            i=math.exp(-1/(4*no_sub_conductors))*math.pow(r,1/no_sub_conductors)*math.pow(sub_cond_spacing,1/no_sub_conductors)
            q=h/i

        elif(no_sub_conductors==3):
            h = math.pow(a1 * b1 * c1, 1 / 3)
            i = math.pow( r * math.pow(sub_cond_spacing, 2) * math.exp(-1 / 4), 1 / no_sub_conductors)
            q = h / i


        elif(no_sub_conductors==4):
            h = math.pow(a1 * b1 * c1, 1 / 3)
            i = math.pow(math.sqrt(2)*r*math.pow(sub_cond_spacing,3)*math.exp(-1/4),1/no_sub_conductors)
            q = h / i

        else:
            h = math.pow(a1 * b1 * c1, 1 / 3)
            i = math.pow(6 * r * math.pow(sub_cond_spacing, 5) * math.exp(-1 / 4), 1 / no_sub_conductors)
            q = h / i

    L=2*math.pow(10,-4)*math.log(q)


#Capacitance---Cap

    if(model==1):
        Cap=0

    else:

        if (system == 1):

            if (no_sub_conductors == 2):
                h = math.pow(no_sub_conductors, 2)
                i = math.pow(phase_spacing_sym, 6) / (math.pow(r, 1)  * math.pow(sub_cond_spacing, 1))
                q = math.pow(i, 1 / h)
                print(q)

            elif (no_sub_conductors == 3):
                h = math.pow(no_sub_conductors, 2)
                i = math.pow(phase_spacing_sym, 6) / math.pow(r * math.pow(sub_cond_spacing, 2) , 3)
                q = math.pow(i, 1 / h)
                print(q)

            elif (no_sub_conductors == 4):
                h = math.pow(no_sub_conductors, 2)
                i = math.pow(phase_spacing_sym, 6) / math.pow(
                    r * math.pow(sub_cond_spacing, 3) * math.sqrt(2), 4)
                q = math.pow(i, 1 / h)

            else:
                h = math.pow(no_sub_conductors, 2)
                i = math.pow(phase_spacing_sym, 6) / math.pow(r * math.pow(sub_cond_spacing, 5) * 6,
                                                              6)
                q = math.pow(i, 1 / h)


        else:
            if (no_sub_conductors == 2):

                h = math.pow(a1 * b1 * c1, 1 / 3)
                i = math.pow(r, 1 / no_sub_conductors) * math.pow(sub_cond_spacing, 1 / no_sub_conductors)
                q = h / i

            elif (no_sub_conductors == 3):
                h = math.pow(a1 * b1 * c1, 1 / 3)
                i = math.pow(r * math.pow(sub_cond_spacing, 2), 1 / no_sub_conductors)
                q = h / i

            elif (no_sub_conductors == 4):
                h = math.pow(a1 * b1 * c1, 1 / 3)
                i = math.pow(math.sqrt(2) * r * math.pow(sub_cond_spacing, 3) , 1 / no_sub_conductors)
                q = h / i

            else:
                h = math.pow(a1 * b1 * c1, 1 / 3)
                i = math.pow(6 * r * math.pow(sub_cond_spacing, 5) , 1 / no_sub_conductors)
                q = h / i

    Cap = 2 * math.pi * 8.84 * math.pow(10, -9) / math.log(q)


# Inductive Reactance---XL

    XL = 2 * math.pi * freq * L * length

#Capacitive Reactance---XC
    if(model==1):
        XC= "infinite"
    else:
        XC=1/(2*math.pi*freq*Cap*length)

# ABCD paramters

    if(model==1):
        Y=0
    else:
        Y=complex(0,1/XC)


    Z=complex(resistance*length,XL)

    if(model == 1):
        A = 1
        B = Z
        C = Y
        D = 1

    elif(model == 2):
        A = 1 + Y * Z / 2
        B = Z
        C = Y * (1 + Y * Z / 4)
        D = 1 + Y * Z / 2

    else:
        A = cmath.cosh(cmath.sqrt(Y * Z))
        B = cmath.sqrt(Z / Y) * cmath.sinh(cmath.sqrt(Y * Z))
        C = cmath.sqrt(Y / Z) * cmath.sinh(cmath.sqrt(Y * Z))
        D = cmath.cosh(cmath.sqrt(Y * Z))


# Sending end line voltage---VS

    # Receving end current---IR

    IR = power_rec * math.pow(10, 6) / (pf * VR * 3)

    VS = (A*VR + IR*B)

#charging current---IC

    if(model==1):
        IC=0

    elif(model==2):
        IC1=VS*Y/2
        IC2=VR*Y/2
        IC=IC1+IC2
    else:
        IC=C*VR

# Sending end line Current----IS

    IS=C*VR + D*IR

#Voltage Regulation---volt_reg

    if(model==1):
        volt_reg= (abs(VS)-abs(VR))*100/abs(VR)

    else:
        volt_reg= (abs(VS/A) -abs(VR))*100/abs(VR)


#Power loss---power_loss

    power_loss = 3 * IR.real * IR.real * resistance * length

#Transmission Efficiency---eff

    eff=power_rec*100/(power_rec+power_loss*math.pow(10,-6))

#Angles

    angle_1=cmath.phase(A) #alpha
    angle_2=cmath.phase(B) #beta

# Sending end circle

    x= -abs(A)*abs(VS/1000)*abs(VS/1000)*math.cos(angle_2-angle_1)/abs(B)
    y=abs(A)*abs(VS/1000)*abs(VS/1000)*math.sin(angle_2-angle_1)/abs(B)
    rad_1=abs(VS/1000)*abs(VR/1000)/abs(B)
    centre_1 = complex(x, y)
    fig, ax = plt.subplots()
    ax.set(xlim=(-5*rad_1, 5*rad_1), ylim=(-5*rad_1, 5*rad_1))
    a_circle = plt.Circle((x, y), rad_1)
    ax.add_artist(a_circle)
    plt.title('Sending End Circle', fontsize=12)
    plt.ylabel("Apparent Power")
    plt.xlabel("Real Power")
    plt.text(100, 200, 'Centre= %s\n Radius = %s\n' % (
        "{:g}".format(centre_1), "{:.3f}".format(rad_1)))


#Receving end circle

    x1 = -abs(A) * abs(VR/1000) * abs(VR/1000) * math.cos(angle_2 - angle_1) / abs(B)
    y1 = -abs(A) * abs(VR/1000) * abs(VR/1000) * math.sin(angle_2 - angle_1) / abs(B)
    rad_2= abs(VS/1000) * abs(VR/1000) / abs(B)
    centre_2 = complex(x1, y1)
    fig, ax_1 = plt.subplots()
    ax_1.set(xlim=(-5*rad_2,5*rad_2), ylim=(-5*rad_2,5*rad_2))
    b_circle = plt.Circle((x1,y1), rad_2,color='r')
    ax_1.add_artist(b_circle)
    plt.title('Receving End Circle', fontsize=12)
    plt.ylabel("Apparent Power")
    plt.xlabel("Real Power")
    plt.text(100, 200, 'Centre= %s\n Radius = %s\n' % (
    "{:g}".format(centre_2),"{:.3f}".format(rad_2)))
    plt.show()

#Dictionary to pass values to output function
    features = {
        "Inductance": L,
        "Capacitance":Cap,
        "XL":XL,
        "XC":XC,
        "IC":IC,
        "A":A,
        "B":B,
        "C":C,
        "D":D,
        "VS":VS,
        "IS":IS,
        "volt_reg":volt_reg,
        "power_loss":power_loss,
        "eff":eff,
        "centre_1":centre_1,
        "rad_1":rad_1,
        "centre_2":centre_2,
        "rad_2":rad_2,
        "model":model
    }
    return features


#Output Function
def output(features):
    L=features["Inductance"]
    Cap=features['Capacitance']
    XL=features["XL"]
    XC=features["XC"]
    IC=features["IC"]
    A=features["A"]
    B=features["B"]
    C=features["C"]
    D=features["D"]
    VS=features["VS"]
    IS=features["IS"]
    volt_reg=features["volt_reg"]
    power_loss=features["power_loss"]
    eff=features["eff"]
    centre_1=features["centre_1"]
    rad_1=features["rad_1"]
    centre_2=features["centre_2"]
    rad_2=features["rad_2"]
    model=features["model"]



    n=datetime.now()
    f = open('output.txt', 'w')
    print('Team Members\n', file=f)
    print('a. Angad Bajwa-107118014', file=f)
    print("b. Mandar Burande-107118056",file=f)
    print('c. Aditya Pethkar-107118072\n',file=f)
    print('Date and time of this submission :\t',n, file=f)
    print('\nInductance per phase per km in H/km :\t', format_e(Decimal(L))," H/km",file=f)
    print('Capacitance per phase per km in F/km :\t', format_e(Decimal(Cap))," F/km",file=f)
    print('Inductive Reactance of the line in Ohm :\t' +"{:.3f}".format(XL),"Ohm", file=f)
    if(model==1):
        print('Capacitive Reactance of the line in Ohm :\t' + XC, "Ohm", file=f)
    else:
        print('Capacitive Reactance of the line in Ohm :\t' + "{:.3f}".format(XC), "Ohm", file=f)
    print('Charging current drawn for the sending end substation :' +"{:.3f}".format(abs(IC))," A", file=f)
    print('ABCD Parameters',file=f)
    if(model==1):
        print("Short Transmission Line",file=f)
    elif(model==2):
        print("Nominal Pi Model Medium Transmission Line",file=f)
    else:
        print("Long Transmission Line",file=f)

    print('A= ',A,"\t\t",file=f)
    print('B= ',B,"Ohm",file=f)
    print('C= ', C,"mho\t\t",file=f)
    print('D= ', D,file=f)
    print('Sending end Voltage :' +"{:.3f}".format(abs(math.sqrt(3)*VS)/1000), "kV", file=f)
    print('Sending end Current :' +"{:.3f}".format(abs(IS)), " A", file=f)
    print('Voltage Regulation :\t' +"{:.3f}".format(volt_reg), " %", file=f)
    print('Power Loss :' +"{:.3f}".format(power_loss*math.pow(10,-6)), " MW", file=f)
    print('Transmission Efficiency :' +"{:.3f}".format(eff), " %", file=f)
    print('Sending End Circle :\t', file=f)
    print('a. Centre :\t', centre_1, file=f)
    print('b. Radius :\t' +"{:.3f}".format(rad_1), file=f)
    print('Receiving End Circle :\t', file=f)
    print('a. Centre :\t', centre_2, file=f)
    print('b. Radius :\t' +"{:.3f}".format(rad_2), file=f)






#Calling functions
parameters=user_input()
features=calculations(parameters)
output(features)
