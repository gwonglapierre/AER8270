import numpy as np
import math as math
import matplotlib.pyplot as plt

def COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S):
    
    # Number of panels
    numPan = len(XC)                                                                # Number of panels/control points
    
    # Initialize arrays
    I = np.zeros([numPan,numPan])                                                   # Initialize I integral matrix
    J = np.zeros([numPan,numPan])                                                   # Initialize J integral matrix
    
    # Compute integral
    for i in range(numPan):                                                         # Loop over i panels
        for j in range(numPan):                                                     # Loop over j panels
            if (j != i):                                                            # If the i and j panels are not the same
                # Compute intermediate values
                A  = -(XC[i]-XB[j])*np.cos(phi[j])-(YC[i]-YB[j])*np.sin(phi[j])     # A term
                B  = (XC[i]-XB[j])**2 + (YC[i]-YB[j])**2                            # B term
                Cn = np.sin(phi[i]-phi[j])                                          # C term (normal)
                Dn = -(XC[i]-XB[j])*np.sin(phi[i])+(YC[i]-YB[j])*np.cos(phi[i])     # D term (normal)
                Ct = -np.cos(phi[i]-phi[j])                                         # C term (tangential)
                Dt = (XC[i]-XB[j])*np.cos(phi[i])+(YC[i]-YB[j])*np.sin(phi[i])      # D term (tangential)
                E  = np.sqrt(B-A**2)                                                # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):       # If E term is 0 or complex or a NAN or an INF
                    I[i,j] = 0                                                      # Set I value equal to zero
                    J[i,j] = 0                                                      # Set J value equal to zero
                else:
                    # Compute I (needed for normal velocity), Ref [1]
                    term1  = 0.5*Cn*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in I equation
                    term2  = ((Dn-A*Cn)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in I equation
                    I[i,j] = term1 + term2                                          # Compute I integral
                    
                    # Compute J (needed for tangential velocity), Ref [2]
                    term1  = 0.5*Ct*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in I equation
                    term2  = ((Dt-A*Ct)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in I equation
                    J[i,j] = term1 + term2                                          # Compute J integral
                
            # Zero out any problem values
            if (np.iscomplex(I[i,j]) or np.isnan(I[i,j]) or np.isinf(I[i,j])):      # If I term is complex or a NAN or an INF
                I[i,j] = 0                                                          # Set I value equal to zero
            if (np.iscomplex(J[i,j]) or np.isnan(J[i,j]) or np.isinf(J[i,j])):      # If J term is complex or a NAN or an INF
                J[i,j] = 0                                                          # Set J value equal to zero
    
    return I, J                                                                     # Return both I and J matrices

# User-defined knowns
Vinf = 1                                                                        # Freestream velocity
AoA  = 0                                                                        # Angle of attack [deg]
numB = 9                                                                        # Number of boundary points (including endpoint)
tO   = (360/(numB-1))/2                                                         # Boundary point angle offset [deg]
AoAR = AoA*(np.pi/180)                                                          # Convert AoA to radians [rad]

# Plotting flags
flagPlot = [1,      # Shape polygon with panel normal vectors
            1,      # Geometry boundary pts, control pts, first panel, second panel
            1,      # Analytical and SPM pressure coefficient plot
            1,      # Streamline plot
            1]      # Pressure coefficient contour plot

# %% CREATE CIRCLE BOUNDARY POINTS

# Angles used to compute boundary points
theta = np.linspace(0,360,numB)                                                 # Create angles for computing boundary point locations [deg]
theta = theta + tO                                                              # Add panel angle offset [deg]
theta = theta*(np.pi/180)                                                       # Convert from degrees to radians [rad]

# Boundary points
XB = np.cos(theta)                                                              # Compute boundary point X-coordinate [radius of 1]
YB = np.sin(theta)                                                              # Compute boundary point Y-coordinate [radius of 1]

# Number of panels
numPan = len(XB)-1                                                              # Number of panels (control points)

# %% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

# Check for direction of points
edge = np.zeros(numPan)                                                         # Initialize edge value array
for i in range(numPan):                                                         # Loop over all panels
    edge[i] = (XB[i+1]-XB[i])*(YB[i+1]+YB[i])                                   # Compute edge values

sumEdge = np.sum(edge)                                                          # Sum of all edge values

# If panels are CCW, flip them (don't if CW)
if (sumEdge < 0):                                                               # If panels are CCW
    XB = np.flipud(XB)                                                          # Flip the X-data array
    YB = np.flipud(YB)                                                          # Flip the Y-data array

# %% PANEL METHOD GEOMETRY - REF [1]

# Initialize variables
XC  = np.zeros(numPan)                                                          # Initialize control point X-coordinate
YC  = np.zeros(numPan)                                                          # Initialize control point Y-coordinate
S   = np.zeros(numPan)                                                          # Initialize panel length array
phi = np.zeros(numPan)                                                          # Initialize panel orientation angle array

# Find geometric quantities of the airfoil
for i in range(numPan):                                                         # Loop over all panels
    XC[i]   = 0.5*(XB[i]+XB[i+1])                                               # X-value of control point
    YC[i]   = 0.5*(YB[i]+YB[i+1])                                               # Y-value of control point
    dx      = XB[i+1]-XB[i]                                                     # Change in X between boundary points
    dy      = YB[i+1]-YB[i]                                                     # Change in Y between boundary points
    S[i]    = (dx**2 + dy**2)**0.5                                              # Length of the panel
    phi[i]  = math.atan2(dy,dx)                                                 # Angle of panel (positive X-axis to inside face)
    if (phi[i] < 0):                                                            # Make all panel angles positive [rad]
        phi[i] = phi[i] + 2*np.pi

# Compute angle of panel normal w.r.t. horizontal and include AoA
delta                = phi + (np.pi/2)                                          # Angle of panel normal [rad]
beta                 = delta - AoAR                                             # Angle of panel normal and AoA [rad]
beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi                           # Make all panel angles between 0 and 2pi [rad]

# %% COMPUTE SOURCE PANEL STRENGTHS - REF [5]

# Geometric integral (normal [I] and tangential [J])
# - Refs [2] and [3]
I, J = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S)                                        # Compute geometric integrals

# Populate A matrix
# - Simpler option: A = I + np.pi*np.eye(numPan,numPan)
A = np.zeros([numPan,numPan])                                                   # Initialize the A matrix
for i in range(numPan):                                                         # Loop over all i panels
    for j in range(numPan):                                                     # Loop over all j panels
        if (i == j):                                                            # If the panels are the same
            A[i,j] = np.pi                                                      # Set A equal to pi
        else:                                                                   # If panels are not the same
            A[i,j] = I[i,j]                                                     # Set A equal to geometric integral

# Populate b array
# - Simpler option: b = -Vinf*2*np.pi*np.cos(beta)
b = np.zeros(numPan)                                                            # Initialize the b array
for i in range(numPan):                                                         # Loop over all panels
    b[i] = -Vinf*2*np.pi*np.cos(beta[i])                                        # Compute RHS array

# Compute source panel strengths (lam) from system of equations
lam = np.linalg.solve(A,b) 
Adim=lam/(2*np.pi)                                                     # Compute all source strength values

# Check the sum of the source strengths
# - This should be very close to zero for a closed polygon
print("Sum of L: ",sum(lam*S))                                                  # Print sum of all source strengths

# %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

# Compute velocities
# - Simpler method: Vt = Vinf*np.sin(beta) + np.dot(J,lam)/(2*np.pi)
#                   Cp = 1 - (Vt/Vinf)**2
Vt = np.zeros(numPan)                                                           # Initialize tangential velocity array
Cp = np.zeros(numPan)                                                           # Initialize pressure coefficient array
for i in range(numPan):                                                         # Loop over all i panels
    addVal = 0                                                                  # Reset the summation value to zero
    for j in range(numPan):                                                     # Loop over all j panels
        addVal = addVal + (lam[j]/(2*np.pi))*J[i,j]                             # Sum all tangential source panel terms
    
    Vt[i] = Vinf*np.sin(beta[i]) + addVal                                       # Compute tangential velocity by adding uniform flow term
    Cp[i] = 1 - (Vt[i]/Vinf)**2                                                 # Compute pressure coefficient

# Analytical angles and pressure coefficients
analyticTheta = np.linspace(0,2*np.pi,200)                                      # Analytical theta angles [rad]
analyticCP    = 1 - 4*np.sin(analyticTheta)**2                                  # Analytical pressure coefficient []

# %% COMPUTE LIFT AND DRAG

# Compute normal and axial force coefficients
CN = -Cp*S*np.sin(beta)                                                         # Normal force coefficient []
CA = -Cp*S*np.cos(beta)                                                         # Axial force coefficient []

#%% Plotting

fig = plt.figure(1)                                                         # Create figure
plt.cla()                                                                   # Get ready for plotting
plt.plot(analyticTheta*(180/np.pi),analyticCP,'b-',label='Analytical')      # Plot analytical pressure coefficient
plt.plot(beta*(180/np.pi),Cp,'ks',markerfacecolor='r',label='SPM')          # Plot panel method pressure coefficient
plt.xlabel('Angle [deg]')                                                   # Set X-label
plt.ylabel('Pressure Coefficient [Cp]')                                          # Set Y-label
plt.title('Pressure Coefficient Comparison')                                # Set title
plt.xlim(0, 360)                                                            # Set X-limits
plt.ylim(-3.5, 1.5)                                                         # Set Y-limits
plt.legend()                                                                # Show legend
plt.show()                                                                  # Display plot

fig = plt.figure(2)                                                         # Create figure
plt.cla()                                                                   # Get ready for plotting
plt.plot(XB,YB,'k-',label='Panels')                                         # Plot polygon
plt.plot([XB[0], XB[1]],[YB[0], YB[1]],'b-',label='First Panel')            # Plot first panel
plt.plot([XB[1], XB[2]],[YB[1], YB[2]],'g-',label='Second Panel')           # Plot second panel
plt.plot(XB,YB,'ko',markerfacecolor='k',label='Boundary Points')            # Plot boundary points
plt.plot(XC,YC,'ko',markerfacecolor='r',label='Control Points')             # Plot control points
plt.xlabel('X-Axis')                                                        # Set X-label
plt.ylabel('Y-Axis')                                                        # Set Y-label
plt.axis('equal')                                                           # Set axes equal
plt.legend()                                                                # Show legend
plt.show()                   