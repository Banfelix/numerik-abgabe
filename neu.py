"""
###############################################################################
                          FV-Solver Shallow Water
###############################################################################
"""
# -------------------------------------------------------------------------
# DO NOT TOUCH START !!!!!!! 
# -------------------------------------------------------------------------
import numpy as np
#from numpy import *
import params
import matplotlib.pyplot as plt
# -------------------------------------------------------------------------
# DO NOT TOUCH ENDE !!!!!!! 
# -------------------------------------------------------------------------


def MAIN():
    # -------------------------------------------------------------------------
    # DO NOT TOUCH START !!!!!!! 
    # -------------------------------------------------------------------------

    # Testcase:----------------------------------------------------------------
    # const water height = 1
    # normal dam  Break  = 2
    #testcase = params.testcase
    testcase = 2

    # constant at 0m     = 1
    # bump ground        = 2
    ground   = params.ground 

    # homogenes system   = 1
    # inhomogenes system = 2
    system   = params.system

    # Ende der Simulation:-----------------------------------------------------
    t_end    = 1

    #Hinweis: Startzeit t=0

    # CFL-Bedingung:----------------------------------------------------------- 
    cfl      = params.cfl

    # Gravitationskonstante:---------------------------------------------------
    g        = params.g
    
    # Gitter:------------------------------------------------------------------ 
    N        = params.N
    x_0      = params.x_0       
    x_end    = params.x_end

    # -------------------------------------------------------------------------
    # DO NOT TOUCH ENDE !!!!!!! 
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Initialisierung Gitter: 
    # -------------------------------------------------------------------------

    #####################################
    dx = (x_end - x_0) / N
    #####################################

    # -------------------------------------------------------------------------
    # Initialisierung vom Testcase 
    # -------------------------------------------------------------------------

    # Allokieren:--------------------------------------------------------------

    H_loc  = np.zeros(N+2)
    vx_loc = np.zeros(N+2)
    bx_loc = np.zeros(N+2)
    
    # primitive Variablen
    w      = np.zeros((N+2 ,3))
    # konservative Variablen
    u      = np.zeros((N+2, 2))

    # Setzen der Randbedingungen:----------------------------------------------

    if testcase  ==1 :    # Const Water Height   
        #####################################
        H_loc[:] = 0.5
        #####################################
    elif testcase == 2:     # normal Dam Break
        #####################################
        for j in range(N+2):
            x_j = x_0 + (j - 0.5) * dx
            if x_j >= 0:
                H_loc[j] = 4.0
            else:
                H_loc[j] = 1.2
    print(H_loc)
        #####################################

    if ground == 1:       # const at 0m
        #####################################
        bx_loc[:] = 0
        #####################################
    elif ground == 2:       # bump ground
        #####################################
        None
        #hier muss etwas implementiert werden
        #####################################

    
	#homogen
    if system == 1:       # homogen
        #####################################
        None
        #hier muss etwas implementiert werden
        #####################################
    elif system == 2:       # inhomogen
        #####################################
        None
        #hier muss etwas implementiert werden
        #####################################

	# Befüllen von w und u:
    #####################################
    w[:, 0] = H_loc[:] - bx_loc[:]
    w[:, 1] = vx_loc[:]
    w[:, 2] = bx_loc[:]

    u[:, 0] = H_loc[:]
    u[:, 1] = (H_loc[:] - bx_loc[:]) * vx_loc[:]
    #####################################

    '''
    # Gitterpunkte berechnen (Zellzentren)
    x_vals = np.array([x_0 + (i - 0.5) * dx for i in range(N+2)])
    # Figure & Achsen vorbereiten
    plt.ion()  # interaktiver Modus an
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    line_H, = ax1.plot(x_vals, u[:, 0], 'b-', label='H(x)')
    line_u, = ax2.plot(x_vals, w[:, 1], 'r-', label='u(x)')
    ax1.set_ylabel('Wasserhöhe H')
    ax2.set_ylabel('Geschwindigkeit u')
    ax2.set_xlabel('x')
    ax1.grid(True)
    ax2.grid(True)
    ax1.legend()
    ax2.legend()
    '''



    f_flux = np.zeros((N+2, 2))   
    g_flux = np.zeros((N+2, 2))
    g_flux_rechts = np.zeros((N+2, 2))
    g_flux_links =np.zeros((N+2, 2))
    u_rechts = np.zeros((N+2, 2))
    u_links = np.zeros((N+2, 2))
    slope = np.zeros(2)
    
    f_flux_links = np.zeros((N+2, 2))
    f_flux_rechts = np.zeros((N+2, 2))

    t=0
    while t<t_end:

		# Bestimmung der CFL Zahl: ----------------------------------------------

        ##################################### 
        a = np.zeros(N+2)
        a[:] = np.abs(w[:, 1]) + np.sqrt(g * w[:, 0])
        a_max = np.max(a[1:-1])
        dt = cfl * dx / a_max

        #####################################
        
        #Zeitoperator und Raumoperator: -----------------------------------------

        #####################################
        
        # Berechnung von den Zellrandwerten wobei rechts/links [i] der Rechte und linke Wert
        # Der i-ten Zelle sind 
        u[0, :] = u[1, :]
        u[N+1] = u[N, :]
        delta_rechts = np.zeros((N+2, 2))
        delta_links = np.zeros((N+2, 2))
        for i in range(1, N+1): # von 1 bis N 
            for k in range(2):
                delta_links = u[i, k] - u [i-1, k]
                delta_rechts = u[i+1, k] - u [i, k]
            
                if delta_links * delta_rechts <= 0:
                    slope[k] = 0
                else:
                    slope[k] = delta_links if abs(delta_links) < abs(delta_rechts) else delta_rechts

            u_links[i, :] = u[i, :] - 0.5 * dx * slope[:]
            u_rechts[i, :] = u[i, :] + 0.5 * dx * slope[:]
        
        #Beide Ränder der linken Ghost Zelle
        u_rechts[0, :] = u_links[1, :]   
        u_links[0, :] = u_rechts[0, :]

        #Beide Ränder der rechten Ghost Zelle
        u_links[N+1, :] = u_rechts[N, :]
        u_rechts[N+1, :] = u_links[N+1, :]
    


        # Zeitrekonstruktion
        for i in range(1, N+1):
            
            # f(u) der ersten Gleichung 
            f_flux_links[i, 0] = u_links[i, 1]
            f_flux_rechts[i, 0] = u_rechts[i,1]

            # f(u) der zweiten Gleichung 
            f_flux_links[i, 1] = (u_links[i, 1]**2) / u_links[i, 0] # -b(x) fehlt, falls nicht = 0
            f_flux_rechts[i, 1] = (u_rechts[i, 1]**2) / u_rechts[i, 0]

            # Neue Randwerte mit halbem Zeitschritt (Rekonstruktion)
            u_links[i, :] = u_links[i, :] - (dt / (2 * dx)) *  f_flux_rechts[i, :] - f_flux_links[i, :]                 
            u_rechts[i, :] = u_rechts[i, :] - (dt / (2 * dx)) *  f_flux_rechts[i, :] - f_flux_links[i, :]               
        
        #Beide Ränder der linken Ghost Zelle
        u_rechts[0, :] = u_links[1, :]   
        u_links[0, :] = u_rechts[0, :]

        #Beide Ränder der rechten Ghost Zelle
        u_links[N+1, :] = u_rechts[N, :]
        u_rechts[N+1, :] = u_links[N+1, :]                                                                                                                 



        alpha = np.zeros(2)
        alpha = np.ones(2)
        print(u_links[:,0])
        for i in range(1, N+1):                                               #u_links[i+1] "Wert der von rechts kommt" ????
            alpha[0] = max(np.abs(u_rechts[i, 0] + np.sqrt(g*u_rechts[i, 0])) , np.abs(u_links[i+1, 0]) + np.sqrt(g*u_links[i+1, 0])) # -b(x) Falls Notwendig ToDo um von H auf h zu kommen (h = u[:, 0])
            alpha[1] = max(np.abs(u_rechts[i, 1] + np.sqrt(g*u_rechts[i, 0])) , np.abs(u_links[i+1, 1]) + np.sqrt(g*u_links[i+1, 0])) # -b(x) Falls Notwendig ToDo
            
            g_flux[i, 0] = 0.5 * (u_rechts[i, 1] + u_links[i+1, 1]) - 0.5 * alpha[0] * (u_links[i+1, 0] - u_rechts[i, 0])
            g_flux[i, 1] = 0.5 * ((u_rechts[i, 1]**2) / u_rechts[i, 0] + (u_links[i+1, 1]**2) / u_rechts[i, 0]) - 0.5 * alpha[0] * (u_links[i+1, 0] - u_rechts[i, 0]) # -b(x) Falls Notwendig ToDo
        
        for i in range(1, N+1):
            u[i, :] = u[i, :] - (dt/dx) * (g_flux[i, :] - g_flux[i-1, :])
        
        u[0, :] = u[1, :]
        u[N+1] = u[N, :]
        
        '''
            # --- Live-Plot aktualisieren ---
        if int(t/dt) % 10 == 0:  # alle ~10 Zeitschritte
            line_H.set_ydata(u[:, 0])
            line_u.set_ydata(w[:, 1])
            ax1.set_ylim(0, np.max(u[:, 0]) * 1.1)
            ax2.set_ylim(np.min(w[:, 1]) * 1.1, np.max(w[:, 1]) * 1.1)
            plt.pause(0.01)
        '''
        #####################################

        #Update der Variablen: -------------------------------------------------

        #####################################
        w[:, 0] = u[:, 0] - bx_loc[:]
        w[:, 1] = u[:, 1] / w[:, 0]
        w[:, 2] = bx_loc[:]
        #####################################

        # Update Zeit: ----------------------------------------------------------
        t=t+dt

    plt.ioff()
    plt.show()

    # ---------------------------------------------------------------------------
    # DO NOT TOUCH START !!!!!!!
    # ---------------------------------------------------------------------------
    print('OUTPUT: DO NOT TOUCH!')
    H_output = u[int(round((params.x_output - (x_0 - dx/2)) / dx)), 0]
    u_output = w[int(round((params.x_output - (x_0 - dx/2)) / dx)), 1]
    print('======================================================================')
    print('Finaler Zeitschritt                     : dt =',dt)
    print('Gitterweite                             : dx =',dx)
    print('Absolute Wasserhöhe         bei x=',params.x_output,': H  =',H_output)
    print('Horizontale Geschwindigkeit bei x=',params.x_output,': u  =',u_output)
    #print('Absolute Wasserhöhe         bei x=',params.x_output,': H  =',H_output)
    #print('Horizontale Geschwindigkeit bei x=',params.x_output,': u  =',u_output)
    print('======================================================================')

    #np.savetxt(__file__[:-3]+'.csv',np.array(dt,dx,H_output,u_output),delimiter=',')
    # ---------------------------------------------------------------------------
    # DO NOT TOUCH ENDE !!!!!!!
    # ---------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# DO NOT TOUCH START !!!!!!!
# -------------------------------------------------------------------------------
if __name__ == '__main__':
    MAIN()
# -------------------------------------------------------------------------------
# DO NOT TOUCH ENDE !!!!!!!
# -------------------------------------------------------------------------------

