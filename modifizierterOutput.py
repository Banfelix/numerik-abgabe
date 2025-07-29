"""
###############################################################################
                          FV-Solver Shallow Water
###############################################################################
"""
# -------------------------------------------------------------------------
# DO NOT TOUCH START !!!!!!! 
# -------------------------------------------------------------------------
import numpy as np
from numpy import *
import testParams
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
    testcase = 2

    # constant at 0m     = 1
    # bump ground        = 2
    ground   = testParams.ground 

    # homogenes system   = 1
    # inhomogenes system = 2
    system   = testParams.system

    # Ende der Simulation:-----------------------------------------------------
    t_end    = testParams.t_end

    #Hinweis: Startzeit t=0

    # CFL-Bedingung:----------------------------------------------------------- 
    cfl      = testParams.cfl

    # Gravitationskonstante:---------------------------------------------------
    g        = testParams.g

    
    # Gitter:------------------------------------------------------------------ 
    N        = testParams.N
    x_0      = testParams.x_0       
    x_end    = testParams.x_end

    # -------------------------------------------------------------------------
    # DO NOT TOUCH ENDE !!!!!!! 
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Initialisierung Gitter: 
    # -------------------------------------------------------------------------

    #####################################
    delta_x = (x_end - x_0) / (N+2)                 
    dx = delta_x       
    x_i_einhalb = linspace(-1, 1, (N+3))                              # 103 Wände   # 101 + 2 Ghost     # Wand bei -1 und 1 
    x_i = linspace((-1+ 0.5 * delta_x), (1-0.5 * delta_x), (N+2))     # 102 Zellen  # 100 + 2 Ghost

    #####################################

    # -------------------------------------------------------------------------
    # Initialisierung vom Testcase 
    # -------------------------------------------------------------------------

    # Allokieren:--------------------------------------------------------------

    H_loc  = np.zeros(N+2)
    vx_loc = np.zeros(N+2)
    bx_loc = np.zeros(N+2)
    
    # primitive Variablen
    w = np.zeros((3, N+2))
    # konservative Variablen
    u = np.zeros((2, N+2))

    # Setzen der Randbedingungen:----------------------------------------------

    if testcase == 1:    # Const Water Height   
        #####################################
        H_loc[:] = 0.5  
        vx_loc[:] = 0.0
        bx_loc[:] = 0.0

        #####################################
    elif testcase == 2:     # normal Dam Break
        #####################################
        for j in range(N+2):
            if x_i[j] <= 0:
                H_loc[j] = 1.2
            else:
                H_loc[j] = 4.0
        vx_loc[:] = 0.0
        bx_loc[:] = 0.0
        #####################################

    if ground == 1:       # const at 0m
        #####################################
        None
        #hier muss etwas implementiert werden
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
    w[0, :] = H_loc - bx_loc
    w[1, :] = vx_loc
    w[2, :] = bx_loc

    u[0, :] = H_loc
    u[1, :] = (H_loc - bx_loc) * vx_loc
    #####################################
    t=0.
    while t<t_end:

		# Bestimmung der CFL Zahl: ----------------------------------------------

        #####################################

        i = 0
        dt = 0.01
        a_i = abs(vx_loc[i]) + sqrt(g * (H_loc[i] - bx_loc[i]))
        if dt > (cfl * dx / a_i):
            dt = (cfl * dx / a_i)

        #####################################
        
        #Zeitoperator und Raumoperator: -----------------------------------------

        #####################################
        flux = np.zeros((2, N+2))
        g_flux = np.zeros((2, N+1))
        flux[0, :] = u[1, :]            # hu            # Physikalischer Fluss f(U)
        flux[1, :] = u[1, :] * w[1, :]  # hu^2
        # flux[1, :] = u[1, :] * w[1, :] + g * u[0, :] * w[0, :]  # hu^2 + ghH 

        for j in range(N+1):            # Numerischer Fluss g_{j+1/2}
            fL = flux[:, j]
            fR = flux[:, j+1]

            hL = w[0, j]
            hR = w[0, j+1]
            uL = w[1, j]
            uR = w[1, j+1]

            val1 = abs(uL) + sqrt(g * hL)
            val2 = abs(uR) + sqrt(g * hR)
            a = val1 if val1 > val2 else val2 
            
            g_flux[:, j] = 0.5 * (fL + fR) - 0.5 * a * (u[:, j+1] - u[:, j])
        
        #####################################

        #Update der Variablen: -------------------------------------------------

        #####################################
        for j in range(1, N+1):     
            u[:, j] = u[:, j] - (dt / dx) * (g_flux[:, j] - g_flux[:, j-1])         # i= 0 und i = N+1 sind die Ghost Zellen  #
        u[:, 0] = u[:, 1]       # Ghost Zellen
        u[:, N+1] = u[:, N]
        
        w[0, :] = u[0, :] - bx_loc
        w[1, :] = u[1, :] / w[0, :]
        w[2, :] = bx_loc

        #####################################

        # Update Zeit: ----------------------------------------------------------
        t=t+dt

    # ---------------------------------------------------------------------------
    # DO NOT TOUCH START !!!!!!!
    # ---------------------------------------------------------------------------
    print('OUTPUT: DO NOT TOUCH!')
    H_output = u[0, int(round((testParams.x_output - (x_0 - dx/2)) / dx))]                  # Indexe umgedreht bei H output
    u_output = w[1, int(round((testParams.x_output - (x_0 - dx/2)) / dx))]
    print('======================================================================')
    print('Finaler Zeitschritt                     : dt =',dt)
    print('Gitterweite                             : dx =',dx)
    print('Absolute Wasserhöhe         bei x=',testParams.x_output,': H  =',H_output)
    print('Horizontale Geschwindigkeit bei x=',testParams.x_output,': u  =',u_output)
    print('======================================================================')

    #np.savetxt(__file__[:-3]+'.csv',np.array(dt,dx,H_output,u_output),delimiter=',')    erstmal auskommentiert
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
