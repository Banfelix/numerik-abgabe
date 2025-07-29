import numpy as np
from numpy import *
import params


def MAIN():
    # -------------------------------------------------------------------------
    # DO NOT TOUCH START !!!!!!! 
    # -------------------------------------------------------------------------

    # Testcase:----------------------------------------------------------------
    # const water height = 1
    # normal dam  Break  = 2
    testcase = params.testcase

    # constant at 0m     = 1
    # bump ground        = 2
    ground   = params.ground 

    # homogenes system   = 1
    # inhomogenes system = 2
    system   = params.system

    # Ende der Simulation:-----------------------------------------------------
    t_end    = params.t_end

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
    delta_x = (x_end - x_0) / 102                         
    x_i_einhalb = linspace(-1, 1, 103)                              # 103 Wände
    x_i = linspace((-1+ 0.5 * delta_x), (1-0.5 * delta_x), 102)     # 102 Zellen
    #####################################   

    # -------------------------------------------------------------------------
    # Initialisierung vom Testcase 
    # -------------------------------------------------------------------------

    # Allokieren:--------------------------------------------------------------
    H_loc  = [0] * 102
    print(H_loc)
    print(x_i_einhalb)
    print(len(x_i_einhalb))
    print(x_i)
    print(len(x_i))
    print(x_i[0])
    print((x_i_einhalb[0]+ x_i_einhalb[1])/2)

# -------------------------------------------------------------------------------
# DO NOT TOUCH START !!!!!!!
# -------------------------------------------------------------------------------
if __name__ == '__main__':
    MAIN()
# -------------------------------------------------------------------------------
# DO NOT TOUCH ENDE !!!!!!!
# -------------------------------------------------------------------------------


        flux = np.zeros((2, N+2))
        g_flux = np.zeros((2, N+1))
        flux[0, :] = u[1, :]            # hu            # Physikalischer Fluss f(U)
        flux[1, :] = u[1, :] * w[1, :]  # hu^2

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