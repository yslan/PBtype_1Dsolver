# PBtype_1Dsolver

## cases
Every running cases located in cases/
They should include .rea file and userf.m

rea collect every input parameter we use in code
userf is a user specified function f = userf(C_1,C_2)

 - Our code solve the following 1D equation
 
   - Eqs:
 
      Psi_1 = k_B*T * ln(C_1) + e*z_1*phi + f_1(C_1,C_2) ... (1)
   
      Psi_2 = k_B*T * ln(C_2) + e*z_2*phi + f_2(C_1,C_2) ... (2)
   
      \- (epsilon * Phi')' = e*(z_1*C_1+z_2*C_2) + e*Q    ... (3)
 
   - BC:
 
      C1, C2: Left noflux + Right Dirichlet
    
      Phi: both Dirichlet
