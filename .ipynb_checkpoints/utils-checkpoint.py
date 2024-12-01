import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import linalg
from tqdm import tqdm
import itertools



def integrate_system_control(A, B, u, z0, tk):
    T = np.linspace(0, tk, int(tk*10))
    
    def right_side(z, t, A, B, u):
        return A@z - B@u.get_control(t)
        # return A@z - B*u.get_control(t)
    
    z = integrate.odeint(right_side, z0, T, args=(A, B, u))
    
    return z


def integrate_system_disturbance(A, C, v, y0, tk):
    T = np.linspace(0, tk, int(tk*10))
    
    def right_side(y, t, A, C, v):
        return A@y + C@v.get_disturbance(t)
        # return A@y + C*v.get_disturbance(t)
    
    y = integrate.odeint(right_side, y0, T, args=(A, C, v))
    
    return y


class Psi:
    def __init__(self, A, c,
                 tk
                ):
        self.A = A
        self.c = c
        self.tk = tk
    
    def exponential_matrix(self, matrix, t):
        
        res = linalg.expm(matrix*t)

        return res
    
    def func(self, t):
        exp_mat = self.exponential_matrix(self.A.T, self.tk - t)
        
        return np.dot(exp_mat, self.c)

    
class Control:
    def __init__(self, Psi, power_min, power_max, sigma_min,     sigma_max, matrix_B):
        self.Psi = Psi
        self.power_min = power_min
        self.power_max = power_max
        self.sigma_min = sigma_min
        self.sigma_max = sigma_max
        self.matrix_B = matrix_B
      
    
    def get_control(self, t):
        psi_t = self.Psi.func(t) 
        
        if psi_t@self.matrix_B[:, 0] < 0:
            power = self.power_max
        else:
            power = self.power_min
        
        if psi_t@self.matrix_B[:, 1] < 0:
            sigma = self.sigma_max
        else:
            sigma = self.sigma_min
            
        return np.array([power, sigma])

    
class Disturbance:
    def __init__(self, Psi, v1_min, v1_max, v2_min, v2_max, martirx_C):
        self.Psi = Psi
        self.v1_min = v1_min
        self.v1_max = v1_max
        self.v2_min = v2_min
        self.v2_max = v2_max
        self.matrix_C = matrix_C

    
    def get_disturbance(self, t):
        psi = self.Psi.func(t)
        
        if psi@matrix_C[:, 0] < 0:
            v1 = self.v1_max
        else:
            v1 = self.v1_min
        
        if psi@matrix_C[:, 1] < 0:
            v2 = self.v2_max
        else:
            v2 = self.v2_min
            
        return np.array([v1, v2])

    
def get_attainability_domains(tk, z0, y0, grid_c, matrix_A, matrix_B, matrix_C,
                              u1_min, u1_max, u2_min, u2_max, v1_min, v1_max, v2_min, v2_max):
    
    control_domain = []

    for vector_c in tqdm(grid_c):
        function_psi = Psi(matrix_A,
                       vector_c, 
                       tk
                      )

        function_control = Control(function_psi, power_min=u1_min, power_max=u1_max,
                                   sigma_min= u2_min, sigma_max=u2_max,
                                   matrix_B = matrix_B)
        z = integrate_system_control(matrix_A, matrix_B, function_control, z0, tk)
        control_domain.append(z[-1])

    control_domain = np.array(control_domain)


    disturbance_domain = []

    for vector_c in tqdm(grid_c):
        function_psi = Psi(matrix_A,
                       vector_c, 
                       tk
                      )

        function_disturbance = Disturbance(function_psi, v1_min= v1_min, v1_max=v1_max,
                                           v2_min= v2_min, v2_max=v2_max,
                                           martirx_C=matrix_C)
        y = integrate_system_disturbance(matrix_A, matrix_C, function_disturbance, y0, tk)
        disturbance_domain.append(y[-1])

    disturbance_domain = np.array(disturbance_domain)
    
    return control_domain, disturbance_domain


def find_minmax(control_domain, disturbance_domain):
    # поиск минимакса 

    list_y_max_distance = []
    list_max_distance = []

    for z in control_domain:

        max_distance = -np.inf
        y_max_distance = disturbance_domain[0]

        for y in disturbance_domain:
            distance = np.linalg.norm(y - z)

            if distance > max_distance:
                max_distance = distance
                y_max_distance = y

        list_y_max_distance.append(y_max_distance)
        list_max_distance.append(max_distance)

    idx_minmax = np.argmin(list_max_distance)

    minmax = list_max_distance[idx_minmax]

    saddle_point_minmax = np.array([control_domain[idx_minmax], list_y_max_distance[idx_minmax]])
    
    return minmax, saddle_point_minmax


def find_maxmin(control_domain, disturbance_domain):
    # поиск максимина

    list_z_min_distance = []
    list_min_distance = []

    for y in disturbance_domain:

        min_distance = np.inf
        z_min_distance = control_domain[0]

        for z in control_domain:
            distance = np.linalg.norm(y - z)

            if distance < min_distance:
                min_distance = distance
                z_min_distance = z

        list_z_min_distance.append(z_min_distance)
        list_min_distance.append(min_distance)

    idx_maxmin = np.argmax(list_min_distance)

    maxmin = list_min_distance[idx_maxmin]

    saddle_point_maxmin = np.array([list_z_min_distance[idx_maxmin], disturbance_domain[idx_maxmin]])
    
    return maxmin, saddle_point_maxmin





