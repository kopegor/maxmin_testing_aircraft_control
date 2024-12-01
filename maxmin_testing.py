import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import linalg
import itertools
from tqdm import tqdm
from scipy.integrate import solve_ivp


# Як - 52
# Cy0 = 4.81
# Cya = 4.81
# Cx0 = 0.0375
# B = 0.062
# S = 15
# M = 1000
# b = 1.64

# Як - 55
Cy0 = 0
Cya = 4.3
Cx0 = 0.035
B = 0.07

rho = 1.2
S = 14.805
M = 1000
b = 1.746


# # Jz = 1550  #!!!???
# # mza = 600  #!!!???
# # mzs = 900  #!!!???

Jz = 5000
mza = 0.005
mzs = 0.003

g = 9.8


# Определение параметров самолета и условий полета
M = 1000  # Масса ЛА (кг)
S = 14.805  # Площадь поверхности крыла (кв. м)
b = 1.746  # Расстояние от ц.м. до центра давления ЛА (м)
g = 9.81  # Ускорение свободного падения (м/с^2)
Jz = 5000  # Момент инерции корпуса ЛА относительно оси z (кг*м^2)
rho = 1.225  # Плотность воздуха (кг/м^3)

# Определение аэродинамических коэффициентов
c0x = 0.035
c0y = 0
c_alpha_y = 4.3
B = 0.07
m_alpha_z = 0.3
m_sigma_z = 0.7

# Определение времени моделирования и шага интегрирования
t_start = 0
t_end = 100
dt = 0.01
t_span = (t_start, t_end)

# Определение начальных условий
V0 = 50  # Начальная скорость (м/с)
theta0 = 0  # Начальный траекторный угол (рад)
phi0 = 0  # Начальный угол тангажа (рад)
Omega0 = 0  # Начальная абсолютная угловая скорость (рад/с)
initial_state = [V0, theta0, phi0, Omega0]

# Определение времени моделирования и шага интегрирования
t_start = 0
t_end = 30
dt = 0.01
t_span = (t_start, t_end)

# Определение функций изменения параметров со временем
def delta_V(t):
    # изменение скорости со временем
    return 0.8 * np.sin(0.1 * t)

def delta_alpha(t):
    # изменение угла атаки со временем
    return 0.15 * np.sin(0.2 * t)

def sigma(t):
    # изменение отклонения руля высоты со временем
    return 0.1 * np.sin(0.5 * t)


# Определение изменения силы тяги со временем
def P(t):
    # сила тяги линейно 
    return 0.7 * np.sin(0.7 * t)


# Определение функции правых частей системы дифференциальных уравнений
def equations_of_motion(t, state):
    V, theta, phi, Omega = state
    alpha = phi - theta  # Угол атаки
    delta_alpha_val = delta_alpha(t)
    delta_V_val = delta_V(t)
    sigma_val = sigma(t)
    X = 0.5 * rho * (V + delta_V_val)**2 * S * (c0x + B * (c0y + c_alpha_y * (alpha + delta_alpha_val))**2)
    Y = 0.5 * rho * (V + delta_V_val)**2 * S * (c0y + c_alpha_y * (alpha + delta_alpha_val))
    Mz = -0.5 * rho * (V + delta_V_val)**2 * S * b * (m_alpha_z * (alpha + delta_alpha_val) + m_sigma_z * sigma_val)
    
    dVdt = (-M * g * np.sin(theta) + P(t) * np.cos(alpha) - X) / M
    dthetadt = (P(t) * np.sin(alpha) + Y) / (M * V)
    dphidt = Omega
    dOmegadt = -Mz / Jz
    
    return [dVdt, dthetadt, dphidt, dOmegadt]


# Решение системы дифференциальных уравнений
solution = solve_ivp(equations_of_motion, t_span, initial_state, t_eval=np.arange(t_start, t_end, dt))

# Визуализация результатов
t_values = solution.t
V_values, theta_values, phi_values, Omega_values = solution.y

plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(t_values, V_values)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Velocity vs Time')

plt.subplot(2, 2, 2)
plt.plot(t_values, theta_values)
plt.xlabel('Time (s)')
plt.ylabel('Trajectory angle (rad)')
plt.title('Trajectory angle vs Time')

plt.subplot(2, 2, 3)
plt.plot(t_values, phi_values)
plt.xlabel('Time (s)')
plt.ylabel('Pitch angle (rad)')
plt.title('Pitch angle vs Time')

plt.subplot(2, 2, 4)
plt.plot(t_values, Omega_values)
plt.xlabel('Time (s)')
plt.ylabel('Angular velocity (rad/s)')
plt.title('Angular velocity vs Time')

plt.tight_layout()
plt.show()



# Определение функции правых частей линеаризованной системы
def linearized_equations_of_motion(t, state):
    return matrix_A @ state + matrix_B @ np.array([Delta_P(t), Delta_sigma(t)]) + matrix_C @ np.array([delta_V(t), delta_alpha(t)])

# Определение изменения управления со временем
def Delta_P(t):
    # изменение силы тяги со временем
    return 0.1 * np.sin(0.1 * t)

def Delta_sigma(t):
    # изменение отклонения руля высоты со временем
    return 0.05 * np.sin(0.5 * t)

# Определение изменения возмущений со временем
def delta_V(t):
    #  изменение скорости со временем
    return 0.01 * np.sin(0.1 * t)

def delta_alpha(t):
    #  изменение угла атаки со временем
    return 0.01 * np.sin(0.2 * t)

# Заполнение матриц коэффициентов
P_star = 100.0  # Пример: сила тяги в условиях стандартной атмосферы
V_star = 100.0  # Пример: скорость в условиях стандартной атмосферы
alpha_star = 0.1  # Пример: угол атаки в условиях стандартной атмосферы
theta_star = 0.2  # Пример: траекторный угол в условиях стандартной атмосферы

# Определение параметров самолета и условий полета
M = 1000  # Масса ЛА (кг)
S = 14.805  # Площадь поверхности крыла (кв. м)
b = 1.746  # Расстояние от ц.м. до центра давления ЛА (м)
g = 9.81  # Ускорение свободного падения (м/с^2)
Jz = 5000  # Момент инерции корпуса ЛА относительно оси z (кг*м^2)
# Jz = 250
rho = 1.225  # Плотность воздуха (кг/м^3)

# Определение аэродинамических коэффициентов
c0x = 0.035
c0y = 0
c_alpha_y = 4.3
B = 0.07
m_alpha_z = 0.3
m_sigma_z = 0.7



#================================
V_prog = 150 # 300
alpha_prog = -5
theta_prog = -15
P_prog = 0
sigma_prog = 0


a11 = -rho*S*V_prog*(Cx0 + B*(Cy0**2) +2*B*Cy0*Cya*alpha_prog + B*(Cya**2)*(alpha_prog**2)) / M
a12 = (P_prog*np.cos(alpha_prog) + rho*B*S*Cya*(V_prog**2)*(Cy0 + Cya*alpha_prog)) / M
a13 = -a12

a21 = (rho*S*V_prog*(Cy0 + Cya*alpha_prog) + M*g*np.cos(theta_prog) / (V_prog**2))
a22 = -(rho*S*Cya*(V_prog**2) + P_prog*np.cos(alpha_prog) / V_prog) / M
a23 = -a22

a41 = -rho*b*S*V_prog*(mza*alpha_prog + mzs*sigma_prog) / Jz
a42 = 0.5*rho*S*b*mza*(V_prog**2) / Jz
a43 = -a42

matrix_A = np.array([[a11, a12, a13, 0.],
              [a21, a22, a23, 0.],
              [0., 0., 0., 1.],
              [a41, a42, a43, 0.]])

b11 = np.cos(alpha_prog) / M
b21 = np.sin(alpha_prog) / (M * V_prog)
b42 = -0.5*rho*S*mzs*(V_prog**2) / Jz


matrix_B = np.array([[b11, 0.],
              [b21, 0.],
              [0., 0.],
              [0., b42]])



c11 = rho*S*V_prog*(Cx0 + B*(Cy0**2) + 2*B*Cy0*Cya*alpha_prog + B*(Cya**2)*(alpha_prog**2)) / M
    
c12 = rho*B*S*(V_prog**2)*Cya*(Cy0 + Cya*alpha_prog) / M

c21 = -rho*S*V_prog*(Cy0 + Cya*alpha_prog) / M

c22 = -0.5*rho*S*Cya*(V_prog**2) / M

c41 = -rho*b*S*V_prog*(mza*alpha_prog + mzs*sigma_prog)

c42 = -0.5*rho*b*S*mza*(V_prog**2) / Jz

matrix_C = np.array([[c11, c12],
              [c21, c22],
              [0., 0.],
              [c41, c42]])



# Определение начальных условий
# x0 = np.array([0, 0, 0, 0])  # Начальное состояние системы в отклонениях
x0 = np.random.random(4)


# Определение времени моделирования
t_start = 0
t_end = 10
dt = 0.01
t_span = (t_start, t_end)

# Решение системы дифференциальных уравнений
solution_linearized = solve_ivp(linearized_equations_of_motion, t_span, x0, t_eval=np.arange(t_start, t_end, dt))

# Визуализация результатов
t_values_linearized = solution_linearized.t
x_values, theta_values, phi_values, Omega_values = solution_linearized.y

plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(t_values_linearized, x_values)
plt.xlabel('Time (s)')
plt.ylabel('Velocity Deviation (m/s)')
plt.title('Velocity Deviation vs Time')

# plt.ylim(-1000, 1000)



plt.subplot(2, 2, 2)
plt.plot(t_values_linearized, theta_values)
plt.xlabel('Time (s)')
plt.ylabel('Theta Deviation (rad)')
plt.title('Theta Deviation vs Time')

# plt.ylim(-1000, 1000)


plt.subplot(2, 2, 3)
plt.plot(t_values_linearized, phi_values)
plt.xlabel('Time (s)')
plt.ylabel('Phi Deviation (rad)')
plt.title('Phi Deviation vs Time')

# plt.ylim(-1000, 1000)
# plt.xlim(0, 0.01)

plt.subplot(2, 2, 4)
plt.plot(t_values_linearized, Omega_values)
plt.xlabel('Time (s)')
plt.ylabel('Omega Deviation (rad/s)')
plt.title('Omega Deviation vs Time')

# plt.ylim(-1000, 1000)


plt.tight_layout()
plt.show()

from scipy import linalg
import itertools
from tqdm import tqdm

np.random.seed(42)
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
        
        # res = linalg.expm(matrix*t)
        # res = np.eye(matrix.shape[0]) + t*matrix + 0.5*(t**2)*(matrix@matrix) + (t**3)*(matrix@matrix)@matrix/6
        
        res = np.eye(matrix.shape[0])
        factorial = 1
        matrix_power_i = matrix
        
        for i in range(1, 7):
            res += (t**i) * matrix_power_i / factorial
            matrix_power_i = matrix_power_i @ matrix
            factorial *= i+1
            
        return res
    
    def func(self, t):
        exp_mat = self.exponential_matrix(self.A.T, self.tk - t)
        
        return np.dot(exp_mat, self.c)
    
    
    
class Control:
    def __init__(self, Psi, power_min, power_max, sigma_min, sigma_max, matrix_B):
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
    def __init__(self, Psi, v1_min, v1_max, v2_min, v2_max, matrix_C):
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
                                           matrix_C=matrix_C)
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


tk = 20
z0 = np.array([0, 0, 0, 0])
y0 = np.array([5, -3, -7, 4])#*100

grid_c = np.array([x / np.linalg.norm(x) for x in itertools.product(np.linspace(-1, 1, 4),
                                                                    np.linspace(-1, 1, 4),
                                                                    np.linspace(-1, 1, 4),
                                                                    np.linspace(-1, 1, 4))])    

list_diff_disturbance = []
list_diff_trajectories_disturbance = []
list_diff_vector_c_disturbance = []
disturbance_domain = []

for vector_c in tqdm(grid_c):
    function_psi = Psi(matrix_A,
                   vector_c, 
                   tk
                  )

    function_disturbance = Disturbance(function_psi, v1_min= -10, v1_max=10,
                                       v2_min= -3, v2_max=3,
                                       matrix_C=matrix_C)
    values_disturbance = [function_disturbance.get_disturbance(t_cur).tolist() for t_cur in np.linspace(0, tk, int(tk*10))]

    if values_disturbance not in list_diff_disturbance:
        list_diff_disturbance.append(values_disturbance)
        list_diff_vector_c_disturbance.append(vector_c)
    
    
    y = integrate_system_disturbance(matrix_A, matrix_C, function_disturbance, y0, tk)
    disturbance_domain.append(y[-1])
    y = y.tolist()
    if y not in list_diff_trajectories_disturbance:
        list_diff_trajectories_disturbance.append(y)
    

disturbance_domain = np.array(disturbance_domain)
list_diff_disturbance = np.array(list_diff_disturbance)
list_diff_trajectories_disturbance = np.array(list_diff_trajectories_disturbance)
list_diff_vector_c_disturbance = np.array(list_diff_vector_c_disturbance)
# print(len(list_diff_disturbance))


list_diff_control = []
list_diff_trajectories_control = []
list_diff_vector_c_control = []
control_domain = []

for vector_c in tqdm(grid_c):
    function_psi = Psi(matrix_A,
                   vector_c, 
                   tk
                  )

    function_control = Control(function_psi, power_min=0, power_max=15, sigma_min= -5, sigma_max=5,
                                       matrix_B=matrix_B)
    values_control = [function_control.get_control(t_cur).tolist() for t_cur in np.linspace(0, tk, int(tk*10))]

    if values_control not in list_diff_control:
        list_diff_control.append(values_control)
        list_diff_vector_c_control.append(vector_c)
    
    
    z = integrate_system_control(matrix_A, matrix_B, function_control, z0, tk)
    control_domain.append(z[-1])
    z = z.tolist()
    if z not in list_diff_trajectories_control:
        list_diff_trajectories_control.append(z)
    

control_domain = np.array(control_domain)
list_diff_control = np.array(list_diff_control)
list_diff_trajectories_control = np.array(list_diff_trajectories_control)
list_diff_vector_c_control = np.array(list_diff_vector_c_control)
# print(len(list_diff_control))

fig = plt.figure(figsize=(12, 8))

fig.add_subplot(3, 3, 1)
plt.scatter(control_domain[:, 0], control_domain[:, 1], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 1], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()

fig.add_subplot(3, 3, 2)
plt.scatter(control_domain[:, 0], control_domain[:, 2], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 2], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()

fig.add_subplot(3, 3, 3)
plt.scatter(control_domain[:, 0], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()

fig.add_subplot(3, 3, 5)
plt.scatter(control_domain[:, 1], control_domain[:, 2], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 1], disturbance_domain[:, 2], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()

fig.add_subplot(3, 3, 6)
plt.scatter(control_domain[:, 1], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 1], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()

fig.add_subplot(3, 3, 9)
plt.scatter(control_domain[:, 2], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 2], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")
plt.grid()
plt.legend()
plt.show()


# поиск максимина и минимакса для определения наличия седловой точки
minmax, saddle_point_minmax = find_minmax(control_domain, disturbance_domain)
maxmin, saddle_point_maxmin = find_maxmin(control_domain, disturbance_domain)


fig = plt.figure(figsize=(12, 8))

fig.add_subplot(3, 3, 1)
plt.scatter(control_domain[:, 0], control_domain[:, 1], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 1], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][0], saddle_point_minmax[0][1], color="red")
plt.scatter(saddle_point_minmax[1][0], saddle_point_minmax[1][1], color="red")
plt.plot([saddle_point_minmax[0][0], saddle_point_minmax[1][0]],
         [saddle_point_minmax[0][1], saddle_point_minmax[1][1]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][0], saddle_point_maxmin[0][1], color="blue")
plt.scatter(saddle_point_maxmin[1][0], saddle_point_maxmin[1][1], color="blue")
plt.plot([saddle_point_maxmin[0][0], saddle_point_maxmin[1][0]],
         [saddle_point_maxmin[0][1], saddle_point_maxmin[1][1]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()



fig.add_subplot(3, 3, 2)
plt.scatter(control_domain[:, 0], control_domain[:, 2], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 2], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][0], saddle_point_minmax[0][2], color="red")
plt.scatter(saddle_point_minmax[1][0], saddle_point_minmax[1][2], color="red")
plt.plot([saddle_point_minmax[0][0], saddle_point_minmax[1][0]],
         [saddle_point_minmax[0][2], saddle_point_minmax[1][2]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][0], saddle_point_maxmin[0][2], color="blue")
plt.scatter(saddle_point_maxmin[1][0], saddle_point_maxmin[1][2], color="blue")
plt.plot([saddle_point_maxmin[0][0], saddle_point_maxmin[1][0]],
         [saddle_point_maxmin[0][2], saddle_point_maxmin[1][2]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()




fig.add_subplot(3, 3, 3)
plt.scatter(control_domain[:, 0], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 0], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][0], saddle_point_minmax[0][3], color="red")
plt.scatter(saddle_point_minmax[1][0], saddle_point_minmax[1][3], color="red")
plt.plot([saddle_point_minmax[0][0], saddle_point_minmax[1][0]],
         [saddle_point_minmax[0][3], saddle_point_minmax[1][3]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][0], saddle_point_maxmin[0][3], color="blue")
plt.scatter(saddle_point_maxmin[1][0], saddle_point_maxmin[1][3], color="blue")
plt.plot([saddle_point_maxmin[0][0], saddle_point_maxmin[1][0]],
         [saddle_point_maxmin[0][3], saddle_point_maxmin[1][3]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()



fig.add_subplot(3, 3, 5)
plt.scatter(control_domain[:, 1], control_domain[:, 2], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 1], disturbance_domain[:, 2], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][1], saddle_point_minmax[0][2], color="red")
plt.scatter(saddle_point_minmax[1][1], saddle_point_minmax[1][2], color="red")
plt.plot([saddle_point_minmax[0][1], saddle_point_minmax[1][1]],
         [saddle_point_minmax[0][2], saddle_point_minmax[1][2]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][1], saddle_point_maxmin[0][2], color="blue")
plt.scatter(saddle_point_maxmin[1][1], saddle_point_maxmin[1][2], color="blue")
plt.plot([saddle_point_maxmin[0][1], saddle_point_maxmin[1][1]],
         [saddle_point_maxmin[0][2], saddle_point_maxmin[1][2]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()




fig.add_subplot(3, 3, 6)
plt.scatter(control_domain[:, 1], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 1], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][1], saddle_point_minmax[0][3], color="red")
plt.scatter(saddle_point_minmax[1][1], saddle_point_minmax[1][3], color="red")
plt.plot([saddle_point_minmax[0][1], saddle_point_minmax[1][1]],
         [saddle_point_minmax[0][3], saddle_point_minmax[1][3]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][1], saddle_point_maxmin[0][3], color="blue")
plt.scatter(saddle_point_maxmin[1][1], saddle_point_maxmin[1][3], color="blue")
plt.plot([saddle_point_maxmin[0][1], saddle_point_maxmin[1][1]],
         [saddle_point_maxmin[0][3], saddle_point_maxmin[1][3]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()



fig.add_subplot(3, 3, 9)
plt.scatter(control_domain[:, 2], control_domain[:, 3], label="control domain", color="green", marker="x")
plt.scatter(disturbance_domain[:, 2], disturbance_domain[:, 3], label="disturbance domain", color="orange", marker="o")

plt.scatter(saddle_point_minmax[0][2], saddle_point_minmax[0][3], color="red")
plt.scatter(saddle_point_minmax[1][2], saddle_point_minmax[1][3], color="red")
plt.plot([saddle_point_minmax[0][2], saddle_point_minmax[1][3]],
         [saddle_point_minmax[0][3], saddle_point_minmax[1][3]],
         color="red", label="minmax", linewidth=2)

plt.scatter(saddle_point_maxmin[0][2], saddle_point_maxmin[0][3], color="blue")
plt.scatter(saddle_point_maxmin[1][2], saddle_point_maxmin[1][3], color="blue")
plt.plot([saddle_point_maxmin[0][2], saddle_point_maxmin[1][2]],
         [saddle_point_maxmin[0][3], saddle_point_maxmin[1][3]],
         color="blue", label="maxmin", linewidth=2)

plt.grid()
plt.legend()
plt.show()


