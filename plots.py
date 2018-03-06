import numpy as np
import matplotlib.pyplot as plt
from Input_variables import *
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from internal_forces import *
from shear_force import *
from matplotlib import cm
'''
Generates plot of forces, moments, stresses, displacements, curvatures, slopes (Also in 3D)
'''
def get_simple_plot1(function, n_points, C1, C2, int_constants_z, forces_y):
    x_points = np.linspace(0, l_a, num=n_points)
    y_points = np.zeros((n_points))
    for i in range(np.size(x_points)):
        y_points[i] = function(x_points[i], C1, C2, int_constants_z, forces_y)

    plt.plot(x_points,y_points)

    plt.show()

def get_simple_plot(function, n_points, C1, C2, int_constants_z):
    x_points = np.linspace(0, l_a, num=n_points)
    y_points = np.zeros((n_points))
    for i in range(np.size(x_points)):
        y_points[i] = function(x_points[i], C1, C2, int_constants_z)

    plt.plot(x_points,y_points)

    plt.show()

def plot_3D(n_points):
    x_points = np.linspace(0, l_a, num = n_points)
    y_points = np.zeros((n_points))
    z_points = np.zeros((n_points))
    for i in range(np.size(x_points)):
        y_points[i] = get_displacement(x_points[i], C1, -C2, int_constants_z, forces_y)
        z_points[i] = get_displacement(x_points[i], -C3, C1, RF_y[3:], forces_y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(y_points,x_points, z_points)

    ax.set_xlabel('Aileron deflection in z axis [m]', fontsize=8)
    ax.set_ylabel('Span-wise location[m], from hinge 3 to 1', fontsize=8)
    ax.set_zlabel('Aileron deflection in y axis [m]', fontsize=8)
    plt.title("3D deflection of aileron's hinge")

    plt.show()



def von_mises(x_position):
    normal_stress = get_normal_stress(x_position)
    normal_stress_tot = np.vstack((normal_stress, normal_stress[0]))
    q_flow_origin = compute_q(x_position)
    q_flow = q_flow_origin[:-1]
    q_flow_revers = np.flip(q_flow[:-1], 0)
    q_flow_tot = np.append(q_flow, q_flow_revers)
    mises = np.zeros((18,3))
    for i in range(np.shape(normal_stress)[0]):
        normal = (normal_stress_tot[i][2] + normal_stress_tot[i+1][2])/2
        z = (normal_stress_tot[i][0] + normal_stress_tot[i+1][0])/2
        y = (normal_stress_tot[i][1] + normal_stress_tot[i+1][1])/2
        q = q_flow_tot[i]/t_sk
        avg_mises = np.sqrt((normal**2.)+3.*(q**2.))
        mises[i] = np.array([z,y,avg_mises])
    normal = (normal_stress_tot[2][2] + normal_stress_tot[15][2])/2
    z = (normal_stress_tot[2][0] + normal_stress_tot[15][0])/2
    y = (normal_stress_tot[2][1] + normal_stress_tot[15][1])/2
    q = q_flow_origin[-1] / t_sp
    avg_mises = np.sqrt((normal ** 2.) + 3. * (q ** 2.))
    mises[-1] = np.array([z,y,avg_mises])
    return mises


''' #get von mises at Ribs
n_points = 1000
x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((n_points,4))
for i in range(np.size(x_points)):
    run = von_mises(x_points[i])
    idx = np.argmax(run[:,2])

    y_points[i] = np.append(x_points[i], run[int(idx)])
idx = np.argmax(y_points[:,3])
massimo = y_points[idx]
print(massimo)
x_points= np.array([x_R_1y, x_F_I, x_P, x_R_2y])
n_points = 1000
#x_points = np.linspace(0, l_a, num=n_points)
y_points = np.zeros((len(x_points),4))
for i in range(np.size(x_points)):
    run = von_mises(x_points[i])
    idx = np.argmax(run[:,2])
    y_points[i] = np.append(x_points[i], run[idx])
print(y_points)
'''



''' Different plots'''
number_points = 1000
#get_simple_plot(get_moment_Mx, number_points, 1, 1, 1)
get_simple_plot1(get_displacement, number_points, C1, -C2, int_constants_z, forces_y) # displacement plot z
#get_simple_plot1(get_displacement, number_points, -C3, C1, RF_y[3:]) # displacemnent plot y

#get_simple_plot(get_slope, number_points, C1, -C2, int_constants_z) # slope plot z
#get_simple_plot(get_slope, number_points, -C3, C1, RF_y[3:]) # slope plot y
#get_simple_plot(get_curvature, number_points, C1, -C2, int_constants_z) # curvature plot z
#get_simple_plot(get_curvature, number_points, -C3, C1, RF_y[3:]) # curvature plot y
#get_simple_plot(get_moment_My, number_points, C1, C1, RF_y[3:]) # moment plot My
#get_simple_plot(get_moment_Mz, number_points, C1, C1, RF_y[3:]) # moment plot Mz
#get_simple_plot(get_shear_Fz, number_points, C1, C1, RF_y[3:]) # shear plot Fz
#get_simple_plot(get_shear_Fy, number_points, C1, C1, RF_y[3:]) # shear plot Fy
##get_simple_plot(get_normal_stress, number_points, C1, C1, RF_y[3:])

plot_3D(number_points)

''' plot 3D von mises'''
n_points = 100 #or 1000
x_points = np.linspace(0, l_a, num = n_points)
mises = np.array([0.,0.,0.])
x_points_array = np.array([0])
for i in range(np.size(x_points)):
    mises = np.vstack((mises,von_mises(x_points[i])))
    point = np.ones((18,1))*x_points[i]
    x_points_array = np.vstack((x_points_array,point))
mises = mises[1:]
x_points_array = x_points_array[1:]

import plotly.plotly as py
import plotly.graph_objs as go
import plotly
import numpy as np
import matplotlib as mpl
#x, y, z = np.random.multivariate_normal(np.array([0,0,0]), np.eye(3), 400).transpose()
norm = mpl.colors.Normalize(vmin=np.min(mises[:,2]), vmax = np.max(mises[:,2]))
c = [norm(i) for i in mises[:,2]]
trace1 = go.Scatter3d(
    x=x_points_array[:].reshape((len(x_points_array),)),
    y=mises[:,1].reshape((len(x_points_array),)),
    z=mises[:,0].reshape((len(x_points_array),)),
    mode='markers',
    marker=dict(
        size=12,
        color=mises[:,2],       # set color to an array/list of desired values
        colorscale='Viridis',   # choose a colorscale
        showscale = True,
        opacity=1.,

    ),

)

data = [trace1]
layout = go.Layout(
                    scene = dict(
                    xaxis = dict(
                         range = [0.,2.7],
                         nticks=4,
                         title='span-wise location [m]' ),
                    yaxis = dict(
                         range = [-0.3,0.3],
                         nticks=4,
                        title='y location [m]'),
                    zaxis = dict(
                         range = [-0.3,0.3],
                         nticks=4,
                        title='chord-wise location [m]'),),
                    width=700,
                    margin=dict(
                    r=20, l=10,
                    b=10, t=10)
                  )
fig = go.Figure(data=data,layout=layout)
plotly.offline.plot(fig, filename='von_mises.html')

