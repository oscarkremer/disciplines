import numpy as np
import pandas as pd
from numpy.linalg import inv
import matplotlib.pyplot as plt

if __name__=='__main__':
    m1 = 1 
    m2 = 1 
    L1 = 1
    L2 = 1
    delta_time = 0.002
    final_time = 5
    t = np.linspace(0,final_time,int(final_time/delta_time)+1)
    t_plot = t[:-2]
    sizes = t.shape[0]
    y = 0.3*np.sin(5*t) + 0.4*np.ones(sizes)
    x = 0.3*np.cos(5*t) + 0.2*np.ones(sizes)
    cos_theta2 = (np.power(x,2) + np.power(y,2) - (L1**2)*np.ones(sizes) - (L2**2)*np.ones(sizes))/(2*L1*L2);
    sin_theta2 = np.sqrt((1 - np.power(cos_theta2,2)))
    theta2_d = np.arctan2(sin_theta2, cos_theta2)
    theta1_d = np.arctan2(y, x) - np.arctan2(L2*np.sin(theta2_d), L1*np.ones(sizes) + L2*np.cos(theta2_d));
    theta1_dotd = (1/delta_time)*np.diff(theta1_d)
    theta2_dotd = (1/delta_time)*np.diff(theta2_d)
    theta1_dot2d = (1/delta_time)*np.diff(theta1_dotd)
    theta2_dot2d = (1/delta_time)*np.diff(theta2_dotd)
    q_des = np.array([theta1_d, theta2_d])
    q_dot_des = np.array([theta1_dotd, theta2_dotd])
    q_dot2_des = np.array([theta1_dot2d, theta2_dot2d])
    q0 = [0.1, 0]
    q_dot0 = [0, 0]
    k_p = 2500
    k_v = 2*k_p**0.5
    k_i = 1000
    g = 9.8
    q_plot1 = np.zeros(sizes-2)
    q_plot2 = np.zeros(sizes-2)
    computed_torque1 = np.zeros((sizes-2,2))
    computed_torque2 = np.zeros((sizes-2,2))
    perturb = 5
    error_old = q_des[:,0] - q0
    error_plot = np.zeros((sizes-2,2))
    int_error = [0, 0]
    for i in range(sizes-2):
        M = np.zeros((2,2))
        M[0][0] = (m1+m2)*L1**2 + m2*L2**2 + 2*m2*L1*L2*np.cos(q0[1])
        M[0][1] = m2*L2**2+m2*L1*L2*np.cos(q0[1])
        M[1][0] = m2*L2**2 + m2*L1*L2*np.cos(q0[1])
        M[1][1] = m2*L2**2
        V = np.array([[-m2*L1*L2*(2*q_dot0[0]*q_dot0[1]+q_dot0[1]**2)*np.sin(q0[1])],
            [m2*L1*L2*(q_dot0[0]**2)*np.sin(q0[1])]])
        G = np.array([[(m1+m2)*g*L1*np.cos(q0[0])+m2*g*L2*np.cos(q0[0]+q0[1])],[m2*g*L2*np.cos(q0[0]+q0[1])]])   

        tau = M.dot(q_dot2_des[:,i]+k_v*(q_dot_des[:,i]-q_dot0) + k_p*(q_des[:,i]-q0)).reshape((2,1))+V+G
        
        q_dot2 = inv(M).dot(tau-V-G)
        q_dot = np.array(q_dot0).reshape((2,1))+ delta_time*q_dot2
        q = np.array(q0).reshape((2,1)) + delta_time*q_dot + 0.5*q_dot2*(delta_time)**2;
        error_plot[i,:] = q_des[:,i]-q0
        q0 = list(q.reshape(2))

        q_dot0 = list(q_dot.reshape(2))
        computed_torque1[i] = tau[0]
        computed_torque2[i] = tau[1]
        q_plot1[i] = q[0]
        q_plot2[i] = q[1]



    plt.plot(t_plot, q_plot1)
    plt.title('Trajectory Following - First Rotational Joint')
    plt.show()


'''
    subplot(2,1,2);
    plot(t_plot, q_plot(:,2), t, theta2_d);
    title('Trajectory Following - Second Rotational Joint')


    figure(4)
    subplot(2,1,1);
    plot(t_plot, computed_torque1);
    title('Computed Torque')

    subplot(2,1,2);
    plot(t_plot, computed_torque2);
    title('Computed Force')


    x_real1  = L1*cos(q_plot(:,1));
    y_real1  = L1*sin(q_plot(:,1));
    x_real2  = L1*cos(q_plot(:,1)) + L2*cos(q_plot(:,1)+q_plot(:,2));
    y_real2  = L1*sin(q_plot(:,1)) + L2*sin(q_plot(:,1)+q_plot(:,2));


    figure(5)
    subplot(2,1,1);
    plot(t_plot, x_real2);
    title('X - Axis')

    subplot(2,1,2);
    plot(t_plot, y_real2);
    title('Y - Axis')


    figure(6)
    subplot(2,1,1);
    plot(t_plot, error_plot(:,1));
    title('Error - Joint 1')

    subplot(2,1,2);
    plot(t_plot, error_plot(:,2));
    title('Error - Joint 2')
'''