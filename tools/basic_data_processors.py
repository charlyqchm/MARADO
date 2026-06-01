import numpy as np
import matplotlib.pyplot as plt

au_to_eV = 27.2114

def read_file(filename, col_x, col_y):
    x_array = []
    y_array = []

    infile = open(filename, "r")

    for val in infile:
        read = not any(char == "#" for char in val)
        if read:
            x_array.append(float(val.split()[col_x]))
            y_array.append(float(val.split()[col_y]))

    infile.close()

    x_array = np.array(x_array)
    y_array = np.array(y_array)

    return x_array, y_array

def get_fourier_transform(data, time_step):
    n = len(data)

    # Calculate the frequencies corresponding to the Fourier transform
    freq = np.fft.fftfreq(10*n, d=time_step)

    # Compute the Fourier transform of the data
    ft_data = np.fft.fft(data, n=10*n)

    # Return the frequencies and the Fourier-transformed data
    return freq*2*np.pi, ft_data

def plot_abs_spectrum_1D(only_pulse_dir=None, sample_dir=None,
                         Ex_detector=1, Hy_detector=2,  
                         x_lim=[0, 20], y_lim=None, help=False):

    if help:
        print("This function plots the absorption spectrum of a 1D simulation using OMxRTA.\n"
                "For this, previously run simulations printing the Ex field and Hy field are required.\n"
                "The following parameters can be specified:\n"
            "       only_pulse_dir -> Directory of the simulation with only the pulse (reference simulation). \n"
            "       sample_dir     -> Directory of the simulation with the sample (test simulation). \n"
            "       Ex_detector    -> Index of the Ex detector. Default: 1.\n"
            "       Hy_detector    -> Index of the Hy detector. Default: 2.\n"
            "       x_lim          -> Limits of the x-axis in eV. Default: [0, 20].\n"
            "       y_lim          -> Limits of the y-axis.\n")
        return

    #Check the existence of the file
    if only_pulse_dir is None:
        print("Error: only_pulse_dir must be specified.")
        return
    if sample_dir is None:
        print("Error: sample_dir must be specified.")
        return
    
    filename = only_pulse_dir+"/output_detector_"+str(Ex_detector).zfill(7)+"/point.dat"

    time, E0xL_t = read_file(filename, 0, 1)

    filename = only_pulse_dir+"/output_detector_"+str(Hy_detector).zfill(7)+"/point.dat"
    time, H0yL_t = read_file(filename, 0, 1)

    dt = time[2]-time[1]

    freq, E0xL_w = get_fourier_transform(E0xL_t, dt)
    freq, H0yL_w = get_fourier_transform(H0yL_t, dt)

    filename = sample_dir+"/output_detector_"+str(Ex_detector).zfill(7)+"/point.dat"

    time, ExL_t = read_file(filename, 0, 1)

    filename = sample_dir+"/output_detector_"+str(Hy_detector).zfill(7)+"/point.dat"
    time, HyL_t = read_file(filename, 0, 1)

    dt = time[2]-time[1]

    freq, ExL_w = get_fourier_transform(ExL_t, dt)
    freq, HyL_w = get_fourier_transform(HyL_t, dt)

    Pw = - np.conj(ExL_w-E0xL_w) * (HyL_w-H0yL_w)

    freq = freq*au_to_eV

    fig, ax1 = plt.subplots()

    N = len(freq)
    spec = np.real(Pw[:N//2])
    ax1.plot(freq[:N//2], spec)

    ax1.set_xlim(x_lim[0], x_lim[1])

    if y_lim is not None:
        ax1.set_ylim(y_lim[0], y_lim[1])

    plt.tight_layout()
    plt.show()

    return

def plot_trans_spectrum_1D(directory=None, Ex_detector=1, Hy_detector=2, 
                           x_lim=[0, 20], y_lim=None, help=False):

    if help:
        print("This function plots the transmission spectrum of a 1D simulation using OMxRTA.\n"
                "For this, a previously run simulation printing the Ex field and Hy field is required.\n"
                "The following parameters can be specified:\n"
            "       directory      -> Directory of the simulation. \n"
            "       Ex_detector    -> Index of the Ex detector. Default: 1.\n"
            "       Hy_detector    -> Index of the Hy detector. Default: 2.\n"
            "       x_lim          -> Limits of the x-axis in eV. Default: [0, 20].\n"
            "       y_lim          -> Limits of the y-axis.\n")
        return


    filename = directory+"/output_detector_"+str(Ex_detector).zfill(7)+"/point.dat"

    time, ExL_t = read_file(filename, 0, 1)

    filename = directory+"/output_detector_"+str(Hy_detector).zfill(7)+"/point.dat"
    time, HyL_t = read_file(filename, 0, 1)

    dt = time[2]-time[1]

    freq, ExL_w = get_fourier_transform(ExL_t, dt)
    freq, HyL_w = get_fourier_transform(HyL_t, dt)

    Pw = - np.conj(ExL_w) * (HyL_w)

    freq = freq*au_to_eV

    fig, ax1 = plt.subplots()

    N = len(freq)
    spec = np.real(Pw[:N//2])
    ax1.plot(freq[:N//2], spec)

    ax1.set_xlim(x_lim[0], x_lim[1])

    if y_lim is not None:
        ax1.set_ylim(y_lim[0], y_lim[1])

    plt.tight_layout()
    plt.show()

    return