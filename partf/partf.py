import subprocess
import numpy as np
import matplotlib.pyplot as plt


delts = [3e-6, 2e-6, 1e-6, 5e-7,1e-7]
delts = [5e-8, 1e-8]
#delts = [0.000001]

def partf():
    print("Running Part f")
    t_end = 0.1

    # temporal convergence
    xnum = 50
    #delts = [0.0001, 0.00001]
    for delt in delts:
        output_file = f'./partf/files/{delt:.0e}'.replace('e-', 'em').replace('e+', 'ep') + '.txt'
        print(output_file)
        process = subprocess.Popen(
            ["./mp2", str(xnum), str(delt), str(t_end), str(output_file)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            print(line, end='')

        process.wait()

        #with open(output_file, 'r') as f:
        #    lines = f.readlines()
        #    x_vec1 = (lines[0].strip()).split()
        #    vel_vec1 = (lines[-1].strip()).split()
        #x_vec1 = np.array(x_vec1[1:], dtype=float)
        #vel_vec1 = np.array(vel_vec1[1:], dtype=float)


def read_field(filename):
    fields = {}
    current = None
    buffer = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Detect new field
            if line.endswith(":"):
                if current is not None:
                    fields[current] = np.array(buffer, dtype=float)
                current = line[:-1]
                buffer = []
                continue

            # Skip empty lines
            if not line:
                continue

            # Read numeric data
            if current is not None:
                buffer.append([float(x) for x in line.split()])

    # Save last field
    if current is not None:
        fields[current] = np.array(buffer, dtype=float)

    return fields

def plot_functions():

    #for delt in delts:
    output_file = f'./partf/files/{delts[0]:.0e}'.replace('e-', 'em').replace('e+', 'ep') + '.txt'
    
    # -------- Load data --------
    data = read_field(output_file)

    dens = data["dens"]
    xvel = data["xvel"]
    yvel = data["yvel"]
    ener = data["ener"]

    ny, nx = dens.shape
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))

    # -------- Density --------
    plt.figure()
    plt.contourf(X, Y, dens, 50)
    plt.colorbar(label="Density")
    plt.title("Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()

    # -------- Energy --------
    plt.figure()
    plt.contourf(X, Y, ener, 50)
    plt.colorbar(label="Energy")
    plt.title("Energy")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()

    plt.figure()
    plt.contourf(X, Y, xvel, 50)
    plt.colorbar(label="|Velocity|")
    plt.title("X Velocity")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()

    plt.figure()
    plt.contourf(X, Y, yvel, 50)
    plt.colorbar(label="|Velocity|")
    plt.title("Y Velocity")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()
    # -------- Velocity Magnitude --------
    velmag = np.sqrt(xvel**2 + yvel**2)

    plt.figure()
    plt.contourf(X, Y, velmag, 50)
    plt.colorbar(label="|Velocity|")
    plt.title("Velocity Magnitude")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.show()

    # -------- Velocity Vectors --------
    # Normalize vectors so all arrows have the same length
    u = xvel / velmag
    v = yvel / velmag

    skip = max(nx // 25, 1)

    plt.figure(figsize=(10, 8))
    Q = plt.quiver(
        X[::skip, ::skip],
        Y[::skip, ::skip],
        u[::skip, ::skip],
        v[::skip, ::skip],
        velmag[::skip, ::skip],  # color = magnitude
        cmap='viridis',
        scale=25,           # HIGHER number = smaller arrows
        scale_units='width', # keeps arrows consistent across different plot sizes
        width=0.003,        # arrow shaft thickness
        headwidth=3,        # arrow head width
        headlength=4,       # arrow head length
        minlength=0.1       # don't draw arrows that are too small
    )
    plt.colorbar(Q, label="Velocity magnitude")
    plt.gca().set_aspect('equal')  # equal aspect ratio
    plt.gca().invert_yaxis()  # top-left = (0,0)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Velocity Field (uniform arrow length, color = magnitude)")
    plt.tight_layout()
    plt.show()


def convergence_plot():

    max_vel = []
    for delt in delts:
        output_file = f'./partf/files/{delt:.0e}'.replace('e-', 'em').replace('e+', 'ep') + '.txt'
    
        # -------- Load data --------
        data = read_field(output_file)

        dens = data["dens"]
        xvel = data["xvel"]
        yvel = data["yvel"]
        ener = data["ener"]
        velmag = np.sqrt(xvel**2 + yvel**2)
        max_vel.append(np.max(velmag))
        # Take the "exact" solution as the one with the smallest delt
    exact_vel = max_vel[np.argmin(delts)]

    # Compute relative error
    rel_error = np.abs(max_vel - exact_vel) / exact_vel

    # Log-log plot
    plt.figure()
    plt.loglog(delts, rel_error, marker='o', linestyle='-')
    plt.loglog(delts, delts, 'r--', label=f'O(Δt)')
    plt.xlabel("Δt")
    plt.ylabel("Relative Error")
    plt.title("Relative Error of Max Velocity vs Δt")
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    plt.show()

partf()
convergence_plot()
#plot_functions()
