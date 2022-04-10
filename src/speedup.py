import matplotlib.pyplot as plt



if __name__ == "__main__":
    # strong scaling N=3200, NT=400, L=1.0, T=1.0e3, u=5.0e-7, v=2.85e-7
    datax = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    datay = []
    times = [34.6037, 23.4254, 17.9786, 15.3101, 13.4425, 12.5666, 12.9395, 12.9384, 13.1408, 12.6289, 12.5282, 12.5314]
    for i in range(len(times)):
        datay.append(times[0] / times[i])
    plt.plot(datax, datay, "-o")
    plt.title("Strong Scaling - N=3200")
    plt.xlabel("Number of Cores (n)")
    plt.ylabel("Speedup (S[n])")
    plt.savefig("../out/speedup1.png")
    plt.clf()

    # strong scaling N=200, NT=400, L=1.0, T=1.0e3, u=5.0e-7, v=2.85e-7
    datax = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    datay = []
    times = [0.156394, 0.117354, 0.097103, 0.088419, 0.085739, 0.087316, 0.091249, 0.092673, 0.094445, 0.09261, 0.09621, 0.09923]
    for i in range(len(times)):
        datay.append(times[0] / times[i])
    plt.plot(datax, datay, "-o")
    plt.title("Strong Scaling - N=200")
    plt.xlabel("Number of Cores (n)")
    plt.ylabel("Speedup (S[n])")
    plt.savefig("../out/speedup2.png")
    plt.clf()

    # weak scaling
    datax = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    datay = []
    times = [2.18952, 2.91594, 3.24746, 3.74945, 4.17908, 4.6611, 5.81468, 6.63324, 7.53734, 8.07801, 8.7695, 9.53799]
    for i in range(len(times)):
        datay.append(times[0] / times[i])
    plt.plot(datax, datay, "-o")
    plt.title("Weak Scaling Analysis")
    plt.xlabel("Number of Cores (n)")
    plt.ylabel("Speedup (S[n])")
    plt.savefig("../out/speedup3.png")
    plt.clf()