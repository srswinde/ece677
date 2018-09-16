import pandas as pd
import matplotlib.pyplot as plt

def amdahl(P, N, S):
    return 1/((P/N)+S)


def main(upperlimit=10000):
    speedups=[]
    for N in range(1, upperlimit):
        P=0.50
        S=1-P
        speedups.append( amdahl(P, N, S) )

    dt = pd.Series(speedups)
    ax=plt.gca()
    ax.set_xlabel("Number of Computers")
    ax.set_ylabel("Speed Up Factor")
    dt.plot(ax=ax, grid=True, title="Number of Computers VS Speed Up")
    return dt




