# importing libraries
from audioop import reverse
from cmath import pi
from tkinter.tix import MAX
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtp
import math as m
import csv
# plotly

import sys
import os

from sqlalchemy import true

from functions.levy import *
from genetic.genetic import *



def make_graphic(
    folder: str,
    FUNCTION=None,
    N=2,
    MIN=None, MAX=None,
    DIM=400,
    title='NONE function',
    REAL_MINIMA_X = None, #[np.array([1, 1])],
    limits=None,
    **kwargs):

    if (not os.path.isdir(folder)):
        os.makedirs(folder)

    # N-dimensional X values
    x = []
    if limits is None:
        for i in range(N):
            x.append(np.linspace(MIN, MAX, DIM))
    else:
        for i in range(N):
            x.append(np.linspace(limits[i][0], limits[i][1], DIM))

    X = np.array(np.meshgrid(x[0], x[1]))
    F = FUNCTION(X)
    # input("Here")

    # print(X, np.prod(X))

    # Building graph
    fig = plt.figure(figsize=(20, 10))
    # syntax for 3-D plotting
    ax = fig.add_subplot(1, 2, 1,projection ='3d')

    surf = ax.plot_surface(X[0], X[1], F, cmap='coolwarm', linewidth=0, alpha=0.7)
    #, cmap ='viridis', edgecolor ='green')
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    #ax.set_zlim(-20, 100)
    ax_projection = fig.add_subplot(1, 2, 2)
    # print(REAL_MINIMA_X)
    cset = ax_projection.contourf(X[0], X[1], F, 256, cmap='coolwarm', origin='lower')
    fig.colorbar(cset, shrink=0.5, aspect=5)

    if REAL_MINIMA_X is not None:
        for i in range(len(REAL_MINIMA_X)):
            REAL_MINIMA = FUNCTION(REAL_MINIMA_X[i])
            ax.scatter(REAL_MINIMA_X[i][0], REAL_MINIMA_X[i][1], REAL_MINIMA, c="red", marker="x")
            ax_projection.scatter(REAL_MINIMA_X[i][0], REAL_MINIMA_X[i][1], c="red", marker="x")
        if FUNCTION == G3:
            c = plt.Circle((0, 0), 1, fill=False)
            ax_projection.add_artist(c)
    #fig.tight_layout()
    #plt.tight_layout()
    plt.savefig(f"{folder}/raw.png", dpi=300)
    #plt.show()
    return (fig, ax, ax_projection)


def generate(
    stepname: int=0,
    GENERATIONS: int=200,
    POPULATION: int=100, #DIM
    CHROMO_DIM_COOF: int=128,
    need_to_see=[
        0, 1, 10, 25, 50,
        100, 150, 200, 250,
        300, 350, 400, 450, 500],
    FUNCTION=None,
    resdir:str="res",
    N:int=2, DIM:int=400,
    fig=None, ax=None, ax_projection=None,
    title:str="NONE function",
    MIN:int=None, MAX:int=None, 
    limits:tuple=None,
    KEYFUNCTION:Callable=lambda x: 0,
    REAL_MINIMA_X=None,
    **kwargs):

    reverse = False
    if "reverse" in kwargs:
        reverse = kwargs["reverse"]

    REZDIR = f"{resdir}/result_{stepname}"
    if (not os.path.isdir(REZDIR)):
        os.makedirs(REZDIR)
    
    if fig is None or ax is None or ax_projection is None:
        fig, ax, ax_projection = make_graphic(
            resdir, FUNCTION=FUNCTION, N=N,
            MIN=MIN, MAX=MAX, DIM=DIM, title=title, limits=limits, REAL_MINIMA_X=REAL_MINIMA_X)

    mx = (MIN, MAX)
    if MIN is None or MAX is None:
        mx = None

    pop = Population(POPULATION, N, CHROMO_DIM_COOF*DIM, FUNCTION, maxes=mx, limits=limits, KEYFUNCTION=KEYFUNCTION, reverse=reverse)

    CHROMO_SIZE = pop.get_accuracy()

    drawn = []
    
    figExits = plt.figure()
    axeExit = figExits.add_subplot()

    BestInidivid = []
    
    GENDIR = REZDIR + "/generations"
    if (not os.path.isdir(GENDIR)):
        os.makedirs(GENDIR)

    TABLESF = REZDIR + "/tables"
    if (not os.path.isdir(TABLESF)):
        os.makedirs(TABLESF)

    for _ in range(GENERATIONS+1):
        #print(f">{_}")
        rez = pop.count()
        
        #while len(ax.collections) >= 3 : ax.collections.remove(ax.collections[2])
        #print(_)
        if _ % 5 == 0 or _ in need_to_see:#in need_to_see:
            exits = []
            
            if _ in need_to_see:

                for i in drawn:
                    i.set_visible(False)
                    i.remove()
                drawn.clear()

                for i in rez:
                    drawn.append(ax.scatter3D(i[0], i[1], i[2]))
                    drawn.append(ax_projection.scatter(i[0], i[1]))
                    ax.set_title(f'GENERATION {_}')
                    exits.append(i[2])
                fig.canvas.draw()
                fig.canvas.flush_events()
                
                fig.savefig(f"{GENDIR}/generation_{_}.png", dpi=150)
                #plt.show()
                #input()
            else:   
                for i in rez:
                    exits.append(i[2])
            
            BestInidivid.append( tuple((_, pop.get_best())) )
            axeExit.scatter(list([_]*POPULATION), exits)

            #print(rez)
            # plt.pause(0.00000001)

        pop.next()

    figExits.savefig(f"{REZDIR}/generations.png", dpi=300)
    axeExit.set_ylim(-5, 20)
    figExits.savefig(f"{REZDIR}/generations-better.png", dpi=300)

    # save Arguments
    with open(TABLESF + "/args.csv", "w", encoding="utf8") as argfile:
        writer = csv.writer(argfile, delimiter=";")
        writer.writerow(["Аргумент", "Значение"])

        writer.writerow(["Размерность вектора агрументов", N])
        writer.writerow(["Размер хромосомы (бит)", CHROMO_SIZE])
        writer.writerow(["Размер популяции", POPULATION])
        if limits is not None:
            for i in range(N):
                writer.writerow([f"Минимум аргумента {i+1}", limits[i][0]])
                writer.writerow([f"Максимум аргумента {i+1}", limits[i][1]])
        else:
            writer.writerow([f"Минимум аргумента", MIN])
            writer.writerow([f"Максимум аргумента", MAX])

        writer.writerow(["Достижимая погрешность", round((20) / (2**(CHROMO_SIZE) - 1), 4) ])

        argfile.close()

    # save Values
    with open(TABLESF + "/generations.csv", "w", encoding="utf8", newline='') as argfile:
        writer = csv.writer(argfile, delimiter=";")
        writer.writerow(["Итеррация"] + [f"$x_{i}$" for i in range(1, N+1)] + ["F(x)"])
        # print(BestInidivid)
        
        for i in range(len(BestInidivid)):
            writer.writerow(
                [BestInidivid[i][0]] + \
                ([] + list(BestInidivid[i][1][0])) + \
                [BestInidivid[i][1][1][-1]]
            )
        

        argfile.close()

    return list(BestInidivid[-1])


def RunGenetic(
    title: str,
    FUNCTION: Callable,
    folder:str,
    STEPS:int=10,
    POPULATION=75,
    GENERATIONS=250,
    MIN=None, MAX=None, limits=None, **kwargs):

    if not os.path.isdir(folder): os.makedirs(folder)

    fig, ax3d, subx = make_graphic(
        folder, FUNCTION=FUNCTION, title=title, MIN=MIN, MAX=MAX, limits=limits, **kwargs
    ) 
    bests = []
    N = 2
    for i in range(STEPS):
        bests.append([i] + generate(
                resdir=folder,
                FUNCTION=FUNCTION,
                title=title,
                stepname=i,
                POPULATION=POPULATION,
                GENERATIONS=GENERATIONS,
                fig=fig,
                ax=ax3d,
                ax_projection=subx,
                MAX=MAX,
                MIN=MIN,
                limits=limits,
                **kwargs
            ) 
        )

    with open(folder + "/repeats.csv", "w", encoding="utf8", newline='') as argfile:
            writer = csv.writer(argfile, delimiter=";")
            writer.writerow(["Шаг", "Поколение"] + [f"$x_{i}$" for i in range(1, N+1)] + ["F(x)"])

            for i in range(len(bests)):
                #print(bests[i])
                writer.writerow(
                    [bests[i][0], bests[i][1]] + \
                    ([] + list(bests[i][2][0])) + \
                    [bests[i][2][1]]
                )
            

            argfile.close()

#a = np.array([np.uint8(123), np.uint8(123), np.uint8(123)])
#print( np.unpackbits(a) )

from functions.zakharov import *

from functions.g3 import *

EPS = 0.01

def kf(x):
    # print(x)
    ax = round(abs(x[0]**2 + x[1]**2 - 1), 4)
    return -ax if ax <= EPS else ax*-10000 # *1000

RunGenetic(
    "G3",
    G3,
    "res/g3",
    STEPS=5,
    GENERATIONS=500,
    POPULATION=250,
    limits=[(0, 1), (0, 1)],
    REAL_MINIMA_X=[ (1/np.sqrt(2), 1/np.sqrt(2)) ],
    KEYFUNCTION=kf,
    reverse=True # 1000 if abs(x[0]**2 + x[1]**2 - 1) >= EPS else 0
)

RunGenetic(
    "Levy",
    LevyFull,
    "res/levy",
    STEPS=5,
    limits=( (-10, 10), (-10, 10) ),
    REAL_MINIMA_X=[np.array([1, 1])]
)


RunGenetic(
    "Zakharov",
    Zakharov,
    "res/zakharov",
    STEPS=5,
    limits=( (-10, 10), (-10, 10) ),
    REAL_MINIMA_X=[np.array([0, 0])]
)#'''