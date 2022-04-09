# importing libraries
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

from functions.levy import *
from genetic.genetic import *



def make_graphic(folder: str, FUNCTION=LevyFull, N=2, MIN=None, MAX=None, DIM=400, title='Levy function', REAL_MINIMA_X = np.array([1, 1]), **kwargs):

    if (not os.path.isdir(folder)):
        os.makedirs(folder)

    # N-dimensional X values
    x = []
    for i in range(N):
        x.append(np.linspace(MIN, MAX, DIM))

    X = np.array(np.meshgrid(x[0], x[1]))
    F = FUNCTION(X)

    REAL_MINIMA   = FUNCTION(REAL_MINIMA_X)

    # Building graph
    fig = plt.figure(figsize=(20, 10))
    # syntax for 3-D plotting
    ax = fig.add_subplot(1, 2, 1,projection ='3d')

    surf = ax.plot_surface(X[0], X[1], F, cmap='coolwarm', linewidth=0, alpha=0.7)#, cmap ='viridis', edgecolor ='green')

    ax.scatter(REAL_MINIMA_X[0], REAL_MINIMA_X[1], REAL_MINIMA,c="red", marker="x")

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    #ax.set_zlim(-20, 100)
    ax_projection = fig.add_subplot(1, 2, 2)
    cset = ax_projection.contourf(X[0], X[1], F, 256, cmap='coolwarm', origin='lower')
    ax_projection.scatter(REAL_MINIMA_X[0], REAL_MINIMA_X[1], c="red", marker="x")

    fig.colorbar(cset, shrink=0.5, aspect=5)

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
    need_to_see=[0, 1, 10, 25, 50, 100, 150, 200],
    FUNCTION=LevyFull,
    resdir:str="res",
    N:int=2, MIN:int=None, MAX:int=None, DIM:int=400,
    fig=None, ax=None, ax_projection=None, title:str="Levy function", **kwargs):

    REZDIR = f"{resdir}/result_{stepname}"
    if (not os.path.isdir(REZDIR)):
        os.makedirs(REZDIR)
    
    if fig is None or ax is None or ax_projection is None:
        fig, ax, ax_projection = make_graphic(resdir, FUNCTION=FUNCTION, N=N, MIN=MIN, MAX=MAX, DIM=DIM, title=title)

    pop = Population(POPULATION, N, CHROMO_DIM_COOF*DIM, (MIN, MAX), FUNCTION)

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
        writer.writerow(["Достижимая погрешность", round((20) / (2**(CHROMO_SIZE) - 1), 4) ])

        argfile.close()

    # save Values
    with open(TABLESF + "/generations.csv", "w", encoding="utf8", newline='') as argfile:
        writer = csv.writer(argfile, delimiter=";")
        writer.writerow(["Итеррация"] + [f"$x_{i}$" for i in range(1, N+1)] + ["F(x)"])

        for i in range(len(BestInidivid)):
            writer.writerow(
                [BestInidivid[i][0]] + \
                ([] + list(BestInidivid[i][1][0])) + \
                [BestInidivid[i][1][1]]
            )
        

        argfile.close()

    return list(BestInidivid[-1])


def RunGenetic(title: str, FUNCTION: Callable, folder:str, STEPS:int=10, POPULATION=50, GENERATIONS=250,MIN=-10, MAX=10, **kwargs):

    if not os.path.isdir(folder): os.makedirs(folder)

    fig, ax3d, subx = make_graphic(folder, FUNCTION=FUNCTION, title=title, MIN=MIN, MAX=MAX, **kwargs) 
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

from functions.zakharov import Zakharov


RunGenetic(
    "Levy",
    LevyFull,
    "res/levy",
    MIN=-10,
    MAX=10,
    STEPS=5,
    REAL_MINIMA_X=np.array([1, 1])
)

RunGenetic(
    "Zakharov",
    Zakharov,
    "res/zakharov",
    MIN=-10,
    MAX=10,
    STEPS=5,
    REAL_MINIMA_X=np.array([0, 0])
)