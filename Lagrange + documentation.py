#!/usr/bin/env python
# coding: utf-8

# In[1]:


def inputForLagrange():
    
    """
    Функция для ввода данных функции, для которой будут находиться локальные экстремумы.
    
    Parameters
    -----------
    params: str
        Строка на вход со списком переменных через пробел
    F: str
        Функция в виде строки
    L: str
        Обозначение для наличия ограничений функции('1' - есть, '0' - нет)
    lim1, lim2: str
        Ограничения функции
    Z: str
        Ограничивающая функция
        
    Returns
    -----------
    Final: str
        Словарь со всеми переменными и функциями, которые ранее были введены
    
    """
    
    params = input('Список переменных через пробел : ')
    params = params.split(' ')
    if len(params) > 2:
        print('Введено более двух переменных')
    elif len(params) < 2 :
        print('Введено менее двух переменных')
    F = input('Функция в аналитическом виде :  ')
    L = input('Есть ли ограничения, 1 - да; 0 нет : ')
    if L == '1': 
        lim1 = input(f'Ограничения для {params[0]} через пробел :  ')
        lim2 = input(f'Ограничения для {params[1]} через пробел :  ')
    elif L =='0':
        lim1 = ('999')
        lim2 = ('999')
    else :
        print('Ошибка ввода наличия ограничений')
    lim1 = list(map(float, (lim1.split(' '))))
    lim2 = list(map(float, (lim2.split(' '))))

    Z = input('Ограничивающая функция :  ')
    Final = {'p1': params[0],
             'p2': params[1],
             'func': F,
             'lims1': lim1,
             'lims2': lim2,
             'Z' : Z}
    return Final


from sympy import *

def Lagrange(dictionary):
    
    """
    Функция для поиска экстремумов функции с пом. метода Лагранжа
    
    Parameters
    -----------
    dictionary: dict
        На вход подаётся return функции inputForLagrange() - словарь со всеми переменными и функциями,
        которые ранее были введены
        
    Returns
    -----------
    : dict
        Словарь с функцией и её переменными, координатами точек экстремумов и значением функции в них, типов экстремумов
    
    """
    
    import numpy as np
    # преобразование данных для символьного вычислнения
    from sympy.parsing.sympy_parser import parse_expr

    data = dictionary
    func = data['func'] + ' + l *' + '(' + data['Z'] + ')'
    func = parse_expr(func)
    z = parse_expr(data['func'])
    p1 = data['p1']
    x = Symbol(p1)
    p2 = data['p2']
    y = Symbol(p2)
    l = Symbol('l')
    lim1 = data['lims1']
    lim2 = data['lims2']
    ## реализация метода
    
    dx = func.diff(p1)
    dy = func.diff(p2)
    dl = func.diff(l)
    points = solve((dx,dy,dl), [x,y,l],dict = True)
    
    M = Matrix([[0, z.diff(x), z.diff(y)],
                [z.diff(x), func.diff(x,2),func.diff(x,y)],
                [z.diff(y) , func.diff(x,y), func.diff(y,2)]])
    determinant = M.det()
    
    G = Matrix([[func.diff(x,2), func.diff(x,y)],
                [func.diff(x,y), func.diff(y,2)]])

    for i in points :
        if len(lim1) == len(lim2) == 2 :

            if (float(i[x]) >= lim1[0] and float(i[x]) <= lim1[1]) and (float(i[y]) >= lim2[0] and float(i[y]) <= lim2[1]):

                if determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0:
                    i['тип'] ='условный максимум'
                    i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
                    print(i)
            
                elif determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) < 0:
                    i['тип'] ='условный минимум'
                    i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
                    print(i)
    
                elif (((func.diff(x,2).subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0) * (G.det().subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0)) == 0 or
                ((func.diff(x,2).subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) < 0) * (G.det().subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0))) == 0:
                    i['тип'] ='седловая точка'
                    i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
                    print(i)

                else:
                    i['тип'] = 'требуется дополнительное исследование'
                    i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
                    print(i)    
    else :
        if determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0:

            i['тип'] ='условный максимум'
            i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
            print(i)

        elif determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) < 0:
            i['тип'] ='условный минимум'
            i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
            print(i)
    
        elif (((func.diff(x,2).subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0) * (G.det().subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0)) == 0 or
        ((func.diff(x,2).subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) < 0) * (G.det().subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0))) == 0:
            i['тип'] ='седловая точка'
            i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
            print(i)

        else:
            i['тип'] = 'требуется дополнительное исследование'
            i['Func'] = func.subs(x,i[x]).subs(y,i[y]).subs(l,i[l])
            print(i)   

            
    if len(lim1) == 2 :
        xx = np.linspace(lim1[0] - 5, lim1[1] + 5, 1000)
        yy = np.linspace(lim2[0] - 5, lim2[1] + 5, 1000)
    else:
        xx = np.linspace(-20, 20, 1000)
        yy = np.linspace(-20, 20, 1000)
    
    f = lambdify([x,y], z)
  

    coord_x_p = np.array([float(i[x]) for i in points if 'тип' in i]) 
    coord_y_p = np.array([float(i[y]) for i in points if 'тип' in i])
    coord_z_p = f(coord_x_p, coord_y_p)
    M = input('Масштабируем график, 1-да , 0- нет :')
    if M=='1':
        limx1=(float(sum(coord_x_p)/len(coord_x_p))+max(coord_x_p))
        limx2=(float(sum(coord_x_p)/len(coord_x_p))-max(coord_x_p))
        limy1=(float(sum(coord_y_p)/len(coord_y_p))+max(coord_y_p))
        limy2=(float(sum(coord_y_p)/len(coord_y_p))-max(coord_y_p))
        xx = np.linspace(limx2, limx1, 1000)
        yy = np.linspace(limy2, limx1 , 1000)
        
    elif M =='0':
        limx1 = ('999')
        limx2 = ('999')
        limy1 = ('999')
        limy2 = ('999')
        
    else :
        return 'Ошибка ввода наличия ограничений'
    
    t = [i['тип'] for i in points if i['тип']]
    X, Y = np.meshgrid(xx, yy)

    
    Z = f(X,Y)
        

    
    return {'x' : x, # переменная 1
            'y' : y, # переменная 2
            'z': z, # начальная функция
            'func' : func, # функция с лямбдой
            'x_p' : coord_x_p, # координаты точек экстремумов по 1 пер (после ограничений)
            'y_p' : coord_y_p, # координаты точек экстремумов по 2 пер (после ограничений)
            'z_p' : coord_z_p, # значения функции в экстремумах,
            'X' : X,
            'Y' : Y,
            'Z' : Z,
            'type': t,
           'limx1' : limx1,
           'limx2' : limx2,
           'limy1' : limy1,
           'limy2' : limy2
           }

def graph(d):
    
    """
    Функция для построения трёхмерного графика.
    
    Parameters
    -----------
    d: dict
        Словарь с функцией и её переменными, координатами точек экстремумов и значением функции в них, типов экстремумов
        
    Returns
    -----------
    plot
        График с отмеченным точками разных цветов в зависимости от типа экстремума с возможностью масштабирования
    
    """   
    
    #%matplotlib notebook
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.figsize'] = (8,6)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    

    surf = ax.plot_surface(d['X'], d['Y'], d['Z'], cmap='Purples',
                           alpha=0.5,linewidth=1, antialiased=False)
    
    for i in range(len(d['type'])):
        if d['type'][i] == 'седловая точка':
            ax.scatter3D(d['x_p'][i], d['y_p'][i], d['z_p'][i], c = 'red', s=50, alpha=1)
            
        elif d['type'][i] == 'условный минимум':
            ax.scatter3D(d['x_p'][i], d['y_p'][i], d['z_p'][i], c = 'yellow', s=50, alpha=1)
            
        elif d['type'][i] == 'условный максимум':
            ax.scatter3D(d['x_p'][i], d['y_p'][i], d['z_p'][i], c = 'green', s=50, alpha=1)
        
        elif d['type'][i] == 'требуется дополнительное исследование':
            ax.scatter3D(d['x_p'][i], d['y_p'][i], d['z_p'][i], c = 'orange', s=50, alpha=1)

    ax.set_xlabel(f"{d['x']}")
    ax.set_ylabel(f"{d['y']}")
    ax.set_zlabel(f"{d['z']}")
    plt.show()

def graph_lines(d):
    
    """
    Функция для построения графика проекций с линиями уровня.
    
    Parameters
    -----------
    d: dict
        Словарь с функцией и её переменными, координатами точек экстремумов и значением функции в них, типов экстремумов
        
    Returns
    -----------
    plot
        График с линиями уровня функции
        
    """   
    
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.figsize'] = (6,6)
    fig, ax = plt.subplots()
    
    ax.contour(d['Z'])
    ax.clabel(ax.contour(d['Z']))
    
    ax.set_title(f'Линии уровня функции {d["z"]}')
    plt.show()


# In[2]:


help(inputForLagrange)


# In[3]:


help(Lagrange)


# In[4]:


help(graph)


# In[5]:


help(graph_lines)

