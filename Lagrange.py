def inputForLagrange():
    from numpy import inf
    params = input('Список переменных через пробел : ')
    params = params.split(' ')
    if len(params) > 2:
        print('Введено более двух переменных')
    elif len(params) < 2 :
        print('Введено менее двух переменных')
    F = input('Функция в аналитическом виде :  ')
    if L == '1':
        lim1 = input(f'Ограничения для {params[0]} через пробел :  ')
        lim2 = input(f'Ограничения для {params[1]} через пробел :  ')
    elif L =='0':
        lim1 = (f'{-inf} {inf}')
        lim2 = (f'{-inf} {inf}')
    else :
        return 'Ошибка ввода наличия ограничений'
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
        if determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) > 0:
            print(i, 'условный максимум')
        elif determinant.subs(x,i[x]).subs(y,i[y]).subs(l,i[l]) < 0:
            print(i , 'условный минимум')
        elif ((func.diff(x,2) > 0) * (G.det() > 0)) == 0 or ((func.diff(x,2) < 0) * (G.det() > 0)) == 0:
            print(i , 'седловая точка')
        else:
            print(i , 'требуется дополнительное исследование')    
        
        