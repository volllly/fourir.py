# -*- coding: utf-8 -*-
"""Fourier Calculator by Paul Volavsek"""

import argparse
import sympy
    
#Verwendungsbeispiel Parameter punkte:      -p "0, 0" "H, T/2" "0, T" -g -u
#Verwendungsbeispiel Parameter expressions: -e "sin(t * 2 * pi / T)" "T" -g -u


parser = argparse.ArgumentParser(description='Calculates fourier coefficients and series either from an expression or points (use "pi" as pi and "E" as e)')

group = parser.add_mutually_exclusive_group(required=True)

group.add_argument('-p',
                   '--points',
                    metavar="POINT",
                    nargs='+',
                    dest='points',
                    help='points in the format "x, y"')

group.add_argument('-e',
                   '--expression',
                    metavar="EXPRESSION",
                    nargs='+',
                    dest='expression',
                    help='enter expressions and intervall "EXPRESSION; START, END"')

parser.add_argument('-u',
                    '--unicode',
                    action='store_true',
                    default=False,
                    dest='unicode',
                    help='use unicode printing')

parser.add_argument('-i',
                    '--print-input-graph',
                    action='store_true',
                    default=False,
                    dest='ingraph',
                    help='print input graph')

parser.add_argument('-o',
                    '--print-output-graph',
                    metavar='N',
                    nargs=1,
                    dest='outgraph',
                    type=int,
                    help='print fourier graph with N terms')

parser.add_argument('-t',
                    '--textplot',
                    metavar='DIMENSIONS',
                    nargs=2,
                    dest='textplot',
                    help='use ascii to show graphs giving dimensions as X Y')



args = parser.parse_args()

t, k, d = sympy.symbols('t, k, d')
n, terms = sympy.symbols('n, terms', integer=True, positive= True)
w = sympy.symbols('w', positive= True)

points = []
expressions = []
symbols = set()
T = None
if args.expression:
    for expression in args.expression:
        expressions.append([sympy.sympify(expression.split(';')[0]), [expression.split(';')[1], None]])
        [symbols.add(symbol) for symbol in expressions[-1][0].free_symbols]
        for (i, interval) in enumerate(expressions[-1][1][0].split(',')):
            expressions[-1][1][i] = sympy.sympify(interval)
            [symbols.add(symbol) for symbol in expressions[-1][1][i].free_symbols]
            
    symbols.remove(t)
    T = expressions[-1][1][1]
    
    s = sympy.Piecewise(*((expression[0], sympy.And(expression[1][0] <= t, t < expression[1][1])) for expression in expressions))
else:
    for point in args.points:
        points.append((sympy.sympify(point.split(',')[0]), sympy.sympify(point.split(',')[1])))
        for coordinate in points[-1]:
            [symbols.add(symbol) for symbol in coordinate.free_symbols]
    
    T = points[-1][1]

    pieces = []
    i = 0
    while i < len(points) - 1:
        j = i
        if i > 0:
            if points[i][1] == points[i + 1][1]:
                j += 1
        piece = sympy.solve([points[j][1] * k + d - points[j][0], points[j + 1][1] * k + d - points[j + 1][0]], (k, d))
        pieces.append(piece[k] * t + piece[d])
        i += 1
        
    s = sympy.Piecewise(*((pieces[i], sympy.And(points[i][1] <= t, t < points[i + 1][1])) for i in range(0, len(pieces), 1)))
    
print('\ns(t):\n')

sympy.pprint(s, use_unicode=args.unicode)
print()

if args.ingraph:
    if args.textplot:
        sympy.plotting.textplot(s.subs({symbol:1 for symbol in symbols}).subs({terms: args.outgraph[0]}), 0, 1, int(args.textplot[0]), int(args.textplot[1]))
    else:
        import pylab
        import matplotlib
        pylab.plot(matplotlib.numpy.linspace(0, 1, 1000), sympy.lambdify(t, s.subs({symbol:1 for symbol in symbols}), modules=['numpy'])(matplotlib.numpy.linspace(0, 1, 1000)))
        pylab.show()

if args.expression:
    a0 = sympy.simplify(sympy.trigsimp(  ((2 / T) * sum([sympy.integrate(s,     (t, expression[1][0], expression[1][1])) for expression in expressions])        )))
else:
    a0 = sympy.simplify(sympy.trigsimp(  ((2 / T) * sum([sympy.integrate(piece, (t, points[i][1],     points[i + 1][1])) for (i, piece) in enumerate(pieces)])  )))

print('\n\n\n\na0:\n')
sympy.pprint(a0, use_unicode=args.unicode)

if args.expression:
    an = sympy.simplify(sympy.trigsimp((  (2 / T) * sum([sympy.integrate(sympy.cos(n*((2 * sympy.pi) / T)*t)*s, (t, expression[1][0], expression[1][1])) for expression in expressions])  ).subs({T: (2 * sympy.pi) / w})))
else:
    an = sympy.simplify(sympy.trigsimp((  (2 / T) * sum([sympy.integrate(sympy.cos(n*((2 * sympy.pi) / T)*t)*piece, (t, points[i][1], points[i + 1][1])) for (i, piece) in enumerate(pieces)])  ).subs({T: (2 * sympy.pi) / w})))

print('\n\nan:\n')
sympy.pprint(an, use_unicode=args.unicode)

if args.expression:
    bn = sympy.simplify(sympy.trigsimp((  (2 / T) * sum([sympy.integrate(sympy.sin(n*((2 * sympy.pi) / T)*t)*s, (t, expression[1][0], expression[1][1])) for expression in expressions])  ).subs({T: (2 * sympy.pi) / w})))
else:
    bn = sympy.simplify(sympy.trigsimp((  (2 / T) * sum([sympy.integrate(sympy.sin(n*((2 * sympy.pi) / T)*t)*piece, (t, points[i][1], points[i + 1][1])) for (i, piece) in enumerate(pieces)])  ).subs({T: (2 * sympy.pi) / w})))

print('\n\nbn:\n')
sympy.pprint(bn, use_unicode=args.unicode)


sfourier = sympy.simplify(a0/2 + sympy.simplify(sympy.trigsimp(sympy.Sum(an * sympy.cos(n*((2 * sympy.pi) / T)*t) + bn * sympy.sin(n*((2 * sympy.pi) / T)*t), (n, 1, terms)))))

sfouriereven = sympy.simplify(a0/2 + sympy.simplify(sympy.trigsimp(sympy.Sum((an * sympy.cos(n*((2 * sympy.pi) / T)*t) + bn * sympy.sin(n*((2 * sympy.pi) / T)*t)).subs({n: (2 * n)}), (n, 1, terms)))))
sfourierodd = sympy.simplify(a0/2 + sympy.simplify(sympy.trigsimp(sympy.Sum((an * sympy.cos(n*((2 * sympy.pi) / T)*t) + bn * sympy.sin(n*((2 * sympy.pi) / T)*t)).subs({n: 2 * n + 1}), (n, 1, terms)))))

print('\n\ns_fourier:\n')
sympy.pprint(sfourier.subs({terms: sympy.oo}), use_unicode=args.unicode)

print('\n\ns_fourier_even:\n')
sympy.pprint(sfouriereven.subs({terms: sympy.oo}), use_unicode=args.unicode)

print('\n\ns_fourier_odd:\n')
sympy.pprint(sfourierodd.subs({terms: sympy.oo}), use_unicode=args.unicode)

if args.outgraph:
    print('\n\nfourier graph with', args.outgraph[0], ' terms over 3 Periods:\n')
    if args.textplot:
        sympy.plotting.textplot(sfourier.subs({symbol:1 for symbol in symbols}).subs({terms: args.outgraph[0]}), -1, 2, int(args.textplot[0]), int(args.textplot[1]))
    else:
        import pylab
        import matplotlib
        pylab.plot(matplotlib.numpy.linspace(-1, 2, 3000), sympy.lambdify(t, sfourier.subs({symbol:1 for symbol in symbols}).subs({terms: args.outgraph[0]}), modules=['numpy'])(matplotlib.numpy.linspace(-1, 2, 3000)))
        pylab.show()
    print()
