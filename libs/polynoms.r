REBOL [
	title: {Functions for making polynoms and fitting data to them}
	copyright: {Bioservo Technologies AB}
	revision: {$ Revision: $}
	author: {Johan Ingvast}
]

make-orth-poly: func [
    {Fits a polynom to the data.  Orthogonal bases are used in rising order of the polynom.
     To see the result from one of the bases, set the corresponding k value to 0 and 
     make a calc-estimate which recalculates the estimate and the residual.
     The bases are scaled such that the sum the square of each value in the list is 1.
     }
    x [block!] {x-data}
    o [integer!] {Order of the estimate}
    /local b
] [
    make object! [
	type: 'make-orth-poly
	bases: copy []
	k: copy []
	residual: none
	estimate: none
	order: o
	x-vals: x
	calc-estimate: func [
	    {From x-vals (in objet) and k calculate a new estimate}
	][
	    clear estimate
	    insert/dup estimate 0 length? first bases
	    bases: head bases
	    forall k [ change estimate add-vector-vector estimate mult-vector-scalar first+ bases first k ]
	    bases: head bases
	    estimate: head estimate
	]
	fit: func [ y /local b kk] [
	    clear k
	    insert/dup estimate: copy [] 0 length? y ; make a row of zeroes
	    foreach b bases [
		append k kk: vector-mult b y
		estimate: add-mult-vector-vector estimate 1 b kk
	    ]
	    residual: add-mult-vector-vector y 1 estimate -1
	]
	make-bases: func [
	    /local fun j
	][
	    repeat i 1 + to-integer order [
		j: i - 1
		either j == 0 [
		    fun: func [ x ] [
			either x = 0 [ 1 ][ x ** j ] ; fix of bug in **
		    ]
		] [
		    fun: func [x][ x ** j]
		]
		;b: map-each t x-vals [ t ** ( i - 1 )  ]
		b: map-each t x-vals [ fun t ]
		repeat j i - 1 [ 
		     b-w: bases/:j
		     w: vector-mult b-w b
		     b: map-each t b [ t - ( w * first+ b-w ) ]
		]
		b: mult-vector-scalar b 1 / square-root vector-mult b b
		append/only bases b
	    ]
	]
	integrate-estimate: func [
	    {Integrate the function}
	    /fun f [function!] {Integrate f(estimate) instead}
	    /local
		sum
		y
	][
	    sum: 0
	    estimate: head estimate
	    y: either fun [ map-each t estimate [ f t ] ][ estimate ]
	    y: next y
	    x-vals: head x-vals
	    forall y [
		sum: sum + (y/-1 + y/1 / 2 * ( x-vals/2 - x-vals/1))
		x-vals: next x-vals 
	    ]
	    sum
	]
	extremas: func [ ][
	    reduce [ first minimum-of estimate first maximum-of estimate ]
	]
	mean: does [
	    k/1 * bases/1/1
	]
	slope: does [
	    k/2 * ( ( last bases/2 ) - first bases/2 ) / ( ( last x-vals ) - first x-vals )
	]
	make-bases
    ]
]

pol-diff: func [ 
    {A polynom k is differentiated once and a new polynom returned}
    k
    /local kp
][
    kp: make block! length? k
    repeat i -1 + length? k [
	append kp k/(i + 1) * i
    ]
    kp
]

pol-to-y: func [ k x /local A ] [
    if number? x [ x: reduce [x] ]
    A: copy []
    for i 0 1 + length? k 1 [
	append/only A map-each t x [ t ** i ]
    ]
    first mult-matrix-matrix reduce [k] A 
]

pol-fit: func [ 
    {
    Adapts the data y-est = k/1 + k/2 * x + k/3 * x^2 so that
    sum ( y-est - y )^2 is minimized
    It returns a block with k constants
    }
    x [block!] {x-data}
    y [block!] {y-data}
    order [integer!] { The order of adaptation, 0 - mean 1-slope 2-square}
    /object {Returns a object:
		k block
		A matrix
		est-y function for calculating y for an x
		residual
		}

    /local
	Ai ki
	obj orderi
][
    {
	y = A k  where A are columns of [ 1 x x**2 ... ]
	A'*y = A'*A k which can be solved.
    }
    Ai: copy []
    for i 0 order 1 [
	append/only Ai map-each t x [ t ** i ]
    ]
    AA: mult-matrix-matrix/B-transposed Ai Ai
    Ay: mult-matrix-matrix/B-transposed Ai reduce [ y ]
    if object [
	orderi: order
	return make object! [
	    type: 'pol-fit
	    k: linear-eq-solve AA Ay
	    order: orderi
	    est-y: func [ x /local A ] [
		if number? x [ x: reduce [x] ]
		A: copy []
		for i 0 order 1 [
		    append/only A map-each t x [ t ** i ]
		]
		first mult-matrix-matrix reduce [k] A 
	    ]
	    residual: mult-add-vectors est-y x 1   y -1
	    root-min-sq: 0
	    foreach t residual [ root-min-sq: root-min-sq + ( t * t) ]
	    root-min-sq: square-root root-min-sq
	    to-string: reform [
		    "k:" mold k newline
		    "order:" order newline
		    "root-min-sq:" root-min-sq
	    ]
	]
    ]
    linear-eq-solve AA Ay
]
; vim: ai sw=4 sts=4
