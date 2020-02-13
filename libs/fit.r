REBOL [
	title: {Functions for fitting data to them}
	copyright: {Bioservo Technologies AB}
	revision: {$ Revision: $}
	author: {Johan Ingvast}
]

prebol-used?: not (#do ['false])

either prebol-used? [
    #include %linalg.r
][
    do %linalg.r
]


fit!: make object! [
    type: 'fit!
    bases: copy []
    k: copy []
    residual: none
    estimate: none
    order: 3
    x: none
    x-range: [ 0 1 ]
    calc-estimate: func [
	{From x (in objet) and k calculate a new estimate}
    ][
	clear estimate
	insert/dup estimate 0 length? first bases
	bases: head bases
	forall k [ change estimate add-vector-vector estimate mult-vector-scalar first+ bases first k ]
	bases: head bases
	estimate: head estimate
    ]

    fit: func [ y-samples /local b kk] [
	y: y-samples
	clear k
	insert/dup estimate: copy [] 0 length? y ; make a row of zeroes
	foreach b bases [
	    append k kk: vector-mult b y
	    estimate: add-mult-vector-vector estimate 1 b kk
	]
	residual: add-mult-vector-vector y 1 estimate -1
    ]
    create-new-base: func [
	{Removes all the information from the seed found in bases such that
	 a new orthogonal base remains which is added.}
	b [block!] {The seed}
	/local b-w w
	
    ][
	    repeat j length? bases [ 
		 b-w: bases/:j
		 w: vector-mult b-w b
		 b: map-each t b [ t - ( w * first+ b-w ) ]
	    ]
	    b: mult-vector-scalar b 1 / square-root vector-mult b b
	    append/only bases b
	    b
    ]
    make-bases-sine: func [
	{Creates pairs of bases with cosine and sine
	    n = 1  cos( x * 2 * pi )	   sin( x * 2 * pi )
	    n = 2  cos( x * 2 * 2 *  pi )  sin( x * 2 *2 *  pi )
	           cos( x * n * 2 * pi )   sin( x * n * 2 * pi )
	    ...
	}
	/ord {Order other than order}
	    n
	/local x-start x-scale s c
    ][
	;make-bases-poly/ord 1
	x-start: x-range/1
	x-scale: 360 / ( x-range/2 - x-range/1 )
	append/only bases insert/dup copy [] 1 / square-root length? x length? x
	create-new-base x
	repeat i to-integer any [ n order ] [
	    c: map-each t x [ cosine t - x-start * x-scale * i ]
	    s: map-each t x [ sine   t - x-start * x-scale * i ]
	    create-new-base c
	    create-new-base s
	]
	new-line/all bases on
    ]
    make-bases-poly: func [
	/ord {Order other than order}
	    n
	/local fun j
    ][
	repeat i 1 + to-integer any [ n order ] [
	    j: i - 1
	    either j == 0 [
		fun: func [ x ] [
		    either x = 0 [ 1 ][ x ** j ] ; fix of bug in **
		]
	    ] [
		fun: func [x][ x ** j]
	    ]

	    b: map-each t x [ fun t ]
	    create-new-base b
	]
    ]
    integrate-estimate: func [
	{Integrate the function}
	/fun f [function!] {Integrate f(estimate) instead}
	/local
	    sum
	    y
    ][ sum: 0
	estimate: head estimate
	y: either fun [ map-each t estimate [ f t ] ][ estimate ]
	y: next y
	x: head x
	forall y [
	    sum: sum + (y/-1 + y/1 / 2 * ( x/2 - x/1))
	    x: next x 
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
	make error! {Slope of fit! not yet implemented}
    ]
]

test-fit: does [
    make fit! [
	order: 5
	x: copy []
	repeat t 101 [ append x t - 1 / 100 ]
	x-range: [0 1]
	y: map-each t x [ t ** 2]
	make-bases-sine
	fit y
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
	    root-mean-sq: 0
	    foreach t residual [ root-mean-sq: root-mean-sq + ( t * t) ]
	    root-mean-sq: square-root root-mean-sq
	    to-string: reform [
		    "k:" mold k newline
		    "order:" order newline
		    "root-mean-sq:" root-mean-sq
	    ]
	]
    ]
    linear-eq-solve AA Ay
]

trig-fit-old: func [ 
    {
    Adapts the data y-est = k/1 + k/2 * ( x - (range/1 + range/2) / 2 ) /(range/2 - range/1) * 2 +
			    k/3 * sin( x * omega) + k/4 * cos( x * omega ) ...
    so that sum ( y-est - y )^2 is minimized
    Returns a object
		k block
		A matrix
		est-y function for calculating y for an x
		residual
}
    x [block!] {x-data}
    y [block!] {y-data}
    orderi [integer!] { The number of sine cos pairs  }
    /range r [ block!] {Minimum and maximum for the lowest mode, if not given max and min of x is used}
    /half {Includes a base sine 0 - pi}
    /local
	Ai ki
	obj 
	x-vals y-vals
][
    make object! [
	A: copy []
	x-vals: x
	y-vals: y
	range: none
	x-scale: none
	x-start: none
	set-range: func [ r ] [
	    range: any [ r reduce [ first minimum-of x first maximum-of x ]]
	    x-start: range/1
	    x-scale: 180 / ( range/2 - range/1 )
	]
	x-scale: none
	set-range r
	k: none
	order: orderi
	bases-at: func [ x /local r] [
	    r: reduce [
		1
		2 * x - x-start - range/2 / 2
	    ]
	    ;if half [
		;append r sine  x - x-start * x-scale
	    ;]
	    repeat i order [
		repend r [
		    cosine x - x-start * i * x-scale
		    sine   x - x-start * i * x-scale
		]
	    ]
	    r
	]
	deriv-bases-at: func [ x /local r] [
	    r: reduce [
		0
		1
	    ]
	    ;if half [
		;append r pi / 180 * x-scale * cosine  x - x-start * x-scale
	    ;]
	    repeat i order [
		repend r [
		    negate i * pi / 180 * x-scale * sine x - x-start * i * x-scale
		    i * pi / 180 * x-scale * cosine x - x-start * i * x-scale
		]
	    ]
	    r
	]
	val-at: func [ x ] [ vector-mult bases-at x k ]
	deriv-at: func [ x ] [ vector-mult deriv-bases-at x k ]

	create-bases: does [
	    clear A
	    append/only A head insert/dup copy [] 1 length? x-vals ; the constant added
	    append/only A map-each t x-vals [ 2 * t - x-start - range/2 / 2 ]
	    repeat i order [
		append/only A map-each t x-vals [ cosine t - x-start * i * x-scale ]
		append/only A map-each t x-vals [ sine   t - x-start * i * x-scale ]
	    ]
	    ;print [ "First-A" length? A/1 "x-vals" length? x-vals ]
	]
	create-bases
	fit: does [
	    k: linear-eq-solve
		    mult-matrix-matrix/B-transposed A A
		    mult-matrix-matrix/B-transposed A reduce [ y-vals ] ] fit
	estimate: func [  ] [
	    first mult-matrix-matrix reduce [k] A
	]
	residual: does [ mult-add-vectors estimate  -1   y-vals 1 ]

	root-mean-sq: func [ /local s ] [
	    s: 0
	    foreach t residual [ s: s + ( t * t) ]
	    square-root s
	]
	to-string: does [
	    reform [
		"k:" mold k newline
		"order:" order newline
		"root-mean-sq:" root-mean-sq
	    ]
	]
	extremas: func [ /local m ]  [ reduce [first minimum-of m: estimate first maximum-of m ] ]
	mean: does [ ( vector-sum estimate ) / length? estimate ]
	integrate-estimate: func [
	    {Integrate the function}
	    /fun f [function!] {Integrate f(estimate) instead}
	    /local
		sum
		x y
		est
	][
	    sum: 0
	    est: head estimate
	    y: either fun [ map-each t est [ f t ] ][ est ]
	    y: next y
	    x: head x-vals
	    forall y [
		sum: sum + (y/-1 + y/1 / 2 * ( x/2 - x/1))
		x: next x 
	    ]
	    sum
	]
    ]
]

trig-fit: func [ 
    {
    Adapts the data y-est = k/1 + k/2 * ( x - (range/1 + range/2) / 2 ) /(range/2 - range/1) * 2 +
			    k/3 * sin( x * omega) + k/4 * cos( x * omega ) ...
    so that sum ( y-est - y )^2 is minimized
    Returns a object
		k block
		A matrix
		est-y function for calculating y for an x
		residual
}
    x [block!] {x-data}
    y [block!] {y-data}
    orderi [integer!] { The number of sine cos pairs  }
    /range r [ block!] {Minimum and maximum for the lowest mode, if not given max and min of x is used}
    /local
	Ai ki
	obj 
	x-vals y-vals
][
    make object! [
	type: 'trig-fit
	A: copy []
	x-vals: x
	y-vals: y
	range: none
	x-scale: none
	x-start: none
	set-range: func [ r ] [
	    range: any [ r reduce [ first minimum-of x first maximum-of x ]]
	    x-start: range/1
	    x-scale: 180 / ( range/2 - range/1 )
	]
	x-scale: none
	set-range r
	k: none
	order: orderi
	bases-at: func [ x /local r] [
	    r: reduce [
		1
		2 * x - x-start - range/2 / 2
	    ]
	    ;if half [
		;append r sine  x - x-start * x-scale
	    ;]
	    repeat i order [
		repend r [
		    cosine x - x-start * i * x-scale
		    sine   x - x-start * i * x-scale
		]
	    ]
	    r
	]
	deriv-bases-at: func [ x /local r] [
	    r: reduce [
		0
		1
	    ]
	    ;if half [
		;append r pi / 180 * x-scale * cosine  x - x-start * x-scale
	    ;]
	    repeat i order [
		repend r [
		    negate i * pi / 180 * x-scale * sine x - x-start * i * x-scale
		    i * pi / 180 * x-scale * cosine x - x-start * i * x-scale
		]
	    ]
	    r
	]
	val-at: func [ x ] [ vector-mult bases-at x k ]
	deriv-at: func [ x ] [ vector-mult deriv-bases-at x k ]

	create-bases: does [
	    clear A
	    append/only A head insert/dup copy [] 1 length? x-vals ; the constant added
	    append/only A map-each t x-vals [ 2 * t - x-start - range/2 / 2 ]
	    repeat i order [
		append/only A map-each t x-vals [ cosine t - x-start * i * x-scale ]
		append/only A map-each t x-vals [ sine   t - x-start * i * x-scale ]
	    ]
	    ;print [ "First-A" length? A/1 "x-vals" length? x-vals ]
	]
	create-bases
	;A-lin: A-trig: kt: kl: none
	fit: func [
	    /local A-lin A-trig
	    kl kt
	][
	    dbg: self
	    ; do the fitting in two steps, first do it by removing mean and linear.
	    ; next the trigonometric.
	    A-lin: copy/part A 2
	    A-trig: copy skip A 2
	    kl: linear-eq-solve
		    mult-matrix-matrix/B-transposed A-lin A-lin
		    mult-matrix-matrix/B-transposed A-lin reduce [ y-vals ]
	    if order > 0 [
		kt: linear-eq-solve
			mult-matrix-matrix/B-transposed A-trig A-trig
			mult-matrix-matrix/B-transposed A-trig
			    reduce [
				mult-add-vectors
				    y-vals					1
				    first mult-matrix-matrix reduce [ kl ] A-lin    -1
			    ]
	    ]
	    k: append kl any [ kt [] ]
	]
	fit
	estimate: func [  ] [
	    first mult-matrix-matrix reduce [k] A
	]
	estimate-deriv: func [
	    /local ret
	][
	    ret: copy []
	    repeat i length? x-vals [
		append ret deriv-at pick x-vals i
	    ]
	]
	residual: does [ mult-add-vectors estimate  -1   y-vals 1 ]

	root-mean-sq: func [ /local s ] [
	    s: 0
	    foreach t residual [ s: s + ( t * t) ]
	    square-root s / length? y-vals
	]
	to-string: does [
	    reform [
		"k:" mold k newline
		"order:" order newline
		"root-mean-sq:" root-mean-sq
	    ]
	]
	extremas: func [ /deriv {Extreme derivative} /local m
	]  [
	    reduce [
		first minimum-of m: either deriv [ estimate-deriv] [estimate]
		first maximum-of m
	    ]
	]
	mean: does [ ( vector-sum estimate ) / length? estimate ]
	angle: does [ arctangent k/2 ]
	integrate-estimate: func [
	    {Integrate the function}
	    /fun f [function!] {Integrate f(estimate) instead}
	    /local
		sum
		x y
		est
	][
	    sum: 0
	    est: head estimate
	    y: either fun [ map-each t est [ f t ] ][ est ]
	    y: next y
	    x: head x-vals
	    forall y [
		sum: sum + (y/-1 + y/1 / 2 * ( x/2 - x/1))
		x: next x 
	    ]
	    sum
	]
	cut-part: func [
	    {Cuts out the part of the data where
		the first x data >= start-x and
		the last x data <= end-x
	     returns a new object
	     }
	     start-x [number!]
	     end-x [number!]
	     /local start end
	][
	    while [ x-vals/1 < start-x ] [ x-vals: next x-vals ] start: index? x-vals
	    while [ all [ x-vals/1 x-vals/1 <= end-x ] ] [ x-vals: next x-vals] end: index? x-vals
	    x-vals: head x-vals
	    x-vals: copy/part at x-vals start at x-vals end
	    y-vals: copy/part at y-vals start at y-vals end
	]
    ]
]

test-trig-fit: func [ 
][
    x: copy [] repeat i 101 [ append x i - 1 ]
    y: map-each t x [ ( sine  3 * a: t *  360 / 100)  + ( 1.5 * cosine a) + cosine 3 * a ]

    trig-fit/half x y 3
]

lin-interp: func [
    x-data [ block! ] {Data vector of x-values, must be continously growing}
    y-data [ block! ] {Data vector of y-values}
    x [ number! ] {The x-value to find }
    /local 
][
    xi: back tail x-data
    forall x-data [
	if x < (first x-data) [
	    xi: x-data first xi break
	]
    ]
    x1: back xi x2: next x1
    y1: pick y-data index? x1
    y2: pick y-data index? x2
    x1: first x1 x2: first x2
    ( y2 - y1 ) / ( x2 - x1 ) * ( x - x1 ) + y1
]


; vim: syntax=rebol ai sw=4 sts=4
