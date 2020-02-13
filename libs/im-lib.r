REBOL [
	title: {Basic functions for the image analysis}
	author: {Johan Ingvast}
]

scale-image: func [
    {Scales the image with the factor scl}
    img
    scl
    /local 
] [
    to-image make face [
	edge: font: para: none
	image: none
	color: white
	size: img/size * scl
	effect: [
	    draw [
		scale scl scl image img
	    ]
	]
    ]
]
crop-image: func [
    image
    pos
    /local p1 p2
] [
    p1: as-pair max 0 min pos/1/image-pos/1 pos/2/image-pos/1
		max 0 min pos/1/image-pos/2 pos/2/image-pos/2
    p2: as-pair min image/size/x max pos/1/image-pos/1 pos/2/image-pos/1
		min image/size/y max pos/1/image-pos/2 pos/2/image-pos/2
    image-bounding-box/1: image-bounding-box/1 + p1
    image-bounding-box/2: image-bounding-box/2 + p2
    copy/part at image p1 p2 - p1
]

rot90-image: func [
    img
    /rot angle
] [
    angle: any [ angle 90 ]
    if angle < 0 [ angle: 360 + angle ]
    while [ angle > 360 ] [ angle: angle - 360 ]
    to-image make face [
	edge: font: para: none
	image: img
	size: either find [ 0 180 ] angle [
	    img/size
	][
	    reverse img/size
	]
	effect: [rotate angle ]
    ]
]

rot-image: func [
    img
    angle
    /local tr1 tr2
] [
    tr1: img/size / 2
    tr2: tr1 * -1
    to-image make face [
	edge: font: para: none
	image: none
	color: white
	size: img/size
	effect: [
	    draw [
		translate tr1 rotate angle translate tr2 image img
	    ]
	]
    ]
]

lowpass1: func [
    {Filters a one dimensional array with the function y_n+1 = y_n * ( 1 - factor ) + factor * u_n+1
    Start value is u_1}
    u [block!] {Dataset}
    factor [number!] {The factor}
    /local y ret
][
    ret: make block! length? u
    y: u/1
    fact*: 1 - factor
    append ret y
    foreach uu next u [
	append ret y: y * fact* + ( factor * uu )
    ]
]

filt-filt: func [ 
    {Runs a fiter from start to end and end to start of the dataset}
    u [block!] {Data to filter}
    factor [ number!] {Factor to filter with}
    /filter filt {The filter to use if other than first order, in this case factor is ignored}
    /local
][
    filt: any [ 
	filt
	func [ x ][ lowpass1 x factor ]
    ]
    reverse filt reverse filt u 
]

; vim: ai sw=4 sts=4
