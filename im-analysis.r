#!/usr/bin/rebview  -sv
REBOL [
    title: {
	    Analysis of threaded bar straightness.
	    Starts with an image of the bar.
    }
    encap: [ title "Screw analysis - " ]
    copyright: {Bioservo Technologies AB}
    revision: {$Revision: 5999 $}

    TODO: {Make it possible to put limits on the diameters}
    TODO: {Settings controls the scale on diameter}
    TODO: {Change zoom with scroll button}
    TODO: {State window in table}

    DONE: {Parameter interval where measurements are done}
    DONE: {Give +- g�ngstigning}
    DONE: {Remove first and last thread data}
    DONE: {Analyze pitch and give variations}
    DONE: {Visa diametrar}
    DONE: {Rotate 90 degrees}
    DONE: {Make possible to switch between different tools}
    DONE: {Filter away small dots in the image}
    DONE: {Make possible to filter the image}
    DONE: {Make possible to make fine adjustment in scale (from settings.r)}
    DONE: {Include a checkbox for automatic scale estimation}
    DONE: {In case small scale adjustment, show it}
    DONE: {Put scales into the diameter plot}
    DONE: {Move position window into the image view}
    DONE: {Fix problem with non-visible rectangle when cropping}
]


; In linux, switch to arial as default font
font-name: "arial"
if system/version/4 = 4 [ 
    face/font/name: font-name
    foreach [n f ] svv/vid-styles
    [
	if all[ f/font f/font/name = "helvetica" ] [ f/font/name: font-name ]
    ]
]

prebol-used?: not (#do ['false])

either prebol-used? [
    #include %../../source/view.r
    #include %libs/linalg.r
    #include %libs/polynoms.r
    #include %libs/fit.r
    #include %libs/im-lib.r
    #include %libs/menue.r
][
    lib-path: %libs
    do lib-path/linalg.r
    do lib-path/polynoms.r
    do lib-path/fit.r
    do lib-path/im-lib.r
    do lib-path/menue.r
    do/args %~/misc-proj/pdf-export/face-to-pdf-lib.r  'face-to-pdf
]


; Default settings

do default-settings: {
; Here we set variables for im-analysis.r

; A list of possible resolutions (dpi) to select from
possible-resolutions: [ 300 600 1200 2400 3200 4800 ] ; dpi

; Corrects deviations from the exact dpi setting. 
; A value larger than 1 means the resolution is larger than the dpi setting
resolution-correction-factor-x: 0.998
resolution-correction-factor-y: 0.998

; The image that is opened at startup.  Set to none for fast startup
default-image: none

; The expected value of the pitch
ideal-pitch: 0.25 ; mm

; The expected lenght of the screw, used to estimate the resolution
approx-screw-length: 100 ; mm

gray-threshold: 0.8 * 255

;The order of the polynom to fit. The higher, the closer to the acutal calculated tops and ridges
polynom-order: 9

; The radial distance one thread top can differ from next without being treated as erreous.
max-thread-diff: ideal-pitch / 4

; The distance from left/right end that is ignored in calculating statistical values
skip-end-statistical-average: 5 ; mm

effective-diameter: 1 - ( ideal-pitch * cosine 30 )
}

settings-file: %settings.r

unless info? settings-file [
    write settings-file default-settings
]

do either (get in info? settings-file 'type ) = 'file [
    read settings-file ] [
    print "Warning! No settings file found, using default!"
    default-settings
] 

if  all [ not value? 'im value? 'default-image file? default-image ] [
    im: load image-name: default-image
    image-bounding-box: reduce [ 0x0 image-name/size ]
]
if any [ not value? 'im  not image? im ] [ im: make image! [ 1x1 ] ]

pix-per-mm-x: 2400 / 25.4
pix-per-mm-y: 2400 / 25.4

make-screw-image: func [
    L {Length of screw [mm]}
    diam [ number! function!] {diameter of screw. Function or constant [mm]}
    pitch {pitch of screw [mm]}
    angle {Orientation}
    deflection [number! function!]  {Function or constant of how screw deflections over the length [mm]}
    /local
	size-y im ym grove-depth x-p diam-c r-effective dy p 
] [
    size-y: 10 * pix-per-mm-y * diam L / 2
    im: make image! reduce [ as-pair L + 10 * pix-per-mm-x size-y white ]

    ; Make blobs at the ends
    change/dup
	im
	black
	as-pair 5 * pix-per-mm-x size-y
    change/dup
	at im as-pair L + 5 * pix-per-mm-x 0
	black 
	as-pair 5 * pix-per-mm-x size-y
    
    ym: size-y / 2
    grove-depth: pitch / 2 * cosine 30
    for i 5 * pix-per-mm-x L + 5 * pix-per-mm-x 1 [
	x-p: i / pix-per-mm-x - 5 

	diam-c: diam  x-p
	r-effective: diam-c - grove-depth / 2 * pix-per-mm-y 
	dy: deflection x-p
	p: pix-per-mm-y *  ( 
	    dy
	    + ( pitch / 4 * sine x-p / pitch * 360 )
	    + ( L / 2 - x-p * sine angle )
	)
	change/dup at im as-pair i ym - r-effective + p  black as-pair 1 2 * r-effective
    ]
]

image-threshold: func [ im t 
    /local i
][
    i: make face [
	image: im
	size: im/size
	color: green
	effect: compose [
	    grayscale
	    luma ( 1 - t ) 
	    luma 254
	    key 254.254.254
	]
	edge: font: para: none
    ]
    to-image i
]


thread-filter: func [
    {Convolutes the data with a tent function with mean value zero.
     it us used to filter out the expected threads
     It returns a new vector of the same lenth as the data.
     }
    data [block!] {of data}
    width [integer!] {The width of the tent}
    /local base filtered
] [

    if even? width [ width: width - 1]

    base: reduce [ -1 + (1 / width) ]
    loop (to-integer width / 2) [ append base 2 / ( to-integer width / 2) + last base ]
    append base next reverse copy base ; mirror it

    filtered: copy []
    loop ( length? data ) - width + 1 [
	sum: 0
	base: head base
	repeat i width [ sum: sum + ( data/:i * first+ base ) ]
	append filtered sum
	data: next data
    ]
    insert/dup filtered 0 to-integer width / 2
    insert/dup tail filtered 0 to-integer width / 2
    filtered
]

tools: stylize [
    grouped-panel: panel
	edge [ size: 1x1 effect: 'ibevel ]
	font [ valign: 'top align: 'left name: "arial" style: 'bold shadow: none color: black]
	"Panel"
    rec-button: toggle "rec" "wait" [
	if word? face/draw-box [ face/draw-box: get face/draw-box ]
	use [ button-face ] [
	    button-face: face
	    clear button-face/picked-list
	    face/draw-box/feel/insert-engage 
		any [ face/var 'rec-button ]
		func [ face action event ] [
		    switch action [
			down [
			    append button-face/picked-list event/offset
			    append button-face/picked-list event/offset
			    append clear button-face/draw [
				    fill-pen button-face/fill-pen pen button-face/pen
				    box button-face/picked-list/1 button-face/picked-list/2 ]
			    button-face/down face event button-face/picked-list
			]
			up [
			    clear button-face/draw
			    face/feel/remove-engage any [ face/var 'rec-button ]
			    button-face/state: false 
			    show [ face button-face ]

			    button-face/up face event button-face/picked-list
			]
			over [
			    unless empty? button-face/picked-list  [
				button-face/picked-list/2: event/offset
				show face
			    ]
			    button-face/over face event button-face/picked-list
			]
		    ]
		    false
		]
	]
    ]
    [
	clear face/draw
	face/feel/remove-engage any [ face/var 'rec-button ]
	show face/draw-box
    ]
    with [
	draw-box: {Must be set to the face where to draw the rectangle}
	picked-list: copy []
	draw: copy []
	up: none
	down: none
	over: none
	fill-pen: 0.0.255.100
	pen: none
	line-width: none
    ]

    select-image-area: toggle "Select-area" "wait ..." [
	face/butt-press face ; Pre hook, can be used to set info bar
	if word? face/draw-box [ face/draw-box: get face/draw-box ]
	use [ button-face ] [
	    button-face: face
	    clear face/picked-list
	    clear face/pos-list
	    face/draw-box/feel/insert-engage any [ button-face/var 'select-image-area] ; name 
		func [ face action event ] [   ; function
		    switch action [
			down [
			    append button-face/picked-list event/offset
			    append button-face/picked-list event/offset
			    append clear button-face/draw [
				    fill-pen button-face/fill-pen
				    pen button-face/pen
				    line-width button-face/line-width
				    box button-face/picked-list/1 button-face/picked-list/2 ]
			    append button-face/pos-list make button-face/transformation [ set-window-pos event/offset ]
			    button-face/down face event button-face/pos-list
			]
			up [
			    clear button-face/draw
			    remove/part find face/feel/engage-list any [button-face/var 'select-image-area ] 2
			    button-face/state: false 
			    show [ face button-face ]
			    append button-face/pos-list make button-face/transformation [ set-window-pos event/offset ]
			    button-face/up face event button-face/pos-list
			]
			over [
			    unless empty? button-face/picked-list  [
				button-face/picked-list/2: event/offset
				show face
			    ]
			    button-face/over face event button-face/pos-list
			]
		    ]
		    false
		]
	    ]
	] [
	    clear face/draw
	    face/feel/remove-engage any [ button-face/var 'select-image-area]
	    show face/draw-box
	]
	with [
	    draw-box: {Must be set to the face where to draw the rectangle}
	    picked-list: copy []
	    pos-list: []
	    transformation: none
	    draw: []
	    up: down: over: none
	    butt-press: none
	    fill-pen: 0.0.255.100
	    pen: none
	    line-width: none
	]

    event-box: box 
	feel [ 
	    info: {Give the funtion an argument and add it to the right list with a given name.
	    A function that returns a value continues to evaluate the list actions.}
	    over-list: []
	    engage-list: []
	    insert-action: func [ name [word!] f [function!] list ] [ insert list reduce [name :f] ]
	    append-action: func [ name [word!] f [function!] list ]  [append list reduce [name :f] ]
	    remove-action: func [ name list ] [ remove/part find list name 2 ] 

	    insert-engage: func [ name f ][ insert-action name :f engage-list ]
	    insert-over: func [ name f ][ insert-action name :f over-list ]

	    append-engage: func [ name f ][ append-action name :f engage-list ]
	    append-over: func [ name f ][ append-action name :f over-list ]

	    remove-engage: func [ name ] [ remove-action name engage-list ]
	    remove-over: func [ name ] [ remove-action name over-list ]

	    engage: func [ face action event ][
	       foreach [name e] engage-list [ unless e face action event [break] ]
	    ]
	    over: func [ face action event ][
	       foreach [name e] over-list [ unless e face action event [break] ]
	    ]
	]
    leave-key: key #"^q" [unview]
]


data-coll: make svv/vid-face [
    size: 200x100
    row-offset: 0
    x-offset: 0 
    edge: none
    num-cols:  0
    num-rows: 0
    row-spaceing: 20
    grid: on
    effect: none
    data:  []
    pane: func [ 
	face
	index
	/local i data-row
    ][
	if integer? index [ 
	    data-row: index + row-offset
	    ;if probe any [ ( probe row-offset * ( row/size/y + 1) ) > size/y  data-row > num-rows ][ return none ]
	    if data-row > num-rows [ return none ]
	    row/offset: as-pair x-offset  index - 1 * ( row/size/y + 1)
	    repeat i num-cols  [
		row/pane/:i/text: data/:data-row/:i
		;row/pane/:i/text: as-pair i index
	    ]
	    return row
	]
	i: to-integer index/y / (row/size/y + 1) + 1
    ]
    cell: make svv/vid-face [
	size: 60x19
	feel: make svv/vid-styles/image/feel [
	    engage: func [ f a e ][
		if a = 'up [ f/parent-face/parent-face/action f f/text ]
	    ]
	] 
	font: make face/font [ halign: 'center valign: 'top style: none]
	para: make face/para [ wrap?: true]
    ]

    row: make face [
	pane: []
	edge: none
    ]
    words: reduce [ 'data func [ new args ] [ new/data: args/2 next args] ]

    init: [ 
	row: make row [] 
	row/parent-face: self
	cell: make cell [ parent-face: row ]
	update
    ]
    check-data: does [
	num-rows: length? data
	num-cols: any [ all [ data/1 length? data/1 ] 0 ]
	;forall data [
	 ;   unless num-cols = length? first data [ make error! reform [ {Not consistent row lengths in data} new-line/all data on ] ]
	;]
    ]
    update: does [
	check-data
	row/size: as-pair num-cols * cell/size/x cell/size/y
	clear row/pane
	repeat i num-cols [ append row/pane make cell [ offset: as-pair i - 1 * size/x 0 ] ]
	row/size: as-pair cell/size/x * num-cols cell/size/y
	effect: all [ grid reduce [ 'grid as-pair 0 row-spaceing as-pair 0 row-spaceing - 1 ] ]
    ]
    remove-col: func [ n ] [
	forall data [ remove at first data n ]
    ]
    remove-row: func [ n ] [
	remove at data n
    ]
]

table: make svv/vid-face [
    size: 200x100
    edge: none
    font: make face/font []
    row-spaceing: 20
    effect: reduce [ 'grid as-pair 0 row-spaceing as-pair 0 row-spaceing - 1 ]

    data: [[]]

    h-scroll: v-scroll: header: left-coll: data-table: none
    scroller-factor-x: none
    scroller-factor-y: none
    pane: []
    init: [
	data-table: make data-coll [
	    grid: off
	    color: none
	]
	left-coll: make data-coll [
	    type: 'left-coll
	    size/x: 50
	    offset: 0x0
	    data: []
	    cell: make cell [
		font: make font [ style: 'bold ]
	    ]
	    color: none
	    para: make para [ tabs: 10 ]
	]
	header: make data-coll [
	    type: 'header
	    size/y: row-spaceing 
	    data: [[]]
	    cell: make cell [
		font: make font [ style: 'bold ]
	    ]
	    color: none
	]
	h-scroll: make svv/vid-styles/scroller [
	    size: 100x15 
	    action: func [ face value ] [
		data-table/x-offset: header/x-offset: negate value * scroller-factor-x
		show [ data-table header ]
	    ]
	]
	v-scroll: make svv/vid-styles/scroller [
	    size: 15x100
	    step: 0.2
	    ratio: 0.4
	    action: func [ face value ] [
		data-table/row-offset: left-coll/row-offset: round value * scroller-factor-y 
		show [ data-table left-coll ]
	    ]
	]

	do h-scroll/init h-scroll/init: none
	do v-scroll/init v-scroll/init: none

	make object! [
	]
	update

	repend pane [ data-table left-coll header ]
	foreach p pane [ do p/init p/init: none ]
	update
    ]

    append-row: func [ row /header txt ][
	append/only data row
	all [ header repend/only left-coll/data [ txt ] ]
	update
	show self
    ]

    append-col: func [ col /header txt ][
	forall data [
	    append first data first+ col
	]
	all [ header append self/header/data/1 txt ]
	update
	show self
    ]
    remove-col: func [ index /local ][
	remove at header/data/1 index
	forall data [ remove at first data index ]
	update
	show self
    ]

    last-appended-row: 0
    append-data: func [ d /below {default} /across ] [
	unless across [ below: true ]
	if block? d [
	    forall d [
		either below [
		    append-data first d
		][
		    append-data/across first d
		]
	    ]
	    exit
	]
	data: any [ data copy []]
	if empty? data [ append/only [] last-appended-row: 0 ]
	if below [ last-appended-row: last-appended-row + 1 ]
	if last-appended-row > length? data [ append/only data copy [] ]
	append data/(last-appended-row) d
    ]
    data-return: does [ last-appended-row: 0 ]

    update-size: does [

	either 0 > scroller-factor-x: data-table/row/size/x - data-table/size/x [
	    ; don't show
	    remove find pane h-scroll 
	    data-table/size/y:  size/y - header/size/y

	][ ; show
	    unless find pane h-scroll [ append pane h-scroll ]
	    data-table/size/y: size/y - header/size/y - h-scroll/size/y 
	]

	either 0 >= scroller-factor-y: data-table/num-rows - ( data-table/size/y / data-table/row-spaceing )  [
	    remove find pane v-scroll
	    data-table/size/x: size/x - left-coll/size/x
	][
	    unless find pane v-scroll [ append pane v-scroll ]
	    data-table/size/x: size/x - left-coll/size/x - v-scroll/size/x   
	]

	left-coll/size/y: data-table/size/y
	header/size/x: data-table/size/x
	h-scroll/size/x: data-table/size/x
	v-scroll/size/y: data-table/size/y

	data-table/offset: as-pair left-coll/size/x row-spaceing
	left-coll/offset: as-pair 0 header/size/y
	header/offset: as-pair left-coll/size/x 0
	h-scroll/offset: as-pair left-coll/size/x size/y - h-scroll/size/y
	v-scroll/offset: as-pair size/x - v-scroll/size/x header/size/y

	h-scroll/page:  data-table/cell/size/x / scroller-factor-x 

	v-scroll/resize none
	h-scroll/resize none

    ]

    update: does [
	for i 1 + length? left-coll/data length? data 1 [ append/only left-coll/data reduce [ i ] ]
	header/data
	for i 1 + length? header/data/1 length? data/1 1 [ append header/data/1 to-string now/time ]

	data-table/data: data
	data-table/cell/font: font
	foreach p reduce [ data-table left-coll header ] [ p/update ]
	update-size
    ]

    words: reduce [
	'data func [ new args ] [ new/data: args/2 next args]
	'width func [ new args ] [
	    append new/init/object! compose [
		new/data-table/cell/size/x: new/header/cell/size/x: (args/2)
	    ]
	    next args
	]
	'left-col 'header func [
	    new args
	    /local fc
	] [
	    append new/init/object! compose/deep [
		fc: get in new select [ left-col left-coll header header ] (to-lit-word first args)
	    ]

	    foreach x second args [
		case [
		    find [ left center right ] x [
			append new/init/object! new-line compose/deep [ set-font fc/cell align (to-lit-word x) ] on
		    ]
		    find [ top middle bottom ] x [
			append new/init/object! new-line compose/deep [ set-font fc/cell valign (to-lit-word x) ] on
		    ]
		    find [ bold italic underline ] x [
			append new/init/object! new-line compose/deep [ set-font fc/cell style (to-lit-word x) ] on
		    ]
		    integer? x [ ; column-width 
			append new/init/object! new-line compose [ fc/cell/size/x: fc/size/x: (x) ] on
		    ]
		    block? x [
			make error! {Function for left column not implemented}
		    ]
		    true [
			switch first args [
			    left-col	[ append new/init/object! new-line compose/deep [ repend/only fc/data [ (x) ] ] on ]
			    header	[ append new/init/object! new-line compose [ append fc/data/1 (x) ] on ]
			]
		    ]
		]
	    ]
	    next args
	]
    ]
    to-array: has [ row arr ][
	arr: copy []
	append/only arr row: copy []
	append row copy ""

	append row header/data/1

	row-head: left-coll/data
	foreach r data [
	    row: copy []
	    append/only arr row

	    append row first first+ row-head
	    append row r
	]
	new-line/all arr on
    ]

    to-string*: func [ /transpose /local arr str ] [
	str: copy {}
	arr: to-array
	if transpose [ arr: transpose-matrix arr ]

	foreach r arr [
	    append str mold first+ r
	    foreach c r [
		append str tab
		append str mold c
	    ]
	    append str newline
	]
	str
    ]
]

repend tools [ 'data-coll data-coll 'table table ]

test-table: does [
    view layout [
	styles tools
	t: table 200x200
	    data [
		[12 343 33]
		["lkj" #"j" 333]
		["Kalle" "Johan" "Fredrik" ]
	    ]
	    underline
	    header [ left "first" "second" "third" ]
	    left-col [ 70 right "johan" "bosse" "lalal" "wolla" underline ]
	button "Append row" [
	    t/append-row/header head insert/dup copy [] "newrow" length? t/data/1 "new head"
	]
	button "Append column" [
	    t/append-col/header head insert/dup copy [] "newcol" length? t/data "new head"
	]
	;em: table 400x150 left-col [ center "one" "two" "three" ]
	;button "Set data" [ append em/data [ [ "ole" "dole" "doff"]["kinke" "lane" "doff"]["binke" "bane" "poff" ]] em/update show em ]
	key #"q" [unview]
    ]
]

draw-zoom: make object! [
    info: {
	A object for handling zooming of the content in draw.
	It is controlled with mouse or buttons. See the example "test" method.
	}
    translate: func [
	offset [pair! block!]
    ] [
	T/5: T/5 + offset/1
	T/6: T/6 + offset/2
	update-Txy
    ]

    update-Txy: does [
	T-x/1: T/1
	T-x/5: T/5
	T-y/4: T/4
	T-y/6: T/6
    ]

    zoom: func [
	offset [pair! block!]
	factor [number!]
    ][
    {
       [s1 0 x1 ]   [s2 0 x2 ]  [s1*s2  0 s1*x2+x1 ]
       [0 s1 y1 ]   [0 s2 y2 ] =[0  s1*s2 s1*y2+y1 ]
       [0 0   1 ]   [0 0   1 ]  [0 0   1 ]
    }

	T/5: factor * T/5 + ( offset/1 * ( 1 - factor ) )
	T/6: factor * T/6 + ( offset/2 * ( 1 - factor ) )
	T/1: factor * T/1
	T/4: factor * T/4
	update-Txy
    ]
    get-average-scale: does [ square-root T/1 * T/4 ]
    inverse: func [ x /local ret ] [
	ret: make x x
	ret/1: x/1 - T/5  / T/1
	ret/2: x/2 - T/6  / T/4
	ret
    ]
    ; x is given in face pixels forward returns image pixels (zoomed
    forward: func [ x /local ret ] [
	ret: make x x
	ret/1: T/1 * x/1 + T/5
	ret/2: T/4 * x/2 + T/6 
	ret
    ]

    T: [1 0 0 1 0 0]
    T-x: [1 0 0 1 0 0 ] ; Scales and translates only in x
    T-y: [1 0 0 1 0 0 ] ; Scales and translates only in y
    old-T: T
    reset: does [
	change T [1 0 0 1 0 0]
	update-Txy
    ]
    faces: []

    update: does [ show faces ]

    dist: func [x y ][ square-root x/1 - y/1 ** 2 + ( x/2 - y/2 ** 2 ) ]

    button: none ; normal or alt
    last-mouse-offset: 0x0
    last-clicked-pos: 0x0
    last-alt-clicked-pos: 1000x1000

    feel-over: func [ face action pos ][
	last-mouse-offset: pos - face/offset
	true
    ]
    feel-engage: func [ face action event ][
	switch action [
	    alt-down [
		last-alt-clicked-pos: event/offset
		button: 'alt
		old-T: copy T
	    ]
	    down [
		last-clicked-pos: event/offset
		button: 'normal
	    ]
	    over [
		switch button [
		    normal [
			translate event/offset - last-clicked-pos
			last-clicked-pos: event/offset
			update
		    ]
		    alt [
			change T old-T
			zoom last-clicked-pos ( dist last-clicked-pos  event/offset ) / max 1 dist last-clicked-pos last-alt-clicked-pos 
			update
		    ]
		]
	    ]
	    up [
	    ]
	]
	false
    ]
    set-feel: func [ o ][
	o/feel/insert-over 'track-pos :feel-over
	o/feel/append-engage 'mouse-zoom :feel-engage
    ]

    test: func [
	/local
	tint-setting
	luma-setting
	contrast-setting
	tr
	image-face
	imm
	lns
	update-bb
	bb
	Z
    ][
	Z: make self []
	tint-setting: 0
	luma-setting: 0
	contrast-setting: 0
	tr: 0x0

	image-face: make face [ 
	    image: im
	    size: im/size
	    color: green
	    effect: [
		luma luma-setting
		contrast contrast-setting
		tint tint-setting
	    ]
	    font: para: edge: feel: none
	]
	show image-face

	update-bb: does [
	    show image-face
	    imm: to-image image-face
	    show bb
	]
	imm: to-image image-face
	lns: []

	Z/zoom 0x0 1500 / im/size/x

	view/options  compose layout [
	    bb: box 1500x500 effect [ draw [ matrix Z/T image 0x0 imm push lns ] ]
		edge [ size: 1x1 color: black ]
		feel [
		    over: get in Z 'feel-over
		    engage: get in Z 'feel-engage
		]
		do [ append Z/faces bb ]
	    key #"Z" [ Z/zoom Z/last-mouse-offset 1 / 1.2 Z/update ]
	    key #"z" [ Z/zoom Z/last-mouse-offset 1.2 Z/update ]
	    key #"0" [ Z/reset Z/update ]
	    key #"c" [ repend lns [ 'pen 'green 'line-width 10 'circle inverse Z/last-mouse-offset 100 ] show bb ]
	    key #"q" [unview]

	    across
	    label "Tint" 80 right slider 400x20 (tint-setting / 512 + 0.5) [ tint-setting: to-integer value - 0.5 * 512  update-bb ]
	    return
	    label "luma" 80 right slider 400x20 (luma-setting / 512 + 0.5) [ luma-setting: to-integer value - 0.5 * 512  update-bb ]
	    return
	    label "Contrast" 80 right slider 400x20 (contrast-setting / 255 + 0.5)  [contrast-setting: to-integer value - 0.5 * 512 update-bb]
	] [all-over]
	Z
    ]
]

pos-obj!: make object! [
    image-pos: [none none]
    window-pos: does [ Z/forward image-pos ]
    world-pos: does [ reduce [	image-pos/1 - xy-origo/1 / pix-per-mm-x
				image-pos/2 - xy-origo/2 / pix-per-mm-y ] ]
    set-image-pos: func [ x ] [ change image-pos reduce [x/1 x/2] ]
    set-window-pos: func [ x ] [ change image-pos Z/inverse reduce [ x/1 x/2] ]
    set-world-pos: func [ x
    ] [
	change image-pos
	    reduce [
		x/1 * pix-per-mm-x + xy-origo/1
		x/2 * pix-per-mm-y + xy-origo/2
	    ]
    ]
]

make object! [
    time: none
    set 'tic does [ time: now/precise/time ]
    set 'toc func[ /name str]  [if name [ print [str now/precise/time - time ] ] now/precise/time - time ]
]

num-to-string: func [ n l /local str ] [   
    str: copy either positive? n [ "" ][ "-" ]
    n: abs n
    if n >= 1 [
	append str to-integer n
	n: mod n 1
    ]
    if l > length? str [ append str #"." ]
    while [ l > length? str ] [ n: n *  10 append str to-integer n  n: mod n 1 ]
    str
] 

form-num: func [
    b [block! number!] {Block of numbers}
    l [integer!] {Max number of chars of each number}
    /local str num
][
    str: copy ""
    if number? b [ b: to-block b ]
    foreach x b [
	append str num-to-string x l
	append str " "
    ]
    clear tail back str
    head str
]


pix-per-mm-x: 2400 / 25.4 * resolution-correction-factor-x
pix-per-mm-y: 2400 / 25.4 * resolution-correction-factor-y
xy-origo: copy [ 0 0]

pitch-offset: 0

Z: make draw-zoom [
    thin-linewidth: 
    normal-linewidth: 
    fat-linewidth: 3 
    get-average-scale-x: does [ square-root T-x/1 * T-x/4 ]
    update: does [
	thin-linewidth: 0.5 / get-average-scale
	normal-linewidth: 2 / get-average-scale
	fat-linewidth: 3 / get-average-scale

	thin-linewidth-x: 0.5 / get-average-scale-x
	normal-linewidth-x: 2 / get-average-scale-x
	fat-linewidth-x: 3 / get-average-scale-x

	show faces
    ]
]

edge-path: make object! [
    x: []
    y: []
    reset: does [ x: head x y: head y ]
    clear: does [ system/words/clear x system/words/clear y ]
    clone: does [ make self [] ]
    cut-part: func [
	{Cuts out the part of the data where
	    the first x data >= start-x and
	    the last x data <= end-x
	 returns a new object
	 }
	 start-x [number!]
	 end-x [number!]
	 /clone {Returns a new object, leaves the original untouched}
	 /local start end new
    ][
	new: either clone [ self/clone ][ self ]
	reset
	while [ x/1 < start-x ] [ x: next x ] start: index? x
	;until [ any [ not x/2  (first+ x) > end-x ] ] end: -1 + index? x
	while [ all [ x/1 x/1 <= end-x ] ] [ x: next x] end: index? x
	reset
	new/x: copy/part at x start at x end
	new/y: copy/part at y start at y end
	new
    ]

    to-pairs: func [ /y-scale ys /local p] [
	p: copy []
	ys: any [ ys 1 ]
	reset
	forall x [
	    append p as-pair (first x) (first+ y) * ys
	]
	y: head y
	p
    ]
    color: black
    linewidth: 'Z/normal-linewidth
    y-scale: 1
    draw: does [
	head append reduce [
		'line-width linewidth 
		'pen color
		'line
	    ] to-pairs/y-scale y-scale
    ]
    select-indexes: func [ ind ][
	x: map-each t ind [ pick x t ]
	y: map-each t ind [ pick y t ]
    ]
]

pick-indexes: func [ b [block!] indexes [block!] ] [ map-each t indexes [ pick b t ] ]

set-image: func [
    img [ image! file!]
    /local image-face name
] [
    if file? img [
	name: to-string img
	lay/text: any [
	    next find/last name "/"
	    name
	]
	lay/changes: 'text show lay
	if error? try [ img: load img ] [ exit ]
	image-bounding-box: reduce [ 0x0 img/size ]
    ]
    sz: im/size
    ;top-line-x: copy []
    ;top-line-y: copy []
    top: make edge-path []
    bot: make edge-path []

    pitch-plot: copy []
    pitch-fit-plot: copy []
    diameter-plot: copy []
    flange-width-plot: copy []
    edge-plot: copy []
    plot: copy []

    image-face: make face [ 
	image: img
	size: img/size
	color: white
	effect: [
	    key white
	]
	font: para: edge: feel: none
    ]
    im: img
    shown-image-gray: to-image image-face
    update-image

    Z/reset 
    Z/zoom 0x0 view-sz/x / im/size/x
    Z/translate as-pair 0 view-sz/y - ( Z/T/1 * im/size/y ) / 2

    hide fix-angle-button

    show lay
]

update-image: func [ ] [
    shown-image-bw: image-threshold shown-image-gray to-integer gray-threshold
    if all [ value? 'bw-image bw-image ] [ show bw-image]
]

find-passes: func [
    im [image!] {The image to search in}
    along [ object!] {Positions can be found by val-at x-pos}
    range [ block! ] {Between these pixels}
    data [object!] {The path to add the top edge to}
    color [tuple!] {The color that makes up the screw}
    /local
	;x last-pix x-to-pair p this-pix
	

][
    x-to-pair: func [ x ] [ as-pair x along/val-at x ]

    entry: none
    last-pix: color = pick im x-to-pair range/1

    for x range/1 + 1 range/2 1 [
	this-pix: color = pick im x-to-pair x
	either last-pix [
	    if  all [ not this-pix entry ] [
		; going out of material, ie from low values to high
		; interpolate between previous and last value
		y1: first pick shown-image-gray x-to-pair x - 1
		y2: first pick shown-image-gray x-to-pair x
		x1: x - 1
		x2: x
		xx:  ( gray-threshold - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1
		comment {
		    y = (y2 - y1 ) / ( x2 - x1 ) * ( x - x1 ) + y1 = c
		    x - x1 = ( c - y1 )/ ( dy / dx )
		    x = x1 + (c - y1 ) * dx / dy
		}

		append data/x entry +  xx / 2
		append data/y xx - entry

		;unless all [ y1 <= gray-threshold y2 >= gray-threshold ] [ print [ "error going out:" y1 y2 ] ]
		;unless all [  xx >= x1 xx <= x2 ] [ print [ "error x-out" x xx ] ]
	    ]
	][
	    if this-pix [
		y1: first pick shown-image-gray x-to-pair x - 1
		y2: first pick shown-image-gray x-to-pair x
		x1: x - 1
		x2: x
		entry: ( gray-threshold - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1
		;unless all [ y1 >= gray-threshold y2 <= gray-threshold ] [ print [ "error going in:" y1 y2 ] ]
		;unless all [  entry >= x1 entry <= x2 ] [ print [ "error x-in" x xx ] ]
	    ]
	]
	last-pix: this-pix
    ]
]

find-edges: func [
    im [image!] {The image to search in}
    top [object!] {The path to add the top edge to}
    bot [object!] {The path to add the bot edge to}
    color [tuple!] {The color that makes up the screw}
    ;/local slice b x
][

    y-size: im/size/y 
    x-size: im/size/x 

    find-non-empty: does [
	c: x + 1
	while [ c < im/size/x ] [
	    ; cut out a slice of the image at a time and search in it.
	    b: find slice: copy/part at im c + 1 as-pair 1 im/size/y  color
	    if b [ 
		top-y: index? b
		bot-y: (1 + length? head b ) - index? find reverse slice color
		break
	    ]
	    c: c + 1
	]
	x: c 
    ]

    search: func [
	{Searches from y (default up)
	 Returns the position of the new color}
	y {start at this pos, it is not checked}
	:cmp [function!] {Search until this is true}
	/up {the default}
	/down {Search downwards}
	/local p
    ] [
	step: either down [ x-size ] [ negate x-size ]
	p:  y + 1 * x-size + x + 1 
	loop y  [
	    if cmp pick im p  [ return to-integer p - 1 / x-size ]
	    p: p + step
	    if p > length? im [ break ]
	]
	none
    ]

    x: 0
    top-y: none
    bot-y: none
    cmp-screw: func [ col ] [ col = color ]
    cmp-noscrew: func [ col ] [ col != color ]
    while [ x < im/size/x ] [
	if  until [
		unless all [ top-y bot-y ] [find-non-empty]
		; ? bot-y
		if x + 1 >= im/size/x [ break/return false ]
		all [ 
		    either cmp-screw pick im top-y * x-size + x + 1
			[ top-y: search/up top-y cmp-noscrew top-y: all [ top-y top-y + 1 ] ]
			[ top-y: search/down top-y cmp-screw ] 

		    either cmp-screw pick im bot-y * x-size + x + 1
			[ 
			    bot-y: search/down bot-y cmp-noscrew 
			    ;probe reduce [ "down"  x bot-y ]
			    bot-y
			] [
			    bot-y: search/up bot-y cmp-screw
			    bot-y: all [ bot-y bot-y + 1 ]
			    ;probe reduce [ "Up  " x bot-y ]
			    bot-y
			] 
		]
	    ]
	[
	    append top/x x
	    append top/y top-y
	    append bot/x x
	    append bot/y bot-y
	]
	x: x + 1
    ]
    not empty? top/x
]

evaluate-screw: func [
    {Initiates all the calculations to evaluate the set screw}
][
    calc-edge-paths
    top/color: magenta  top/linewidth: 'Z/fat-linewidth
    bot/color: magenta bot/linewidth: 'Z/fat-linewidth
    clear plot
    append plot top/draw
    append plot bot/draw
    show [ bw-image fix-angle-button ]
    handle-fit

    inf-area/append-data reduce [
	    last parse/all image-name "/" 
	    image-bounding-box/1 
	    image-bounding-box/2 - image-bounding-box/1
    ]
    inf-area/update
    show inf-area

]

calc-edge-paths: func [
    /local from to  x 
][
    estimate-used-resolution: func [ pix ][
	either auto-scale-checkbox/data [
	    res: possible-resolutions
	    len-err: map-each x res [ abs ( approx-screw-length - (pix / x * 25.4) ) ]
	    resolution-butt/data: at head resolution-butt/data index? minimum-of len-err
	    show resolution-butt
	    pick res index? minimum-of len-err
	][
	    to-integer first parse resolution-butt/text none
	]
    ]

    top/clear
    bot/clear

    ;find-edges shown-image-gray top bot green
    find-edges shown-image-bw top bot green

    ;identifgy the edges of the bar by searching from the middle
    ; use the top-line
    x: top/x
    x: at head x (length? x) / 2
    while [ all [ x/-1 x/1 - 1 = x/-1 ] ] [ x: back x ]
    from: x/1
    clear range-selection/pos-list
    append range-selection/pos-list make pos-obj! [ set-image-pos reduce [ from 0 ] ]
    xy-origo: reduce [ from pick top/y index? x ]
    x: at head x (length? x) / 2
    while [ all [ x/2 x/1 + 1 = x/2 ] ] [ x: next x ]
    to: x/1
    append range-selection/pos-list make pos-obj! [ set-image-pos reduce [ to 0 ] ]

    pix-per-mm-x: ( image-dpi: estimate-used-resolution to - from ) / 25.4 * resolution-correction-factor-x
    pix-per-mm-y:  image-dpi / 25.4 * resolution-correction-factor-y
    
]

evaluate-profile: func [
    {Calculates the widht of the flanges at a fix distance from the centerline.  }
    /local e 
	s ;top-flange-width bot-flange-width
	faults
][
    top-middle-est: e: make mean-est [ x-vals: y-vals: none ]
    e/k/1: e/k/1 - ( effective-diameter / 2 * pix-per-mm-y )

    bot-middle-est: e: make mean-est [ x-vals: y-vals: none ]
    e/k/1: e/k/1 + ( effective-diameter / 2 * pix-per-mm-y )
    
    find-passes 
	shown-image-bw
	top-middle-est
	reduce [ analysis-range/1/image-pos/1  analysis-range/2/image-pos/1 ]
	top-flange-width: make edge-path [ ]
	green

    find-passes
	shown-image-bw
	bot-middle-est
	reduce [ analysis-range/1/image-pos/1  analysis-range/2/image-pos/1 ]
	bot-flange-width: make edge-path [ ]
	green

    top-flange-width: make top-flange-width [ y: map-each t y [ t / pix-per-mm-x ] ]
    bot-flange-width: make bot-flange-width [ y: map-each t y [ t / pix-per-mm-x ] ]

    faults: remove-outliers top-flange-width 0.5 0.02
    faults: faults + remove-outliers bot-flange-width 0.5 0.02

    flange-width: append copy top-flange-width/y bot-flange-width/y

    either 1 < length? flange-width [
	flange-width-mean: ( vector-sum flange-width ) / (length? flange-width ) 
	flange-width-std: if 1 < length? flange-width [ 
	    s: 0
	    foreach t flange-width [ s: t - flange-width-mean ** 2 + s ]
	    square-root s  / ( -1 + length? flange-width )
	]
	flange-width-mean-str: to-dec form-num flange-width-mean  5
	flange-width-std-str:  to-dec form-num flange-width-std 5 
    ] [
	flange-width-mean: "N/A"
	flange-width-std: "N/A"
	flange-width-mean-str: "N/A"
	flange-width-std-str: "N/A"
    ]

    ; Transform the data to scale to the window
    ; min 0.10 max 0.20
    box-height: diagram-box/size/y
    top-flange-width: make top-flange-width [
	color: red
	y-scale: 1
	y: map-each t y [ box-height - ( t - 0.10 / 0.15 * box-height ) ]
	;repeat i length? y [
	    ;poke y i  box-height - ( w - 0.10 / 0.10 * box-height )
	;]
    ]
    bot-flange-width: make bot-flange-width [
	color: blue
	y-scale: 1
	y: map-each t y [ box-height - ( t - 0.1 / 0.15 * box-height ) ]
    ]
    flange-width-plot: top-flange-width/draw
    append flange-width-plot bot-flange-width/draw

    show diagram-box

    either faults * 50 > ( ( length? top-flange-width/y ) + length? bot-flange-width/y ) [
	unless running-batch [
	    alert reform [ "Found" faults "outliers found in profile. Please do manual checking!" ]
	]
	std-dev-prefix: "-"
    ] [ std-dev-prefix: "" ]

    inf-area/append-data  reduce [
	flange-width-mean-str
	rejoin [  std-dev-prefix flange-width-std-str ]
    ]
]

to-dec: func [x][ to-decimal trim x]

remove-outliers: func [
    {Removes data points that are unexpected. Selects the faults by comparing the data to a low pass filtered value.
    retruns the number of faults removed.}
    o [object!] {A edge-path object}
    liteness [number!] {How much of the actual value is used when comparing}
    limit [ number! ] {The limit of acceptable outlayers}
    /local filt-y
	i 
	faults
][
    ; Remove thread items differing more than max-thread-diff from its low passed filtereed
    ; value. Run filt-filt
    filt-y: filt-filt o/y liteness
    i: 1
    faults: 0
    while [ i <= length? o/y ] [ 
	;either ( abs o/y/:i - first+ filt-y ) > ( max-thread-diff * pix-per-mm-y )[
	either ( abs o/y/:i - first+ filt-y ) > limit [
	    remove at o/y i
	    remove at o/x i
	    faults: faults + 1
	][
	    i: i + 1
	]
    ]
    faults
]

handle-fit: func [
    /local x y
	p
	est-path resid-path curv-path
	o
	sel-top sel-bot
	top-filt bot-filt
	top-extemas bot-extremas
	est-obj
	;top-inner-est top-outer-est bot-inner-est bot-outer-est
	;mean-est
	curvyness
	d1 d2

] [
    p: range-selection/pos-list
    if all [ shown-image-bw/size/x > 10 shown-image-bw/size/y > 10 2 = length? p ] [
	clear edge-plot

	sel-top: top/cut-part/clone p/1/image-pos/1 p/2/image-pos/1
	sel-bot: bot/cut-part/clone p/1/image-pos/1 p/2/image-pos/1


	inf-area/data-return
	inf-area/append-data reduce [
	    to-dec form-num (pick p/2/world-pos 1 ) - (pick p/1/world-pos 1) 6
	]
	inf-area/update
	show [ edge-plot-box inf-area pitch-box ]

	sel-bot/color: gray sel-bot/linewidth: 'Z/normal-linewidth
	sel-top/color: gray sel-top/linewidth: 'Z/normal-linewidth

	top-filt: make sel-top [ y: thread-filter sel-top/y to-integer pix-per-mm-x * ideal-pitch ]
	bot-filt: make sel-bot [ y: thread-filter sel-bot/y to-integer pix-per-mm-x * ideal-pitch ]

	top-extremas: calc-extrema-index top-filt/y
	bot-extremas: calc-extrema-index bot-filt/y

	top-outer: make sel-top [ select-indexes top-extremas/minima-index ]
	top-inner: make sel-top [ select-indexes top-extremas/maxima-index ]
	bot-outer: make sel-bot [ select-indexes bot-extremas/maxima-index ]
	bot-inner: make sel-bot [ select-indexes bot-extremas/minima-index ]

	foreach o [ top-outer top-inner bot-outer bot-inner] [
	    remove-outliers get o 0.4 max-thread-diff * pix-per-mm-y
	]

	est-obj: make object! [
	    range: reduce [
		first minimum-of reduce [ top-outer/x/1 bot-outer/x/1 top-inner/x/1 bot-inner/x/1 ]
		first maximum-of reduce [ last top-outer/x last bot-outer/x last top-inner/x last bot-inner/x ]
	    ]
	]

	top-inner-est: trig-fit/range top-inner/x top-inner/y polynom-order est-obj/range
	top-outer-est: trig-fit/range top-outer/x top-outer/y polynom-order est-obj/range
	bot-inner-est: trig-fit/range bot-inner/x bot-inner/y polynom-order est-obj/range
	bot-outer-est: trig-fit/range bot-outer/x bot-outer/y polynom-order est-obj/range

	make sel-top [ y: top-inner-est/estimate x: top-inner/x linewidth: 'Z/normal-linewidth color: forest  append edge-plot draw]
	make sel-top [ y: top-outer-est/estimate x: top-outer/x linewidth: 'Z/normal-linewidth color: forest  append edge-plot draw]
	make sel-bot [ y: bot-inner-est/estimate x: bot-inner/x linewidth: 'Z/normal-linewidth color: forest  append edge-plot draw]
	make sel-bot [ y: bot-outer-est/estimate x: bot-outer/x linewidth: 'Z/normal-linewidth color: forest  append edge-plot draw]

	show  edge-plot-box

	start-pos: make p/1 [ set-world-pos reduce [ skip-end-statistical-average 0] ]
	end-pos: make p/2 [ set-world-pos reduce [ (pick world-pos 1) - skip-end-statistical-average 0] ]
	analysis-range: reduce [ start-pos end-pos ]

	outer-est: make bot-outer-est [
	    k: add-mult-vector-vector k 1.0 top-outer-est/k -1.0

	    cut-part start-pos/image-pos/1 end-pos/image-pos/1
	    create-bases
	    y-vals: none
	]
	inner-est: make bot-inner-est [
	    k: add-mult-vector-vector k 1.0 top-inner-est/k -1.0
	    cut-part start-pos/image-pos/1 end-pos/image-pos/1
	    create-bases
	    y-vals: none
	]
	clear diameter-plot

	mean-est: make top-outer-est [ 
	    k: add-mult-vector-vector k 0.25 top-inner-est/k +0.25
	    k: add-mult-vector-vector k 1 bot-inner-est/k +0.25
	    k: add-mult-vector-vector k 1 bot-outer-est/k +0.25
	    ;k: mult-vector-scalar k 0.25
	    make edge-path [
		foreach [ x1 x2 x3 x4 ] top-outer/x [ append x x1 ]
		foreach [ y1 y2 y3 y4 ] estimate [ append y y1 ]
		unless ( last x ) = last top-outer/x [
		    append x last top-outer/x 
		    append y last estimate
		]
		color: violet
		linewidth: 'Z/fat-linewidth
		set 'tmp self
		append edge-plot draw 
	    ]
	]
	
	curvyness: make mean-est [
	    k/1: 0 k/2: 0 ; remove offset and slope
	]
	inf-area/append-data reduce [
	    to-dec form-num screw-rotation: mean-est/angle 4

	    to-dec form-num outer-est/mean / pix-per-mm-y 5 
	    to-dec form-num (first tmp: outer-est/extremas ) / pix-per-mm-y 5 
	    to-dec form-num tmp/2  / pix-per-mm-y 5 
	]

	inf-area/append-data reduce [
	    to-dec form-num inner-est/mean / pix-per-mm-y 5
	    to-dec form-num (first tmp: inner-est/extremas )  / pix-per-mm-y 5
	    to-dec form-num tmp/2 / pix-per-mm-y 5
	]

	inf-area/append-data reduce [
	    to-dec  form-num ( curvyness/integrate-estimate/fun func[x][ abs x ] ) / (pix-per-mm-x * pix-per-mm-y) 4
	    to-dec  form-num (first tmp: curvyness/extremas) / pix-per-mm-y  4 
	    to-dec  form-num tmp/2 / pix-per-mm-y 4
	]
	inf-area/update
	show inf-area

	d1: make edge-path [ 
	    x: outer-est/x-vals
	    y: mult-vector-scalar outer-est/estimate 1 / pix-per-mm-y ; to mm
	    y-min: first minimum-of y
	    y-max: first maximum-of y
	    y-diff: y-max - y-min
	    if y-diff = 0 [ y-diff: 0.5 ]
	    y: add-vector-scalar y negate y-min ; -0.5  ; show range 0.5 - 1.0
	    y: mult-vector-scalar y diagram-box/size/y / ( y-min - y-max ) ;-0.5  ; scale 
	    y: add-vector-scalar y diagram-box/size/y  
	    color: magenta
	    line-width: 'Z/normal-linewidth-x
	    append diameter-plot draw
	]

	d2: make edge-path [
	    x: inner-est/x-vals
	    y: mult-vector-scalar inner-est/estimate 1 / pix-per-mm-y ; to mm
	    y-min: first minimum-of y
	    y-max: first maximum-of y
	    y-diff: y-max - y-min
	    if y-diff = 0 [ y-diff: 0.5 ]
	    y: add-vector-scalar y negate y-min ; -0.5  ; show range 0.5 - 1.0
	    y: mult-vector-scalar y diagram-box/size/y / ( y-min - y-max ) ;-0.5  ; scale 
	    y: add-vector-scalar y diagram-box/size/y  
	    color: green
	    line-width: 2
	    append diameter-plot draw
	]
	show diagram-box

	update-pitch-offset

	calc-true-pitch

	evaluate-profile

    ]
]

calc-extrema-index: func [ 
    {Find the (local) extreme points of y and returns the indexes where it has change direction
     first will always be a top.
     Returns  a object with top-index and bot-index
    }
    y [block!]
][
    make object! [
	maxima-index: copy []
	minima-index: copy []
	repeat i -1 + length? y [
	    either ( length? maxima-index ) = length? minima-index [
		if (pick y i + 1) < pick y i [
		    append maxima-index i
		]
	    ][
		if (pick y i + 1) > pick y i [
		    append minima-index i
		]
	    ]
	]
	; Remove first and last hits
	;remove head maxima-index
	;remove back tail maxima-index
	;remove head minima-index
	;remove back tail minima-index
    ]
]

calc-pitch: func [ 
    {Find the extreme points of pth and return a path with x-values taken mod pitch}
    pth [object!]
    pitch  [decimal!]
    offset [number!]
    /y y-vals [block!] {Use these for finding with /top /bot otherwise the y-values in pth}
    /local x-top x-valey xx top bot
][
    ; find peeks
    x-top: copy [] x-valey: copy []
    fun: none
    if any[ top bot ] [
	fun-top: func [ ][
	]
    ]
    y-vals: any [ y-vals pth/y ]

    top: copy []
    bot: copy []
    x-bot: copy []
    repeat i -1 + length? pth/x [
	either ( length? x-top ) = length? x-valey [
	    if pth/y/(i + 1) < pth/y/:i [
		append x-top pth/x/:i
		append top y-vals/:i
	    ]
	][
	    if pth/y/(i + 1) > pth/y/:i [
		append x-valey pth/x/:i + (pitch / 2 * pix-per-mm-x )
		append x-bot pth/x/:i
		append bot y-vals/:i
	    ]
	]
    ]
    make edge-path [
	x: sort append copy x-top x-valey
	y: map-each t x [ (mod t / pix-per-mm-y - offset  pitch) / pitch ]
	y-max: top
	x-max: x-top
	y-min: bot
	x-min: x-bot
    ]
]

unfold-pitch: func [ 
    {A vector that has been built with mod pitch will be be connected at the limits 0 and pitch}
    pitch [number!] {the assumed pitch}
    y [block!] {The values built with mod pitch}
    /local ret low-p new
    {Current solution adds one pitch if difference to previous is larger than half a pitch (hence the difference will be reduced)
     There is a problem with this in the following case with this series of pitch values
     [ 0 0.51 0.1 ] (pitch 1) 
     which should be [ 0 -0.49 0.1 ] but becomes [ 0 -0.49 -0.9 ]
     the reason is that the first step is larger than pitch/2 but the second is not
     A possibility would be not to compare with the last only, but with a low pass version of the calculated values
    }
     

][
    ret: copy [ 0 ]
    low-p: 0
    y: next y
    forall y [
	diff: (first y )- first back y
	new: (last ret)
	    + diff
	    + case [
		diff > (pitch / 2)	[ negate pitch ]
		diff < (pitch / -2 )	[ pitch ]
		true			[ 0 ]
	    ]
	while [ ( low-p - (new - pitch) ) < ( new - low-p) ][ new: new - pitch ]
	while [ ( low-p - (new + pitch) ) > ( new - low-p) ][ new: new + pitch ]
	append ret new
	low-p: 2 * low-p + new / 3
    ]
    y: head y
    ret
]

update-pitch-offset: func [ /local y-size ] [
    y-size: pitch-box/size/y
    top-outer-pitch: make top-outer [
	y: map-each t x [ t / pix-per-mm-x + pitch-offset // ideal-pitch ]
	y-scale: y-size / ideal-pitch 
	color: red
    ]
    top-inner-pitch: make top-inner [ 
	y: map-each t x [ t / pix-per-mm-x + ( ideal-pitch / 2) + pitch-offset // ideal-pitch ]
	y-scale: y-size / ideal-pitch
	color: brick
    ]
    bot-inner-pitch: make bot-inner [
	y: map-each t x [ t / pix-per-mm-x + pitch-offset // ideal-pitch ]
	y-scale: y-size / ideal-pitch
	color: sky
    ]
    bot-outer-pitch: make bot-outer [
	y: map-each t x [ t / pix-per-mm-x + ( ideal-pitch / 2) + pitch-offset // ideal-pitch ]
	y-scale: y-size / ideal-pitch
	color: blue
    ]

    pitch-plot: top-outer-pitch/draw
    append pitch-plot top-inner-pitch/draw
    append pitch-plot bot-inner-pitch/draw
    append pitch-plot bot-outer-pitch/draw

    if value? 'pitch-fit [
	pitch-fit-draw: make edge-path [
	    x: pitch-fit/x-vals
	    y: add-vector-scalar pitch-fit/estimate  (ideal-pitch / 2) + pitch-offset
	    y-scale: pitch-box/size/y / ideal-pitch
	    color: green
	]

	pitch-fit-plot: pitch-fit-draw/draw
    ]

    show [ pitch-box]
]

calc-true-pitch: func [
    /local
	data-x data-y data data2
	start-pos end-pos
][
    ;start-pos: make range-selection/pos-list/1 [ set-world-pos reduce [ skip-end-calc-pitch 0] ]
    ;end-pos:   make range-selection/pos-list/2 [ set-world-pos reduce [ (pick world-pos 1) - skip-end-calc-pitch 0] ]

    top-outer-pitch/cut-part analysis-range/1/image-pos/1 analysis-range/2/image-pos/1 ; start-pos/image-pos/1 end-pos/image-pos/1
    top-inner-pitch/cut-part analysis-range/1/image-pos/1 analysis-range/2/image-pos/1 ; start-pos/image-pos/1 end-pos/image-pos/1
    bot-inner-pitch/cut-part analysis-range/1/image-pos/1 analysis-range/2/image-pos/1 ; start-pos/image-pos/1 end-pos/image-pos/1
    bot-outer-pitch/cut-part analysis-range/1/image-pos/1 analysis-range/2/image-pos/1 ; start-pos/image-pos/1 end-pos/image-pos/1

    data-x: copy top-outer-pitch/x 
    data-y: unfold-pitch ideal-pitch top-outer-pitch/y

    append data-x top-inner-pitch/x
    append data-y unfold-pitch ideal-pitch top-inner-pitch/y

    append data-x bot-inner-pitch/x
    append data-y unfold-pitch ideal-pitch bot-inner-pitch/y

    append data-x bot-outer-pitch/x
    append data-y unfold-pitch ideal-pitch bot-outer-pitch/y 

    data: transpose-matrix reduce [ data-x data-y ]
    data: sort data 
    clear data-x clear data-y
    for i 10 length? data 10 [
	append data-x
	    data/(i - 9)/1 +
	    data/(i - 8)/1 +
	    data/(i - 7)/1 +
	    data/(i - 6)/1 +
	    data/(i - 5)/1 +
	    data/(i - 4)/1 +
	    data/(i - 3)/1 +
	    data/(i - 2)/1 +
	    data/(i - 1)/1 +
	    data/(i)/1 / 10
	append data-y
	    data/(i - 9)/2 +
	    data/(i - 8)/2 +
	    data/(i - 7)/2 +
	    data/(i - 6)/2 +
	    data/(i - 5)/2 +
	    data/(i - 4)/2 +
	    data/(i - 3)/2 +
	    data/(i - 2)/2 +
	    data/(i - 1)/2 +
	    data/(i)/2 / 10
    ]


    pitch-cal: trig-fit data-x data-y 0


    {Vad betyder lutningen?
     Lutningen * L / ideal-pitch är det antal var som inte får plats på stången
     så antalet varv är L / ideal-pitch - L * lutning / ideal-pitch = L * ideal-pitch ( 1 - lutning)
     Och den exakta pitchen ideal-pitch / ( 1 - lutning)
     och det kan väl lika väl vara
     ideal-pitch ( 1 + lutning)
    }
    calculated-pitch: ideal-pitch * ( 1 + ( pitch-cal/k/2 * pix-per-mm-x ) )
    pitch-info/text:  form-num calculated-pitch 7
    inf-area/append-data form-num calculated-pitch 6 


    pitch-fit: trig-fit pitch-cal/x-vals pitch-cal/residual polynom-order
    pitch-fit-draw: make edge-path [
	x: pitch-fit/x-vals
	y: add-vector-scalar pitch-fit/estimate  (ideal-pitch / 2)
	y-scale: pitch-box/size/y / ideal-pitch
	color: green
    ]

    comment[
	use [ dy ][
	    dy: copy []
	    y: next copy/part pitch-fit/estimate ( length? pitch-fit/y-vals ) - 20
	    forall y [ append dy y/1 - y/-1 ]
	    pitch-extreme: (first maximum-of dy ) - (first minimum-of dy )
	]
	dbg: [ pitch-fit pitch-extreme dy ]
    ]
    pitch-extremes: map-each t copy/part
	skip pitch-fit/estimate-deriv 10
	(length? pitch-fit/y-vals ) - 20
	[
	    calculated-pitch * t * pix-per-mm-x
	]
    pitch-extreme: (first maximum-of pitch-extremes ) - first minimum-of pitch-extremes

    ;inf-area/append-data to-dec form-num max negate pitch-extremes/1 pitch-extremes/2 5
    inf-area/append-data to-dec form-num pitch-extreme 5
    inf-area/update

    pitch-fit-plot: pitch-fit-draw/draw

    show [ inf-area pitch-info pitch-box ]
]


handle-pitch: func [ /local x y mean ][
    top-filt: make top [ y: thread-filter top/y to-integer pix-per-mm-x * ideal-pitch ]
    top-filt/cut-part 
		    range-selection/pos-list/1/image-pos/1
		    range-selection/pos-list/2/image-pos/1
    bot-filt: make bot [ y: thread-filter bot/y to-integer pix-per-mm-x * ideal-pitch ]
    bot-filt/cut-part 
		    range-selection/pos-list/1/image-pos/1
		    range-selection/pos-list/2/image-pos/1

    top-extremas: calc-extrema-index top-filt/y
    bot-extremas: calc-extrema-index bot-filt/y

    update-pitch-offset
]

running-av: func [ x y start-x end-x l /local ret-x ret-y ] [
    until [ ( first+ x ) >= start-x ] start-i: ( index? x ) - 1 
    until [ ( first+ x ) >= end-x ] end-i: ( index? x) - 1
    x: head x

    av: func [ x l /local sum] [
	x: skip x l / -2
	sum: 0
	loop l [ sum: sum + first+ x ]
	sum / l 
    ]
    ret-x: copy []
    ret-y: copy []
    for i start-i end-i 1 [ append ret-x x/:i
	append ret-y av at y i l
    ]
    reduce [ ret-x ret-y]
]

shown-image-bw: shown-image-gray: make image! 1x1

view-sz: 1600x200
residual-sz: view-sz / 1x2
update-image

mkfun: func [ str /local fun][
    if  error? try [
	fun: func[x] load to-block str
	1 * fun 0
    ][
	alert "Function des not evaluate correctly"
	return none
    ]
    :fun
]
create-ideal: func [
    /local
	savename-field length-field diam-field pitch-field
	angle-field deflection-field 
][
    view/new layout [
	style tt text 200 right
	across
	tt "Name of image file no-extension" savename-field: field 200 
	return
	tt "Length of screw"  length-field: field 100	"100"
	return
	tt "Outer diameter expression in x" diam-field: field 100 "0.92"
	return
	tt "pitch [mm]" pitch-field: field 100 "0.25"
	return
	tt "Angle"	angle-field: field 100 "0"
	return
	tt "Deflection"	deflection-field: field 100 "0"
	return
	button "Save"  [
	    save/png to-rebol-file join savename-field/text ".png"
		make-screw-image
		    to-decimal length-field/text
		    mkfun diam-field/text
		    to-decimal pitch-field/text
		    to-decimal angle-field/text
		    mkfun deflection-field/text
	]
	button "Close" [unview]
    ]
]

lay: layout [
    styles tools

    space 0 across
    at 0x0 button "Open" #"^o" [
	image-name: request-file/filter/only ["*.png" "*.jpg"]
	all [
	    image-name
	    exists? image-name
	    set-image image-name
	]
	running-batch: false
    ]
    button "Batch analysis" #"^b" [ 
	files: request-file/title/filter "Select files to process" "Done" "*.png"
	files: any [ files [] ]
	running-batch: true
	foreach x files [
	    image-name: x
	    ;print [ "Analysis of file" image-name ]
	    set-image image-name
	    evaluate-screw
	    wait 0.1
	]
    ]
    button "To pdf" #"^p" [ write %test.pdf face-to-pdf lay ]
    ; button "Ideal screw" [ create-ideal ]
    button "Help" [ view/new help-layout ]

    space 20x20 below
    origin 20x30
    orig: at
    im-image: box view-sz
	effect [ draw [ line-cap 1 matrix Z/T image shown-image-gray pen forest box 0x0 shown-image-gray/size ] ]
	white
	do [ append Z/faces im-image ]
    at orig
    bw-image: box view-sz
	effect [ draw [ matrix Z/T image shown-image-bw pen push plot ] key white ]
	edge [ size: 3 color: black ]
	do [ append Z/faces bw-image ]
    at orig
    edge-plot-box: box view-sz effect [ draw [ matrix Z/T push edge-plot ] ]
	do [ append Z/faces edge-plot-box ]
    at orig
    event-draw-box: event-box view-sz
	effect [ grid 100x0 0x0 blue draw [
	    line-cap 1
	    push range-selection/draw
	]
    ] feel [
	    transformer: make pos-obj! []
	    insert-over 'pos-update
		func [ face action pos /local f ][
		    transformer/set-window-pos pos - face/offset
		    pos-field/text: form-num transformer/world-pos 6
		    show pos-field
		    true
		]
	]
	do [ append Z/faces event-draw-box ]
    tmp: at
    at orig
    panel 200.200.200 [ space 0 across pos-field: text 100 right text "mm" ]
    at tmp
    pitch-box: panel white residual-sz effect [ grid 100x50 0x0 draw [
	    line-cap 1
	    matrix Z/T-x
	    push pitch-plot
	    push pitch-fit-plot
    ]]
    feel [
	p: p-offset: none
	engage: func [f a e][
	    switch a [
		down [
		    p: e/offset/y
		    p-offset: pitch-offset
		]
		over [
		    pitch-offset: ( e/offset/y - p ) / f/size/y * ideal-pitch + p-offset
		    update-pitch-offset
		    show  pitch-box
		]
	    ]
	]
	over: func [ f e p /local x][
	    if all [ value? 'pitch-fit pitch-fit ] [
		p: make pos-obj! [ set-window-pos p - f/offset ]
		x: first p/image-pos
		pitch-box-text/text: rejoin [
		    "local pitch " form-num calculated-pitch * ( 1 + ( (pitch-fit/deriv-at x)* pix-per-mm-x)) 8
		]
		show pitch-box-text
	    ]
	]
   ] [
	    at 0x0 pitch-box-text: text "Exact local pitch xx mm dfads" with [ color: 200.200.200 ]
    ]
    do [ append Z/faces pitch-box ]
    diagram-box: panel white residual-sz effect [ grid 100x20 0x0 draw [
		line-cap 1
		matrix Z/T-x
		push flange-width-plot
		push diameter-plot
	]]
	feel [
	    over: func [ f e p /local x][
		p: make pos-obj! [ set-window-pos p - f/offset]
		x: first p/image-pos
		if all [ value? 'outer-est outer-est ] [
		    diameter-text/text: reform [
			;"do=" form-num (lin-interp outer-est/x-vals outer-est/estimate x ) / pix-per-mm-y 4
			"do=" form-num (outer-est/val-at x ) / pix-per-mm-y 4
			"di=" form-num (inner-est/val-at x ) / pix-per-mm-y 4
			"mm"
		    ]
		    show diameter-text
		]
	    ]
	]
	    [   ;at 0x-4 text leaf "1.0" with [ color: white]
		;at 0x30 text leaf "0.80" with [ color: white]
		;at 0x70 text leaf "0.60" with [ color: white]
		at 0x0 diameter-text: text "do=______^/di=______ mm" with [ color: 200.200.200.128 ] wrap 60x30
		at 1565x0 panel [
		    at 0x10 text 35 "0.22" right with [ color: white ]
		    at 0x30 text 35 "0.19" right with [ color: white ]
		    at 0x50 text 35 "0.16" right with [ color: white ]
		    at 0x70 text 35 "0.13" right with [ color: white ]
		] with [ color: none ]
	    ]
	do [ append Z/faces diagram-box ]
    across

    grouped-panel  "Coordinate systems" [
	space 20x10 origin 10x30 across
	zero-button: select-image-area "Set zero"
	    with [
		draw-box: 'event-draw-box
		transformation: pos-obj!
		up: func [ face event pos ][
		    xy-origo: pos/2/image-pos
		]
	    ]
	resolution-butt: rotary "4800 dpi" [
	    pix-per-mm-x: (to-decimal first parse value "") / 25.4 * resolution-correction-factor-x
	    pix-per-mm-y: (to-decimal first parse value "") / 25.4 * resolution-correction-factor-y
	] 
	return
	text "Auto scale detection" auto-scale-checkbox: check on 15x15
	return
	info with [
	    text: unless all [ resolution-correction-factor-x = 1 resolution-correction-factor-y = 1 ] [
		reform [ "Scale correction" resolution-correction-factor-x resolution-correction-factor-y ]
	    ]
	] ;edge [ size: 1x2 effect: 'ibevel ]
    ]
    grouped-panel "Level of material" [
	space 20x10 origin 10x30 across
	;show-grayscale-checkbox: check 15x15 [ show-grayscale?: value update-image ]
	gray-level-slider: slider gray-threshold / 255 250x15 [
	    gray-level-field/text: to-string gray-threshold: round 255 * value 
	    show gray-level-field
	    update-image
	]
	gray-level-field: field to-string round gray-threshold 40 [
	    gray-threshold: to-integer value
	    gray-level-slider/data: gray-threshold / 255
	    show gray-level-slider
	    update-image
	]

	return
	button "Calc edge" #"^t" [
		;material-color: black
		if all [ shown-image-bw/size/x > 10 shown-image-bw/size/y > 10 ] [ evaluate-screw ]
	]
	panel [
	    across
	    range-selection: select-image-area "Select range" "wait..."  #"^r"
		with [
		    draw-box: 'event-draw-box
		    transformation: pos-obj!
		    up: :handle-fit
		]
	    ;ord-field: field 35 [polynom-order: to-integer value handle-fit] with [ text: to-integer polynom-order]
	]
    ]
	;return
    across
    grouped-panel 120 "Manipulate image" [
	origin 10x50 space 20x5
	;crop-button: select-image-area "Crop" #"C" 90
	    ;with [
		;butt-press: does [ set-info "Select the two corners that surround the cropped image" ]
		;draw-box: 'event-draw-box
		;transformation: pos-obj!
		;up: func [ face event pos ][
		    ;im: crop-image im pos
		    ;set-image im
		    ;set-info none
		    ;show [ lay ]
		;]
		;down: does [ set-info "Hold and release at other corner" ]
	    ;]

	rotate90-button: button "Rotate 90" #"R" 90
	    [ set-image rot90-image im show lay ]
	fix-angle-button: button "Fix angle" #"F" 90 [ 
	    set-image rot-image im negate screw-rotation
	    show [lay]
	] 
	origin 10x10
    ]
    grouped-panel "Pitch" [
	origin 20x30 space 20
	panel [ across space 0
	    text "Desired pitch" underline return
	    ideal-pitch-field: field to-string ideal-pitch 65 [
		    ideal-pitch: do face/text
		    handle-pitch
		]
	    panel [ origin 0x0
		at 0x0 arrow  up [ ideal-pitch: ideal-pitch + 0.0005 ideal-pitch-field/text: to-string ideal-pitch 
				    show ideal-pitch-field handle-pitch ]
		with [ size: as-pair 10 ideal-pitch-field/size/y + ideal-pitch-field/edge/size/y / 2 ]
		at 0x10 arrow down [ ideal-pitch: ideal-pitch - 0.0005 ideal-pitch-field/text: to-string ideal-pitch 
				    show ideal-pitch-field handle-pitch ]
		with [ size: as-pair 10 ideal-pitch-field/size/y + ideal-pitch-field/edge/size/y / 2 ]
	    ]
	]
	space 20x5
	text "Calculated pitch" underline
	pitch-info: info "" 80 [
	    ideal-pitch-field/text: value
	    show ideal-pitch-field
	    ideal-pitch: do value handle-pitch
	]
    ]
    inf-area: table  600x180 data [[]]
	left width 80
	header [ left ]
	left-col new-line/all [ left bold  140
	    "Length [mm]"
	    "Angle [deg]"
	    "Outer diam-Mean [mm]"
	    "Outer diam min [mm]"
	    "Outer diam max [mm]"
	    "Inner diam Mean [mm]"
	    "Inner diam min [mm]"
	    "Inner diam max [mm]"
	    "Area to midline [mm2]"
	    "Deflection min [mm]"
	    "Deflection max [mm]"
	    "Calc. pitch [mm]"
	    "Calc. pm [mm]"
	    "Flange mean [mm]"
	    "Flange std [mm]"
	    "Image name"
	    "UL [pix]"
	    "size [pix]"
	] off
	with [
	    append init [
		header/cell/feel: make header/cell/feel [
		    detect: func [ face event ] [
			if  event/type = 'alt-up [
			    ; find out what position in the header this has
			    ; call delete column for that position 
			    i: index? find head face/parent-face/pane face
			    cntxt-menue face [
				item "Delete column" [ remove-col i unview]
				item "Change head text" [
				    header/data/1/:i: any [ request-text/title "New heading name" header/data/1/:i ]
				    show header
				    unview
				]
			    ]
			    ;print [ "bring up the menue for column" i ]
			    ;remove-col i
			    ; later add menue with delete as one option and edit possibilities another
			]
			event
		    ]
		]
	    ]
	]


    key #"^Z" [ Z/zoom Z/last-mouse-offset 1 / 1.2 Z/update ]
    key #"^z" [ Z/zoom Z/last-mouse-offset 1.2 Z/update ]
    key #"^x" [
	Z/reset
	Z/zoom 0x0 view-sz/x / im/size/x
	Z/translate as-pair 0 view-sz/y - ( Z/T/1 * im/size/y ) / 2
	Z/update
    ]
    key #"^C" [
	write clipboard:// inf-area/to-string*/transpose
    ]
    leave-key
    return
    indent -20
    info-bar: text 640 "Press any button" edge [size: 1x2 color: white effect: 'bevel]
	do [ set-info: func [ txt ][ info-bar/text: txt ] ]
]

help-layout: layout [
    title "Help for analysing the screws"
    area 400x500
	{Shortcut keys:
	#"^^t" Start calculating the edges, pitch and diameters
	#"^^o"  Open a new image
	#"^^b"  Open images and start calculation
	#"^^R"  Reopen the last opened image
	#"^^z" Zoom into the image with the cursor as center
	#"^^Z" Zoom out of the image
	#"^^x" Fit zoom
	#"^^c" Copies the data table to clipboard to be pasted to a text file

Draging inte the main window translates the image.
Right-clicking and dragging zooms around the last clicked point.

Vertical drag in the pitch image translates the pitch.

Line coloring scheeme:
    Top	    Red
    Bot	    Blue
    Averaged	thicker, behind
    inner   magenta/lighter
    outer   green/harder
}
    key #"q" [unview]
]

resolution-butt/texts: resolution-butt/data:  map-each x possible-resolutions [ reform [ x "dpi" ] ]

lay/size/y: info-bar/offset/y + info-bar/size/y
info-bar/size/x: lay/size/x 



Z/set-feel event-draw-box
set-image im

view/options/new  lay [all-over]

hide fix-angle-button
do-events

tool-butt-example: does [
    view layout [
	styles tools
	b: rec-button with [
	    draw-box: 'c
	    up: func [ f pos-list ] [
		print mold pos-list
	    ]
	]

	c: event-box 400x400 red  effect[ draw [ push b/draw push e/draw sel/draw ]]

	e: rec-button "make box" "Wait ..." with [ draw-box: 'c up: func[f picked-list][ print picked-list ] ]

	sel: select-image-area "in image" with [
	    draw-box: 'c
	    transformation: pos-obj!
	    pen: green
	    up: func [face e pos-list ][ print mold pos-list ]
	]

	leave-key
    ]
]

; vim: ai sw=4 sts=4
