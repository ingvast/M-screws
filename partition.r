#!/usr/bin/rebview  -sv
REBOL [
    title: {Program for partition an image of screws}
    author: {Johan Ingvast}
    copyright: "Bioservo Technologies AB 2015"
    encap: [ title "Partition - " ]
]

; ************************************
; Debugging code


lin: func [ x /local r ][ r: copy [] forall x [ append r index? x ] r]


;**************************************


prebol-used?: not (#do ['false])

; In linux, switch to arial as default font
if system/version/4 = 4 [ 
    font-sans-serif: "arial"
    face/font/name: font-sans-serif
    foreach [n f ] svv/vid-styles
    [
	if all[ f/font f/font/name = "helvetica" ] [ f/font/name: font-sans-serif ]
    ]
]

; Read libraries

either prebol-used? [
    #include %../../source/view.r
    #include  %libs/im-lib.r
][
    do %libs/im-lib.r
]

; Partition
; Analyses an image and cuts out the individual screws. Each screw must fit into an orthogonal box

level: 0.7 

image-scale: 8

cc-space: charset " ^-"
cc-number: charset [ #"0" - #"9" ]

find-cuts: func [
    histvect
    level
    /local  ;negative-cross positive-cross
][
    positive-cross: copy []
    negative-cross: copy []

    cuts: copy []

    histvect: next head histvect
    forall histvect [
	if all [
	    histvect/1 >= level 
	    histvect/-1 < level
	] [ append positive-cross index? histvect ]
	if all [
	    histvect/1 <= level 
	    histvect/-1 > level
	] [ append negative-cross index? histvect ]
    ]
    histvect: head histvect
    if histvect/1 > level
	[ insert positive-cross 0 ]
    append negative-cross len-across - 1

    for i 1 length? positive-cross 1 [
	append cuts
	    ( positive-cross/:i +
		any [ negative-cross/:i positive-cross/:i ]
	    ) / 2
    ]
    cuts
]

find-cut-length: func [
    vector
    /local ;len
][
    len: length? vector
    val-at-screw: pick vector 0.4 * len
    ; the values starts at low value, so the left edge should be just after when value has nearly reached val-at-screw.
    left-index: len * 0.4
    for i left-index 0. -1 [
	if val-at-screw * 0.9 > pick vector i [
	    left-index: i  - ( len * 0.003 )
	    break
	]
    ]
    right-index: len * 0.4
    for i right-index len * 1.0 1 [
	if val-at-screw * 0.9 > pick vector i [
	    right-index: i + ( len * 0.003 )
	    break
	]
    ]
    reduce [ left-index right-index ]
]

analysis-image: none

make-histogram: func [
    {Sums up the rows of the image, result is 0 - 1. 0 is black and 1 is white}
    im [ image! ]  {The image to be analyzed. Only the red color is analyzed. The image should be gray}
    /vertical {The screws are oriented vertically}

    /local sz step-along step-across ;len-along ;len-across
	histvect
][
    comment [
	analysis-image: im
	orig-image: im

	image-scale: any [ sc 1 ]
	if scale [ analysis-image: scale-image analysis-image 1 / sc ]
    ]

    sz: either vertical [
	reverse analysis-image/size
    ][
	analysis-image/size
    ]

    ; Make a histogram of the image in the direction selected
    len-along: first sz
    len-across: second sz

    either vertical [
	step-along: sz/y
	step-across: 1
    ][
	step-along: 1
	step-across: sz/x
    ]

    histvect: make block! len-across
    start: 1
    loop len-across [
	this: start
	sum: 0
	loop len-along [
	    sum: sum + first pick analysis-image this
	    this: this + step-along
	]
	append histvect  sum / len-along / 256.
	start: start + step-across
    ]
    histvect
]

save-image: func [
    i [ integer! ] {Index in image, starting at 1}
    name [string! file!] {Name of image to save (without postfix}
    /local 
	;p-image
	y
][
    y: sort copy/part at cuts i 2
    p-image: copy/part		at orig-image as-pair
					max 0 image-scale * cut-length/1
					image-scale * y/1
				as-pair
					((min cut-length/2 orig-image/size/x) - cut-length/1 ) * image-scale
					image-scale * ( y/2 + 1 - y/1 )
    save/png rejoin [ directory name ".png" ] p-image 
]

align-right: func [ str "String to put it into" val ] [
    unless string? val [ val: to-string val ]
    head change skip str (length? str) - length? val  val
]

save-all: func [
][
    repeat i -1 + length? cuts [
	save-image
	    i
	    to-file rejoin [ trim basename-field/text align-right "000" i - 1 + to-integer start-number-field/text ]
    ]
]

;pos-draw:   copy []
;neg-draw:   copy []
cuts-draw:  copy []
cuts-length-draw: copy []


view-image: func[ /local basename ] 
[
    basename: to-string copy/part file-name find/last file-name #"."
    view/new/title lay: layout [
	pos: at
	image analysis-image
	at pos
	line-box: box analysis-image/size with [ color: none ]
	    effect [
		draw [
		    push cuts-draw
		]
		draw [
		    push cuts-length-draw
		]
	    ]

	across
	label "Level" 
	slider 500x15 level [ find-cuts histvect level: value show-image ]
	return
	return
	label "Basename" basename-field: field 70 basename
	label "Start-number"
	start-number-field: field "1" 30 [
	    unless parse/all value [ any cc-space some cc-number any cc-space ] [
		face/text: copy "1"
		show face
	    ]
	    show-image
	]
	button "Reverse (r)" #"^r" [ reverse cuts show-image ]
	button "Rotate 90 (R)" #"^R" [
	    unview
	    orig-image: rot90-image orig-image
	    analyse-image orig-image
	    view-image
	    do-events
	]
	button "Save all (S)" #"^S" [ save-all ]
	button "Quit (q)" #"^q" [ unview ]
    ]  to-string file-name
    show-image
]

show-image: does [
    clear cuts-draw
    clear cuts-length-draw

    foreach p cuts [
	repend cuts-draw [
	    'pen 'green 'push reduce [ 'line as-pair 0 p as-pair len-across p ] ]
    ]

    foreach p cut-length [
	repend cuts-length-draw [
	    'pen 'blue 'push reduce [ 'line as-pair p 0 as-pair p len-along ] ]
    ]

    line-box/pane: copy []
    move-line-handle: make face [
	size: 50x20
	number: none
	color: white
	down: none
	dir: 'y ; direction of motion
	cut: none
	cut-draw: none
	change-line: func [ f offset ] [
	    f/cut/(f/number): f/cut-draw/(f/number * 4)/2/(f/dir): f/cut-draw/(f/number * 4)/3/(f/dir): offset
	]
	feel: make feel [
	    engage: func [ f a e][
		switch a [
		    down [ f/down: e/offset ]
		    over away [
			f/offset: f/offset + e/offset - f/down
			f/change-line f  f/offset/(f/dir) + (f/size/(f/dir) / 2) 
			show [ f line-box ]
		    ]
		]
	    ]
       ]
    ]
    n: 1
    for i 1 length? cuts 1 [
	; create boxes that can be put into the line-box/pane
	; moving the box should also move the line
	append  line-box/pane make move-line-handle [
	    dir: 'y
	    cut: cuts 
	    cut-draw: cuts-draw
	    offset: as-pair n * 20 cuts/:i - ( size/y / 2 )
	    text: n - 1 + to-integer start-number-field/text
	    number: n
	]
	n: n + 1
    ]
    for i 1 length? cut-length 1 [
	; create boxes that can be put into the line-box/pane
	; moving the box should also move the line
	append  line-box/pane make move-line-handle [
	    dir: 'x
	    cut: cut-length
	    cut-draw: cuts-length-draw
	    offset: as-pair cut/:i - ( size/x / 2 ) n * 20
	    text: pick [ "left" "right" ] i 
	    number: i
	]
	n: n + 1
    ]
    show lay
]

load-image: func [
    {Loads the image and makes the histogram and a proposal of how to partition it}
    f-name
][
    orig-image: load f-name
    analyse-image orig-image
]

analyse-image: func [ im ][
    analysis-image: scale-image im 1 / image-scale
    histvect: make-histogram analysis-image
    histvect: filt-filt histvect 0.15
    lengthvect: make-histogram/vertical analysis-image
    cuts: find-cuts histvect level
    cut-length: find-cut-length lengthvect
]

until [
    file-names: request-file/filter/title/path [ "*.jpg" "*.png" ] "Select image to partition" "Open"
]

directory: file-names/1
first+ file-names

either 1 < length? file-names  [
    auto: on
] [
    auto: off
]

foreach f-name file-names [
    file-name: f-name
    num: none
    parse/all file-name [ any [ [ 4 cc-number any cc-number ] | [ copy num 1 4 cc-number to end ] | skip ] ]
	
    load-image join directory file-name
    view-image

    comment [
	make-histogram/scale load join directory file-name image-scale
	find-cuts level
    ]

    if num [ start-number-field/text: num show-image ]
    if auto [
	save-all
	unview
    ]
]
either auto [quit][ do-events ]

; vim: ai sw=4 sts=4
