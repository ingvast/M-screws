Hur vill jag att det ska fungera?
ska jag anropa med en funktion och argumenten 

For this to work on linux system make sure that the window handler detects the
name of a window that ends with prop:menue and unsets all the frames and tiltes
on it.

REBOL [
    title: "Context enues for VID"
    copyright: "Bioservo Technologies AB"
    author: "Johan Ingvast"
    version: {$ revision: $}
    example:  none
]


os-name: select [
    4 linux
    3 windows
    2 apple
    7 freeBSD
    9 openBSD
] fourth system/version 
;print ["Check the os. " os-name "is selected" ]


; Adding a function to always get the mouse position within the window 
remove at system/view/screen-face/feel/event-funcs 2

append system/view/screen-face/feel/event-funcs func [ face event ][
    if event/type = 'move [ last-window-cursor-offset: event/offset ]
    event
]

view-menue: func [ 
    lay
    xy
][
    switch os-name [
	linux [
	    lay/text: "prop:menue"
	    lay/offset: xy
	    view/new l
	]
	windows [
	    view/new/offset/options/title lay xy [ no-title no-border]  ""
	]
    ]
]

test-view-menue: does [
    l: layout/tight [ text "lalal" text "Johan"  text "Close" [unview] ]
    view-menue l 800x800
    do-events
]

menue-styles: stylize [
    item: text para [ origin: 5x2 ]
    separator: bar 
    header: text underline
]

layout-menue: func [
    definition
    /local
	lay
	label
	fun
	l
][
    lay: copy [below]
    append lay definition
    
    l: layout/tight/styles lay  menue-styles 
    l/feel: make l/feel [
	detect: func [f e][
	    ;print e/type
	    if e/type = 'inactive [ unview ]
	    e
	]
    ]
    l/edge: make l/edge [ size: 1x1 color: black ]
    l/size/y: l/size/y + 10
    l
]

cntxt-menue: func [ parent-face def /local xy ][
    if any [
	none? last-window-cursor-offset
	last-window-cursor-offset = 32767x32767
    ] [ exit ]
    xy: -10x-10
	+ last-window-cursor-offset
	+ get in find-window parent-face 'offset
    view-menue l: layout-menue def xy
]


test-menue: does [
    view layout [
	button "Menue" [ unview ] [
	    cntxt: cntxt-menue face [
		header "Testmenue"
		separator
		item "close the menue" [ unview ] 
		item "quit" [ unview/all ] 
		item "Introspection" [ introspect face ] 
	    ]
	]
	key #"q" [unview]
    ]
]


; vim: sw=4 sts=4 ai textwidth=0 
