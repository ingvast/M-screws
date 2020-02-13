
; Here we set variables for im-analysis.r

; A list of possible resolutions (dpi) to select from
possible-resolutions: [ 300 600 1200 2400 3200 4800 ] ; dpi

; Corrects deviations from the exact dpi setting. 
; A value larger than 1 means the resolution is larger than the dpi setting
;resolution-correction-factor-x: 0.998
;resolution-correction-factor-y: 0.998
resolution-correction-factor-x: 1
resolution-correction-factor-y: 1

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
