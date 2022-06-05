


(defparameter x (units:seconds 1))
(defparameter y (units:meters 5))

(bformat t "\n")
(bformat t "x = %s\n" x)
(bformat t "y = %s\n" y)
(bformat t "x*y = %s\n" (units:mul x y))
(bformat t "x/y = %s\n" (units:div x y))


(defparameter v1 (units:div (units:meters 5) (units:seconds)))
(defparameter v2 (units:div (units:meters 10) (units:seconds)))
(defparameter v3 (units:mul v2 v2))

(bformat t "v1 = %s\n" v1)
(bformat t "v2 = %s\n" v2)
(bformat t "v3 = %s\n" v3)
(bformat t "reciprocal v3 = %s\n" (units:reciprocal v3))
(bformat t "sqrt v3 = %s\n" (units:sqrt v3))
(bformat t "v1+v2 = %s\n" (units:+ v1 v2))
(bformat t "v1*v2 = %s\n" (units:* v1 v2))
(bformat t "v1/v2 = %s\n" (units:/ v1 v2))
(bformat t "v1-v2 = %s\n" (units:- v1 v2))
(bformat t "v1^2 = %s\n" (units:pow v1 2))
(bformat t "v1*2.0 = %s\n" (units:mul-by-scalar v1 2.0))
(bformat t "v1/2.0 = %s\n" (units:div-by-scalar v1 2.0))
(bformat t "(amount v1/v2) = %s\n" (units:amount (units:div v1 v2)))

(bformat t "(amount v1*(m^-1*s)) = %s\n" (units:amount (units:mul v1 (units:Quantity 1.0 :meters -1 :seconds))))
(bformat t "kB = %s\n" units:kB)
(bformat t "R = %s\n" units:R)
(bformat t "R in kcal = %s\n" (units:energy-to-kcal units:R))
