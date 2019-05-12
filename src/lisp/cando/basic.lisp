(in-package :cando)

(defmacro |:=| (a b)
  `(defparameter ,a ,b))

