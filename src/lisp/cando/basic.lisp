(in-package :cando)

(defmacro |:=| (a b)
  `(defparameter ,a ,b))

(defmacro |<-| (a b)
  `(defparameter ,a ,b))

