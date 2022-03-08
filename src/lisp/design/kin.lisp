(in-package :kin)

(defclass plug ()
  ((plug-name :initarg :plug-name :accessor plug-name)
   (smarts :initarg :smarts :initform nil :accessor smarts)
   (cap :initarg :cap :initform nil :accessor cap)))

(defclass out-plug (plug) ())

(defclass in-plug (plug) ())


(defun make-out-plug (name &key smarts cap)
  (make-instance 'out-plug :plug-name name :smarts smarts :cap cap))

(defun make-in-plug (name &key smarts cap)
  (make-instance 'in-plug :plug-name name :smarts smarts :cap cap))

(eval-when (:compile-toplevel :execute :load-toplevel)
  (export '(out-plug
            plug-name
            smarts
            cap
            in-plug
            make-out-plug make-in-plug)))
