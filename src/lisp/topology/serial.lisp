
(in-package #:serial)

(defclass serializable () ())

(defmethod print-object ((obj serializable) stream)
  (if *print-readably*
      (cando.serialize:print-object-readably-with-slots obj stream) ; in subclass: (call-next-method)
      (print-unreadable-object (obj stream :type t)
        ;; in subclass: whatever you want here
                                )))

(defmethod conspack:encode-object append
    ((obj serializable) &key &allow-other-keys)
  (loop for slotd in (clos:class-slots (class-of obj))
        for slot-name = (clos:slot-definition-name slotd)
        when (slot-boundp obj slot-name)
          collect (cons slot-name (slot-value obj slot-name))))

(defmethod conspack:decode-object-initialize progn
    ((object serializable) class alist &key &allow-other-keys)
  (declare (ignore class))
  (loop for (name . value) in alist
        do (setf (slot-value object name) value)))

