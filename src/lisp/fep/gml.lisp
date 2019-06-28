;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package #:fep-gml)

(defun render-node (stream id name &key (color "white") (shape "rectangle"))
  (format stream "  node~%")
  (format stream "  [~%")
  (format stream "    id ~a~%" id)
  (format stream "    label ~s~%" (format nil "~a:~a" name id))
  (format stream "    graphics~%")
  (format stream "    [~%")
  (format stream "      type ~s~%" (string-downcase (string shape)))
  (format stream "      fill ~s~%" (string-downcase (string color)))
  (format stream "    ]~%")
  (format stream "  ]~%"))

(defun render-morph (stream source target label &key type)
  (when source
    (format stream "  edge~%")
    (format stream "  [~%")
    (format stream "    source ~a~%" source)
    (format stream "    target ~a~%" target)
    (format stream "    label ~s~%" label)
    (format stream "    graphics~%")
    (format stream "    [~%")
    (format stream "      width 2~%")
    (format stream "      type ~s~%" "line")
    (format stream "    ]~%")
    (format stream "  ]~%")))



