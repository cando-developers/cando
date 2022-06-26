(defpackage #:cando-user-install
  (:use #:common-lisp)
  (:export #:install
           #:*snapshot-timestamp*
           #:update
           #:updatep))