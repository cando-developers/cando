(defpackage #:cando-widgets
  (:use #:common-lisp)
  (:nicknames :cw)
  (:local-nicknames (:j :jupyter)
                    (:jw :jupyter-widgets))
  (:export
    #:add-page
    #:container
    #:make-file-task-page
    #:make-molecule-select
    #:make-simple-task-page
    #:molecules
    #:molecule-select
    #:page
    #:selected
    #:sketch-molecule
    #:sketch-molecules))

