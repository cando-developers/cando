
(require :asdf)
(push (translate-logical-pathname #P"cando:cando;") asdf:*central-registry*)
(asdf:load-system :cando)

(defparameter *s* #(1 2 3 4))
(nreverse *s*)
*s*
