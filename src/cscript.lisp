(k:recurse #P"adapt/"
           #P"geom/"
           #P"units/"
           #P"chem/"
           #P"kinematics/"
           #P"lisp/"
           #P"main/"
           #P"data/"
           #P"examples/")

(k:library "expat" :required t :min-version "2.2.0")

(k:sources :install-code
           #~"lisp/"
           #~"data/")
