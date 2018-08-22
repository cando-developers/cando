(defsystem "language.smiles"
  :version    (:read-file-form "version-string.sexp")
  :depends-on ("alexandria"

               (:version "esrap"                         "0.16")

               (:version "parser.common-rules"           "0.4")

               (:version "architecture.builder-protocol" "0.9"))

  :components ((:module    "parser"
                :pathname  "src/smiles/parser"
                :serial    t
                :components ((:file       "package")
                             (:file       "variables")
                             (:file       "protocol")
                             (:file       "grammar"))))

  :in-order-to ((test-op (test-op "language.smiles/test"))))

(defsystem "language.smiles/playground"
  :depends-on ("language.smiles"

               "architecture.builder-protocol.print-tree"))

(defsystem "language.smiles/test"
  :depends-on ("alexandria"

               "fiveam"

               "language.smiles")

  :components ((:module     "test"
                :pathname   "test"
                :serial     t
                :components ((:file       "package")))

               (:module     "parser"
                :depends-on ("test")
                :pathname   "test/smiles/parser"
                :serial     t
                :components ((:file       "package")
                             (:file       "grammar" ))))

  :perform    (test-op (operation component)
                (uiop:symbol-call '#:language.smiles.test '#:run-tests)))
